import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from lab.discovery import DiscoveryCombo, DiscoveryResult
from lab.exp_disco_degrade import run_disco_degrade

FAKE_METRICS = {'weight_coverage': 0.5, 'skipprob': 0.1,
                 'duration_coverage': 0.6, 'alignment_coverage': 0.7}


def degrade_stub(log, level):
    if level == 0.0:
        return log, set()
    return f'{log}_degraded', {f'dropped_at_{level}'}


class ZeroLevelDedupTest(unittest.TestCase):
    # Level 0.0 drops nothing regardless of dimension, so it's the same
    # (log, tree) computation under every dim - run_disco_degrade should
    # call compute_metrics for it once, not once per dimension.

    def setUp(self):
        self.tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        self.combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult('FAKE_TREE'))}
        self.degradations = {'activity': degrade_stub, 'trace': degrade_stub}

    def test_compute_metrics_called_once_for_level_zero(self):
        with patch('lab.exp_disco_degrade.compute_metrics',
                   return_value=dict(FAKE_METRICS)) as mock_cm:
            with patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'):
                run_disco_degrade(['fake_log.xes'], combos=self.combos,
                                   degradations=self.degradations, levels=[0.0, 0.5],
                                   out_csv=str(self.tmp_out))
            # 1 call for the shared level-0.0 point + 2 dims x 1 non-zero level
            self.assertEqual(mock_cm.call_count, 3)

    def test_both_dims_get_a_row_at_level_zero_with_matching_metrics(self):
        with patch('lab.exp_disco_degrade.compute_metrics',
                   return_value=dict(FAKE_METRICS)):
            with patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'):
                df = run_disco_degrade(['fake_log.xes'], combos=self.combos,
                                        degradations=self.degradations, levels=[0.0],
                                        out_csv=str(self.tmp_out))
            self.assertEqual(len(df), 2)
            self.assertEqual(set(df['degradation_dim']), {'activity', 'trace'})
            for _, row in df.iterrows():
                self.assertEqual(row['status'], 'ok')
                self.assertEqual(row['weight_coverage'], FAKE_METRICS['weight_coverage'])
                self.assertEqual(row['alignment_coverage'], FAKE_METRICS['alignment_coverage'])

    def test_nonzero_levels_still_computed_per_dim_per_level(self):
        with patch('lab.exp_disco_degrade.compute_metrics',
                   return_value=dict(FAKE_METRICS)) as mock_cm:
            with patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'):
                df = run_disco_degrade(['fake_log.xes'], combos=self.combos,
                                        degradations=self.degradations, levels=[0.5, 1.0],
                                        out_csv=str(self.tmp_out))
            # no level-0.0 shared call here, so it's a plain 2 dims x 2 levels
            self.assertEqual(mock_cm.call_count, 4)
            self.assertEqual(len(df), 4)


class NotImplementedComboTest(unittest.TestCase):
    # A combo whose discover() raises NotImplementedError should be
    # skipped entirely - no compute_metrics call at all, including at
    # level 0.0, and every row reports status='not_implemented'.

    def test_skipped_without_computing_anything(self):
        def not_implemented(log):
            raise NotImplementedError('no discovery yet')

        combos = {'stub': DiscoveryCombo('stub', not_implemented)}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        with patch('lab.exp_disco_degrade.compute_metrics') as mock_cm:
            with patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'):
                df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                        degradations=degradations, levels=[0.0, 0.5],
                                        out_csv=str(tmp_out))
            mock_cm.assert_not_called()
            self.assertTrue((df['status'] == 'not_implemented').all())
            self.assertTrue(df['weight_coverage'].isna().all())


class ComputeMetricsErrorTest(unittest.TestCase):
    # A compute_metrics failure (for the shared level-0.0 call or a
    # regular per-level call) should be recorded as an error row rather
    # than crashing the sweep, and should not be mistaken for 'ok'.

    def test_level_zero_error_is_recorded_and_not_crashing(self):
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult('FAKE_TREE'))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        with patch('lab.exp_disco_degrade.compute_metrics',
                   side_effect=RuntimeError('boom')):
            with patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'):
                df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                        degradations=degradations, levels=[0.0],
                                        out_csv=str(tmp_out))
            self.assertEqual(len(df), 1)
            self.assertIn('error: boom', df.iloc[0]['status'])
            self.assertIsNone(df.iloc[0]['weight_coverage'])

    def test_nonzero_level_error_does_not_block_other_levels(self):
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult('FAKE_TREE'))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        def flaky(log, tree, slpn_path, ppt_weights=None):
            if log == 'FAKE_LOG':
                return dict(FAKE_METRICS)
            raise RuntimeError('boom')

        with patch('lab.exp_disco_degrade.compute_metrics', side_effect=flaky):
            with patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'):
                df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                        degradations=degradations, levels=[0.0, 0.5],
                                        out_csv=str(tmp_out))
            self.assertEqual(len(df), 2)
            zero_row = df[df['degradation_level'] == 0.0].iloc[0]
            nonzero_row = df[df['degradation_level'] == 0.5].iloc[0]
            self.assertEqual(zero_row['status'], 'ok')
            self.assertIn('error: boom', nonzero_row['status'])


if __name__ == '__main__':
    unittest.main()
