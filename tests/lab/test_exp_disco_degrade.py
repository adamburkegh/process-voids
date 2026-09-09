import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

from skipalignments import Activity, Tau, Xor

from lab.discovery import DiscoveryCombo, DiscoveryResult
from lab.exp_disco_degrade import (
    run_disco_degrade, main, _node_rows, CLASSICAL_METRIC_KEYS, PER_NODE_METRIC_KEYS,
)
from process_voids.coveragemass import TREE_METRIC_KEYS

FAKE_METRICS = {'weight_coverage': 0.5, 'weight_voidage': 0.5, 'skipprob': 0.1,
                 'salign_coverage': 0.7}
FAKE_CLASSICAL_METRICS = {k: 0.42 for k in CLASSICAL_METRIC_KEYS}
FAKE_DV = object()  # compute_metrics(..., return_dv=True)'s second value - opaque here,
                    # since _classical_metrics itself is mocked in every test below


def degrade_stub(log, level):
    if level == 0.0:
        return log, set()
    return f'{log}_degraded', {f'dropped_at_{level}'}


# Every test below also mocks build_id_net (a cheap conversion in
# reality, but 'FAKE_TREE' isn't a real ProcessTree so it can't run for
# real here) and _classical_metrics itself (bypasses voidmass_table_pn/
# coverage_by_alignment_pn - the expensive classical-alignment
# computation, out of scope for these orchestration tests - see
# tests/process_voids/test_voidmass_pn_prototype.py for that).


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
                   return_value=(dict(FAKE_METRICS), FAKE_DV)) as mock_cm, \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {})), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            run_disco_degrade(['fake_log.xes'], combos=self.combos,
                               degradations=self.degradations, levels=[0.0, 0.5],
                               out_csv=str(self.tmp_out))
            # 1 call for the shared level-0.0 point + 2 dims x 1 non-zero level
            self.assertEqual(mock_cm.call_count, 3)

    def test_both_dims_get_a_row_at_level_zero_with_matching_metrics(self):
        with patch('lab.exp_disco_degrade.compute_metrics',
                   return_value=(dict(FAKE_METRICS), FAKE_DV)), \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {})), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=self.combos,
                                    degradations=self.degradations, levels=[0.0],
                                    out_csv=str(self.tmp_out))
            self.assertEqual(len(df), 2)
            self.assertEqual(set(df['degradation_dim']), {'activity', 'trace'})
            for _, row in df.iterrows():
                self.assertEqual(row['status'], 'ok')
                self.assertEqual(row['weight_coverage'], FAKE_METRICS['weight_coverage'])
                self.assertEqual(row['salign_coverage'], FAKE_METRICS['salign_coverage'])
                self.assertEqual(row['voidmass_deficit'], FAKE_CLASSICAL_METRICS['voidmass_deficit'])

    def test_nonzero_levels_still_computed_per_dim_per_level(self):
        with patch('lab.exp_disco_degrade.compute_metrics',
                   return_value=(dict(FAKE_METRICS), FAKE_DV)) as mock_cm, \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {})), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=self.combos,
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

        with patch('lab.exp_disco_degrade.compute_metrics') as mock_cm, \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {})), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                    degradations=degradations, levels=[0.0, 0.5],
                                    out_csv=str(tmp_out))
            mock_cm.assert_not_called()
            self.assertTrue((df['status'] == 'not_implemented').all())
            self.assertTrue(df['weight_coverage'].isna().all())
            self.assertTrue(df['voidmass_deficit'].isna().all())


class ComputeMetricsErrorTest(unittest.TestCase):
    # A compute_metrics failure (for the shared level-0.0 call or a
    # regular per-level call) should be recorded as an error row rather
    # than crashing the sweep, and should not be mistaken for 'ok'.

    def test_level_zero_error_is_recorded_and_not_crashing(self):
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult('FAKE_TREE'))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        with patch('lab.exp_disco_degrade.compute_metrics',
                   side_effect=RuntimeError('boom')), \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {})), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                    degradations=degradations, levels=[0.0],
                                    out_csv=str(tmp_out))
            self.assertEqual(len(df), 1)
            self.assertIn('RuntimeError: boom', df.iloc[0]['status'])
            self.assertIsNone(df.iloc[0]['weight_coverage'])
            self.assertIsNone(df.iloc[0]['voidmass_deficit'])

    def test_nonzero_level_error_does_not_block_other_levels(self):
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult('FAKE_TREE'))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        def flaky(log, tree, slpn_path, ppt_weights=None, return_dv=False):
            if log == 'FAKE_LOG':
                return dict(FAKE_METRICS), FAKE_DV
            raise RuntimeError('boom')

        with patch('lab.exp_disco_degrade.compute_metrics', side_effect=flaky), \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {})), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                    degradations=degradations, levels=[0.0, 0.5],
                                    out_csv=str(tmp_out))
            self.assertEqual(len(df), 2)
            zero_row = df[df['degradation_level'] == 0.0].iloc[0]
            nonzero_row = df[df['degradation_level'] == 0.5].iloc[0]
            self.assertEqual(zero_row['status'], 'ok')
            self.assertIn('RuntimeError: boom', nonzero_row['status'])


class NodeRowsTest(unittest.TestCase):
    """
    _node_rows against a small REAL tree (Xor(Tau, a)) - mass_by_weight/
    voidage_by_weight/mandatory_node_count/total_node_count run for
    real; coverage_by_alignment/coverage_by_alignment_pn are mocked
    since their own correctness is covered elsewhere (test_coveragemass,
    test_voidmass_pn_prototype) - this test is about the row-assembly
    schema (one row per vm_table node, right ids, right columns), not
    re-deriving the metric math.
    """

    def setUp(self):
        self.a = Activity(None, 'a', 100000)
        self.a.id = '1'
        self.tau = Tau(None, 'tau', 0)
        self.tau.id = '2'
        self.choice = Xor(None, [self.tau, self.a])
        self.choice.id = '3'
        self.a.set_parent(self.choice)
        self.tau.set_parent(self.choice)
        # mass_by_weight/infer_operator_weights need .weight already set
        # on every node (including Tau) - normally done by
        # transfer_pt_weights inside compute_metrics before _node_rows
        # ever sees the tree; set directly here since this test bypasses
        # that pipeline.
        self.a.weight = 1
        self.tau.weight = 0
        self.choice.weight = 1

        class FakeDv:
            def __init__(self, skip_probs):
                self.skip_probs = skip_probs

        self.dv = FakeDv({self.choice: 0.1, self.a: 0.2, self.tau: 0.0})
        self.vm_table = {
            self.choice: {'deficit': 1.0, 'movecount': 2.0,
                           'voidmass_subprocess': 0.5, 'voidmass_process': 0.5},
            self.a: {'deficit': 0.5, 'movecount': 1.0,
                     'voidmass_subprocess': 0.5, 'voidmass_process': 0.5},
            self.tau: {'deficit': 0.0, 'movecount': 0.0,
                       'voidmass_subprocess': 0.0, 'voidmass_process': 1.0},
        }

    def test_one_row_per_vm_table_node(self):
        with patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=0.77), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.88):
            rows = _node_rows('mylog', 'mycombo', 'trace', 0.5, self.dv, self.vm_table)
        self.assertEqual({r['node_id'] for r in rows}, {'1', '2', '3'})

    def test_row_carries_cell_identity_and_node_identity(self):
        with patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=0.77), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.88):
            rows = _node_rows('mylog', 'mycombo', 'trace', 0.5, self.dv, self.vm_table)
        row = next(r for r in rows if r['node_id'] == '1')
        self.assertEqual(row['log'], 'mylog')
        self.assertEqual(row['combo'], 'mycombo')
        self.assertEqual(row['degradation_dim'], 'trace')
        self.assertEqual(row['degradation_level'], 0.5)
        self.assertEqual(row['node_type'], 'Activity')
        self.assertEqual(row['alphabet'], 'a')

    def test_row_has_every_metric_column_with_correct_values(self):
        with patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=0.77), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.88):
            rows = _node_rows('mylog', 'mycombo', 'trace', 0.5, self.dv, self.vm_table)
        row = next(r for r in rows if r['node_id'] == '1')

        for key in PER_NODE_METRIC_KEYS + CLASSICAL_METRIC_KEYS + TREE_METRIC_KEYS:
            self.assertIn(key, row)

        self.assertEqual(row['skipprob'], 0.2)
        self.assertEqual(row['voidmass_deficit'], 0.5)
        self.assertEqual(row['voidmass_movecount'], 1.0)
        self.assertEqual(row['salign_coverage'], 0.77)
        self.assertEqual(row['alignment_coverage_pn'], 0.88)
        # a's own subtree is just itself, no silent alternative from its
        # own perspective (mandatory_node_count/total_node_count are
        # evaluated AT that node, not from an outside ancestor's view).
        self.assertEqual(row['total_node_count'], 1)

    def test_tau_node_gets_a_row_too(self):
        with patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=0.0), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.0):
            rows = _node_rows('mylog', 'mycombo', 'trace', 0.5, self.dv, self.vm_table)
        tau_row = next(r for r in rows if r['node_id'] == '2')
        self.assertEqual(tau_row['node_type'], 'Tau')
        # Tau is excluded from mandatory/total counts by definition.
        self.assertEqual(tau_row['mandatory_node_count'], 0)
        self.assertEqual(tau_row['total_node_count'], 0)


class DryRunTest(unittest.TestCase):
    """--dry-run prints the resolved Experiment and exits without
    calling compute_metrics or run_disco_degrade at all - main()'s
    describe-then-return path must never fall through into computing
    anything."""

    def test_dry_run_with_run_name_prints_and_computes_nothing(self):
        with patch('sys.argv', ['exp_disco_degrade', '--run', 'smoke', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('lab.exp_disco_degrade.compute_metrics') as mock_cm, \
             patch('lab.exp_disco_degrade.build_id_net') as mock_build_id_net:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_cm.assert_not_called()
        mock_build_id_net.assert_not_called()
        self.assertIn('Experiment: smoke', output)
        self.assertIn('rtfm', output)
        self.assertIn('inductive', output)
        self.assertIn('cells:', output)

    def test_dry_run_with_ad_hoc_logs_prints_and_computes_nothing(self):
        with patch('sys.argv', ['exp_disco_degrade', 'fake_log.xes',
                                 '--combos', 'inductive', '--levels', '0.0', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('lab.exp_disco_degrade.compute_metrics') as mock_cm, \
             patch('lab.exp_disco_degrade.build_id_net') as mock_build_id_net:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_cm.assert_not_called()
        mock_build_id_net.assert_not_called()
        self.assertIn('Experiment: ad hoc', output)
        self.assertIn('fake_log', output)

    def test_dry_run_accepts_claims_fixture_combo_and_degradation_names(self):
        """--combos/--degradations must also resolve claims_known/
        appeal_seq etc - these aren't in ALL_COMBOS/ALL_DEGRADATIONS
        (real-discovery only), they come from lab.claims_fixture's
        registrations, merged in only for this ad hoc lookup."""
        with patch('sys.argv', ['exp_disco_degrade', 'data/claims.xes',
                                 '--combos', 'claims_known',
                                 '--degradations', 'appeal_seq', 'loop_block',
                                 '--levels', '0.0', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('lab.exp_disco_degrade.compute_metrics') as mock_cm, \
             patch('lab.exp_disco_degrade.build_id_net') as mock_build_id_net:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_cm.assert_not_called()
        mock_build_id_net.assert_not_called()
        self.assertIn('claims_known', output)
        self.assertIn('appeal_seq', output)
        self.assertIn('loop_block', output)

    def test_unknown_combo_name_errors_out(self):
        with patch('sys.argv', ['exp_disco_degrade', 'fake_log.xes',
                                 '--combos', 'not_a_real_combo', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'):
            with self.assertRaises(SystemExit):
                main()


if __name__ == '__main__':
    unittest.main()
