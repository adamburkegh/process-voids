import io
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

import pandas as pd

from lab.exp_disco_degrade import run_disco_degrade, main, ALL_METRICS
from lab.params import ALL_COMBOS, ALL_DEGRADATIONS, ALL_LEVELS


class RunDiscoDegradeDelegatesToLabRunTest(unittest.TestCase):
    """
    run_disco_degrade is a thin wrapper over lab.run.run(), preserving
    exp_disco_degrade's original signature, defaults and output
    filenames - the cell loop/metric roster is lab.run's own, already
    covered by tests/lab/test_run.py (including a parity test proving
    the two produce identical output for identical inputs, and the
    metrics-exclusion behaviour --exclude-metric relies on). This file
    only needs to check the wrapping itself, plus this module's own CLI.
    """

    def test_delegates_to_lab_run_with_the_given_arguments(self):
        fake_df = pd.DataFrame([{'status': 'ok'}])
        fake_node_df = pd.DataFrame([{'node_id': '1'}])
        fake_timings_df = pd.DataFrame([{'metric_or_stage': 'dv'}])

        with patch('lab.exp_disco_degrade.lab_run.run',
                   return_value=(fake_df, fake_node_df, fake_timings_df)) as mock_run:
            result = run_disco_degrade(
                ['fake_log.xes'], combos={'c': 'COMBO'}, degradations={'d': 'DEGRADE'},
                levels=[0.0, 0.5], metrics=['SOME_METRIC'], out_csv='out.csv',
                node_out_csv='out_nodes.csv', timings_out_csv='out_timings.csv')

        mock_run.assert_called_once_with(
            ['fake_log.xes'], combos={'c': 'COMBO'}, degradations={'d': 'DEGRADE'},
            levels=[0.0, 0.5], metrics=['SOME_METRIC'], out_csv='out.csv',
            node_out_csv='out_nodes.csv', timings_out_csv='out_timings.csv')
        self.assertEqual(result, (fake_df, fake_node_df, fake_timings_df))

    def test_defaults_match_the_full_real_discovery_catalog(self):
        fake = pd.DataFrame()
        with patch('lab.exp_disco_degrade.lab_run.run', return_value=(fake, fake, fake)) as mock_run:
            run_disco_degrade(['fake_log.xes'])

        _args, kwargs = mock_run.call_args
        self.assertEqual(kwargs['combos'], ALL_COMBOS)
        self.assertEqual(kwargs['degradations'], ALL_DEGRADATIONS)
        self.assertEqual(kwargs['levels'], ALL_LEVELS)
        self.assertEqual(kwargs['metrics'], ALL_METRICS)
        self.assertEqual(kwargs['out_csv'], 'var/lab/results/exp_disco_degrade.csv')


class DryRunTest(unittest.TestCase):
    """--dry-run prints the resolved Experiment and exits without calling
    run_disco_degrade at all - main()'s describe-then-return path must
    never fall through into computing anything."""

    def test_dry_run_with_run_name_prints_and_computes_nothing(self):
        with patch('sys.argv', ['exp_disco_degrade', '--run', 'smoke', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('lab.exp_disco_degrade.run_disco_degrade') as mock_run:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_run.assert_not_called()
        self.assertIn('Experiment: smoke', output)
        self.assertIn('payment_approval', output)
        self.assertIn('inductive_noise20', output)
        self.assertIn('cells:', output)

    def test_dry_run_with_ad_hoc_logs_prints_and_computes_nothing(self):
        with patch('sys.argv', ['exp_disco_degrade', 'fake_log.xes',
                                 '--combos', 'inductive_noise20', '--levels', '0.0', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('lab.exp_disco_degrade.run_disco_degrade') as mock_run:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_run.assert_not_called()
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
             patch('lab.exp_disco_degrade.run_disco_degrade') as mock_run:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_run.assert_not_called()
        self.assertIn('claims_known', output)
        self.assertIn('appeal_seq', output)
        self.assertIn('loop_block', output)

    def test_unknown_combo_name_errors_out(self):
        with patch('sys.argv', ['exp_disco_degrade', 'fake_log.xes',
                                 '--combos', 'not_a_real_combo', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'):
            with self.assertRaises(SystemExit):
                main()


class ExcludeMetricCliTest(unittest.TestCase):
    """--exclude-metric computes ALL_METRICS minus the given ids and
    passes that inclusive complement to run_disco_degrade's own
    `metrics=` - the actual exclusion behaviour (null-fallback matching,
    no-column-for-excluded-metric) is lab.run's, covered in
    tests/lab/test_run.py."""

    def test_excluded_metric_is_removed_from_what_gets_passed_through(self):
        with patch('sys.argv', ['exp_disco_degrade', '--run', 'smoke',
                                 '--exclude-metric', 'voidsat2', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('lab.exp_disco_degrade.run_disco_degrade') as mock_run:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        # --dry-run still prints and computes nothing, but confirms the
        # exclusion was parsed and would have applied.
        mock_run.assert_not_called()
        self.assertIn('excluding metrics', output)
        self.assertIn('voidsat2', output)

    def test_excluded_metric_reaches_run_disco_degrade_as_the_complement(self):
        fake = pd.DataFrame()
        with patch('sys.argv', ['exp_disco_degrade', '--run', 'smoke',
                                 '--exclude-metric', 'voidsat2']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('lab.exp_disco_degrade.run_disco_degrade',
                   return_value=(fake, fake, fake)) as mock_run:
            main()

        _args, kwargs = mock_run.call_args
        scored_ids = {m.id for m in kwargs['metrics']}
        self.assertNotIn('voidsat2', scored_ids)
        self.assertIn('skipprob', scored_ids)

    def test_unknown_exclude_metric_name_errors_out(self):
        with patch('sys.argv', ['exp_disco_degrade', '--run', 'smoke',
                                 '--exclude-metric', 'not_a_real_metric', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'):
            with self.assertRaises(SystemExit):
                main()


if __name__ == '__main__':
    unittest.main()
