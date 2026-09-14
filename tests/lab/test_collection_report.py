'''
lab.collection_report answers "has this metric ever been collected for
this log", not "what's its latest value" - a real result CSV's 'log'
column (not the filename) says which log a row belongs to, and a
metric counts as collected if some row has a real (non-null) value for
it, anywhere across every root-level result CSV ever written.
'''

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from lab.collection_report import collected_metrics, to_markdown_table


FAKE_METRIC_IDS = ['weight_coverage', 'voidsat']


def _write_csv(dir_path, name, rows):
    pd.DataFrame(rows).to_csv(Path(dir_path) / name, index=False)


class CollectedMetricsTest(unittest.TestCase):
    def setUp(self):
        self.results_dir = tempfile.mkdtemp()
        # patch at the point of use - collected_metrics reads the ids
        # off lab.collection_report's own imported name, not the
        # registry module directly
        patcher = patch('lab.collection_report.METRIC_IDS', FAKE_METRIC_IDS)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_marks_a_metric_done_when_a_csv_has_a_real_value(self):
        _write_csv(self.results_dir, 'payment_approval_20260101-000000.csv', [
            {'log': 'payment_approval', 'weight_coverage': 0.8, 'voidsat': None},
        ])
        df = collected_metrics(self.results_dir)
        self.assertTrue(df.loc['weight_coverage', 'payment_approval'])
        self.assertFalse(df.loc['voidsat', 'payment_approval'])

    def test_a_column_thats_present_but_all_null_is_not_done(self):
        _write_csv(self.results_dir, 'rtfm_20260101-000000.csv', [
            {'log': 'rtfm', 'weight_coverage': None},
            {'log': 'rtfm', 'weight_coverage': None},
        ])
        df = collected_metrics(self.results_dir)
        self.assertFalse(df.loc['weight_coverage', 'rtfm'])

    def test_nodes_and_timings_companions_are_not_scanned(self):
        # a genuine value only in the *_nodes.csv companion must not
        # count - the root-level CSV is the collection record
        _write_csv(self.results_dir, 'rtfm_20260101-000000_nodes.csv', [
            {'log': 'rtfm', 'voidsat': 0.5},
        ])
        _write_csv(self.results_dir, 'rtfm_20260101-000000_timings.csv', [
            {'log': 'rtfm', 'voidsat': 0.5},
        ])
        df = collected_metrics(self.results_dir)
        self.assertNotIn('rtfm', df.columns)

    def test_a_csv_with_no_log_column_is_skipped_without_error(self):
        _write_csv(self.results_dir, 'mass_term_probe.csv', [
            {'combo': 'inductive_noise20', 'rep': 0, 'seconds': 1.2},
        ])
        df = collected_metrics(self.results_dir)  # must not raise
        self.assertEqual(df.shape[1], 0)

    def test_two_csvs_for_the_same_log_are_combined(self):
        _write_csv(self.results_dir, 'rtfm_20260101-000000.csv', [
            {'log': 'rtfm', 'weight_coverage': 0.9, 'voidsat': None},
        ])
        _write_csv(self.results_dir, 'rtfm_20260102-000000.csv', [
            {'log': 'rtfm', 'weight_coverage': None, 'voidsat': 0.3},
        ])
        df = collected_metrics(self.results_dir)
        self.assertTrue(df.loc['weight_coverage', 'rtfm'])
        self.assertTrue(df.loc['voidsat', 'rtfm'])

    def test_every_registry_metric_id_gets_a_row_even_if_never_seen(self):
        _write_csv(self.results_dir, 'rtfm_20260101-000000.csv', [
            {'log': 'rtfm', 'weight_coverage': 0.9},
        ])
        df = collected_metrics(self.results_dir)
        self.assertEqual(list(df.index), FAKE_METRIC_IDS)
        self.assertFalse(df.loc['voidsat', 'rtfm'])


class ToMarkdownTableTest(unittest.TestCase):
    def test_renders_done_and_blank_cells(self):
        df = pd.DataFrame(
            {'payment_approval': [True, False], 'rtfm': [False, True]},
            index=['weight_coverage', 'voidsat'],
        )
        table = to_markdown_table(df)
        lines = table.splitlines()
        self.assertEqual(lines[0], '| metric | payment_approval | rtfm |')
        self.assertEqual(lines[1], '| --- | --- | --- |')
        self.assertEqual(lines[2], '| weight_coverage | Done |  |')
        self.assertEqual(lines[3], '| voidsat |  | Done |')


if __name__ == '__main__':
    unittest.main()
