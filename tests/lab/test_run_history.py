"""
The run-history CSV: one row per run, genuinely appended.

The point of the file is that a result CSV stays readable a year later
without timestamp-matching a text log, so the tests here are about
durability - two runs leave two rows, an unset seed is still recorded,
and a failure to write history never costs a completed run its results.
"""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from lab.run_history import RUN_HISTORY_COLUMNS, append_run_history, run_history_row


class AppendSemanticsTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / 'history' / 'run_history.csv'

    def _row(self, **overrides):
        row = {k: '' for k in RUN_HISTORY_COLUMNS}
        row.update(overrides)
        return row

    def test_a_first_write_creates_the_file_with_the_declared_columns(self):
        append_run_history(self.path, self._row(out_csv='a.csv'))
        df = pd.read_csv(self.path)
        self.assertEqual(list(df.columns), list(RUN_HISTORY_COLUMNS))

    def test_two_runs_leave_two_rows(self):
        append_run_history(self.path, self._row(out_csv='a.csv'))
        append_run_history(self.path, self._row(out_csv='b.csv'))
        df = pd.read_csv(self.path)
        self.assertEqual(len(df), 2)
        self.assertEqual(list(df['out_csv']), ['a.csv', 'b.csv'])

    def test_a_rerun_to_the_same_out_csv_appends_rather_than_replacing(self):
        """The upsert helper the other scripts use would purge the
        earlier row here. History is the one file that must not lose
        what it already holds."""
        append_run_history(self.path, self._row(out_csv='same.csv', run_timestamp='t1'))
        append_run_history(self.path, self._row(out_csv='same.csv', run_timestamp='t2'))
        df = pd.read_csv(self.path)
        self.assertEqual(len(df), 2)
        self.assertEqual(list(df['run_timestamp']), ['t1', 't2'])

    def test_column_order_is_stable_across_appends(self):
        append_run_history(self.path, self._row(out_csv='a.csv'))
        append_run_history(self.path, self._row(out_csv='b.csv'))
        self.assertEqual(list(pd.read_csv(self.path).columns), list(RUN_HISTORY_COLUMNS))


class RowContentTest(unittest.TestCase):

    def _row(self, **kwargs):
        kwargs.setdefault('out_csv', 'var/lab/results/out.csv')
        kwargs.setdefault('run_name', 'smoke')
        kwargs.setdefault('log_paths', ['logs/rtfm.xes.gz'])
        kwargs.setdefault('combos', ['inductive_noise20'])
        kwargs.setdefault('degradations', ['activity_gradual', 'trace'])
        kwargs.setdefault('levels', [0.0, 0.5])
        kwargs.setdefault('metric_ids', ['skipprob'])
        return run_history_row(**kwargs)

    def test_row_has_exactly_the_declared_columns(self):
        self.assertEqual(set(self._row()), set(RUN_HISTORY_COLUMNS))

    def test_config_lists_are_joined_readably(self):
        row = self._row()
        self.assertEqual(row['combos'], 'inductive_noise20')
        self.assertEqual(row['degradations'], 'activity_gradual;trace')
        self.assertEqual(row['levels'], '0.0;0.5')

    def test_logs_are_recorded_by_the_same_name_the_result_rows_use(self):
        """Joinability is the requirement, not prettiness: lab.run names
        a log Path(log_path).stem in every result row, so history has to
        spell it identically - including for a .xes.gz path, where that
        stem keeps a trailing '.xes'."""
        for log_path in ('data/payment_approval.xes',
                         'data/logs/rtfm_fine_appeal.xes.gz'):
            with self.subTest(log_path=log_path):
                self.assertEqual(self._row(log_paths=[log_path])['logs'],
                                 Path(log_path).stem)

    def test_out_csv_is_recorded_with_forward_slashes(self):
        """Recorded as a portable path rather than this machine's
        separator, so a history file stays legible - and joinable
        against a path typed by hand - away from Windows."""
        row = self._row(out_csv=Path('var') / 'lab' / 'results' / 'out.csv')
        self.assertEqual(row['out_csv'], 'var/lab/results/out.csv')

    def test_excluded_metrics_are_the_complement_of_what_was_scored(self):
        row = self._row(metric_ids=['skipprob'])
        excluded = row['excluded_metrics'].split(';')
        self.assertIn('voidsat2', excluded)
        self.assertNotIn('skipprob', excluded)

    def test_seed_records_what_the_interpreter_actually_saw(self):
        with patch.dict('os.environ', {'PYTHONHASHSEED': '0'}):
            self.assertEqual(self._row()['pythonhashseed'], '0')

    def test_an_unset_seed_is_recorded_as_unset_not_left_blank(self):
        """Blank reads back from CSV as NaN, indistinguishable from a
        column that was never written - the whole question this file
        exists to answer."""
        import os
        environ = {k: v for k, v in os.environ.items() if k != 'PYTHONHASHSEED'}
        with patch.dict('os.environ', environ, clear=True):
            self.assertEqual(self._row()['pythonhashseed'], 'unset')

    def test_versions_carry_git_state_not_just_the_release_number(self):
        with patch('lab.run_history._version_line', return_value='0.5.0 (git abc1234 dirty)'), \
             patch('lab.run_history.dependency_version_line', return_value='0.3.0 (git unknown)'):
            row = self._row()
        self.assertEqual(row['process_voids_version'], '0.5.0 (git abc1234 dirty)')
        self.assertEqual(row['skipalignments_version'], '0.3.0 (git unknown)')


class WriteFailureTest(unittest.TestCase):

    def test_a_failed_history_write_warns_and_does_not_raise(self):
        """A completed sweep is expensive; losing its results because
        the bookkeeping file could not be written would be the worse
        failure by far."""
        with patch('lab.run_history.pd.DataFrame.to_csv', side_effect=OSError('read-only')), \
             self.assertLogs('lab.run_history', level='WARNING') as logs:
            append_run_history(Path(tempfile.mkdtemp()) / 'h.csv',
                               {k: '' for k in RUN_HISTORY_COLUMNS})
        self.assertTrue(any('run history' in m.lower() for m in logs.output))


if __name__ == '__main__':
    unittest.main()
