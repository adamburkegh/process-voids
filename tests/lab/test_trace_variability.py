'''
Tests for lab.trace_variability.

The trace-degradation dose-response curves are useful to look at but
too flat/uneventful to be worth a figure in the paper - this reduces
each curve (one per metric x log x combo, restricted to
degradation_dim='trace') to a single number, its standard deviation
across levels, so Adam can quote a min/max across every curve instead
of showing all of them.
'''

import unittest

import pandas as pd

from lab.trace_variability import (
    std_dev_for, summarize, to_latex_table, to_markdown_table, variability_table)


def _row(combo, dim, level, **metrics):
    return {'combo': combo, 'degradation_dim': dim, 'degradation_level': level,
            'status': 'ok', **metrics}


def fake_df():
    return pd.DataFrame([
        _row('inductive_noise20', 'trace', 0.0, voidsalign3=0.1, voidsat2=0.2,
             voidmass_process_lower=0.05, voidmass_process_upper=0.15),
        _row('inductive_noise20', 'trace', 0.5, voidsalign3=0.1, voidsat2=0.2,
             voidmass_process_lower=0.05, voidmass_process_upper=0.15),
        _row('inductive_noise20', 'trace', 1.0, voidsalign3=1.0, voidsat2=1.0,
             voidmass_process_lower=1.0, voidmass_process_upper=1.0),
        # a different dimension present in the same file - must be ignored
        _row('inductive_noise20', 'activity_gradual', 0.5, voidsalign3=0.9, voidsat2=0.9,
             voidmass_process_lower=0.9, voidmass_process_upper=0.9),
    ])


def fake_df_missing_voidsalign3():
    return pd.DataFrame([
        _row('inductive_noise20', 'trace', 0.0, voidsat2=0.3,
             voidmass_process_lower=0.1, voidmass_process_upper=0.2),
        _row('inductive_noise20', 'trace', 0.5, voidsat2=0.6,
             voidmass_process_lower=0.25, voidmass_process_upper=0.35),
    ])


class StdDevForTest(unittest.TestCase):
    def test_zero_std_dev_when_flat(self):
        df = fake_df()
        # level 1.0 is excluded by _exclude_degenerate in load_logs, not
        # here - std_dev_for operates on whatever it's given, so pass an
        # already-filtered frame to isolate this from that behaviour.
        flat = df[df['degradation_level'] != 1.0]
        self.assertAlmostEqual(std_dev_for(flat, 'inductive_noise20', 'voidsalign3'), 0.0)

    def test_nonzero_std_dev_reflects_spread(self):
        df = fake_df()
        std = std_dev_for(df, 'inductive_noise20', 'voidsalign3')
        self.assertGreater(std, 0.0)

    def test_banded_metric_uses_midpoint(self):
        df = fake_df()
        flat = df[df['degradation_level'] != 1.0]
        self.assertAlmostEqual(std_dev_for(flat, 'inductive_noise20', 'voidmass_process'), 0.0)

    def test_missing_metric_column_returns_none(self):
        df = fake_df_missing_voidsalign3()
        self.assertIsNone(std_dev_for(df, 'inductive_noise20', 'voidsalign3'))

    def test_ignores_other_degradation_dimensions(self):
        '''The activity_gradual row in fake_df would blow the std dev up
        if it leaked in - confirms _series_for's own dim filter is doing
        its job here.'''
        df = fake_df()
        flat = df[df['degradation_level'] != 1.0]
        self.assertAlmostEqual(std_dev_for(flat, 'inductive_noise20', 'voidsalign3'), 0.0)


class VariabilityTableTest(unittest.TestCase):
    def test_shape_is_metric_combo_by_log(self):
        log_dfs = {'log_a': fake_df(), 'log_b': fake_df()}
        df = variability_table(log_dfs, combos=['inductive_noise20'], metric_ids=['voidsalign3'])
        self.assertEqual(list(df.columns), ['log_a', 'log_b'])
        self.assertEqual(list(df.index), [('Void by Skip Alignment', 'inductive_noise20')])

    def test_missing_metric_reads_as_nan_not_error(self):
        log_dfs = {'complete': fake_df(), 'partial': fake_df_missing_voidsalign3()}
        df = variability_table(log_dfs, combos=['inductive_noise20'], metric_ids=['voidsalign3'])
        self.assertTrue(pd.isna(df.loc[('Void by Skip Alignment', 'inductive_noise20'), 'partial']))


class SummarizeTest(unittest.TestCase):
    def test_reports_min_and_max_ignoring_nan(self):
        df = pd.DataFrame(
            {'log_a': [0.1, float('nan')], 'log_b': [0.4, 0.2]},
            index=pd.MultiIndex.from_tuples([('m1', 'c1'), ('m2', 'c1')], names=['metric', 'combo']),
        )
        min_val, min_where, max_val, max_where = summarize(df)
        self.assertAlmostEqual(min_val, 0.1)
        self.assertEqual(min_where, ('m1', 'c1', 'log_a'))
        self.assertAlmostEqual(max_val, 0.4)
        self.assertEqual(max_where, ('m1', 'c1', 'log_b'))


class TableFormattingTest(unittest.TestCase):
    def _sample_df(self):
        return pd.DataFrame(
            {'rtfm': [0.0, float('nan')], 'bpic2020_rfp': [0.00854, 0.0006]},
            index=pd.MultiIndex.from_tuples(
                [('Void by Aligned Durations', 'inductive_noise20'),
                 ('Void by Skip Alignment', 'inductive_noise20')],
                names=['metric', 'combo']),
        )

    def test_markdown_table_formats_values_and_missing(self):
        text = to_markdown_table(self._sample_df())
        self.assertIn('Void by Aligned Durations', text)
        self.assertIn('0.0085', text)
        self.assertIn('inductive_noise20', text)
        # a NaN cell must not render as the literal string 'nan'
        self.assertNotIn('nan', text)

    def test_latex_table_formats_values_and_missing(self):
        text = to_latex_table(self._sample_df())
        self.assertIn(r'\begin{tabular}', text)
        self.assertIn(r'\end{tabular}', text)
        self.assertIn('0.0085', text)
        self.assertNotIn('nan', text)


if __name__ == '__main__':
    unittest.main()
