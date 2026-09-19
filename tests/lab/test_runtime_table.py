'''
Tests for lab.runtime_table.

A metric's own row in a *_timings.csv is usually near-zero - the real
cost sits in shared stages (classical, dv, ...) its ProcessMetric.needs
declares. total_seconds sums a metric's own row(s) plus its declared
stage rows, so the reported number reads as "what this metric alone
would have cost", not something additive across metrics without
double-counting a shared stage.
'''

import unittest

import pandas as pd

from lab.runtime_table import (
    format_duration, load_timings, runtime_table, to_latex_table, to_markdown_table,
    total_seconds)


def _timings_df(rows):
    return pd.DataFrame([
        {'combo': combo, 'metric_or_stage': stage, 'seconds': seconds}
        for combo, stage, seconds in rows
    ])


class TotalSecondsTest(unittest.TestCase):
    def test_sums_own_row_and_declared_stages(self):
        df = _timings_df([
            ('inductive_noise20', 'dv', 10.0),
            ('inductive_noise20', 'executions_cache', 0.5),
            ('inductive_noise20', 'voidsalign3', 0.01),
            ('inductive_noise20', 'classical', 800.0),  # unrelated stage, must not be counted
        ])
        self.assertAlmostEqual(total_seconds(df, 'inductive_noise20', 'voidsalign3'), 10.51)

    def test_falls_back_to_voidsalign2_when_voidsalign3_absent(self):
        '''A pre-refresh timings CSV has voidsalign2's own row, not
        voidsalign3's - same version gap as lab.cross_log_plots.'''
        df = _timings_df([
            ('inductive_noise20', 'dv', 5.0),
            ('inductive_noise20', 'executions_cache', 0.2),
            ('inductive_noise20', 'voidsalign2', 0.02),
        ])
        self.assertAlmostEqual(total_seconds(df, 'inductive_noise20', 'voidsalign3'), 5.22)

    def test_counts_voidmass_process2_own_rows(self):
        df = _timings_df([
            ('inductive_noise20', 'classical', 100.0),
            ('inductive_noise20', 'voidmass_process2_lower', 0.25),
            ('inductive_noise20', 'voidmass_process2_upper', 0.5),
        ])
        self.assertAlmostEqual(total_seconds(df, 'inductive_noise20', 'voidmass_process'), 100.75)

    def test_different_combo_is_isolated(self):
        df = _timings_df([
            ('inductive_noise20', 'classical', 100.0),
            ('inductive_noise20', 'voidmass_process_lower', 0.0),
            ('toothpaste_noise10', 'classical', 9000.0),
            ('toothpaste_noise10', 'voidmass_process_upper', 0.0),
        ])
        self.assertAlmostEqual(total_seconds(df, 'inductive_noise20', 'voidmass_process'), 100.0)
        self.assertAlmostEqual(total_seconds(df, 'toothpaste_noise10', 'voidmass_process'), 9000.0)


class FormatDurationTest(unittest.TestCase):
    def test_seconds_only(self):
        self.assertEqual(format_duration(45), '45s')

    def test_minutes_and_seconds(self):
        self.assertEqual(format_duration(125), '2m 5s')

    def test_hours_and_minutes(self):
        self.assertEqual(format_duration(3725), '1h 2m')

    def test_zero(self):
        self.assertEqual(format_duration(0), '0s')


class RuntimeTableTest(unittest.TestCase):
    def test_shape_is_metric_combo_by_log(self):
        log_timings = {
            'rtfm': _timings_df([
                ('inductive_noise20', 'dv', 10.0),
                ('inductive_noise20', 'executions_cache', 1.0),
                ('inductive_noise20', 'voidsalign3', 0.0),
            ]),
            'bpic2020_rfp': _timings_df([
                ('inductive_noise20', 'dv', 20.0),
                ('inductive_noise20', 'executions_cache', 2.0),
                ('inductive_noise20', 'voidsalign3', 0.0),
            ]),
        }
        df = runtime_table(log_timings, combos=['inductive_noise20'], metric_ids=['voidsalign3'])
        self.assertEqual(list(df.columns), ['rtfm', 'bpic2020_rfp'])
        self.assertEqual(list(df.index), [('Void by Skip Alignment', 'inductive_noise20')])
        self.assertAlmostEqual(df.loc[('Void by Skip Alignment', 'inductive_noise20'), 'rtfm'], 11.0)
        self.assertAlmostEqual(
            df.loc[('Void by Skip Alignment', 'inductive_noise20'), 'bpic2020_rfp'], 22.0)


class TableFormattingTest(unittest.TestCase):
    def _sample_df(self):
        return pd.DataFrame(
            {'rtfm': [125.0], 'bpic2020_rfp': [3725.0]},
            index=pd.MultiIndex.from_tuples([('Void by Skip Alignment', 'inductive_noise20')],
                                             names=['metric', 'combo']),
        )

    def test_markdown_table_has_formatted_durations(self):
        text = to_markdown_table(self._sample_df())
        self.assertIn('Void by Skip Alignment', text)
        self.assertIn('inductive_noise20', text)
        self.assertIn('2m 5s', text)
        self.assertIn('1h 2m', text)

    def test_latex_table_has_formatted_durations(self):
        text = to_latex_table(self._sample_df())
        self.assertIn(r'\begin{tabular}', text)
        self.assertIn(r'\end{tabular}', text)
        self.assertIn('Void by Skip Alignment', text)
        self.assertIn('2m 5s', text)


class LoadTimingsTest(unittest.TestCase):
    def test_reads_one_df_per_log(self):
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'a_timings.csv'
            _timings_df([('inductive_noise20', 'dv', 1.0)]).to_csv(path, index=False)
            loaded = load_timings({'rtfm': str(path)})
            self.assertEqual(list(loaded['rtfm']['metric_or_stage']), ['dv'])


if __name__ == '__main__':
    unittest.main()
