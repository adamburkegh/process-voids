'''
Tests for lab.cross_log_plots.

lab.plots compares combos or dimensions within ONE log's own result
CSV. This module compares the SAME combo/dimension across DIFFERENT
logs' result CSVs - each log may live in a differently-versioned CSV
(eg one log's sweep predates a metric another log's includes), so the
central behaviour under test is: a log missing a metric's column is
skipped for that line, not an error - same convention as
lab.check_run's missing-column handling.

Smoke tests for plot_cross_log check that it runs and writes the
expected files, not the rendered pixels - same approach as
tests.lab.test_plots.
'''

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
import pandas as pd

from lab.cross_log_plots import METRIC_SPECS, _series_for, load_logs, plot_cross_log


def _row(combo, dim, level, **metrics):
    return {'combo': combo, 'degradation_dim': dim, 'degradation_level': level,
            'status': 'ok', **metrics}


def fake_df_with_all_metrics():
    return pd.DataFrame([
        _row('inductive_noise20', 'trace', 0.0, voidsalign3=0.1, voidsat2=0.2,
             voidmass_process_lower=0.05, voidmass_process_upper=0.15),
        _row('inductive_noise20', 'trace', 0.5, voidsalign3=0.4, voidsat2=0.5,
             voidmass_process_lower=0.2, voidmass_process_upper=0.3),
        _row('inductive_noise20', 'trace', 1.0, voidsalign3=1.0, voidsat2=1.0,
             voidmass_process_lower=1.0, voidmass_process_upper=1.0),
    ])


def fake_df_missing_voidsalign3():
    '''A log whose sweep predates voidsalign3 - has voidsat2/voidmass_process
    but no voidsalign3 column at all, the shape a pre-refresh result CSV
    actually has.'''
    return pd.DataFrame([
        _row('inductive_noise20', 'trace', 0.0, voidsat2=0.3,
             voidmass_process_lower=0.1, voidmass_process_upper=0.2),
        _row('inductive_noise20', 'trace', 0.5, voidsat2=0.6,
             voidmass_process_lower=0.25, voidmass_process_upper=0.35),
    ])


class SeriesForTest(unittest.TestCase):
    def test_missing_column_returns_none(self):
        df = fake_df_missing_voidsalign3()
        self.assertIsNone(_series_for(df, 'inductive_noise20', 'trace', ('voidsalign3',)))

    def test_present_column_returns_sorted_rows(self):
        df = fake_df_with_all_metrics()
        series = _series_for(df, 'inductive_noise20', 'trace', ('voidsat2',))
        self.assertEqual(list(series['degradation_level']), [0.0, 0.5, 1.0])

    def test_unmatched_combo_or_dim_returns_none(self):
        df = fake_df_with_all_metrics()
        self.assertIsNone(_series_for(df, 'toothpaste_noise10', 'trace', ('voidsat2',)))
        self.assertIsNone(_series_for(df, 'inductive_noise20', 'activity_frequency_gradual',
                                       ('voidsat2',)))

    def test_banded_metric_needs_both_columns_present(self):
        df = fake_df_missing_voidsalign3()
        # voidmass_process_lower/upper are both present here
        series = _series_for(df, 'inductive_noise20', 'trace',
                              ('voidmass_process_lower', 'voidmass_process_upper'))
        self.assertIsNotNone(series)


class LoadLogsTest(unittest.TestCase):
    def test_drops_degenerate_rows(self):
        path = Path(tempfile.mkdtemp()) / 'result.csv'
        df = pd.DataFrame([
            _row('inductive_noise20', 'trace', 0.0, voidsat2=0.1),
            _row('inductive_noise20', 'trace', 1.0, voidsat2=1.0),
            {**_row('inductive_noise20', 'trace', 0.5, voidsat2=0.5), 'status': 'error: boom'},
        ])
        df.to_csv(path, index=False)
        loaded = load_logs({'fake_log': str(path)})
        self.assertEqual(list(loaded['fake_log']['degradation_level']), [0.0])


class PlotCrossLogTest(unittest.TestCase):
    def test_writes_one_file_per_combo_dim_pair(self):
        log_dfs = {'log_a': fake_df_with_all_metrics(), 'log_b': fake_df_with_all_metrics()}
        with tempfile.TemporaryDirectory() as out_dir:
            written = plot_cross_log(log_dfs, combos=['inductive_noise20'],
                                      degradation_dims=['trace'], out_dir=out_dir)
            self.assertEqual(len(written), 1)
            self.assertTrue(written[0].exists())
            self.assertEqual(written[0].name, 'inductive_noise20_trace.png')

    def test_multiple_combos_and_dims_multiply_out(self):
        log_dfs = {'log_a': fake_df_with_all_metrics()}
        with tempfile.TemporaryDirectory() as out_dir:
            written = plot_cross_log(log_dfs, combos=['inductive_noise20', 'toothpaste_noise10'],
                                      degradation_dims=['trace', 'activity_frequency_gradual'],
                                      out_dir=out_dir)
            self.assertEqual(len(written), 4)

    def test_log_missing_a_metric_column_does_not_raise(self):
        log_dfs = {'complete_log': fake_df_with_all_metrics(),
                   'partial_log': fake_df_missing_voidsalign3()}
        with tempfile.TemporaryDirectory() as out_dir:
            written = plot_cross_log(log_dfs, combos=['inductive_noise20'],
                                      degradation_dims=['trace'],
                                      metric_ids=['voidsalign3', 'voidsat2'], out_dir=out_dir)
            self.assertEqual(len(written), 1)
            self.assertTrue(written[0].exists())

    def test_ylim_applied_to_every_panel(self):
        '''fake_df_with_all_metrics' values (0.05-1.0) autoscale to a range
        distinguishable from an explicit (0.0, 1.0), so capturing the real
        axes and reading back get_ylim() is a genuine check, not a tautology.'''
        import matplotlib.pyplot as plt
        real_subplots = plt.subplots
        created_axes = []

        def spy_subplots(*args, **kwargs):
            fig, axes = real_subplots(*args, **kwargs)
            created_axes.extend(axes if hasattr(axes, '__iter__') else [axes])
            return fig, axes

        log_dfs = {'log_a': fake_df_with_all_metrics()}
        with patch('lab.cross_log_plots.plt.subplots', side_effect=spy_subplots):
            with tempfile.TemporaryDirectory() as out_dir:
                plot_cross_log(log_dfs, combos=['inductive_noise20'], degradation_dims=['trace'],
                                out_dir=out_dir, ylim=(0.0, 1.0))
        self.assertEqual(len(created_axes), len(METRIC_SPECS))
        for ax in created_axes:
            self.assertEqual(ax.get_ylim(), (0.0, 1.0))

    def test_panel_by_log_writes_by_log_suffixed_file(self):
        log_dfs = {'log_a': fake_df_with_all_metrics(), 'log_b': fake_df_with_all_metrics()}
        with tempfile.TemporaryDirectory() as out_dir:
            written = plot_cross_log(log_dfs, combos=['inductive_noise20'],
                                      degradation_dims=['trace'], out_dir=out_dir, panel_by='log')
            self.assertEqual(len(written), 1)
            self.assertEqual(written[0].name, 'inductive_noise20_trace_by_log.png')

    def test_panel_by_log_metric_ignores_combo_as_figure_axis(self):
        '''panel_by='log_metric' makes combo a LINE, not a figure axis - two
        combos and two dims should still write one file per dim, not per
        (combo, dim), since every combo's line lands in the same panel.'''
        log_dfs = {'log_a': fake_df_with_all_metrics(), 'log_b': fake_df_with_all_metrics()}
        with tempfile.TemporaryDirectory() as out_dir:
            written = plot_cross_log(log_dfs, combos=['inductive_noise20', 'toothpaste_noise10'],
                                      degradation_dims=['trace', 'activity_frequency_gradual'],
                                      out_dir=out_dir, panel_by='log_metric')
            self.assertEqual(len(written), 2)
            names = sorted(p.name for p in written)
            self.assertEqual(names, ['activity_frequency_gradual_by_log_metric.png',
                                      'trace_by_log_metric.png'])

    def test_panel_by_log_metric_has_one_panel_per_log_metric_pair(self):
        log_dfs = {'log_a': fake_df_with_all_metrics(), 'log_b': fake_df_with_all_metrics()}
        with patch('lab.cross_log_plots.plt.subplots', wraps=plt.subplots) as spy:
            with tempfile.TemporaryDirectory() as out_dir:
                plot_cross_log(log_dfs, combos=['inductive_noise20'], degradation_dims=['trace'],
                                out_dir=out_dir, panel_by='log_metric',
                                metric_ids=['voidsalign3', 'voidsat2'])
        # 2 logs x 2 metrics = 4 panels
        spy.assert_called_once_with(1, 4, figsize=(20, 4))

    def test_invalid_panel_by_raises(self):
        log_dfs = {'log_a': fake_df_with_all_metrics()}
        with tempfile.TemporaryDirectory() as out_dir:
            with self.assertRaises(ValueError):
                plot_cross_log(log_dfs, combos=['inductive_noise20'], degradation_dims=['trace'],
                                out_dir=out_dir, panel_by='nonsense')

    def test_metric_specs_cover_default_ids(self):
        for metric_id in ('voidsalign3', 'voidsat2', 'voidmass_process'):
            self.assertIn(metric_id, METRIC_SPECS)


if __name__ == '__main__':
    unittest.main()
