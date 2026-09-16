'''
Tests for lab.plots.

The panels come from lab.metric_registry, and so do this file's fixture
columns: every fake result frame here is built from metric_panels()
rather than from hand-typed column names. That is the point - four times
a metric was renamed or retired without the plots following, and the
suite never noticed, because the fixtures had been hand-typed with the
same stale names as the code.

Smoke tests for plot_dose_response check that it runs and writes the
expected file, not the rendered pixels. _exclude_degenerate is tested
directly as the plain DataFrame filter it is.
'''

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from lab.metric_registry import Metric
from lab.plots import (
    _exclude_degenerate, average_over_nodes, metric_panels, plot_dose_response)
from lab.run import NODE_ROW_COLUMNS


def _panel_columns():
    return [col for _label, cols in metric_panels() for col in cols]


def _banded_panel():
    return next(cols for _label, cols in metric_panels() if len(cols) == 2)


def _values(plain, lower, upper):
    '''One value per panel column: `plain` for a single-column panel,
    (lower, upper) for a banded one.'''
    values = {}
    for _label, cols in metric_panels():
        if len(cols) == 1:
            values[cols[0]] = plain
        else:
            values[cols[0]], values[cols[1]] = lower, upper
    return values


def fake_disco_df():
    rows = []
    for level in (0.0, 0.5, 1.0):
        rows.append({
            'log': 'fake_log', 'combo': 'inductive', 'degradation_dim': 'activity',
            'degradation_level': level, 'status': 'ok',
            **_values(0.5, 0.4, 0.6),
        })
    return pd.DataFrame(rows)


def fake_claims_shaped_df():
    '''Same schema as fake_disco_df, but shaped the way the claims
    fixture's registration (lab.claims_fixture's CLAIMS_COMBOS/
    CLAIMS_DEGRADATIONS) produces it: a single fixed combo
    ('claims_known'), degradation_dim holding the ablation target name
    instead of a generic dimension.'''
    rows = []
    for target in ('assess', 'loop_block'):
        for level in (0.0, 0.5):
            rows.append({
                'log': 'claims', 'combo': 'claims_known', 'degradation_dim': target,
                'degradation_level': level, 'status': 'ok',
                **_values(0.5, 0.01 * level, 0.01 * level),
            })
    rows.append({'log': 'claims', 'combo': 'claims_known', 'degradation_dim': 'assess',
                 'degradation_level': 0.9, 'status': 'error: boom'})
    return pd.DataFrame(rows)


def fake_node_df():
    '''The per-node CSV's shape: one row per (log, combo, degradation_dim,
    degradation_level, node_id), no status column. Every panel column of
    an Activity row reads `value`; a Tau row reads 1.0, far outside the
    Activity rows' range, so an average it contaminated would show it.'''
    def row(node_id, node_type, value, level=0.0):
        return {
            'log': 'fake_log', 'combo': 'inductive', 'degradation_dim': 'activity',
            'degradation_level': level, 'node_id': node_id, 'node_type': node_type,
            **_values(value, value, value),
        }
    return pd.DataFrame([
        row('1', 'Activity', 0.0),
        row('2', 'Activity', 0.5),
        row('3', 'Tau', 1.0),
        # a second cell (level=1.0) that must be dropped entirely
        row('1', 'Activity', 1.0, level=1.0),
    ])


def _metric(id, **overrides):
    fields = dict(id=id, description='d', source='s', scripts=('exp_disco_degrade',),
                  scale='coverage')
    fields.update(overrides)
    return Metric(**fields)


class MetricPanelsQueryTest(unittest.TestCase):
    '''
    The rule metric_panels applies, against a synthetic registry so the
    tests name no real metric - which is what lets a real rename or
    retirement leave them untouched.
    '''

    def _panels(self, *metrics):
        with patch('lab.metric_registry.METRICS', {m.id: m for m in metrics}):
            return metric_panels()

    def test_a_live_scaled_plotted_metric_emitted_by_the_runner_is_a_panel(self):
        self.assertEqual(self._panels(_metric('m')), [('m', ('m',))])

    def test_a_retired_metric_is_not_a_panel(self):
        self.assertEqual(self._panels(_metric('m', status='retired')), [])

    def test_an_unscaled_metric_is_not_a_panel(self):
        self.assertEqual(self._panels(_metric('m', scale=None)), [])

    def test_an_unplotted_metric_is_not_a_panel(self):
        self.assertEqual(self._panels(_metric('m', plotted=False)), [])

    def test_a_metric_the_runner_does_not_emit_is_not_a_panel(self):
        self.assertEqual(self._panels(_metric('m', scripts=('exp_voidmass',))), [])

    def test_lower_and_upper_fold_into_one_banded_panel_lower_first(self):
        panels = self._panels(_metric('m_upper'), _metric('m_lower'))
        self.assertEqual(panels, [('m', ('m_lower', 'm_upper'))])

    def test_panels_follow_registry_order(self):
        panels = self._panels(_metric('b'), _metric('a'))
        self.assertEqual([label for label, _cols in panels], ['b', 'a'])


class RenameAndRetireTest(unittest.TestCase):
    '''
    The acceptance test for the whole change: renaming or retiring a live
    metric in the registry needs no edit to lab.plots or to these tests.
    Before, each of those was a KeyError at plot time.
    '''

    def test_a_renamed_metric_is_plotted_under_its_new_name(self):
        with patch('lab.metric_registry.METRICS', {'renamed': _metric('renamed')}):
            self.assertEqual(metric_panels(), [('renamed', ('renamed',))])
            written = plot_dose_response(fake_disco_df(), out_dir=tempfile.mkdtemp())
        self.assertTrue(Path(written[0]).exists())

    def test_a_retired_metric_leaves_the_plots_and_they_still_draw(self):
        registry = {'kept': _metric('kept'), 'gone': _metric('gone', status='retired')}
        with patch('lab.metric_registry.METRICS', registry):
            self.assertEqual(metric_panels(), [('kept', ('kept',))])
            # a result frame without the retired column at all
            written = plot_dose_response(fake_disco_df(), out_dir=tempfile.mkdtemp())
        self.assertTrue(Path(written[0]).exists())


class RealRegistryPanelsTest(unittest.TestCase):
    '''The real registry's panels, checked by shape rather than by name.'''

    def test_every_panel_column_is_one_the_runner_writes(self):
        '''The KeyError this change exists to prevent: a panel reading a
        column that no current result CSV carries.'''
        missing = [col for col in _panel_columns() if col not in NODE_ROW_COLUMNS]
        self.assertEqual(missing, [])

    def test_every_panel_is_a_one_or_two_column_tuple(self):
        for label, cols in metric_panels():
            with self.subTest(label=label):
                self.assertIn(len(cols), (1, 2))

    def test_no_column_is_plotted_twice(self):
        columns = _panel_columns()
        self.assertEqual(len(columns), len(set(columns)))


class AverageOverNodesTest(unittest.TestCase):
    def _cell(self):
        averaged = average_over_nodes(fake_node_df())
        return averaged[(averaged['degradation_dim'] == 'activity')
                        & (averaged['degradation_level'] == 0.0)].iloc[0]

    def test_excludes_tau_nodes_from_the_average(self):
        '''Mean of the two Activity rows (0.0, 0.5), not pulled toward the
        Tau row's 1.0 - for every panel column.'''
        row = self._cell()
        for col in _panel_columns():
            with self.subTest(column=col):
                self.assertAlmostEqual(row[col], 0.25)

    def test_averages_both_bound_columns_for_a_banded_metric(self):
        row = self._cell()
        lower, upper = _banded_panel()
        self.assertAlmostEqual(row[lower], 0.25)
        self.assertAlmostEqual(row[upper], 0.25)

    def test_adds_an_ok_status_column(self):
        averaged = average_over_nodes(fake_node_df())
        self.assertTrue((averaged['status'] == 'ok').all())

    def test_drops_degradation_level_one_before_averaging(self):
        averaged = average_over_nodes(fake_node_df())
        self.assertNotIn(1.0, set(averaged['degradation_level']))

    def test_result_feeds_directly_into_plot_dose_response(self):
        out_dir = tempfile.mkdtemp()
        written = plot_dose_response(average_over_nodes(fake_node_df()), out_dir=out_dir)
        self.assertEqual(len(written), 1)
        self.assertTrue(Path(written[0]).exists())


class ExcludeDegenerateTest(unittest.TestCase):
    def test_level_one_is_excluded(self):
        filtered = _exclude_degenerate(fake_disco_df())
        self.assertNotIn(1.0, set(filtered['degradation_level']))
        self.assertEqual(set(filtered['degradation_level']), {0.0, 0.5})

    def test_error_status_is_excluded(self):
        df = fake_disco_df()
        df.loc[df['degradation_level'] == 0.5, 'status'] = 'error: boom'
        filtered = _exclude_degenerate(df)
        self.assertEqual(set(filtered['degradation_level']), {0.0})


class PlotDoseResponseTest(unittest.TestCase):
    def test_writes_one_png_per_log_and_dim_with_default_line_by(self):
        out_dir = tempfile.mkdtemp()
        written = plot_dose_response(fake_disco_df(), out_dir=out_dir)
        self.assertEqual(len(written), 1)  # one (log, dim) pair
        self.assertTrue(Path(written[0]).exists())
        self.assertTrue(str(written[0]).endswith('fake_log_activity.png'))

    def test_line_by_degradation_dim_facets_by_combo_instead(self):
        # the claims-shaped case: one combo (claims_known), multiple
        # degradation_dim values (ablation targets) - line_by should
        # put those on the SAME plot as separate lines, not split them
        # into separate figures the way line_by='combo' (default) would.
        out_dir = tempfile.mkdtemp()
        written = plot_dose_response(fake_claims_shaped_df(), out_dir=out_dir,
                                      line_by='degradation_dim')
        self.assertEqual(len(written), 1)  # one (log, combo) pair
        self.assertTrue(Path(written[0]).exists())
        self.assertTrue(str(written[0]).endswith('claims_claims_known.png'))

    def test_line_by_rejects_an_unknown_value(self):
        with self.assertRaises(ValueError):
            plot_dose_response(fake_disco_df(), line_by='nonsense')

    def test_error_rows_are_excluded_without_crashing(self):
        out_dir = tempfile.mkdtemp()
        written = plot_dose_response(fake_claims_shaped_df(), out_dir=out_dir,
                                      line_by='degradation_dim')
        self.assertTrue(Path(written[0]).exists())

    def test_excluded_level_one_row_can_hold_unplottable_data(self):
        # a non-numeric value in the (excluded) level=1.0 row would break
        # matplotlib's ax.plot if it ever reached it - if this test fails
        # (an exception, not just a wrong plot), _exclude_degenerate has
        # stopped being applied inside plot_dose_response
        df = fake_disco_df()
        lower, upper = _banded_panel()
        df[lower] = df[lower].astype(object)
        df[upper] = df[upper].astype(object)
        df.loc[df['degradation_level'] == 1.0, [lower, upper]] = 'not a number'
        out_dir = tempfile.mkdtemp()
        written = plot_dose_response(df, out_dir=out_dir)
        self.assertTrue(Path(written[0]).exists())

    def test_one_axis_per_panel_labelled_with_its_registry_id(self):
        with patch('matplotlib.axes.Axes.set_ylabel') as mock_ylabel:
            plot_dose_response(fake_disco_df(), out_dir=tempfile.mkdtemp())
        self.assertEqual([call.args[0] for call in mock_ylabel.call_args_list],
                         [label for label, _cols in metric_panels()])

    def test_title_suffix_is_appended_to_the_figure_title(self):
        out_dir = tempfile.mkdtemp()
        with patch('matplotlib.figure.Figure.suptitle') as mock_suptitle:
            plot_dose_response(fake_disco_df(), out_dir=out_dir, title_suffix=' (per-node average)')
        mock_suptitle.assert_called_once_with('fake_log - activity (per-node average)')

    def test_title_suffix_defaults_to_empty(self):
        out_dir = tempfile.mkdtemp()
        with patch('matplotlib.figure.Figure.suptitle') as mock_suptitle:
            plot_dose_response(fake_disco_df(), out_dir=out_dir)
        mock_suptitle.assert_called_once_with('fake_log - activity')

    def test_a_bound_pair_with_no_divergence_still_plots(self):
        # the common case (no timed-out variants, lower == upper
        # everywhere) - fake_claims_shaped_df gives every banded panel
        # equal bounds, and the banded path must not require real
        # divergence to work.
        out_dir = tempfile.mkdtemp()
        written = plot_dose_response(fake_claims_shaped_df(), out_dir=out_dir,
                                      line_by='degradation_dim')
        self.assertTrue(Path(written[0]).exists())


if __name__ == '__main__':
    unittest.main()
