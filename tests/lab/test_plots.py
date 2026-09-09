'''
Smoke tests for lab.plots - that plot_dose_response runs and writes the
expected file, for both line_by modes; not pixel-checking the rendered
plot. _exclude_degenerate is tested directly (a plain DataFrame filter,
no matplotlib involved) rather than trying to inspect what actually got
plotted.
'''

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from lab.plots import plot_dose_response, average_over_nodes, _exclude_degenerate


def fake_disco_df():
    rows = []
    for level in (0.0, 0.5, 1.0):
        rows.append({
            'log': 'fake_log', 'combo': 'inductive', 'degradation_dim': 'activity',
            'degradation_level': level, 'status': 'ok',
            'weight_coverage': 0.8, 'skipprob': 0.2, 'salign_coverage': 0.98,
            'alignment_coverage_pn': 0.95, 'voidmass_subprocess': 0.1,
            'voidmass_process': 0.05,
        })
    return pd.DataFrame(rows)


def fake_claims_shaped_df():
    '''Same schema as fake_disco_df, but shaped the way the claims
    fixture's registration (lab.claims_fixture's CLAIMS_COMBOS/
    CLAIMS_DEGRADATIONS) produces it via exp_disco_degrade: a single
    fixed combo ('claims_known'), degradation_dim holding the ablation
    target name instead of a generic dimension.'''
    rows = []
    for target in ('assess', 'loop_block'):
        for level in (0.0, 0.5):
            rows.append({
                'log': 'claims', 'combo': 'claims_known', 'degradation_dim': target,
                'degradation_level': level, 'status': 'ok',
                'weight_coverage': 0.8, 'skipprob': 0.2, 'salign_coverage': 0.98,
                'alignment_coverage_pn': 0.95, 'voidmass_subprocess': 0.05 * level,
                'voidmass_process': 0.01 * level,
            })
    rows.append({'log': 'claims', 'combo': 'claims_known', 'degradation_dim': 'assess',
                 'degradation_level': 0.9, 'status': 'error: boom'})
    return pd.DataFrame(rows)


def fake_node_df():
    '''Shape of exp_disco_degrade's *_nodes.csv (see lab.exp_disco_degrade
    _node_rows/PER_NODE_METRIC_KEYS): one row per (log, combo,
    degradation_dim, degradation_level, node_id), node_skip_prob not
    skipprob, no status column at all (a node CSV only ever holds
    successful cells - see run_disco_degrade).'''
    def row(node_id, node_type, weight_coverage, node_skip_prob,
             voidmass_subprocess, voidmass_process, level=0.0):
        return {
            'log': 'fake_log', 'combo': 'inductive', 'degradation_dim': 'activity',
            'degradation_level': level, 'node_id': node_id, 'node_type': node_type,
            'alphabet': 'a', 'weight_coverage': weight_coverage,
            'weight_voidage': 1 - weight_coverage, 'node_skip_prob': node_skip_prob,
            'salign_coverage': 0.9, 'voidmass_deficit': 0.1, 'voidmass_movecount': 1.0,
            'voidmass_subprocess': voidmass_subprocess, 'voidmass_process': voidmass_process,
            'alignment_coverage_pn': 0.9, 'mandatory_node_count': 1, 'total_node_count': 1,
        }
    return pd.DataFrame([
        row('1', 'Activity', 1.0, 0.0, 0.0, 1.0),
        row('2', 'Activity', 0.5, 0.5, 0.5, 0.5),
        # a Tau row at the same cell, with values far outside the two
        # Activity rows' range - must not pull the average toward it
        row('3', 'Tau', 0.0, 1.0, 1.0, 0.0),
        # a second cell (level=1.0) that must be dropped entirely
        row('1', 'Activity', 0.0, 1.0, 1.0, 0.0, level=1.0),
    ])


class AverageOverNodesTest(unittest.TestCase):
    def test_excludes_tau_nodes_from_the_average(self):
        averaged = average_over_nodes(fake_node_df())
        row = averaged[(averaged['degradation_dim'] == 'activity')
                        & (averaged['degradation_level'] == 0.0)].iloc[0]
        # mean of the two Activity rows only (1.0, 0.5) -> 0.75, not
        # pulled toward the Tau row's 0.0
        self.assertAlmostEqual(row['weight_coverage'], 0.75)
        self.assertAlmostEqual(row['skipprob'], 0.25)

    def test_renames_node_skip_prob_to_skipprob(self):
        averaged = average_over_nodes(fake_node_df())
        self.assertIn('skipprob', averaged.columns)
        self.assertNotIn('node_skip_prob', averaged.columns)

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
        df['voidmass_subprocess'] = df['voidmass_subprocess'].astype(object)
        df.loc[df['degradation_level'] == 1.0, 'voidmass_subprocess'] = 'not a number'
        out_dir = tempfile.mkdtemp()
        written = plot_dose_response(df, out_dir=out_dir)
        self.assertTrue(Path(written[0]).exists())


if __name__ == '__main__':
    unittest.main()
