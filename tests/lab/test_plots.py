'''
Smoke tests for lab.plots - that plotting functions run and write the
expected file for their experiment's shape; not pixel-checking the
rendered plot. _exclude_degenerate is tested directly (a plain
DataFrame filter, no matplotlib involved) rather than trying to inspect
what actually got plotted.
'''

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from lab.plots import plot_claims_degrade, plot_dose_response, _exclude_degenerate


def fake_df():
    rows = []
    for target in ('assess', 'loop_block'):
        for n in (0, 2, 4):
            rows.append({
                'target': target, 'n_drop_cases': n, 'status': 'ok',
                'weight_coverage': 0.8, 'skipprob': 0.2,
                'salign_coverage': 0.98, 'alignment_coverage_pn': 0.95,
                'voidmass_subprocess': 0.05 * n, 'voidmass_process': 0.01 * n,
            })
    rows.append({'target': 'assess', 'n_drop_cases': 6, 'status': 'error: boom'})
    return pd.DataFrame(rows)


class PlotClaimsDegradeTest(unittest.TestCase):
    def test_writes_one_png(self):
        out_dir = tempfile.mkdtemp()
        written = plot_claims_degrade(fake_df(), out_dir=out_dir)
        self.assertEqual(len(written), 1)
        self.assertTrue(Path(written[0]).exists())
        self.assertTrue(str(written[0]).endswith('claims_degrade.png'))

    def test_error_rows_are_excluded_without_crashing(self):
        # the error row (n_drop_cases=6, no metric columns) must not
        # blow up plotting the other, valid rows
        out_dir = tempfile.mkdtemp()
        written = plot_claims_degrade(fake_df(), out_dir=out_dir)
        self.assertTrue(Path(written[0]).exists())


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
    def test_writes_one_png_per_log_and_dim(self):
        out_dir = tempfile.mkdtemp()
        written = plot_dose_response(fake_disco_df(), out_dir=out_dir)
        self.assertEqual(len(written), 1)  # one (log, dim) pair
        self.assertTrue(Path(written[0]).exists())
        self.assertTrue(str(written[0]).endswith('fake_log_activity.png'))

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
