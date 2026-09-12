import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from lab.discovery import DiscoveryCombo, DiscoveryResult
from lab.exp_surprise import _node_rows, run_surprise


class FakeNode:
    """Minimal stand-in for a ProcessTree node - _node_rows only needs
    .id and .get_leaf_labels()."""
    def __init__(self, node_id, leaf_labels):
        self.id = node_id
        self._leaf_labels = leaf_labels

    def get_leaf_labels(self):
        return self._leaf_labels


class NodeRowsTest(unittest.TestCase):
    """self and baseline land in separate columns on the same row - not
    a shared id plus a distribution column, which is what this was
    refactored away from (see lab.exp_surprise's module docstring)."""

    def test_self_and_baseline_are_separate_columns_on_one_row(self):
        n1 = FakeNode('n1', ['a'])
        self_totals = ({n1: 1.0}, {n1: 2.0})
        baseline_totals = ({n1: 3.0}, {n1: 4.0})
        rows = _node_rows('log', 'combo', 'dim', 0.0, self_totals, baseline_totals)
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertNotIn('distribution', row)
        self.assertEqual(row['containment_bits'], 1.0)
        self.assertEqual(row['predecessor_bits'], 2.0)
        self.assertEqual(row['containment_bits_baseline'], 3.0)
        self.assertEqual(row['predecessor_bits_baseline'], 4.0)

    def test_node_present_on_only_one_side_still_gets_a_row(self):
        n1 = FakeNode('n1', ['a'])
        n2 = FakeNode('n2', ['b'])
        self_totals = ({n1: 1.0}, {})
        baseline_totals = ({n2: 5.0}, {})
        rows = _node_rows('log', 'combo', 'dim', 0.0, self_totals, baseline_totals)
        by_id = {row['node_id']: row for row in rows}
        self.assertEqual(set(by_id), {'n1', 'n2'})
        self.assertEqual(by_id['n1']['containment_bits'], 1.0)
        self.assertEqual(by_id['n1']['containment_bits_baseline'], 0.0)
        self.assertEqual(by_id['n2']['containment_bits'], 0.0)
        self.assertEqual(by_id['n2']['containment_bits_baseline'], 5.0)


FAKE_SELF_FIELDS = {'n_events': 10, 'headline_bits': 3.0, 'bits_per_event': 0.3,
                     'ambiguous_event_count': 1, 'unattributable_event_count': 0,
                     'out_of_alphabet_event_count': 0}
FAKE_BASELINE_FIELDS = {'n_events': 10, 'headline_bits': 5.0, 'bits_per_event': 0.5,
                         'ambiguous_event_count': 1, 'unattributable_event_count': 0,
                         'out_of_alphabet_event_count': 0}
FAKE_CELL = (({}, {}), FAKE_SELF_FIELDS, ({}, {}), FAKE_BASELINE_FIELDS)


def degrade_stub(log, level):
    if level == 0.0:
        return log, set()
    return f'{log}_degraded', {f'dropped_at_{level}'}


class RunSurpriseColumnsTest(unittest.TestCase):
    """run_surprise's summary output uses separate self/baseline columns,
    not a shared metric plus a 'distribution' column - see NodeRowsTest's
    docstring for why."""

    def setUp(self):
        self.tmp_out = Path(tempfile.mkdtemp()) / 'nodes.csv'
        self.tmp_summary = Path(tempfile.mkdtemp()) / 'summary.csv'
        self.combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult('FAKE_TREE'))}
        self.degradations = {'trace': degrade_stub}

    def test_summary_row_has_baseline_columns_not_a_distribution_column(self):
        with patch('lab.exp_surprise.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_surprise.observed_intervals', return_value={}), \
             patch('lab.exp_surprise.discover_cached', return_value=('FAKE_TREE', None)), \
             patch('lab.exp_surprise.compute_predecessors', return_value={}), \
             patch('lab.exp_surprise._compute_cell', return_value=FAKE_CELL):
            node_df, summary_df = run_surprise(
                ['fake_log.xes'], combos=self.combos, degradations=self.degradations,
                levels=[0.0], out_csv=str(self.tmp_out), summary_csv=str(self.tmp_summary))

            self.assertNotIn('distribution', summary_df.columns)
            self.assertNotIn('distribution', node_df.columns)
            self.assertIn('headline_bits', summary_df.columns)
            self.assertIn('headline_bits_baseline', summary_df.columns)

            row = summary_df.iloc[0]
            self.assertEqual(row['status'], 'ok')
            self.assertEqual(row['headline_bits'], 3.0)
            self.assertEqual(row['bits_per_event'], 0.3)
            self.assertEqual(row['headline_bits_baseline'], 5.0)
            self.assertEqual(row['bits_per_event_baseline'], 0.5)

    def test_one_summary_row_per_cell_not_two(self):
        """One summary row per cell, carrying both the self and baseline
        columns - not one row per distribution."""
        with patch('lab.exp_surprise.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_surprise.observed_intervals', return_value={}), \
             patch('lab.exp_surprise.discover_cached', return_value=('FAKE_TREE', None)), \
             patch('lab.exp_surprise.compute_predecessors', return_value={}), \
             patch('lab.exp_surprise._compute_cell', return_value=FAKE_CELL):
            _, summary_df = run_surprise(
                ['fake_log.xes'], combos=self.combos, degradations=self.degradations,
                levels=[0.0], out_csv=str(self.tmp_out), summary_csv=str(self.tmp_summary))
            self.assertEqual(len(summary_df), 1)


if __name__ == '__main__':
    unittest.main()
