'''
lab.fixture_metric_table's comparison step: which (fixture, node,
metric) values differ between two tables.

The scoring step itself is the real pipeline on the real fixtures, and
is covered by the per-metric test modules; what needs its own test is
the diff, since a comparison that silently drops a column would report
a change as no change.
'''

import unittest

import pandas as pd

from lab.fixture_metric_table import FIXTURES, diff_rows


def table(**metrics):
    rows = [{'fixture': 'payment_partial', 'node_id': '9', 'node_type': 'Sequence'},
            {'fixture': 'payment_partial', 'node_id': '1', 'node_type': 'Activity'}]
    for name, values in metrics.items():
        for row, value in zip(rows, values):
            row[name] = value
    return pd.DataFrame(rows)


class DiffRowsTest(unittest.TestCase):

    def test_identical_tables_differ_nowhere(self):
        before = table(voidsalign3=[0.25, 0.0])
        self.assertTrue(diff_rows(before, before.copy()).empty)

    def test_a_changed_value_is_reported_with_its_size(self):
        moved = diff_rows(table(voidsalign3=[0.25, 0.0]),
                          table(voidsalign3=[0.5, 0.0]))
        self.assertEqual(len(moved), 1)
        row = moved.iloc[0]
        self.assertEqual((row.node_id, row.metric), ('9', 'voidsalign3'))
        self.assertAlmostEqual(row.before, 0.25)
        self.assertAlmostEqual(row.after, 0.5)
        self.assertAlmostEqual(row.change, 0.25)

    def test_floating_point_noise_is_not_a_change(self):
        moved = diff_rows(table(voidsalign3=[1 / 3, 0.0]),
                          table(voidsalign3=[1 / 3 + 1e-15, 0.0]))
        self.assertTrue(moved.empty)

    def test_every_metric_column_is_compared(self):
        moved = diff_rows(table(voidsalign3=[0.25, 0.0], voidsat2=[0.25, 0.0]),
                          table(voidsalign3=[0.25, 0.0], voidsat2=[0.75, 0.0]))
        self.assertEqual(list(moved['metric']), ['voidsat2'])

    def test_a_metric_the_earlier_table_lacks_is_called_out(self):
        '''Not silently skipped: a new metric column is a difference
        between the two tables, even though no value can be compared.'''
        moved = diff_rows(table(voidsalign3=[0.25, 0.0]),
                          table(voidsalign3=[0.25, 0.0], voidsat2=[0.25, 0.0]))
        self.assertEqual(list(moved['metric']), ['voidsat2'])
        self.assertEqual(list(moved['note']), ['absent before'])


class FixtureCatalogueTest(unittest.TestCase):
    '''The named fixtures build, and build a log and a tree - checked
    without scoring anything, so this costs no alignment search.'''

    def test_every_fixture_builds_a_log_and_a_tree(self):
        for name, (build_log, build_tree) in FIXTURES.items():
            with self.subTest(fixture=name):
                log, tree = build_log(), build_tree()
                self.assertTrue(len(log) > 0)
                self.assertTrue(tree.children)


if __name__ == '__main__':
    unittest.main()
