"""
lab.normal_form_compare: rescoring the void metrics under the reverse
normal-form precedence.

Driven by hand-built alignments on the §2 witness of the normal-form
order brief, so the comparison machinery is tested without an ebi run.
The values themselves are pinned in
tests/process_voids/test_voidsat2.py's NormalFormPrecedencePinTest.
"""

import unittest
from datetime import datetime, timedelta

from skipalignments import Activity, And, Sequence, Skip

from lab.normal_form_compare import comparison_rows, independent_pairs, reordered
from process_voids.coveragemass import _variant_key
from lab.normal_form import SKIP_FIRST, move_kind

ACT_COST = 100000
T0 = datetime(2026, 1, 1, 12, 0, 0)

TRACE = [{'concept:name': 'a', 'time:timestamp': T0},
         {'concept:name': 'c', 'time:timestamp': T0 + timedelta(hours=10)}]
KEY = _variant_key(('a', 'c'))


def build_witness():
    """seq( and(a, b), c ), and its aligner-shaped normal form for <a, c>."""
    a, b, c = (Activity(None, name, ACT_COST) for name in 'abc')
    a.id, b.id, c.id = '1', '2', '3'
    parallel = And(None, [a, b])
    parallel.id = '4'
    a.set_parent(parallel)
    b.set_parent(parallel)
    tree = Sequence(None, [parallel, c])
    tree.id = '5'
    parallel.set_parent(tree)
    c.set_parent(tree)
    path = [('a0', a), ('>>', Skip(b, ACT_COST)), ('c1', c)]
    return tree, parallel, path


class ComparisonTest(unittest.TestCase):

    def setUp(self):
        self.tree, self.parallel, path = build_witness()
        self.paths = {KEY: [path]}
        self.variant_probs = {('a', 'c'): 1.0}
        self.skip_probs = {node: 0.0 for node in
                           [self.tree, self.parallel, *self.parallel.children,
                            self.tree.children[1]]}

    def test_the_witness_has_one_independent_pair_of_differing_kinds(self):
        pairs = independent_pairs(self.paths)
        self.assertEqual([(key, i, kinds) for key, i, kinds in pairs],
                         [(KEY, 0, ('sync', 'skip'))])

    def test_reordering_rewrites_only_what_the_precedence_demands(self):
        flipped = reordered(self.paths, self.tree, SKIP_FIRST)[KEY][0]
        self.assertEqual([move_kind(m) for m in flipped], ['skip', 'sync', 'sync'])

    def test_voidsat_moves_and_voidsalign_does_not(self):
        rows = comparison_rows(self.tree, [TRACE], self.paths,
                               self.variant_probs, self.skip_probs)
        at_parallel = rows[rows['node_id'] == self.parallel.id].set_index('metric')

        self.assertAlmostEqual(at_parallel.loc['voidsalign3', 'change'], 0.0, places=12)
        self.assertAlmostEqual(at_parallel.loc['voidsat2', 'log_first'], 1.0, places=12)
        self.assertAlmostEqual(at_parallel.loc['voidsat2', 'skip_first'], 0.0, places=12)

    def test_a_row_per_node_and_metric(self):
        rows = comparison_rows(self.tree, [TRACE], self.paths,
                               self.variant_probs, self.skip_probs)
        self.assertEqual(len(rows), 5 * 2)
        self.assertEqual(set(rows['metric']), {'voidsalign3', 'voidsat2'})


if __name__ == '__main__':
    unittest.main()
