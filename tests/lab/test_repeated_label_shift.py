"""
lab.repeated_label_shift: per-node voidmass_process under the
label-based and leaf-based computations, from one pass of alignments.
"""

import unittest

from hypothesis import HealthCheck, given, settings
from skipalignments import Activity, Loop, Sequence, Xor

from lab.repeated_label_shift import LABEL, LEAF, shift_table
from process_voids.voidmass_pn import build_id_net, voidmass_table_pn
from tests.process_voids.test_voidmass_process_additivity import (
    models_logs_and_cuts, models_with_repeated_labels, variant_probs_of)

TOLERANCE = 1e-9


def _leaf(label, node_id):
    node = Activity(None, label, 100000)
    node.id = node_id
    return node


def _join(operator, node_id, *children):
    node = operator(None, list(children))
    node.id = node_id
    for child in children:
        child.set_parent(node)
    return node


class ShiftTest(unittest.TestCase):

    def test_seq_a_a_moves_its_leaves_but_not_its_root(self):
        first, second = _leaf('a', 'first'), _leaf('a', 'second')
        tree = _join(Sequence, 'root', first, second)
        shift = shift_table(tree, {('a',): 1.0})
        self.assertAlmostEqual(shift.voidmass_process(LABEL, tree, tree), 0.5)
        self.assertAlmostEqual(shift.voidmass_process(LEAF, tree, tree), 0.5)
        for leaf in (first, second):
            self.assertAlmostEqual(shift.voidmass_process(LABEL, leaf, tree), 0.5)
            self.assertAlmostEqual(shift.voidmass_process(LEAF, leaf, tree), 0.25)

    def test_the_root_moves_where_deduplication_merged_distinct_alignments(self):
        """
        xor(loop(b, a), b) against <a> has three tied optimal alignments:
        through the loop, b-a-b (deficit 2 of 3 moves); the loop's b
        skipped with a as a log move (1 of 1); and the xor's other b
        skipped the same way (1 of 1). The last two differ only in which
        b fired, so a label signature merges them and averages two stories
        - 3/2 over 2, 0.75 - where the definition averages all three, 4/3
        over 5/3, 0.8.
        """
        b_in_loop, a, b_in_xor = _leaf('b', 'b_loop'), _leaf('a', 'a'), _leaf('b', 'b_xor')
        tree = _join(Xor, 'root', _join(Loop, 'loop', b_in_loop, a), b_in_xor)
        shift = shift_table(tree, {('a',): 1.0})
        self.assertAlmostEqual(shift.voidmass_process(LABEL, tree, tree), 0.75)
        self.assertAlmostEqual(shift.voidmass_process(LEAF, tree, tree), 0.8)

    @settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    @given(models_logs_and_cuts())
    def test_with_distinct_labels_the_two_agree_at_every_node(self, case):
        tree, traces, _cut = case
        shift = shift_table(tree, variant_probs_of(traces))
        for node in shift.label:
            self.assertAlmostEqual(shift.voidmass_process(LABEL, node, tree),
                                   shift.voidmass_process(LEAF, node, tree), delta=TOLERANCE)

    @settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    @given(models_with_repeated_labels())
    def test_leaf_based_is_what_voidmass_table_pn_computes(self, case):
        """The leaf-based side stands for the fixed implementation only if
        it matches it. (Before the fix, the label-based side was checked
        against voidmass_table_pn the same way, and matched; that is what
        makes it a faithful record of the old values.)"""
        tree, traces, _cut = case
        probs = variant_probs_of(traces)
        shift = shift_table(tree, probs)
        net, im, fm, activity_to_id, tau_ids, loops = build_id_net(tree)
        result = voidmass_table_pn(tree, probs, net, im, fm, activity_to_id, tau_ids,
                                   id_loop_list=loops, timeout=30)
        if result.timed_out_count or shift.timed_out:
            return
        for node, row in result.table.items():
            self.assertAlmostEqual(shift.voidmass_process(LEAF, node, tree),
                                   row['voidmass_process_lower'], delta=TOLERANCE,
                                   msg=f'node {node.id} of {tree}, traces {traces}')


if __name__ == '__main__':
    unittest.main()
