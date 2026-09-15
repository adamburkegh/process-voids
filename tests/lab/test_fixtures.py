"""
Tests for lab.fixtures' payment_partial fixture.

Its purpose is a case no other fixture here produces: a subprocess
TRAVERSED BUT RECORDED INCOMPLETELY, so a partial term has something to
fire on. Every absence in payment_approval is a whole traversal - a
subprocess either runs and is fully recorded, or does not run - so
skip probability alone accounts for all of it.

sigma7 is the discriminating trace. A loop's language is a (e a)*, so
<a, e> is not a complete traversal: it needs a closing a. These pin that
the aligner really does produce that traversal, that the executions
machinery attributes its model move to the loop, and that a
product-form metric reads 0 at that node while a complement-form one
does not.
"""

import unittest

from skipalignments import Activity, Aligner, Loop, Sequence, Skip, Tau, Xor

from lab.fixtures import (build_partial_sequence_log, build_partial_sequence_tree,
                          build_payment_partial_log, build_payment_partial_tree,
                          build_running_example_log)
from process_voids import pvoid
from process_voids.coveragemass import (deficit, executions, matchcount, movecount,
                                        voidsat)
from process_voids.voidsalign2 import voidsalign2

ACT_COST = 100000

Aligner.set_level_incentive(0)


def align(tree, trace):
    states, _ = Aligner(tree).align_normal_form(
        list(trace), [ACT_COST] * len(trace), True, timeout=100)
    return states


class TreeShapeTest(unittest.TestCase):

    def setUp(self):
        self.tree = build_payment_partial_tree()

    def test_is_the_discovered_shape(self):
        """seq(o, xor(tau, loop(a, e)), s, p) - approval optional under a
        silent branch, s mandatory. Not build_running_example_tree's
        seq(o, loop(a, e), xor(s, tau), p)."""
        o, choice, s, p = self.tree.children
        self.assertIsInstance(self.tree, Sequence)
        self.assertEqual([n.name for n in (o, s, p)], ['o', 's', 'p'])
        self.assertIsInstance(choice, Xor)
        tau, loop = choice.children
        self.assertIsInstance(tau, Tau)
        self.assertIsInstance(loop, Loop)
        self.assertEqual([n.name for n in loop.children], ['a', 'e'])

    def test_keeps_every_payment_approval_trace_and_adds_two(self):
        existing = build_running_example_log()
        partial = build_payment_partial_log()
        self.assertEqual(partial['case:concept:name'].nunique(), 8)
        for case in existing['case:concept:name'].unique():
            with self.subTest(case=case):
                self.assertEqual(
                    list(partial[partial['case:concept:name'] == case]['concept:name']),
                    list(existing[existing['case:concept:name'] == case]['concept:name']))

    def test_the_two_new_traces_are_one_variant_with_different_timings(self):
        """sigma7 and sigma8 share a variant, so only the duration-based
        metrics can tell them apart - the role sigma3/sigma4 play for the
        existing fixture."""
        log = build_payment_partial_log()
        spans = {}
        for case in ('sigma7', 'sigma8'):
            events = log[log['case:concept:name'] == case]
            spans[case] = (tuple(events['concept:name']),
                           (events['time:timestamp'].max()
                            - events['time:timestamp'].min()).total_seconds() / 3600)
        self.assertEqual(spans['sigma7'][0], spans['sigma8'][0])
        self.assertEqual((spans['sigma7'][1], spans['sigma8'][1]), (11.0, 7.0))


class PartialTraversalTest(unittest.TestCase):
    """
    The aligner's own output for sigma7's variant, which the fixture's
    value depends on - checked rather than assumed.
    """

    def setUp(self):
        self.tree = build_payment_partial_tree()
        self.o, self.choice, self.s, self.p = self.tree.children
        self.tau, self.loop = self.choice.children
        self.a, self.e = self.loop.children
        self.states = align(self.tree, ('o', 'a', 'e', 's', 'p'))

    def _loop_executions(self):
        return [ex for state in self.states for ex in executions(state.path, self.loop)]

    def test_a_partial_traversal_is_among_the_optimal_alignments(self):
        """
        One alignment pays a model move for the loop's closing a: two
        synchronous moves and one model move in a single traversal, so
        the loop ran and a third of its moves went unobserved.
        """
        partial = [ex for ex in self._loop_executions()
                   if matchcount(ex) == 2 and movecount(ex) == 3]
        self.assertEqual(len(partial), 1)
        self.assertEqual(deficit(partial[0]), 1)

    def test_the_model_move_is_a_skip_of_the_loops_own_a(self):
        """If this move were attributed elsewhere, or dropped, the
        partial term would read 0 - that would be a defect in the
        executions machinery, not a property of the fixture."""
        partial = next(ex for ex in self._loop_executions() if movecount(ex) == 3)
        skips = [model_elem for _log_elem, model_elem in partial
                 if isinstance(model_elem, Skip)]
        self.assertEqual(len(skips), 1)
        self.assertIs(skips[0].node, self.a)

    def test_the_partial_traversal_is_tied_with_a_log_move_explanation(self):
        """
        The aligner also finds an alignment that discards e as a log
        move, leaving the loop a single synchronous move. Both cost one
        deviation, so both are optimal - anything averaging over tied
        alignments sees the partial traversal at half weight.
        """
        self.assertEqual(len(self.states), 2)
        shapes = sorted((matchcount(ex), movecount(ex)) for ex in self._loop_executions())
        self.assertEqual(shapes, [(1, 1), (2, 3)])


class PipelineTest(unittest.TestCase):
    """The fixture through the full ebi-backed pipeline: the partial
    term has to survive real skip-probability estimation, not just
    hand-built alignments."""

    @classmethod
    def setUpClass(cls):
        cls.tree = build_payment_partial_tree()
        cls.o, cls.choice, cls.s, cls.p = cls.tree.children
        cls.tau, cls.loop = cls.choice.children
        cls.a, cls.e = cls.loop.children
        cls.log = build_payment_partial_log()
        cls.dv = pvoid.skipprob(cls.log, cls.tree, 'var/lab/test_payment_partial.slpn')

    def test_the_loop_is_never_skipped(self):
        """It is traversed whenever reached, so skip probability alone
        reports nothing about it."""
        self.assertEqual(self.dv.skip_probs[self.loop], 0.0)

    def test_the_loop_still_carries_a_deficit(self):
        """What skip probability misses: pooled over every execution of
        every optimal alignment, the loop is short of moves."""
        total_deficit = 0.0
        for variant, weight in self.dv.pl.items():
            states = self.dv.skip_dict_backup.get(', '.join(variant), [])
            share = weight / len(states) if states else 0.0
            for state in states:
                for ex in executions(state.path, self.loop):
                    total_deficit += share * deficit(ex)
        self.assertGreater(total_deficit, 0.0)

    def test_a_product_form_metric_reads_zero_where_a_complement_form_does_not(self):
        """
        The point of the fixture. voidsat is skip_prob * mass, so a
        subprocess that always ran reads exactly 0 however incompletely
        it was recorded. voidsalign2 is 1 - (1 - skip_prob) * mass, so
        the same incompleteness shows up in the mass.
        """
        self.assertEqual(
            voidsat(self.loop, self.tree, self.dv, self.log), 0.0)
        self.assertGreater(
            voidsalign2(self.loop, self.dv.skip_dict_backup, self.dv.pl,
                        self.dv.skip_probs), 0.0)


class PartialSequenceTest(unittest.TestCase):
    """
    partial_sequence: seq(o, seq(x, y, z), p), the fixture that pins the
    definitional property rather than field behaviour. The inner
    sequence is traversed in every trace, so its skip probability is 0
    throughout and only its completeness varies - and every trace has
    exactly one optimal alignment, so nothing is diluted by ties.
    """

    @classmethod
    def setUpClass(cls):
        cls.tree = build_partial_sequence_tree()
        cls.o, cls.inner, cls.p = cls.tree.children
        cls.x, cls.y, cls.z = cls.inner.children
        cls.log = build_partial_sequence_log()
        cls.dv = pvoid.skipprob(cls.log, cls.tree, 'var/lab/test_partial_sequence.slpn')

    def test_each_trace_has_exactly_one_optimal_alignment(self):
        """
        The property payment_partial's loop cannot offer. Completing this
        subprocess costs one model move per unrecorded step, while
        discarding what WAS recorded to skip the block wholesale costs
        more, so the partial explanation wins outright rather than tying.
        """
        for trace in (('o', 'x', 'y', 'z', 'p'), ('o', 'x', 'z', 'p'), ('o', 'x', 'p')):
            with self.subTest(trace=trace):
                self.assertEqual(len(align(self.tree, trace)), 1)

    def test_the_three_traces_give_partial_zero_a_third_and_two_thirds(self):
        cases = [(('o', 'x', 'y', 'z', 'p'), 3, 0),
                 (('o', 'x', 'z', 'p'), 2, 1),
                 (('o', 'x', 'p'), 1, 2)]
        for trace, expected_match, expected_deficit in cases:
            with self.subTest(trace=trace):
                (state,) = align(self.tree, trace)
                (execution,) = executions(state.path, self.inner)
                self.assertEqual(matchcount(execution), expected_match)
                self.assertEqual(movecount(execution), 3)
                self.assertEqual(deficit(execution), expected_deficit)

    def test_the_subprocess_is_never_skipped(self):
        self.assertEqual(self.dv.skip_probs[self.inner], 0.0)

    def test_a_product_form_metric_reads_zero_where_a_complement_form_reads_a_third(self):
        """
        Undiluted, so the complement form's value is exactly derivable:
        the observed masses are 1, 2/3 and 1/3 over three equally
        weighted variants, so the mass is 2/3 and the void reading 1/3.
        """
        self.assertEqual(voidsat(self.inner, self.tree, self.dv, self.log), 0.0)
        self.assertAlmostEqual(
            voidsalign2(self.inner, self.dv.skip_dict_backup, self.dv.pl,
                        self.dv.skip_probs),
            1 / 3, places=9)


if __name__ == '__main__':
    unittest.main()
