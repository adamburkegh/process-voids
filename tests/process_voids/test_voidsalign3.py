'''
Tests for process_voids.voidsalign3 (Definitions [Coverage by Skip
Alignment Correspondence] / [Void by Skip Alignment Correspondence]).

voidsalign3 differs from voidsalign2 in one term only: smovecount
weights a skip move by the number of labelled activities in the
subprocess it skips, |leaves(msub') \\ {silent}|, rather than by
aligncost(<>, msub'). Everything else - the observation conditioning,
the 1/|Upsilon| and 1/|P| averages, the complement - is voidsalign2's
arithmetic unchanged, via coveragemass.observed_alignment_mass. So
these tests concentrate on where the two weights disagree, and pin the
rest end to end rather than re-deriving it.

The weights disagree wherever a subprocess has a traversal performing
fewer labelled activities than it contains:

    skipped subprocess   aligncost(<>, .)   labelled leaves
    a                            1                 1
    seq(b, c)                    2                 2
    xor(a, tau)                  0                 1
    xor(a, b)                    1                 2
    loop(a, e)                   1                 2
'''

import unittest
from types import SimpleNamespace

from skipalignments import Activity, Aligner, Loop, Sequence, Skip, Tau, Xor

from process_voids.voidsalign3 import (
    coversalign3, observed_skip_weighted_mass, smatchcount, smovecount, voidsalign3,
)

ACT_COST = 100000

Aligner.set_level_incentive(0)


def leaf(cls, name, node_id, cost=ACT_COST):
    node = cls(None, name, cost)
    node.id = node_id
    return node


def composite(cls, children, node_id):
    node = cls(None, children)
    node.id = node_id
    for child in children:
        child.set_parent(node)
    return node


def variant_key(trace):
    return ', '.join(trace)


def align(tree, trace):
    states, _ = Aligner(tree).align_normal_form(
        list(trace), [ACT_COST] * len(trace), True, timeout=100)
    return states


def build(tree, traces_with_weights):
    skip_dict, variant_probs = {}, {}
    for trace, weight in traces_with_weights.items():
        skip_dict[variant_key(trace)] = align(tree, trace)
        variant_probs[trace] = weight
    return skip_dict, variant_probs


def skip_of(node):
    return [('>>', Skip(node, ACT_COST))]


class SkipWeightTest(unittest.TestCase):
    '''smovecount of an execution that is a single skip move over each
    shape in the module docstring's table.'''

    def test_a_single_activity_weighs_one(self):
        a = leaf(Activity, 'a', '1')
        self.assertEqual(smovecount(skip_of(a)), 1)

    def test_a_sequence_weighs_its_activities(self):
        b, c = leaf(Activity, 'b', '1'), leaf(Activity, 'c', '2')
        self.assertEqual(smovecount(skip_of(composite(Sequence, [b, c], '3'))), 2)

    def test_an_optional_block_still_weighs_its_activity(self):
        '''
        The case the definition exists for. xor(a, tau) has a traversal
        performing no labelled activity, so aligncost(<>, .) is 0 and a
        skip over it would carry no weight at all - absence of an
        optional block would be invisible in the mass.
        '''
        a, tau = leaf(Activity, 'a', '1'), leaf(Tau, 'skip', '2', cost=0)
        self.assertEqual(smovecount(skip_of(composite(Xor, [a, tau], '3'))), 1)

    def test_a_choice_weighs_every_branch_activity(self):
        a, b = leaf(Activity, 'a', '1'), leaf(Activity, 'b', '2')
        self.assertEqual(smovecount(skip_of(composite(Xor, [a, b], '3'))), 2)

    def test_a_loop_weighs_its_redo_activity_too(self):
        a, e = leaf(Activity, 'a', '1'), leaf(Activity, 'e', '2')
        self.assertEqual(smovecount(skip_of(composite(Loop, [a, e], '3'))), 2)

    def test_a_wholly_silent_subprocess_weighs_nothing(self):
        tau = leaf(Tau, 'skip', '1', cost=0)
        self.assertEqual(smovecount(skip_of(tau)), 0)

    def test_synchronous_moves_count_once_each(self):
        a, b = leaf(Activity, 'a', '1'), leaf(Activity, 'b', '2')
        execution = [('a', a), ('>>', Skip(b, ACT_COST))]
        self.assertEqual(smatchcount(execution), 1)
        self.assertEqual(smovecount(execution), 2)


class MassTest(unittest.TestCase):

    def test_a_lumped_loop_skip_is_weighted_by_its_leaves(self):
        '''
        seq(x, loop(a, e)) with x recorded and the loop skipped as one
        lumped move, hand-built: the root's execution is one synchronous
        move plus a skip of loop(a, e), which weighs 2 - so the ratio is
        1/3, where aligncost(<>, loop) = 1 would give 1/2.
        '''
        a, e = leaf(Activity, 'a', '1'), leaf(Activity, 'e', '2')
        loop = composite(Loop, [a, e], '3')
        x = leaf(Activity, 'x', '4')
        tree = composite(Sequence, [x, loop], '5')
        path = [('x', x), ('>>', Skip(loop, ACT_COST))]
        skip_dict = {'v': [SimpleNamespace(path=path)]}
        variant_probs = {('v',): 1.0}
        # root execution: 1 sync (x) plus a skip of loop(a, e) weighing 2
        self.assertAlmostEqual(
            observed_skip_weighted_mass(tree, skip_dict, variant_probs), 1 / 3, places=9)

    def test_a_fully_recorded_subprocess_has_mass_one(self):
        b, c = leaf(Activity, 'b', '1'), leaf(Activity, 'c', '2')
        bc = composite(Sequence, [b, c], '3')
        a = leaf(Activity, 'a', '4')
        tree = composite(Sequence, [a, bc], '5')
        skip_dict, variant_probs = build(tree, {('a', 'b', 'c'): 1.0})
        self.assertAlmostEqual(
            observed_skip_weighted_mass(bc, skip_dict, variant_probs), 1.0, places=9)


class CoverageAndVoidTest(unittest.TestCase):

    def setUp(self):
        self.b, self.c = leaf(Activity, 'b', '1'), leaf(Activity, 'c', '2')
        self.bc = composite(Sequence, [self.b, self.c], '3')
        self.a = leaf(Activity, 'a', '4')
        self.tree = composite(Sequence, [self.a, self.bc], '5')

    def _probs(self, value):
        return {node: value for node in (self.tree, self.a, self.bc, self.b, self.c)}

    def test_void_is_one_minus_coverage(self):
        skip_dict, variant_probs = build(self.tree, {('a', 'b'): 1.0})
        probs = self._probs(0.25)
        self.assertAlmostEqual(
            voidsalign3(self.bc, skip_dict, variant_probs, probs),
            1 - coversalign3(self.bc, skip_dict, variant_probs, probs), places=9)

    def test_never_missing_reads_zero_void(self):
        skip_dict, variant_probs = build(self.tree, {('a', 'b', 'c'): 1.0})
        self.assertAlmostEqual(
            voidsalign3(self.bc, skip_dict, variant_probs, self._probs(0.0)), 0.0, places=9)

    def test_always_missing_reads_one_void(self):
        skip_dict, variant_probs = build(self.tree, {('a',): 1.0})
        probs = {self.tree: 0.0, self.a: 0.0, self.bc: 1.0, self.b: 1.0, self.c: 1.0}
        self.assertAlmostEqual(
            voidsalign3(self.bc, skip_dict, variant_probs, probs), 1.0, places=9)

    def test_missing_from_half_the_traces_reads_half_void(self):
        skip_dict, variant_probs = build(self.tree, {('a', 'b', 'c'): 0.5, ('a',): 0.5})
        probs = {self.tree: 0.0, self.a: 0.0, self.bc: 0.5, self.b: 0.5, self.c: 0.5}
        self.assertAlmostEqual(
            voidsalign3(self.bc, skip_dict, variant_probs, probs), 0.5, places=9)


class RunningExampleTest(unittest.TestCase):
    '''
    The payment running example (lab.fixtures) through the full
    ebi-backed pipeline. Only the root differs from voidsalign2: sigma3
    and sigma4 (<o, s, p>) align with a skip of loop(a, e), which
    voidsalign2 weighs 1 (its do-child) and voidsalign3 weighs 2 (a and
    e). The root ratio of that variant goes 3/4 -> 3/5, the others read 1
    either way, so the root mass goes 11/12 -> 13/15 and its void reading
    0.0833 -> 2/15. The loop's own row is unchanged: its lumped skip has
    no synchronous move, so it is excluded from the loop's mass.
    '''

    @classmethod
    def setUpClass(cls):
        from lab.fixtures import build_running_example_log, build_running_example_tree
        from process_voids import pvoid

        cls.tree = build_running_example_tree()
        cls.o, cls.approval, cls.sched, cls.p = cls.tree.children
        cls.a, cls.e = cls.approval.children
        cls.s, cls.tau = cls.sched.children
        cls.dv = pvoid.skipprob(build_running_example_log(), cls.tree,
                                'var/lab/test_voidsalign3.slpn')

    def _void(self, node):
        return voidsalign3(node, self.dv.skip_dict_backup, self.dv.pl, self.dv.skip_probs)

    def test_root_weighs_the_lumped_loop_skip_by_both_its_activities(self):
        mass = observed_skip_weighted_mass(self.tree, self.dv.skip_dict_backup, self.dv.pl)
        self.assertAlmostEqual(mass, 13 / 15, places=9)
        self.assertAlmostEqual(self._void(self.tree), 2 / 15, places=9)

    def test_rows_without_a_lumped_skip_match_voidsalign2(self):
        cases = [('o', 'o', 0.0), ('loop', 'approval', 1 / 3), ('a', 'a', 0.0),
                 ('e', 'e', 0.0), ('xor', 'sched', 1 / 6), ('s', 's', 0.0),
                 ('p', 'p', 0.0)]
        for name, attr, expected in cases:
            with self.subTest(node=name):
                self.assertAlmostEqual(self._void(getattr(self, attr)), expected, places=9)


if __name__ == '__main__':
    unittest.main()
