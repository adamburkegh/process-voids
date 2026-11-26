"""
Coverage and Void by Aligned Duration (Definition [Coverage and Void by
Aligned Duration]) - voidsat2.

Pinned against lab.fixtures' partial_sequence, which is built for
exactly this: seq(o, seq(x,y,z), p) over three traces recording the
inner sequence completely, without its middle step, and without two of
three. Every trace has one optimal alignment and the inner sequence is
traversed in all of them, so its skip probability is 0 and only
completeness varies - the ratio reads its full magnitude rather than
one of two ties.

Alignments come from the real aligner, not hand-picked paths.

NormalFormPrecedencePinTest pins the aligner behaviour this metric's
values rest on, as CostModelPinTest does for voidmass_pn's bound.
"""

import unittest
from datetime import datetime, timedelta

from skipalignments import Activity, Aligner, And, Sequence

from lab.fixtures import build_partial_sequence_tree, build_payment_partial_tree
from process_voids.coveragemass import _variant_key, block, consumes, log_to_traces, mdur, obsdur, misdur
from process_voids.voidsat2 import adratio, admass, covsat2, voidsat2

ACT_COST = 100000

Aligner.set_level_incentive(0)

T0 = datetime(2026, 1, 1, 12, 0, 0)


def event(name, hours):
    return {'concept:name': name, 'time:timestamp': T0 + timedelta(hours=hours)}


def trace(*name_hours):
    return [event(name, hours) for name, hours in name_hours]


def align(tree, activities):
    states, _ = Aligner(tree).align_normal_form(
        list(activities), [ACT_COST] * len(activities), True, timeout=100)
    return [s.path for s in states]


def alignments_for(tree, traces):
    """{variant_key: [path, ...]} over the distinct variants in `traces`."""
    result = {}
    for t in traces:
        activities = tuple(e['concept:name'] for e in t)
        result[_variant_key(activities)] = align(tree, activities)
    return result


def nodes_of(tree):
    """{node_id: node} for every node in `tree`."""
    found = {}

    def walk(node):
        found[node.id] = node
        for child in getattr(node, 'children', []) or []:
            walk(child)
    walk(tree)
    return found


class PartialSequenceTest(unittest.TestCase):
    """The inner seq(x,y,z) is recorded whole, then missing y, then
    missing y and z. Observed duration should read 1, 2/3 and 1/3."""

    def setUp(self):
        self.tree = build_partial_sequence_tree()
        by_id = nodes_of(self.tree)
        self.inner = by_id['6']
        self.o = by_id['1']
        self.complete = trace(('o', 0), ('x', 1), ('y', 2), ('z', 3), ('p', 4))
        self.one_missing = trace(('o', 0), ('x', 1), ('z', 3), ('p', 4))
        self.two_missing = trace(('o', 0), ('x', 1), ('p', 4))
        self.log = [self.complete, self.one_missing, self.two_missing]
        self.alignments = alignments_for(self.tree, self.log)
        # traversed in every trace, so nothing here is explained by skipping
        self.skip_probs = {node: 0.0 for node in by_id.values()}

    def _ratio(self, node, t):
        activities = tuple(e['concept:name'] for e in t)
        paths = self.alignments[_variant_key(activities)]
        self.assertEqual(len(paths), 1, 'fixture promises one optimal alignment')
        return adratio(node, self.tree, paths[0], t)

    def test_a_fully_recorded_traversal_reads_one(self):
        self.assertAlmostEqual(self._ratio(self.inner, self.complete), 1.0)

    def test_one_missing_step_of_three_reads_two_thirds(self):
        self.assertAlmostEqual(self._ratio(self.inner, self.one_missing), 2 / 3)

    def test_two_missing_steps_of_three_reads_one_third(self):
        self.assertAlmostEqual(self._ratio(self.inner, self.two_missing), 1 / 3)

    def test_mass_is_the_mean_over_the_three_traces(self):
        self.assertAlmostEqual(
            admass(self.inner, self.tree, self.log, self.alignments), (1 + 2 / 3 + 1 / 3) / 3)

    def test_void_is_the_complement_of_coverage(self):
        cov = covsat2(self.inner, self.tree, self.log, self.alignments, self.skip_probs)
        void = voidsat2(self.inner, self.tree, self.log, self.alignments, self.skip_probs)
        self.assertAlmostEqual(cov, 2 / 3)
        self.assertAlmostEqual(void, 1 / 3)

    def test_a_node_observed_only_before_the_first_consumed_event_reads_one(self):
        """`o` is the first event of every trace, so no interval bounds
        its move and every mdur is 0. It was observed - it has a
        synchronous move - so the ratio is 1, not an absence."""
        self.assertAlmostEqual(
            obsdur(self.o, self.tree, self.alignments[
                _variant_key(('o', 'x', 'y', 'z', 'p'))][0], self.complete), 0.0)
        self.assertAlmostEqual(self._ratio(self.o, self.complete), 1.0)


class UnobservedSubprocessTest(unittest.TestCase):
    """A subprocess with no synchronous move anywhere contributes to no
    O_sigma, so obscount is 0, the mass is 0 and void reads 1 - the
    reading the retired voidsat could not produce, since skip_prob times
    a zero mass is zero however absent the subprocess is."""

    def setUp(self):
        self.tree = build_partial_sequence_tree()
        by_id = nodes_of(self.tree)
        self.y = by_id['3']
        # y is recorded nowhere
        self.log = [trace(('o', 0), ('x', 1), ('z', 3), ('p', 4))]
        self.alignments = alignments_for(self.tree, self.log)

    def test_mass_is_zero_where_nothing_observes_the_subprocess(self):
        self.assertAlmostEqual(admass(self.y, self.tree, self.log, self.alignments), 0.0)

    def test_void_reads_one_for_a_subprocess_observed_nowhere(self):
        self.assertAlmostEqual(
            voidsat2(self.y, self.tree, self.log, self.alignments, {self.y: 0.0}), 1.0)

    def test_missing_duration_is_still_attributed_to_it(self):
        """It is absent, not unmeasured: the skip move takes its share
        of the gap it sits in."""
        path = self.alignments[_variant_key(('o', 'x', 'z', 'p'))][0]
        self.assertGreater(misdur(self.y, self.tree, path, self.log[0]), 0.0)


def _concurrent_tree():
    """seq( and(a, b), c ) - ids 1..5, root last."""
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
    return tree, parallel


class NormalFormPrecedencePinTest(unittest.TestCase):
    """
    voidsat2's values rest on WHERE the aligner puts a skip move, since
    mdur charges a move for the interval its own position falls in. A
    skip alignment is in normal form when transposable moves follow the
    precedence log, then synchronous, then skip, so for a skip that could
    sit on either side of a synchronous move - inside a parallel block,
    where both orderings are in the skip language - the precedence
    decides, and the aligner returns that one representative.

    Pinned behaviourally rather than by reading the rule: on
    seq(and(a, b), c) against <a@0h, c@10h> the skip of b must come after
    the synchronous move on a, which puts it in the 0-10h block and
    charges it half the gap. Transposing it to the front - still a valid
    skip alignment of the same trace - moves it before the first consumed
    event, where no interval bounds it. adratio at and(a, b) is then 0
    one way and 1 the other, so voidsat2 reads 1 or 0: the metric's whole
    range on a tie-break.

    If the precedence changes, this fails, and every voidsat2 value on a
    concurrent model must be re-read rather than silently changing
    meaning. lab.normal_form and lab.normal_form_compare are the
    research side of the same question.
    """

    @classmethod
    def setUpClass(cls):
        cls.tree, cls.parallel = _concurrent_tree()
        cls.trace = trace(('a', 0), ('c', 10))
        cls.delta = align(cls.tree, ('a', 'c'))[0]
        cls.delta_prime = [cls.delta[1], cls.delta[0], *cls.delta[2:]]

    def test_the_aligner_returns_one_alignment_with_the_skip_after_the_sync(self):
        paths = align(self.tree, ('a', 'c'))
        self.assertEqual(len(paths), 1)
        kinds = [('log' if model == '>>' else
                  'skip' if type(model).__name__ in ('Skip', 'TauPath') else 'sync')
                 for _log, model in self.delta]
        self.assertEqual(kinds, ['sync', 'skip', 'sync'])

    def test_the_skip_shares_the_whole_gap_with_the_following_event(self):
        self.assertEqual(consumes(self.delta), [0, 2])
        self.assertEqual(block(self.delta, 1, None), [1, 2])
        self.assertEqual(mdur(self.delta, self.trace, 1), 5 * 3600)

    def test_transposing_it_ahead_of_the_sync_leaves_it_unmeasurable(self):
        self.assertEqual(consumes(self.delta_prime), [1, 2])
        self.assertEqual(mdur(self.delta_prime, self.trace, 0), 0)

    def test_adratio_and_voidsat_span_their_range_on_the_precedence(self):
        skip_probs = {self.parallel: 0.0}
        key = _variant_key(('a', 'c'))
        for name, path, ratio in (('normal form', self.delta, 0.0),
                                  ('transposed', self.delta_prime, 1.0)):
            with self.subTest(alignment=name):
                self.assertAlmostEqual(
                    adratio(self.parallel, self.tree, path, self.trace), ratio, places=12)
                self.assertAlmostEqual(
                    voidsat2(self.parallel, self.tree, [self.trace], {key: [path]},
                             skip_probs), 1 - ratio, places=12)


class PerInstanceDurationTest(unittest.TestCase):
    """Duration is per-instance: payment_partial's sigma7 and sigma8 are
    the SAME variant with different timings, so they share an alignment
    but must not share a ratio. A mass built over deduplicated variants
    could not tell them apart."""

    def setUp(self):
        self.tree = build_payment_partial_tree()
        by_id = nodes_of(self.tree)
        self.approval = by_id['7']
        self.sigma7 = trace(('o', 0), ('a', 1), ('e', 5), ('s', 10), ('p', 11))
        self.sigma8 = trace(('o', 0), ('a', 1), ('e', 5), ('s', 6), ('p', 7))
        self.alignments = alignments_for(self.tree, [self.sigma7])

    def test_the_same_variant_with_different_timings_reads_differently(self):
        key = _variant_key(('o', 'a', 'e', 's', 'p'))
        path = self.alignments[key][0]
        r7 = adratio(self.approval, self.tree, path, self.sigma7)
        r8 = adratio(self.approval, self.tree, path, self.sigma8)
        self.assertNotAlmostEqual(r7, r8)


if __name__ == '__main__':
    unittest.main()
