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
"""

import unittest
from datetime import datetime, timedelta

from skipalignments import Aligner

from lab.fixtures import build_partial_sequence_tree, build_payment_partial_tree
from process_voids.coveragemass import _variant_key, log_to_traces, obsdur, misdur
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
