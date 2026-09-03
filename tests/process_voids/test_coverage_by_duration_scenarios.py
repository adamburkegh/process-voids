
import datetime
import unittest

from skipalignments import Activity, Sequence, Xor, And, Loop

from process_voids.coveragemass import dur, sdur, coverage_by_duration

ACT_COST = 100000
BASE = datetime.datetime(2026, 1, 1)


def ev(activity, t):
    return {'concept:name': activity, 'time:timestamp': BASE + datetime.timedelta(seconds=t)}


def trace(*pairs):
    return [ev(a, t) for a, t in pairs]


def leaf(name, node_id):
    a = Activity(None, name, ACT_COST)
    a.id = node_id
    return a


def zero_skip_probs(*nodes):
    return {n: 0 for n in nodes}


class SequenceFirstEventArtefactTest(unittest.TestCase):
    """Model: seq(a,b). a is first in every trace, so it owns no interval."""

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.log = [trace(('a', 0), ('b', 1)), trace(('a', 0), ('b', 5))]
        self.total_dur = sum(dur(s) for s in self.log)
        self.skip_probs = zero_skip_probs(self.tree, self.a, self.b)

    def test_total_duration(self):
        self.assertEqual(self.total_dur, 6)

    def test_masses(self):
        cases = [(self.tree, 1.0), (self.a, 0.0), (self.b, 1.0)]
        for node, expected in cases:
            with self.subTest(node=node):
                mass = coverage_by_duration(node, self.log, self.skip_probs, self.total_dur)
                self.assertAlmostEqual(mass, expected, places=6)


class ChoicePartitionHoldsTest(unittest.TestCase):
    """Model: seq(o, xor(a,b)). Choice children partition the choice's mass."""

    def setUp(self):
        self.o = leaf('o', '1')
        self.a = leaf('a', '2')
        self.b = leaf('b', '3')
        self.choice = Xor(None, [self.a, self.b])
        self.choice.id = '4'
        self.a.set_parent(self.choice)
        self.b.set_parent(self.choice)
        self.tree = Sequence(None, [self.o, self.choice])
        self.tree.id = '5'
        self.o.set_parent(self.tree)
        self.choice.set_parent(self.tree)
        self.log = [trace(('o', 0), ('a', 2)), trace(('o', 0), ('b', 5))]
        self.total_dur = sum(dur(s) for s in self.log)
        self.skip_probs = zero_skip_probs(self.tree, self.choice, self.o, self.a, self.b)

    def test_total_duration(self):
        self.assertEqual(self.total_dur, 7)

    def test_masses(self):
        cases = [
            (self.tree, 1.0), (self.choice, 1.0),
            (self.a, 2 / 7), (self.b, 5 / 7), (self.o, 0.0),
        ]
        for node, expected in cases:
            with self.subTest(node=node):
                mass = coverage_by_duration(node, self.log, self.skip_probs, self.total_dur)
                self.assertAlmostEqual(mass, expected, places=6)

    def test_choice_children_partition_choice_mass(self):
        mass_a = coverage_by_duration(self.a, self.log, self.skip_probs, self.total_dur)
        mass_b = coverage_by_duration(self.b, self.log, self.skip_probs, self.total_dur)
        self.assertAlmostEqual(mass_a + mass_b, 1.0, places=6)


class ConcurrencyInterleavedPartitionHoldsTest(unittest.TestCase):
    """
    Model: seq(o, and(a,b)). Order of a/b differs between traces; intervals
    still partition (adjacent-gap attribution, not an overall-span one - a
    span implementation would give and(a,b) = 1 on sigma2, not 4).
    """

    def setUp(self):
        self.o = leaf('o', '1')
        self.a = leaf('a', '2')
        self.b = leaf('b', '3')
        self.conc = And(None, [self.a, self.b])
        self.conc.id = '4'
        self.a.set_parent(self.conc)
        self.b.set_parent(self.conc)
        self.tree = Sequence(None, [self.o, self.conc])
        self.tree.id = '5'
        self.o.set_parent(self.tree)
        self.conc.set_parent(self.tree)
        self.log = [trace(('o', 0), ('a', 1), ('b', 2)),
                    trace(('o', 0), ('b', 3), ('a', 4))]
        self.total_dur = sum(dur(s) for s in self.log)
        self.skip_probs = zero_skip_probs(self.tree, self.conc, self.o, self.a, self.b)

    def test_total_duration(self):
        self.assertEqual(self.total_dur, 6)

    def test_masses(self):
        cases = [(self.conc, 1.0), (self.a, 1 / 3), (self.b, 2 / 3)]
        for node, expected in cases:
            with self.subTest(node=node):
                mass = coverage_by_duration(node, self.log, self.skip_probs, self.total_dur)
                self.assertAlmostEqual(mass, expected, places=6)

    def test_concurrency_children_partition_and_mass(self):
        mass_a = coverage_by_duration(self.a, self.log, self.skip_probs, self.total_dur)
        mass_b = coverage_by_duration(self.b, self.log, self.skip_probs, self.total_dur)
        self.assertAlmostEqual(mass_a + mass_b, 1.0, places=6)


class LoopLeafMeanVsCompositeSumTest(unittest.TestCase):
    """
    Model: seq(o, loop(a,b)). Leaves use the mean rule and deliberately do
    NOT sum to the parent's (sum-rule) mass: 1/5 + 1/5 = 2/5 != 1.0.
    """

    def setUp(self):
        self.o = leaf('o', '1')
        self.a = leaf('a', '2')
        self.b = leaf('b', '3')
        self.loop = Loop(None, [self.a, self.b])
        self.loop.id = '4'
        self.a.set_parent(self.loop)
        self.b.set_parent(self.loop)
        self.tree = Sequence(None, [self.o, self.loop])
        self.tree.id = '5'
        self.o.set_parent(self.tree)
        self.loop.set_parent(self.tree)
        self.log = [trace(('o', 0), ('a', 1), ('b', 2), ('a', 3), ('b', 4), ('a', 5))]
        self.total_dur = sum(dur(s) for s in self.log)
        self.skip_probs = zero_skip_probs(self.tree, self.loop, self.o, self.a, self.b)

    def test_total_duration(self):
        self.assertEqual(self.total_dur, 5)

    def test_masses(self):
        cases = [(self.loop, 1.0), (self.a, 1 / 5), (self.b, 1 / 5)]
        for node, expected in cases:
            with self.subTest(node=node):
                mass = coverage_by_duration(node, self.log, self.skip_probs, self.total_dur)
                self.assertAlmostEqual(mass, expected, places=6)

    def test_leaves_do_not_sum_to_loop_mass(self):
        mass_a = coverage_by_duration(self.a, self.log, self.skip_probs, self.total_dur)
        mass_b = coverage_by_duration(self.b, self.log, self.skip_probs, self.total_dur)
        mass_loop = coverage_by_duration(self.loop, self.log, self.skip_probs, self.total_dur)
        self.assertNotAlmostEqual(mass_a + mass_b, mass_loop, places=2)
        self.assertAlmostEqual(mass_a + mass_b, 2 / 5, places=6)


class LoopUnevenIterationsTest(unittest.TestCase):
    """
    Model: seq(o, loop(a,b)), same tree as LoopLeafMeanVsCompositeSumTest but
    with an uneven-duration iteration. A max-based implementation would give
    a = 8; the mean rule
    gives a = 10/3.
    """

    def setUp(self):
        self.o = leaf('o', '1')
        self.a = leaf('a', '2')
        self.b = leaf('b', '3')
        self.loop = Loop(None, [self.a, self.b])
        self.loop.id = '4'
        self.a.set_parent(self.loop)
        self.b.set_parent(self.loop)
        self.tree = Sequence(None, [self.o, self.loop])
        self.tree.id = '5'
        self.o.set_parent(self.tree)
        self.loop.set_parent(self.tree)
        self.log = [trace(('o', 0), ('a', 1), ('b', 2), ('a', 10), ('b', 11), ('a', 12))]
        self.total_dur = sum(dur(s) for s in self.log)
        self.skip_probs = zero_skip_probs(self.tree, self.loop, self.o, self.a, self.b)

    def test_total_duration(self):
        self.assertEqual(self.total_dur, 12)

    def test_masses(self):
        cases = [(self.loop, 1.0), (self.a, 5 / 18), (self.b, 1 / 12)]
        for node, expected in cases:
            with self.subTest(node=node):
                mass = coverage_by_duration(node, self.log, self.skip_probs, self.total_dur)
                self.assertAlmostEqual(mass, expected, places=6)

    def test_a_mean_is_not_the_max_interval(self):
        activities = set(self.a.get_leaf_labels())
        raw_mean = sdur(activities, self.log[0], has_submodels=False)
        self.assertAlmostEqual(raw_mean, 10 / 3, places=6)
        self.assertNotAlmostEqual(raw_mean, 8, places=2)


class NestedCompositeUnderInterleavingTest(unittest.TestCase):
    """
    Model: seq(o, and(seq(a,b), c)). Composite children of a composite
    parent partition its mass even though execution is interleaved:
    seq(a,b) = 2/3, c = 1/3, summing to and(...)'s own mass of 1.0.
    """

    def setUp(self):
        self.o = leaf('o', '1')
        self.a = leaf('a', '2')
        self.b = leaf('b', '3')
        self.c = leaf('c', '4')
        self.seqab = Sequence(None, [self.a, self.b])
        self.seqab.id = '5'
        self.a.set_parent(self.seqab)
        self.b.set_parent(self.seqab)
        self.conc = And(None, [self.seqab, self.c])
        self.conc.id = '6'
        self.seqab.set_parent(self.conc)
        self.c.set_parent(self.conc)
        self.tree = Sequence(None, [self.o, self.conc])
        self.tree.id = '7'
        self.o.set_parent(self.tree)
        self.conc.set_parent(self.tree)
        self.log = [trace(('o', 0), ('a', 1), ('c', 2), ('b', 3))]
        self.total_dur = sum(dur(s) for s in self.log)
        self.skip_probs = zero_skip_probs(self.tree, self.conc, self.seqab, self.o,
                                           self.a, self.b, self.c)

    def test_total_duration(self):
        self.assertEqual(self.total_dur, 3)

    def test_masses(self):
        cases = [
            (self.conc, 1.0), (self.seqab, 2 / 3),
            (self.c, 1 / 3), (self.a, 1 / 3), (self.b, 1 / 3),
        ]
        for node, expected in cases:
            with self.subTest(node=node):
                mass = coverage_by_duration(node, self.log, self.skip_probs, self.total_dur)
                self.assertAlmostEqual(mass, expected, places=6)

    def test_composite_children_partition_and_mass(self):
        mass_seqab = coverage_by_duration(self.seqab, self.log, self.skip_probs, self.total_dur)
        mass_c = coverage_by_duration(self.c, self.log, self.skip_probs, self.total_dur)
        self.assertAlmostEqual(mass_seqab + mass_c, 1.0, places=6)


class RootMassBelowOneTest(unittest.TestCase):
    """
    Model: seq(a,b). Root mass is only 1.0 when every event's activity is
    in the model alphabet - an out-of-vocabulary event (z) leaves a real gap.
    """

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.log = [trace(('a', 0), ('b', 1), ('z', 5))]
        self.total_dur = sum(dur(s) for s in self.log)
        self.skip_probs = zero_skip_probs(self.tree, self.a, self.b)

    def test_total_duration(self):
        self.assertEqual(self.total_dur, 5)

    def test_masses(self):
        cases = [(self.tree, 0.2), (self.b, 0.2), (self.a, 0.0)]
        for node, expected in cases:
            with self.subTest(node=node):
                mass = coverage_by_duration(node, self.log, self.skip_probs, self.total_dur)
                self.assertAlmostEqual(mass, expected, places=6)


class DegenerateTracesTest(unittest.TestCase):

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.skip_probs = zero_skip_probs(self.tree, self.a, self.b)

    def test_single_event_trace_duration_and_sdur_are_zero(self):
        single = trace(('a', 0))
        self.assertEqual(dur(single), 0)
        self.assertEqual(sdur({'a', 'b'}, single, has_submodels=True), 0)
        self.assertEqual(sdur({'a'}, single, has_submodels=False), 0)

    def test_empty_trace_duration_and_sdur_are_zero(self):
        self.assertEqual(dur([]), 0)
        self.assertEqual(sdur({'a', 'b'}, [], has_submodels=True), 0)

    def test_log_of_only_degenerate_traces_raises(self):
        log = [trace(('a', 0)), []]
        with self.assertRaises(ValueError):
            coverage_by_duration(self.tree, log, self.skip_probs)

    def test_degenerate_trace_mixed_into_normal_log_contributes_zero(self):
        normal = trace(('a', 0), ('b', 4))
        degenerate = trace(('a', 10))
        log_with = [normal, degenerate]
        log_without = [normal]

        total_with = sum(dur(s) for s in log_with)
        total_without = sum(dur(s) for s in log_without)
        self.assertEqual(total_with, total_without)

        mass_with = coverage_by_duration(self.tree, log_with, self.skip_probs, total_with)
        mass_without = coverage_by_duration(self.tree, log_without, self.skip_probs, total_without)
        self.assertAlmostEqual(mass_with, mass_without, places=6)


class SharedAlphabetPartitionBreaksTest(unittest.TestCase):
    """
    Model: seq(xor(a,b), xor(a,c)) - 'a' is a leaf in both choices. Known,
    accepted limitation: both choices claim the same interval, so their
    masses sum to 2.0 across the cut rather than partitioning. Each mass
    individually still stays <= 1 (Lemma 1 holds per subprocess) - this is
    expected behaviour, not a bug.
    """

    def setUp(self):
        self.a1 = leaf('a', '1')
        self.b = leaf('b', '2')
        self.xor1 = Xor(None, [self.a1, self.b])
        self.xor1.id = '3'
        self.a1.set_parent(self.xor1)
        self.b.set_parent(self.xor1)

        self.a2 = leaf('a', '4')
        self.c = leaf('c', '5')
        self.xor2 = Xor(None, [self.a2, self.c])
        self.xor2.id = '6'
        self.a2.set_parent(self.xor2)
        self.c.set_parent(self.xor2)

        self.tree = Sequence(None, [self.xor1, self.xor2])
        self.tree.id = '7'
        self.xor1.set_parent(self.tree)
        self.xor2.set_parent(self.tree)

        self.log = [trace(('a', 0), ('a', 5))]
        self.total_dur = sum(dur(s) for s in self.log)
        self.skip_probs = zero_skip_probs(self.tree, self.xor1, self.xor2,
                                           self.a1, self.b, self.a2, self.c)

    def test_total_duration(self):
        self.assertEqual(self.total_dur, 5)

    def test_masses(self):
        cases = [(self.xor1, 1.0), (self.xor2, 1.0)]
        for node, expected in cases:
            with self.subTest(node=node):
                mass = coverage_by_duration(node, self.log, self.skip_probs, self.total_dur)
                self.assertAlmostEqual(mass, expected, places=6)

    def test_shared_alphabet_children_exceed_partition(self):
        mass_1 = coverage_by_duration(self.xor1, self.log, self.skip_probs, self.total_dur)
        mass_2 = coverage_by_duration(self.xor2, self.log, self.skip_probs, self.total_dur)
        self.assertLessEqual(mass_1, 1.0)
        self.assertLessEqual(mass_2, 1.0)
        self.assertAlmostEqual(mass_1 + mass_2, 2.0, places=6)


if __name__ == '__main__':
    unittest.main()
