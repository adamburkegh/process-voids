import unittest
from datetime import datetime, timedelta
from unittest.mock import patch

from skipalignments import Activity, Sequence, Xor, Tau, Loop, Aligner

from process_voids.coveragemass import (
    _variant_key, consumes, nxt, block, mdur, adur, admass, covat, voidat, voidsat,
    make_aligned_duration_cache, log_to_traces,
)

ACT_COST = 100000

Aligner.set_level_incentive(0)


def leaf(name, node_id, cost=ACT_COST):
    node = Activity(None, name, cost)
    node.id = node_id
    return node


def align(tree, trace):
    states, _ = Aligner(tree).align_normal_form(list(trace), [ACT_COST] * len(trace), True, timeout=100)
    return states


def alignments_by_variant(tree, traces_with_weights):
    """{variant_key: [path, ...]} for each distinct trace variant - the
    alignments_by_variant shape admass/adur expect, built from the real
    aligner (not hand-picked paths)."""
    result = {}
    for trace in traces_with_weights:
        states = align(tree, trace)
        result[_variant_key(trace)] = [s.path for s in states]
    return result


def event(name, t):
    return {'concept:name': name, 'time:timestamp': t}


T0 = datetime(2026, 1, 1, 12, 0, 0)


def trace(*name_offsets):
    """trace(('a', 0), ('b', 10)) -> events at T0, T0+10s."""
    return [event(name, T0 + timedelta(seconds=offset)) for name, offset in name_offsets]


class ConsumesNxtBlockTest(unittest.TestCase):
    """Direct unit tests of the Definition [Move Durations] primitives
    against small hand-built paths - no real aligner needed, these are
    pure functions over (log_elem, model_elem) lists."""

    def setUp(self):
        self.a = leaf('a', '1')
        self.tau = Tau(None, 'tau', 0)

    def test_consumes_is_log_and_sync_positions_only(self):
        path = [('a', self.a), ('>>', self.a), ('b', '>>')]
        self.assertEqual(consumes(path), [0, 2])

    def test_nxt_finds_smallest_consuming_position_at_or_after_j(self):
        path = [('a', self.a), ('>>', self.a), ('b', '>>')]
        self.assertEqual(nxt(path, 0), 0)
        self.assertEqual(nxt(path, 1), 2)
        self.assertEqual(nxt(path, 2), 2)

    def test_nxt_is_none_past_the_last_consuming_position(self):
        path = [('a', self.a), ('>>', self.a)]
        self.assertIsNone(nxt(path, 1))

    def test_block_excludes_silent_moves_but_includes_everything_else(self):
        # positions 1 (skip) and 2 (sync) both resolve to nxt=2; the
        # silent (tau) move at position 0 is excluded even though its
        # own nxt is also 2.
        from skipalignments import TauPath, Skip
        path = [('>>', TauPath(self.tau)), ('>>', Skip(self.a, 1)), ('a', self.a)]
        self.assertEqual(block(path, 1), [1, 2])
        self.assertEqual(block(path, 2), [1, 2])

    def test_block_is_empty_when_nxt_is_undefined(self):
        path = [('a', self.a), ('>>', self.a)]
        self.assertEqual(block(path, 1), [])


class SimpleSequenceTest(unittest.TestCase):
    """M: model seq(a,b), one trace <a@T0, b@T0+10s>, fully synchronous.
    No lumping, no sharing - a clean baseline: the first activity has
    nothing preceding it (zero), the second takes the entire gap alone
    (block of size 1, since nothing else shares its target)."""

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)

        self.trace = trace(('a', 0), ('b', 10))
        self.log = [self.trace]
        self.alignments = alignments_by_variant(self.tree, {('a', 'b'): 1.0})

    def test_first_activity_has_no_attributable_gap(self):
        self.assertAlmostEqual(admass(self.a, self.tree, self.log, self.alignments), 0.0)

    def test_second_activity_takes_the_whole_gap_alone(self):
        self.assertAlmostEqual(admass(self.b, self.tree, self.log, self.alignments), 1.0)

    def test_root_also_gets_the_whole_gap(self):
        self.assertAlmostEqual(admass(self.tree, self.tree, self.log, self.alignments), 1.0)


class LumpedMissingSubprocessSharesBlockWithFollowingEventTest(unittest.TestCase):
    """M: model seq(a, seq(x,y), c), trace <a@T0, c@T0+10s> - seq(x,y)
    is entirely missing and skip-alignments lifts it to ONE skip move.

    Confirmed against the real aligner's output (not assumed): that skip
    move and the immediately-following sync move c share the SAME block
    (both resolve to the same nxt, c's own position - c is non-silent,
    so block()'s tau-only exclusion doesn't drop it). The gap therefore
    splits 50/50 between the missing subprocess and c, not "the whole
    gap" to the missing subprocess alone - lumping only reduces how many
    pieces the missing portion itself is chopped into (1 piece here vs.
    one per leaf - x and y separately - on a classical, unlumped
    alignment, where the block would be 3-way instead of 2-way). A
    trailing real event always keeps its own genuine share of the
    block it lands in, lumped or not.
    """

    def setUp(self):
        self.a = leaf('a', '1')
        self.x = leaf('x', '2')
        self.y = leaf('y', '3')
        self.missing = Sequence(None, [self.x, self.y])
        self.missing.id = '4'
        self.x.set_parent(self.missing)
        self.y.set_parent(self.missing)
        self.c = leaf('c', '5')
        self.tree = Sequence(None, [self.a, self.missing, self.c])
        self.tree.id = '6'
        self.a.set_parent(self.tree)
        self.missing.set_parent(self.tree)
        self.c.set_parent(self.tree)

        self.trace = trace(('a', 0), ('c', 10))
        self.log = [self.trace]
        self.alignments = alignments_by_variant(self.tree, {('a', 'c'): 1.0})

    def test_missing_subprocess_gets_half_the_gap(self):
        self.assertAlmostEqual(admass(self.missing, self.tree, self.log, self.alignments), 0.5)

    def test_c_gets_the_other_half(self):
        self.assertAlmostEqual(admass(self.c, self.tree, self.log, self.alignments), 0.5)

    def test_a_has_no_attributable_gap(self):
        self.assertAlmostEqual(admass(self.a, self.tree, self.log, self.alignments), 0.0)

    def test_root_gets_the_whole_gap_since_both_halves_are_under_it(self):
        self.assertAlmostEqual(admass(self.tree, self.tree, self.log, self.alignments), 1.0)

    def test_shared_cache_gives_identical_results_to_no_cache(self):
        # A shared make_aligned_duration_cache is what lab.exp_disco_
        # degrade's per-node loop passes across every node in a report
        # row (see make_executions_cache's equivalent for alignment_mass
        # /coverage_by_alignment_pn) - must not change any node's number.
        cache = make_aligned_duration_cache(self.tree)
        for node in (self.a, self.missing, self.c, self.tree):
            with self.subTest(node=node):
                uncached = admass(node, self.tree, self.log, self.alignments)
                cached = admass(node, self.tree, self.log, self.alignments, cache=cache)
                self.assertAlmostEqual(uncached, cached, places=9)


class MultipleTracesAveragedTest(unittest.TestCase):
    """M: model seq(a,b). Two traces of the same variant <a,b> with
    different real gaps (10s and 30s) - admass averages adur/duration
    PER TRACE INSTANCE, not per deduplicated variant, so this must NOT
    collapse to a single shared number."""

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)

        self.log = [trace(('a', 0), ('b', 10)), trace(('a', 0), ('b', 30))]
        self.alignments = alignments_by_variant(self.tree, {('a', 'b'): 1.0})

    def test_b_still_gets_full_gap_on_every_trace_so_admass_is_one(self):
        # both traces give b a ratio of 1.0 (its own gap / its own
        # trace duration, since it alone occupies the whole trace) -
        # averaging 1.0 and 1.0 is still 1.0, so this alone wouldn't
        # distinguish per-trace from per-variant averaging.
        self.assertAlmostEqual(admass(self.b, self.tree, self.log, self.alignments), 1.0)

    def test_root_is_also_one_regardless_of_the_differing_absolute_gaps(self):
        self.assertAlmostEqual(admass(self.tree, self.tree, self.log, self.alignments), 1.0)


class ZeroDurationSubmodelExcludedFromAverageTest(unittest.TestCase):
    """M: model seq(a,b,c). One trace where a is immediately followed
    by b (no gap ever attributed to a - see SimpleSequenceTest), so a's
    adur is 0 on every trace - admass(a) must be 0.0 (not NaN/crash)
    since L' (traces with observable duration for a) is empty."""

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.c = leaf('c', '3')
        self.tree = Sequence(None, [self.a, self.b, self.c])
        self.tree.id = '4'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.c.set_parent(self.tree)

        self.log = [trace(('a', 0), ('b', 5), ('c', 15))]
        self.alignments = alignments_by_variant(self.tree, {('a', 'b', 'c'): 1.0})

    def test_a_is_zero_not_a_crash(self):
        self.assertEqual(admass(self.a, self.tree, self.log, self.alignments), 0.0)


class EmptyAlignmentsGiveZeroTest(unittest.TestCase):
    """A variant with no alignments recorded (e.g. it errored/timed out
    upstream) contributes nothing rather than crashing."""

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.log = [trace(('a', 0), ('b', 10))]

    def test_missing_variant_gives_zero_admass(self):
        self.assertEqual(admass(self.b, self.tree, self.log, {}), 0.0)

    def test_adur_of_empty_alignment_list_is_zero(self):
        self.assertEqual(adur(self.b, self.tree, [], self.log[0]), 0.0)


class CovatVoidatTest(unittest.TestCase):
    """covat/voidat are just (1-skip_prob)*admass / skip_prob*admass -
    thin wrappers, tested directly against a fixed skip_prob rather than
    a real DerivationPipeline (that's voidsat's job, tested separately)."""

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.log = [trace(('a', 0), ('b', 10))]
        self.alignments = alignments_by_variant(self.tree, {('a', 'b'): 1.0})

    def test_covat_and_voidat_are_complementary_at_full_mass(self):
        self.assertAlmostEqual(
            covat(self.b, self.tree, 0.3, self.log, self.alignments), 0.7, places=6)
        self.assertAlmostEqual(
            voidat(self.b, self.tree, 0.3, self.log, self.alignments), 0.3, places=6)

    def test_zero_admass_gives_zero_regardless_of_skip_prob(self):
        self.assertAlmostEqual(
            covat(self.a, self.tree, 0.3, self.log, self.alignments), 0.0, places=6)
        self.assertAlmostEqual(
            voidat(self.a, self.tree, 0.3, self.log, self.alignments), 0.0, places=6)


class VoidsatUsesRealSkipProbsTest(unittest.TestCase):
    """voidsat wires dv.skip_dict_backup's own State.path objects
    (already in this module's wrapper shape) and dv.skip_probs straight
    through to voidat - a thin fake standing in for a real
    DerivationPipeline's two relevant attributes is enough to test the
    wiring itself; voidat's own numerics are covered above."""

    class _FakeDv:
        def __init__(self, skip_dict_backup, skip_probs):
            self.skip_dict_backup = skip_dict_backup
            self.skip_probs = skip_probs

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.log = [trace(('a', 0), ('b', 10))]

        states = align(self.tree, ('a', 'b'))
        skip_dict_backup = {_variant_key(('a', 'b')): states}
        skip_probs = {self.tree: 0.0, self.a: 0.0, self.b: 0.4}
        self.dv = self._FakeDv(skip_dict_backup, skip_probs)

    def test_voidsat_multiplies_admass_by_the_real_skip_prob(self):
        # b takes the whole gap (admass=1.0) * skip_prob(b)=0.4
        self.assertAlmostEqual(voidsat(self.b, self.tree, self.dv, self.log), 0.4, places=6)

    def test_voidsat_is_zero_when_skip_prob_is_zero(self):
        self.assertAlmostEqual(voidsat(self.a, self.tree, self.dv, self.log), 0.0, places=6)


class CachedTracesTest(unittest.TestCase):
    """
    make_aligned_duration_cache(tree, log) builds log_to_traces(log) - a
    group-by, sort and conversion over the whole log - once, and admass
    reuses it for every node. Only for that same log: a trace list cached
    for a different log would silently give the wrong answer.

    M: model seq(a,b,c). Two logs with the same variant but different
    timing: b's share of the trace is 10/30 in one and 20/30 in the other.
    """

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.c = leaf('c', '3')
        self.tree = Sequence(None, [self.a, self.b, self.c])
        self.tree.id = '4'
        for node in (self.a, self.b, self.c):
            node.set_parent(self.tree)
        self.log_early = [trace(('a', 0), ('b', 10), ('c', 30))]
        self.log_late = [trace(('a', 0), ('b', 20), ('c', 30))]
        self.alignments = alignments_by_variant(self.tree, {('a', 'b', 'c'): 1.0})

    def test_trace_list_built_once_across_every_node(self):
        with patch('process_voids.coveragemass.log_to_traces', wraps=log_to_traces) as counted:
            cache = make_aligned_duration_cache(self.tree, self.log_early)
            for node in (self.tree, self.a, self.b, self.c):
                admass(node, self.tree, self.log_early, self.alignments, cache)
        self.assertEqual(counted.call_count, 1)

    def test_cache_built_for_another_log_is_not_used(self):
        cache = make_aligned_duration_cache(self.tree, self.log_early)
        self.assertAlmostEqual(
            admass(self.b, self.tree, self.log_late, self.alignments, cache), 2 / 3, places=6)
        self.assertAlmostEqual(
            admass(self.b, self.tree, self.log_early, self.alignments, cache), 1 / 3, places=6)


if __name__ == '__main__':
    unittest.main()
