'''
Unit tests for process_voids.surprise (interval surprise metric).

TestRunningExampleWorkedNumbers pins hand-worked figures for the
running example (headline 10.93 bits/17 events, containment
sched=4.96, predecessor approval=9.60) - a regression baseline against
independently-computed reference numbers.

TestComputePredecessorsLoop documents the loop predecessor rule: both
the do-child and the redo-child allocate to the body (do-child) node
unconditionally - the body precedes both the redo and every subsequent
iteration, so there is no first-iteration special case to track.
'''

import datetime
import unittest

from skipalignments import Activity, Tau, Sequence, Xor, And, Loop

from lab.fixtures import build_running_example_log, build_running_example_tree
from process_voids.coveragemass import log_to_traces
from process_voids.surprise import (
    event_intervals, observed_intervals, tail_probability, event_surprise,
    surprise_totals, compute_predecessors, predecessor_totals,
)

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


def tau(name, node_id):
    t = Tau(None, name, 0)
    t.id = node_id
    return t


def set_parent(nodes, parent):
    for node in nodes:
        node.set_parent(parent)


class EventIntervalsTest(unittest.TestCase):
    def test_first_event_owns_no_interval(self):
        tr = trace(('a', 0), ('b', 5), ('c', 12))
        self.assertEqual(event_intervals(tr),
                          [(1, 'b', 5.0), (2, 'c', 7.0)])

    def test_single_event_trace_has_no_intervals(self):
        self.assertEqual(event_intervals(trace(('a', 0))), [])


class ObservedIntervalsTest(unittest.TestCase):
    def test_aggregates_across_traces_per_activity(self):
        traces = [trace(('a', 0), ('b', 5)), trace(('a', 0), ('b', 9), ('b', 11))]
        obs = observed_intervals(traces)
        # sorted ascending, not insertion order - tail_probability needs this
        self.assertEqual(obs['b'], [2.0, 5.0, 9.0])
        self.assertNotIn('a', obs)


class TailProbabilityTest(unittest.TestCase):
    def test_plus_one_denominator_keeps_max_finite(self):
        vals = [1.0, 2.0, 3.0]
        self.assertEqual(tail_probability(vals, 3.0), 1 / 4)

    def test_monotone_non_increasing_in_d(self):
        vals = [1.0, 5.0, 9.0, 20.0]
        ps = [tail_probability(vals, d) for d in (0.0, 5.0, 9.0, 20.0, 100.0)]
        self.assertEqual(ps, sorted(ps, reverse=True))

    def test_empty_observations_floors_at_one(self):
        # never hit via event_surprise(traces) in practice - any activity
        # iterated there necessarily contributed to its own obs entry - but
        # pins tail_probability's own behavior as a standalone function.
        # count=0 (no observations at all) floors to 1, same as any other
        # value more extreme than the reference distribution has ever seen.
        self.assertEqual(tail_probability([], 5.0), 1.0)

    def test_value_beyond_reference_max_floors_at_min_probability(self):
        # a value more extreme than anything ever observed gets the same
        # probability the in-sample maximum would - not 0.0.
        vals = [1.0, 2.0, 3.0]
        self.assertEqual(tail_probability(vals, 100.0), tail_probability(vals, 3.0))
        self.assertEqual(tail_probability(vals, 100.0), 1 / 4)


class EventSurpriseTest(unittest.TestCase):
    def test_longer_interval_is_more_surprising(self):
        traces = [trace(('a', 0), ('b', 3600)), trace(('a', 0), ('b', 36000))]
        rows = event_surprise(traces)
        bits_by_delta = {delta: bits for _, delta, _, bits in rows}
        self.assertGreater(bits_by_delta[36000.0], bits_by_delta[3600.0])

    def test_reuses_supplied_distribution(self):
        traces = [trace(('a', 0), ('b', 2))]
        obs = {'b': [1.0, 2.0, 3.0]}
        rows = event_surprise(traces, obs=obs)
        self.assertEqual(len(rows), 1)
        _, delta, p, bits = rows[0]
        self.assertEqual(p, tail_probability(obs['b'], delta))


class SurpriseTotalsContainmentTest(unittest.TestCase):
    '''Model: seq(a, xor(b, tau)). Additive: root = sum of leaves.'''

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.t = tau('skip-b', '3')
        self.choice = Xor(None, [self.b, self.t])
        self.choice.id = '4'
        set_parent([self.b, self.t], self.choice)
        self.tree = Sequence(None, [self.a, self.choice])
        self.tree.id = '5'
        set_parent([self.a, self.choice], self.tree)

    def test_additive_root_equals_sum_of_leaves(self):
        rows = [('a', 1.0, 0.5, 1.0), ('b', 1.0, 0.5, 1.0), ('b', 1.0, 0.25, 2.0)]
        totals, _ = surprise_totals(self.tree, rows)
        self.assertEqual(totals[self.a], 1.0)
        self.assertEqual(totals[self.b], 3.0)
        self.assertEqual(totals[self.t], 0.0)
        self.assertEqual(totals[self.choice], 3.0)
        self.assertEqual(totals[self.tree], 4.0)

    def test_no_rows_gives_all_zero(self):
        totals, _ = surprise_totals(self.tree, [])
        self.assertTrue(all(v == 0.0 for v in totals.values()))


class SurpriseTotalsDuplicateLeavesTest(unittest.TestCase):
    '''
    Model: seq(before, seq(dup1, mid, dup2), after) - two leaves both
    named 'dup', nested so their LCA (the inner Sequence) is a proper
    internal node, not the root. Their combined surprise should land on
    the LCA and propagate up, with neither the specific leaves nor the
    LCA's siblings getting any of it.
    '''

    def setUp(self):
        self.before = leaf('before', '1')
        self.dup1 = leaf('dup', '2')
        self.mid = leaf('mid', '3')
        self.dup2 = leaf('dup', '4')
        self.inner = Sequence(None, [self.dup1, self.mid, self.dup2])
        self.inner.id = '5'
        set_parent([self.dup1, self.mid, self.dup2], self.inner)
        self.after = leaf('after', '6')
        self.tree = Sequence(None, [self.before, self.inner, self.after])
        self.tree.id = '7'
        set_parent([self.before, self.inner, self.after], self.tree)

    def test_duplicate_bits_land_on_lca_not_leaves(self):
        rows = [('before', 1.0, 0.5, 1.0), ('dup', 1.0, 0.5, 1.0),
                ('mid', 1.0, 0.5, 1.0), ('dup', 1.0, 0.25, 2.0),
                ('after', 1.0, 0.5, 1.0)]
        totals, _ = surprise_totals(self.tree, rows)

        self.assertEqual(totals[self.dup1], 0.0)
        self.assertEqual(totals[self.dup2], 0.0)
        self.assertEqual(totals[self.mid], 1.0)
        # inner = dup's full total (3.0) + mid's own (1.0)
        self.assertEqual(totals[self.inner], 4.0)
        self.assertEqual(totals[self.before], 1.0)
        self.assertEqual(totals[self.after], 1.0)
        # root = everything, still conserved
        self.assertEqual(totals[self.tree], 6.0)


class SurpriseTotalsOutOfAlphabetTest(unittest.TestCase):
    '''
    Model: seq(a, b) - 'c' never appears as a leaf anywhere (eg pruned
    out of a discovered model by noise filtering). Its bits must still
    reach the root rather than vanishing, since the root is supposed to
    hold the log's full surprise total regardless of the model.
    '''

    def setUp(self):
        self.a = leaf('a', '1')
        self.b = leaf('b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        set_parent([self.a, self.b], self.tree)

    def test_out_of_alphabet_bits_land_on_root_not_dropped(self):
        rows = [('a', 1.0, 0.5, 1.0), ('c', 1.0, 0.5, 1.0), ('c', 1.0, 0.25, 2.0)]
        totals, out_of_alphabet_count = surprise_totals(self.tree, rows)

        self.assertEqual(out_of_alphabet_count, 2)
        self.assertEqual(totals[self.a], 1.0)
        # 'c' isn't a leaf anywhere, so it can't reach any node except the root
        self.assertEqual(totals[self.b], 0.0)
        self.assertEqual(totals[self.tree], 1.0 + 1.0 + 2.0)

    def test_no_out_of_alphabet_activity_gives_zero_count(self):
        rows = [('a', 1.0, 0.5, 1.0), ('b', 1.0, 0.5, 1.0)]
        _, out_of_alphabet_count = surprise_totals(self.tree, rows)
        self.assertEqual(out_of_alphabet_count, 0)


class ComputePredecessorsStructuralTest(unittest.TestCase):
    '''Sequence/Xor/And rules, isolated from the loop first-iteration nuance.'''

    def test_sequence_preceding_sibling(self):
        a, b, c = leaf('a', '1'), leaf('b', '2'), leaf('c', '3')
        tree = Sequence(None, [a, b, c])
        tree.id = '4'
        set_parent([a, b, c], tree)

        preds = compute_predecessors(tree)
        self.assertEqual(preds[a], (None, False))
        self.assertEqual(preds[b], (a, False))
        self.assertEqual(preds[c], (b, False))

    def test_choice_ambiguous_on_exit_inherited_on_entry(self):
        a, b, c = leaf('a', '1'), leaf('b', '2'), leaf('c', '3')
        choice = Xor(None, [b, c])
        choice.id = '4'
        set_parent([b, c], choice)
        tree = Sequence(None, [a, choice])
        tree.id = '5'
        set_parent([a, choice], tree)

        preds = compute_predecessors(tree)
        # entering the choice: every branch inherits the choice's own predecessor
        self.assertEqual(preds[b], (a, False))
        self.assertEqual(preds[c], (a, False))

    def test_choice_charges_exit_to_the_choice_node_itself(self):
        a, b, c, d = leaf('a', '1'), leaf('b', '2'), leaf('c', '3'), leaf('d', '4')
        choice = Xor(None, [b, c])
        choice.id = '5'
        set_parent([b, c], choice)
        tree = Sequence(None, [a, choice, d])
        tree.id = '6'
        set_parent([a, choice, d], tree)

        preds = compute_predecessors(tree)
        self.assertEqual(preds[d], (choice, True))

    def test_parallel_same_treatment_as_choice(self):
        a, b, c, d = leaf('a', '1'), leaf('b', '2'), leaf('c', '3'), leaf('d', '4')
        par = And(None, [b, c])
        par.id = '5'
        set_parent([b, c], par)
        tree = Sequence(None, [a, par, d])
        tree.id = '6'
        set_parent([a, par, d], tree)

        preds = compute_predecessors(tree)
        self.assertEqual(preds[b], (a, False))
        self.assertEqual(preds[c], (a, False))
        self.assertEqual(preds[d], (par, True))


class ComputePredecessorsLoopTest(unittest.TestCase):
    '''
    Current rule: do's predecessor is always the redo-child's exit, redo's
    predecessor is always the do-child's exit. Known-imprecise on a
    trace's first iteration (see module docstring) - under revision.
    '''

    def setUp(self):
        self.before = leaf('o', '1')
        self.do = leaf('a', '2')
        self.redo = leaf('e', '3')
        self.after = leaf('p', '4')
        self.loop = Loop(None, [self.do, self.redo])
        self.loop.id = '5'
        set_parent([self.do, self.redo], self.loop)
        self.tree = Sequence(None, [self.before, self.loop, self.after])
        self.tree.id = '6'
        set_parent([self.before, self.loop, self.after], self.tree)

    def test_do_and_redo_both_allocate_to_the_body(self):
        preds = compute_predecessors(self.tree)
        self.assertEqual(preds[self.do], (self.do, False))
        self.assertEqual(preds[self.redo], (self.do, False))

    def test_follower_predecessor_is_the_do_child(self):
        preds = compute_predecessors(self.tree)
        self.assertEqual(preds[self.after], (self.do, False))


class RunningExampleWorkedNumbersTest(unittest.TestCase):
    '''Pins the exact figures from void-entropy-brief-v2's worked example.'''

    def setUp(self):
        self.tree = build_running_example_tree()
        self.traces = log_to_traces(build_running_example_log())
        self.rows = event_surprise(self.traces)

    def _node(self, node_id):
        def _find(n):
            if n.id == node_id:
                return n
            for c in getattr(n, 'children', []):
                found = _find(c)
                if found is not None:
                    return found
            return None
        return _find(self.tree)

    def test_headline(self):
        total = sum(bits for *_, bits in self.rows)
        self.assertEqual(len(self.rows), 17)
        self.assertAlmostEqual(total, 10.93, places=2)

    def test_containment_puts_signal_on_schedule_choice(self):
        totals, _ = surprise_totals(self.tree, self.rows)
        schedule_choice = self._node('8')  # Xor(s, tau)
        self.assertAlmostEqual(totals[schedule_choice], 4.96, places=2)
        self.assertAlmostEqual(totals[self.tree], 10.93, places=2)

    def test_full_predecessor_map(self):
        o, a, e, s = self._node('1'), self._node('2'), self._node('3'), self._node('4')
        skip_schedule, p = self._node('5'), self._node('6')
        schedule_choice = self._node('8')
        predecessors = compute_predecessors(self.tree)
        self.assertEqual(predecessors[o], (None, False))
        self.assertEqual(predecessors[a], (a, False))
        self.assertEqual(predecessors[e], (a, False))
        self.assertEqual(predecessors[s], (a, False))
        self.assertEqual(predecessors[skip_schedule], (a, False))
        self.assertEqual(predecessors[p], (schedule_choice, True))

    def test_predecessor_moves_signal_onto_approval(self):
        predecessors = compute_predecessors(self.tree)
        totals, ambiguous_event_count, unattributable_event_count = \
            predecessor_totals(self.tree, self.rows, predecessors)
        approval = self._node('7')  # Loop(a, e)
        do = self._node('2')        # 'a', the loop's body
        self.assertAlmostEqual(totals[approval], 9.60, places=2)
        self.assertAlmostEqual(totals[self.tree], 10.93, places=2)
        # both do's and redo's own surprise allocate to the body ('a')
        # unconditionally, so the body carries the loop's whole total
        self.assertAlmostEqual(totals[do], totals[approval], places=9)
        # every trace's terminal 'p' resolves through the ambiguous
        # Xor(s, tau) exit - one per trace, six traces
        self.assertEqual(ambiguous_event_count, 6)
        # 'o' has no predecessor, but it's always the first event of
        # every trace (i == 0), so it never owns an interval/row here
        self.assertEqual(unattributable_event_count, 0)

    def test_root_total_agrees_across_attribution_schemes(self):
        containment, _ = surprise_totals(self.tree, self.rows)
        predecessors = compute_predecessors(self.tree)
        pred_totals, _, _ = predecessor_totals(self.tree, self.rows, predecessors)
        self.assertAlmostEqual(containment[self.tree], pred_totals[self.tree], places=9)


class PredecessorTotalsUnattributableTest(unittest.TestCase):
    '''
    An event whose activity has no predecessor at all (only possible for
    a leaf with no preceding sibling anywhere on its path from the root)
    must not be silently dropped - its bits are charged directly to the
    root. Exercised directly via hand-built rows rather than a real
    trace, since a leaf with pred=None structurally executes at most
    once per trace in a well-formed model (repetition needs a Loop, and
    a Loop's own children never get pred=None - see compute_predecessors)
    - RunningExampleWorkedNumbersTest's 'o' never owns an interval for
    exactly this reason, which is why that test alone couldn't catch a
    regression here.
    '''

    def test_unattributable_bits_land_on_root_not_dropped(self):
        a, b = leaf('a', '1'), leaf('b', '2')
        tree = Sequence(None, [a, b])
        tree.id = '3'
        set_parent([a, b], tree)

        predecessors = compute_predecessors(tree)
        self.assertEqual(predecessors[a], (None, False))

        # two rows for 'a' (the None-predecessor leaf) - not realistic
        # from a real trace of this exact tree, but predecessor_totals
        # must handle it correctly regardless of how rows were produced
        rows = [('a', 10.0, 0.5, 1.0), ('a', 20.0, 0.25, 2.0), ('b', 5.0, 0.5, 1.0)]
        totals, _, unattributable_event_count = predecessor_totals(tree, rows, predecessors)

        self.assertEqual(unattributable_event_count, 2)
        containment, _ = surprise_totals(tree, rows)
        self.assertEqual(totals[tree], containment[tree])


class PredecessorTotalsDuplicateLeavesTest(unittest.TestCase):
    '''
    Model: seq(x, dup1, y, dup2, z) - two leaves both named 'dup', each
    with a genuinely different predecessor (dup1's is x, dup2's is y).
    Unlike containment's LCA punt, there's no single converging ancestor
    to charge - surprise is shared evenly across both predecessor chains.
    '''

    def setUp(self):
        self.x = leaf('x', '1')
        self.dup1 = leaf('dup', '2')
        self.y = leaf('y', '3')
        self.dup2 = leaf('dup', '4')
        self.z = leaf('z', '5')
        self.tree = Sequence(None, [self.x, self.dup1, self.y, self.dup2, self.z])
        self.tree.id = '6'
        set_parent([self.x, self.dup1, self.y, self.dup2, self.z], self.tree)

    def test_predecessors_differ_per_duplicate_leaf(self):
        predecessors = compute_predecessors(self.tree)
        self.assertEqual(predecessors[self.dup1], (self.x, False))
        self.assertEqual(predecessors[self.dup2], (self.y, False))

    def test_bits_shared_evenly_across_predecessor_chains(self):
        predecessors = compute_predecessors(self.tree)
        rows = [('dup', 1.0, 0.5, 2.0), ('dup', 1.0, 0.25, 4.0)]
        totals, ambiguous_event_count, unattributable_event_count = \
            predecessor_totals(self.tree, rows, predecessors)

        # each event's bits split 50/50 between x's and y's chains
        self.assertEqual(totals[self.x], 3.0)
        self.assertEqual(totals[self.y], 3.0)
        self.assertEqual(totals[self.tree], 6.0)
        self.assertEqual(ambiguous_event_count, 0.0)
        self.assertEqual(unattributable_event_count, 0.0)

        containment, _ = surprise_totals(self.tree, rows)
        self.assertEqual(totals[self.tree], containment[self.tree])


if __name__ == '__main__':
    unittest.main()
