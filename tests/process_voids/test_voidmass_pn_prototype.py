'''
Prototype: does the classical/standard Petri-net alignment path
(skipalignments.alignall.align_pn_all, via EbiOccurance.build_petri_net)
give voidmass a deficit count that scales with subtree SIZE on total
ablation, unlike the skip-alignment path (executions()) which lumps an
entirely-unwitnessed subtree into one Skip move regardless of size?

This is a throwaway prototype file: process_voids.voidmass_pn does not
exist in the main package yet. If this confirms the classical path is
viable, its contents get folded into coveragemass.py/merged elsewhere;
if not, both this file and voidmass_pn.py get deleted. Either way this
file is not meant to survive as-is.
'''

import unittest
from types import SimpleNamespace
from unittest.mock import patch

from skipalignments.processtree import Sequence, Activity, Xor, Tau, Loop

from skipalignments.alignall import align_pn_all

from lab.fixtures import build_running_example_tree
from process_voids.coveragemass import min_activity_count
from process_voids.voidmass_pn import (
    build_id_net, align_variant, align_variant_all, deficit_by_node, terms_by_node,
    voidmass_table_pn, coverage_by_alignment_pn, timed_out_movecount_bound,
)


def _sequence_of(labels, cost=100000):
    activities = []
    for i, label in enumerate(labels):
        a = Activity(None, label, cost)
        a.id = f'leaf_{label}_{i}'
        activities.append(a)
    tree = Sequence(None, activities)
    tree.id = 'seq_' + '_'.join(labels)
    for a in activities:
        a.set_parent(tree)
    return tree, activities


class TotalAblationSizeSensitivityTest(unittest.TestCase):
    '''
    Mirrors the E1 fixture that caught skip-alignments' lumping: a
    2-activity and an 8-activity subprocess, each entirely unwitnessed
    (the observed trace has NONE of their activities). Under classical
    alignments there is no Skip(subtree) construct - every missing leaf
    must become its own model-move - so deficit should scale 1:1 with
    subtree size, a straight 4x ratio between the two.
    '''

    def setUp(self):
        small_labels = ['s1', 's2']
        large_labels = [f'l{i}' for i in range(8)]
        witness = Activity(None, 'w', 100000)
        witness.id = 'leaf_w'

        small_tree, self.small_leaves = _sequence_of(small_labels)
        large_tree, self.large_leaves = _sequence_of(large_labels)

        root = Sequence(None, [small_tree, large_tree, witness])
        root.id = 'root'
        small_tree.set_parent(root)
        large_tree.set_parent(root)
        witness.set_parent(root)

        self.root = root
        self.small_tree = small_tree
        self.large_tree = large_tree

        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = build_id_net(root)

    def test_deficit_scales_with_subtree_size(self):
        # Only 'w' is observed - both subprocesses are totally ablated.
        alignment = align_variant(['w'], self.net, self.im, self.fm,
                                   self.activity_to_id, self.tau_ids)
        table = deficit_by_node(alignment, self.root, self.activity_to_id, self.tau_ids)

        self.assertEqual(table[self.small_tree], 2)
        self.assertEqual(table[self.large_tree], 8)
        self.assertEqual(table[self.large_tree], 4 * table[self.small_tree])


class E1KnownLimitIsFixedTest(unittest.TestCase):
    '''
    Reuses the EXACT fixture from test_voidmass.py's SizeSensitivityTest
    (same labels: 2 small 's0'/'s1', 8 big 'b0'..'b7', witness 'other') -
    not just a structurally-similar one - so this is directly traceable
    against that test's test_total_ablation_breaks_size_preservation_via_lumping,
    which pins the OLD (skip-alignments) behaviour: both subprocesses
    entirely missing, both lump to deficit=1, ratio 1.0 - documented
    there as a "known limit", not a bug in voidmass_terms itself.

    On the classical path, deficit correctly comes out 2 and 8 (ratio
    4.0, matching the 50%-ablation case's ratio) - confirming that
    "known limit" was an artifact of skip-alignments' lumped normal
    form, not a fundamental property of voidmass.
    '''

    def setUp(self):
        small_acts = [Activity(None, f's{i}', 100000) for i in range(2)]
        for i, act in enumerate(small_acts):
            act.id = f'small{i}'
        self.small = Sequence(None, small_acts)
        self.small.id = 'small'
        for act in small_acts:
            act.set_parent(self.small)

        big_acts = [Activity(None, f'b{i}', 100000) for i in range(8)]
        for i, act in enumerate(big_acts):
            act.id = f'big{i}'
        self.big = Sequence(None, big_acts)
        self.big.id = 'big'
        for act in big_acts:
            act.set_parent(self.big)

        self.other = Activity(None, 'other', 100000)
        self.other.id = 'other'

        self.tree = Sequence(None, [self.other, self.small, self.big])
        self.tree.id = 'root'
        self.other.set_parent(self.tree)
        self.small.set_parent(self.tree)
        self.big.set_parent(self.tree)

        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = build_id_net(self.tree)

    def test_total_ablation_size_preservation_now_holds(self):
        alignment = align_variant(['other'], self.net, self.im, self.fm,
                                   self.activity_to_id, self.tau_ids)
        table = deficit_by_node(alignment, self.tree, self.activity_to_id, self.tau_ids)

        self.assertEqual(table[self.small], 2)
        self.assertEqual(table[self.big], 8)
        self.assertEqual(table[self.big], 4 * table[self.small])


class TauLeafIsNotADeficitTest(unittest.TestCase):
    '''
    Regression test for a bug found and now fixed upstream in
    skip-alignments: align_pn_all's cost function used to only treat a
    model transition as free when its label was None or the literal
    string started with "TAU" - a convention that didn't survive
    EbiOccurance.build_petri_net's id-substitution (transition labels
    become opaque ids like '5', not "TAU..." text). The running
    example's schedule_choice = Xor(s, tau) has exactly this shape:
    'tau' is a genuine skipalignments.processtree.Tau leaf
    (model_move_cost=0), not pm4py-invisible, so it kept a real id
    after renaming and was miscosted as a 100000 deviation instead of a
    free model alternative. Fixed by build_id_net/align_variant now
    passing tau_ids through to align_pn_all's tau_ids parameter. On the
    variant ('o','a','p') - approval resolved with no redo,
    schedule_choice resolved via the tau alternative since 's' never
    appears - root/schedule_choice deficit should be 0 (choosing the
    model's own free branch is not a deviation), not 1.
    '''

    def setUp(self):
        self.tree = build_running_example_tree()
        self.o, self.approval, self.sched, self.p = self.tree.children
        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = build_id_net(self.tree)

    def test_tau_alternative_is_not_counted_as_deficit(self):
        alignment = align_variant(['o', 'a', 'p'], self.net, self.im, self.fm,
                                   self.activity_to_id, self.tau_ids, timeout=30)
        table = deficit_by_node(alignment, self.tree, self.activity_to_id, self.tau_ids)

        self.assertEqual(table[self.sched], 0)
        self.assertEqual(table[self.tree], 0)


class RunningExampleCrossCheckTest(unittest.TestCase):
    '''
    Cross-check against test_voidmass.py's RunningExampleVoidmassTest,
    which pins voidmass-brief.md's own reference table using the
    skip-alignments path. If the classical path reproduces the same
    numbers here, that's expected - see session notes: this fixture's
    only entirely-missing subtree (approval = Loop(a,e), when the whole
    loop is skipped) never actually exposes the lumping bug, because
    only the do-child 'a' is mandatory when the loop runs zero times -
    'e' is legitimately optional either way. A real divergence would
    need a multi-mandatory-leaf subtree instead (see
    TotalAblationSizeSensitivityTest for that case, hand-built since
    this fixture doesn't have one).
    '''

    def setUp(self):
        self.tree = build_running_example_tree()
        self.o, self.approval, self.sched, self.p = self.tree.children
        self.a, self.e = self.approval.children
        self.variant_probs = {
            ('o', 'a', 's', 'p'): 2 / 6,
            ('o', 'a', 'e', 'a', 's', 'p'): 1 / 6,
            ('o', 's', 'p'): 2 / 6,
            ('o', 'a', 'p'): 1 / 6,
        }
        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = build_id_net(self.tree)
        result = voidmass_table_pn(self.tree, self.variant_probs, self.net, self.im,
                                    self.fm, self.activity_to_id, self.tau_ids, timeout=30)
        self.table, self.skip_dict = result.table, result.skip_dict

    def test_matches_skip_alignments_reference_table(self):
        cases = [
            ('N (root)', self.tree, 0.333, 0.080),
            ('approval', self.approval, 0.333, 0.080),
            ('a', self.a, 0.333, 0.080),
            ('e', self.e, 0.0, 0.0),
            ('o', self.o, 0.0, 0.0),
            ('sched', self.sched, 0.0, 0.0),
            ('p', self.p, 0.0, 0.0),
        ]
        for name, node, deficit_exp, variant2_exp in cases:
            with self.subTest(node=name):
                row = self.table[node]
                # No variant times out in this fixture, so lower==upper -
                # either reading is the real (pre-fix) value.
                self.assertAlmostEqual(row['deficit_lower'], deficit_exp, places=3)
                self.assertAlmostEqual(row['deficit_upper'], deficit_exp, places=3)
                self.assertAlmostEqual(row['voidmass_process_lower'], variant2_exp, places=3)
                self.assertAlmostEqual(row['voidmass_process_upper'], variant2_exp, places=3)


class TiedAlignmentDoubleCountingTest(unittest.TestCase):
    '''
    Reproduces the claims-fixture finding (session notes, claim_25): a
    reordered pair inside an optional Xor(Sequence(a, b), Tau) has more
    than one equal-cost repair, and align_pn_all's all-optimal search
    returns every commuting reordering of the irrelevant log/tau moves
    as a SEPARATE "optimal alignment" - inflating whichever causal
    story happens to have more free slots to permute, well before any
    dedup logic runs. b:a (b observed before a, but the model requires
    a before b) has three genuinely distinct causal repairs: skip the
    whole Xor branch (both events become log moves, 0 deficit), model-
    move a + sync b, or sync a + model-move b - each equally valid, but
    the raw alignment count across them is NOT 1:1:1 (some have more
    commuting-order duplicates than others). voidmass_table_pn must not
    let that raw count skew the weighting - see the dedup fix.
    '''

    def setUp(self):
        a = Activity(None, 'a', 100000)
        a.id = 'a'
        b = Activity(None, 'b', 100000)
        b.id = 'b'
        seq = Sequence(None, [a, b])
        seq.id = 'seq'
        a.set_parent(seq)
        b.set_parent(seq)
        tau = Tau(None, 'skip', 0)
        tau.id = 'tau'
        self.tree = Xor(None, [seq, tau])
        self.tree.id = 'root'
        seq.set_parent(self.tree)
        tau.set_parent(self.tree)

        self.seq = seq
        self.a = a
        self.b = b
        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = build_id_net(self.tree)

    def test_deficit_is_symmetric_between_a_and_b_after_dedup(self):
        variant_probs = {('b', 'a'): 1.0}
        table = voidmass_table_pn(self.tree, variant_probs, self.net, self.im, self.fm,
                                  self.activity_to_id, self.tau_ids, timeout=30).table
        # Three equally-valid causal stories (skip / blame-a / blame-b) -
        # uniform-over-signatures weighting gives 'a' and 'b' the same
        # deficit share (1/3 each), not whatever ratio the raw,
        # commuting-order-inflated alignment count happened to produce.
        self.assertAlmostEqual(table[self.a]['deficit_lower'], table[self.b]['deficit_lower'],
                                places=9)
        self.assertAlmostEqual(table[self.a]['deficit_lower'], 1 / 3, places=9)


class PooledAlignmentMassTest(unittest.TestCase):
    '''
    alignment_mass_pooled/voidmass_process: POOLED matchcount/movecount
    (summed across every variant/execution, not averaged per-execution)
    - the divisor voidmass_process/voidmass_subprocess need, and still a
    legitimate quantity in its own right. NOT what coverage_by_alignment_pn
    uses any more - see NonPooledAlignmentCoverageTest below and
    coverage_by_alignment_pn's own docstring for why an earlier version
    of that function used this pooled quantity and that was wrong:
    \\covermove's formal definition (defn:move-coverage) averages
    per-execution, it doesn't pool.
    '''

    def setUp(self):
        small_labels = ['s1', 's2']
        large_labels = [f'l{i}' for i in range(8)]
        witness = Activity(None, 'w', 100000)
        witness.id = 'leaf_w'

        small_tree, _ = _sequence_of(small_labels)
        large_tree, _ = _sequence_of(large_labels)

        root = Sequence(None, [small_tree, large_tree, witness])
        root.id = 'root'
        small_tree.set_parent(root)
        large_tree.set_parent(root)
        witness.set_parent(root)

        self.root = root
        self.small_tree = small_tree
        self.large_tree = large_tree
        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = build_id_net(root)

    def test_full_conformance_gives_pooled_mass_one(self):
        variant_probs = {('s1', 's2', 'l0', 'l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'w'): 1.0}
        table = voidmass_table_pn(self.root, variant_probs, self.net, self.im, self.fm,
                                  self.activity_to_id, self.tau_ids).table
        for node in (self.root, self.small_tree, self.large_tree):
            with self.subTest(node=node):
                self.assertAlmostEqual(table[node]['alignment_mass_pooled_lower'], 1.0, places=6)
                self.assertAlmostEqual(table[node]['alignment_mass_pooled_upper'], 1.0, places=6)

    def test_total_ablation_gives_pooled_mass_zero(self):
        variant_probs = {('w',): 1.0}
        table = voidmass_table_pn(self.root, variant_probs, self.net, self.im, self.fm,
                                  self.activity_to_id, self.tau_ids).table
        for node in (self.small_tree, self.large_tree):
            with self.subTest(node=node):
                self.assertAlmostEqual(table[node]['alignment_mass_pooled_lower'], 0.0, places=6)
                self.assertAlmostEqual(table[node]['alignment_mass_pooled_upper'], 0.0, places=6)

    def test_partial_ablation_matches_expected_ratio(self):
        # 1 of 2 small activities observed, 4 of 8 large - 50% either way
        variant_probs = {('s1', 'l0', 'l1', 'l2', 'l3', 'w'): 1.0}
        table = voidmass_table_pn(self.root, variant_probs, self.net, self.im, self.fm,
                                  self.activity_to_id, self.tau_ids).table
        # not asserting an exact value here (that's SizeSensitivityTest's
        # job on the skip-alignments path) - just that pooling behaves
        # sanely: strictly between 0 and 1, matching neither extreme.
        for node in (self.small_tree, self.large_tree):
            with self.subTest(node=node):
                for key in ('alignment_mass_pooled_lower', 'alignment_mass_pooled_upper'):
                    mass = table[node][key]
                    self.assertGreater(mass, 0.0)
                    self.assertLess(mass, 1.0)


class NonPooledAlignmentCoverageTest(unittest.TestCase):
    '''
    coverage_by_alignment_pn (\\covermove, defn:move-coverage): per-
    execution match/movecount ratios averaged - NOT pooled - verified
    term-by-term against the formal definition (session notes). Uses
    coverage_by_alignment_pn(node, skip_prob, skip_dict, variant_probs)
    - skip_dict is voidmass_table_pn's result.skip_dict.
    '''

    def setUp(self):
        small_labels = ['s1', 's2']
        large_labels = [f'l{i}' for i in range(8)]
        witness = Activity(None, 'w', 100000)
        witness.id = 'leaf_w'

        small_tree, _ = _sequence_of(small_labels)
        large_tree, _ = _sequence_of(large_labels)

        root = Sequence(None, [small_tree, large_tree, witness])
        root.id = 'root'
        small_tree.set_parent(root)
        large_tree.set_parent(root)
        witness.set_parent(root)

        self.root = root
        self.small_tree = small_tree
        self.large_tree = large_tree
        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = build_id_net(root)

    def test_total_ablation_gives_zero_coverage_regardless_of_skip_prob(self):
        variant_probs = {('w',): 1.0}
        skip_dict = voidmass_table_pn(self.root, variant_probs, self.net, self.im, self.fm,
                                      self.activity_to_id, self.tau_ids).skip_dict
        self.assertAlmostEqual(
            coverage_by_alignment_pn(self.small_tree, 0.3, skip_dict, variant_probs), 0.0)

    def test_full_conformance_reduces_to_one_minus_skip_prob(self):
        variant_probs = {('s1', 's2', 'l0', 'l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'w'): 1.0}
        skip_dict = voidmass_table_pn(self.root, variant_probs, self.net, self.im, self.fm,
                                      self.activity_to_id, self.tau_ids).skip_dict
        self.assertAlmostEqual(
            coverage_by_alignment_pn(self.root, 0.0, skip_dict, variant_probs), 1.0)


class PooledVsNonPooledDivergenceTest(unittest.TestCase):
    '''
    Proves the fix has a real effect, not just a refactor: a case with
    two executions of DIFFERENT size (one loop iteration vs three) and a
    deficit only in the smaller one, where pooling (sum counts, divide
    once) and per-execution averaging (average ratios, equal weight per
    execution) give different, hand-derived numbers.

    Tree: Sequence(a, Loop(x, y)) - do=x, redo=y.
      Variant A (weight 0.5): 'x' only - a MISSING (deficit 1), one loop
        iteration, x matched. movecount=2 (a,x), matchcount=1 (x),
        ratio=0.5.
      Variant B (weight 0.5): 'a','x','y','x' - fully conforming, two
        loop iterations. movecount=4, matchcount=4, ratio=1.0.

    Per-execution average (what coverage_by_alignment_pn now computes):
      0.5*0.5 + 0.5*1.0 = 0.75
    Pooled (what it used to compute, still alignment_mass_pooled):
      deficit_sum = 0.5*1 + 0.5*0 = 0.5
      movecount_sum = 0.5*2 + 0.5*4 = 3.0
      1 - 0.5/3.0 = 5/6 (~0.8333)
    0.75 != 5/6 - the two are genuinely different quantities.
    '''

    def setUp(self):
        a = Activity(None, 'a', 100000)
        a.id = 'a'
        x = Activity(None, 'x', 100000)
        x.id = 'x'
        y = Activity(None, 'y', 100000)
        y.id = 'y'
        loop = Loop(None, [x, y])
        loop.id = 'loop'
        x.set_parent(loop)
        y.set_parent(loop)
        self.tree = Sequence(None, [a, loop])
        self.tree.id = 'root'
        a.set_parent(self.tree)
        loop.set_parent(self.tree)

        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = \
            build_id_net(self.tree)

    def test_pooled_and_non_pooled_diverge_as_hand_derived(self):
        variant_probs = {('x',): 0.5, ('a', 'x', 'y', 'x'): 0.5}
        result = voidmass_table_pn(self.tree, variant_probs, self.net, self.im, self.fm,
                                   self.activity_to_id, self.tau_ids, timeout=30)
        table, skip_dict = result.table, result.skip_dict

        # no variant times out in this fixture, so lower==upper
        pooled = table[self.tree]['alignment_mass_pooled_lower']
        non_pooled = coverage_by_alignment_pn(self.tree, 0.0, skip_dict, variant_probs)

        self.assertAlmostEqual(pooled, 5 / 6, places=6)
        self.assertAlmostEqual(non_pooled, 0.75, places=6)
        self.assertNotAlmostEqual(pooled, non_pooled, places=3)


class CoverageByAlignmentPnTimedOutRatioTest(unittest.TestCase):
    '''coverage_by_alignment_pn threads timed_out_ratio straight through
    to coveragemass.alignment_mass (see that function's own tests for the
    underlying mechanics) - hand-built skip_dict/variant_probs here, no
    real alignment search needed, since this is only checking the
    wiring: one variant with a real (perfect-fit) alignment, one variant
    mapped to an empty list standing in for a timed-out search.'''

    def setUp(self):
        self.tree, (self.a, self.b) = _sequence_of(['a', 'b'])
        fit_path = [('a', self.a), ('b', self.b)]
        self.skip_dict = {
            'a, b': [SimpleNamespace(path=fit_path)],
            'a': [],  # timed out - present, zero alignments
        }
        self.variant_probs = {('a', 'b'): 0.5, ('a',): 0.5}

    def test_default_excludes_the_timed_out_variant(self):
        self.assertAlmostEqual(
            coverage_by_alignment_pn(self.tree, 0.0, self.skip_dict, self.variant_probs),
            0.5, places=6)

    def test_upper_bound_treats_it_as_perfectly_matched(self):
        self.assertAlmostEqual(
            coverage_by_alignment_pn(self.tree, 0.0, self.skip_dict, self.variant_probs,
                                      timed_out_ratio=1.0),
            1.0, places=6)

    def test_lower_bound_treats_it_as_fully_missed(self):
        self.assertAlmostEqual(
            coverage_by_alignment_pn(self.tree, 0.0, self.skip_dict, self.variant_probs,
                                      timed_out_ratio=0.0),
            0.5, places=6)


def _patched_timeout(variants_to_time_out):
    '''patch.object context forcing align_variant_all to return zero
    alignments (a per-variant timeout) for the given variants, and run
    for real on every other one.'''
    import process_voids.voidmass_pn as voidmass_pn_module
    real = voidmass_pn_module.align_variant_all
    timed_out = {tuple(v) for v in variants_to_time_out}

    def fake(variant, *args, **kwargs):
        if tuple(variant) in timed_out:
            return []
        return real(variant, *args, **kwargs)

    return patch.object(voidmass_pn_module, 'align_variant_all', side_effect=fake)


def _loop_tree():
    '''loop(seq(a,b), tau) - the model whose loop lets model moves
    outnumber any leaf count: each observed 'a' can force another
    iteration, each needing its own model move of 'b'.'''
    a = Activity(None, 'a', 100000)
    a.id = 'a'
    b = Activity(None, 'b', 100000)
    b.id = 'b'
    body = Sequence(None, [a, b])
    body.id = 'body'
    a.set_parent(body)
    b.set_parent(body)
    tau = Tau(None, 'tau', 0)
    tau.id = 'tau'
    loop = Loop(None, [body, tau])
    loop.id = 'loop'
    body.set_parent(loop)
    tau.set_parent(loop)
    return loop


def _xor_tau_tree():
    '''seq(xor(seq(c,d), tau), e) - cheapest path takes the free tau
    branch, so min_activity_count(root) = 1.'''
    c = Activity(None, 'c', 100000)
    c.id = 'c'
    d = Activity(None, 'd', 100000)
    d.id = 'd'
    cd = Sequence(None, [c, d])
    cd.id = 'cd'
    c.set_parent(cd)
    d.set_parent(cd)
    tau = Tau(None, 'tau', 0)
    tau.id = 'tau'
    choice = Xor(None, [cd, tau])
    choice.id = 'choice'
    cd.set_parent(choice)
    tau.set_parent(choice)
    e = Activity(None, 'e', 100000)
    e.id = 'e'
    root = Sequence(None, [choice, e])
    root.id = 'root'
    choice.set_parent(root)
    e.set_parent(root)
    return root


class TimedOutVariantDoesNotCrashTest(unittest.TestCase):
    '''
    Regression test for a failure seen in real runs: align_variant_all
    returning zero alignments for a variant (a per-variant timeout) used
    to make voidmass_table_pn raise ZeroDivisionError and discard every
    other variant's completed work.

    M: model seq(a,b), two variants each weight 0.5: <a,b> aligns for
    real (perfect fit), <a> is forced to time out. The timed-out variant
    is bounded by X_max = 2|sigma| + min_activity_count(root) = 2*1 + 2
    = 4, the SAME X at every node, so it contributes w*X_max = 2.
    '''

    def setUp(self):
        self.tree, (self.a, self.b) = _sequence_of(['a', 'b'])
        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = \
            build_id_net(self.tree)
        self.variant_probs = {('a', 'b'): 0.5, ('a',): 0.5}

    def test_does_not_raise_and_reports_the_bounds(self):
        with _patched_timeout([('a',)]):
            result = voidmass_table_pn(self.tree, self.variant_probs, self.net, self.im,
                                       self.fm, self.activity_to_id, self.tau_ids, timeout=30)

        self.assertEqual(result.timed_out_count, 1)
        self.assertAlmostEqual(result.timed_out_weight, 0.5, places=6)

        # voidmass_movecount is the OBSERVED total from completed variants
        # only (<a,b>: movecount 2 at weight 0.5); the bounds' shared
        # denominator is reported separately as movecount_bound.
        root = result.table[self.tree]
        self.assertAlmostEqual(root['movecount'], 1.0, places=6)
        self.assertAlmostEqual(root['movecount_bound'], 1.0 + 2.0, places=6)
        self.assertAlmostEqual(root['deficit_lower'], 0.0, places=6)
        self.assertAlmostEqual(root['deficit_upper'], 2.0, places=6)
        self.assertAlmostEqual(root['voidmass_subprocess_lower'], 0.0, places=6)
        self.assertAlmostEqual(root['voidmass_subprocess_upper'], 2.0 / 3.0, places=6)

        # A leaf gets the same w*X_max as the root. That shared X is what
        # keeps voidmass_process's bounds valid: its numerator is the
        # node's deficit, its denominator the ROOT's movecount, and a
        # bound needs the same X on both.
        leaf = result.table[self.a]
        self.assertAlmostEqual(leaf['movecount'], 0.5, places=6)
        self.assertAlmostEqual(leaf['movecount_bound'], 0.5 + 2.0, places=6)
        self.assertAlmostEqual(leaf['voidmass_subprocess_upper'], 2.0 / 2.5, places=6)
        self.assertAlmostEqual(leaf['voidmass_process_lower'], 0.0, places=6)
        self.assertAlmostEqual(leaf['voidmass_process_upper'], 2.0 / 3.0, places=6)

    def test_no_timeout_means_zero_diagnostics_and_coincident_bounds(self):
        result = voidmass_table_pn(self.tree, {('a', 'b'): 1.0}, self.net, self.im, self.fm,
                                   self.activity_to_id, self.tau_ids, timeout=30)
        self.assertEqual(result.timed_out_count, 0)
        self.assertEqual(result.timed_out_weight, 0.0)
        for row in result.table.values():
            self.assertEqual(row['movecount'], row['movecount_bound'])
            self.assertEqual(row['deficit_lower'], row['deficit_upper'])


class TimedOutMovecountBoundTest(unittest.TestCase):
    '''
    timed_out_movecount_bound(|sigma|, tree) = 2|sigma| + C_root must
    bound the movecount of EVERY optimal alignment of sigma, loops
    included. A bound from the tree's shape alone cannot: loop(seq(a,b),
    tau) on <a,a,a> has an optimal alignment that syncs every 'a' and
    model-moves 'b' once per iteration - movecount 6, above a leaf-count
    bound of |sigma| + 2 = 5.
    '''

    def setUp(self):
        self.loop = _loop_tree()
        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = \
            build_id_net(self.loop)

    def _movecounts(self, trace):
        alignments = align_variant_all(trace, self.net, self.im, self.fm, self.activity_to_id,
                                       self.tau_ids, id_loop_list=self.id_loop_list, timeout=30)
        return [terms_by_node(al, self.loop, self.activity_to_id, self.tau_ids)[self.loop][1]
                for al in alignments]

    def test_bound_holds_for_every_optimal_alignment(self):
        for trace in (['a', 'a', 'a'], ['a'] * 5, ['b'], ['x', 'x'], ['a', 'x', 'b']):
            bound = timed_out_movecount_bound(len(trace), self.loop)
            for m in self._movecounts(trace):
                with self.subTest(trace=trace):
                    self.assertLessEqual(m, bound)

    def test_loop_exceeds_a_leaf_count_bound(self):
        movecounts = self._movecounts(['a', 'a', 'a'])
        self.assertIn(6, movecounts)
        self.assertGreater(max(movecounts), 3 + 2)


class CostModelPinTest(unittest.TestCase):
    '''
    timed_out_movecount_bound's derivation rests on align_pn_all's cost
    model: log move = labelled model move = 100000, sync = 0, every
    silent/tau transition = 0. Given that, no optimal alignment costs
    more than "every event a log move, plus a complete cheapest model
    path of C_root non-silent moves", so log + model <= |sigma| + C_root.

    Pinned behaviourally rather than by reading constants: a trace of a
    single unmodelled activity has exactly that alignment as its
    optimum, so its cost must be exactly 100000 * (1 + C_root). That
    also checks the premise the bound needs - a complete (final-state-
    reaching) model path with min_activity_count(root) non-silent moves
    exists. If the cost model changes, this fails and the bound must be
    re-derived, rather than silently going invalid.
    '''

    def test_unmodelled_single_event_costs_one_log_move_plus_cheapest_path(self):
        for tree in (_loop_tree(), _xor_tau_tree()):
            net, im, fm, _activity_to_id, tau_ids, id_loop_list = build_id_net(tree)
            _t, (alignments, code, _first) = align_pn_all(
                ['x'], net, im, fm, id_loop_list, timeout=30, tau_ids=tau_ids)[0]
            with self.subTest(tree=tree.id):
                self.assertEqual(code, 0)
                for alignment in alignments:
                    self.assertEqual(alignment['cost'], 100000 * (1 + min_activity_count(tree)))


class TimedOutBoundsAreSoundTest(unittest.TestCase):
    '''
    The bounds must bracket the value the cell would have had if the
    variant had NOT timed out. Model loop(seq(a,b), tau): <a,b> at weight
    0.9 fits perfectly, <a,a,a,a,a> at weight 0.1 has tied optimal
    alignments iterating the loop up to five times. Computed once for
    real, then again with <a,a,a,a,a> forced to time out.

    Also records the previous fix's error: substituting the cheapest
    path's length (min_activity_count) as the timed-out variant's
    deficit and movecount lands BELOW the real value. Per-variant
    deficit <= movecount says nothing about a pooled ratio, and the
    cheapest path's length bounds a variant's movecount from below, not
    above.
    '''

    def setUp(self):
        self.loop = _loop_tree()
        self.net, self.im, self.fm, self.activity_to_id, self.tau_ids, self.id_loop_list = \
            build_id_net(self.loop)
        self.slow = ('a',) * 5
        self.variant_probs = {('a', 'b'): 0.9, self.slow: 0.1}

    def _run(self):
        return voidmass_table_pn(self.loop, self.variant_probs, self.net, self.im, self.fm,
                                 self.activity_to_id, self.tau_ids,
                                 id_loop_list=self.id_loop_list, timeout=30)

    def test_bounds_bracket_the_real_value_at_every_node(self):
        real = self._run()
        with _patched_timeout([self.slow]):
            bounded = self._run()

        self.assertEqual(real.timed_out_count, 0)
        self.assertEqual(bounded.timed_out_count, 1)
        for node, real_row in real.table.items():
            for metric in ('deficit', 'voidmass_subprocess', 'voidmass_process'):
                truth = real_row[f'{metric}_lower']
                with self.subTest(node=node.id, metric=metric):
                    self.assertEqual(real_row[f'{metric}_lower'], real_row[f'{metric}_upper'])
                    self.assertLessEqual(bounded.table[node][f'{metric}_lower'], truth + 1e-9)
                    self.assertGreaterEqual(bounded.table[node][f'{metric}_upper'], truth - 1e-9)

    def test_cheapest_path_substitution_undershoots_the_real_value(self):
        real = self._run()
        with _patched_timeout([self.slow]):
            bounded = self._run()

        truth = real.table[self.loop]['voidmass_subprocess_lower']
        d0 = bounded.table[self.loop]['deficit_lower']
        m0 = bounded.table[self.loop]['movecount']
        wc = bounded.timed_out_weight * min_activity_count(self.loop)
        previous_upper = (d0 + wc) / (m0 + wc)
        self.assertGreater(truth, previous_upper)


if __name__ == '__main__':
    unittest.main()
