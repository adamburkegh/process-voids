'''
Hand-derived regression tests for process_voids.voidsalign (Definitions
[Skip-Weighted Move Counts] / [Voidage by Skip-Weighted Alignment
Moves], corrected form), following test_coverage_by_alignment.py's
fixture pattern: the real (ebi-free) skip-aligner builds each
skip_dict, and expected values are hand-derived from the alignment it
is known to find for these small models. skip_probs is always supplied
directly (as every other coveragemass-family test file already does),
independent of any real skip-probability estimation - the point under
test is the size-share arithmetic, not skipprob itself.

The previous (superseded) form of this metric was skipprob times a
match/movecount COMPLETENESS ratio - an inverted U that collapses to
zero exactly where a subprocess is wholly missing, since a skipped
execution's ratio is always 0/w. This form is skipprob times a SIZE
share (smovetotal(msub)/smovetotal(whole model)) instead, which does
not collapse under skipping - see the module docstring.
'''

import unittest

from skipalignments import Activity, Sequence, Tau, Xor, Skip, Aligner

from process_voids.coveragemass import make_executions_cache, min_activity_count
from process_voids.voidsalign import smatchcount, smovecount, smovetotal, voidsalign

ACT_COST = 100000

Aligner.set_level_incentive(0)


def leaf(cls, name, node_id, cost=ACT_COST):
    node = cls(None, name, cost)
    node.id = node_id
    return node


def variant_key(trace):
    return ', '.join(trace)


def align(tree, trace):
    states, _ = Aligner(tree).align_normal_form(list(trace), [ACT_COST] * len(trace), True, timeout=100)
    return states


def build(tree, traces_with_weights):
    """traces_with_weights: {trace tuple: weight}. Runs the real (ebi-free)
    alignment search for each distinct variant."""
    skip_dict = {}
    variant_probs = {}
    for trace, weight in traces_with_weights.items():
        skip_dict[variant_key(trace)] = align(tree, trace)
        variant_probs[trace] = weight
    return skip_dict, variant_probs


class SmovecountWeightingTest(unittest.TestCase):
    '''
    smovecount, unlike coveragemass.movecount, weights a skip move by
    the minimum executable length of the subprocess it stands for - a
    lumped composite skip carries the size of its whole subtree, not 1.
    '''

    def test_lumped_composite_skip_carries_subprocess_size(self):
        x = leaf(Activity, 'x', '1')
        y = leaf(Activity, 'y', '2')
        xy = Sequence(None, [x, y])
        xy.id = '3'
        x.set_parent(xy)
        y.set_parent(xy)
        execution = [('>>', Skip(xy, 1))]
        self.assertEqual(min_activity_count(xy), 2)
        self.assertEqual(smovecount(execution), 2)

    def test_leaf_skip_carries_weight_one(self):
        z = leaf(Activity, 'z', '1')
        execution = [('>>', Skip(z, 1))]
        self.assertEqual(smovecount(execution), 1)

    def test_smatchcount_counts_only_synchronous_moves(self):
        a = leaf(Activity, 'a', '1')
        b = leaf(Activity, 'b', '2')
        execution = [('a', a), ('>>', Skip(b, 1))]
        self.assertEqual(smatchcount(execution), 1)
        self.assertEqual(smovecount(execution), 2)


class LinearInSkipprobConstantSizeShareTest(unittest.TestCase):
    '''
    The central property under test: on seq(a, b), b's smovecount
    contribution is exactly 1 whether b is a synchronous move (matched)
    or a leaf skip move (min_activity_count(b) = 1) - so b's SIZE share
    of the whole model's skip-weighted moves stays constant (0.5) as b
    goes from always recorded, to half, to never recorded, even though
    skipprob varies from 0 to 0.5 to 1. voidsalign(b) is then exactly
    linear in skipprob: 0, 0.25, 0.5 - never collapsing back to zero at
    full absence, unlike the superseded ratio-based form.
    '''

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)

    def _voidsalign_b(self, traces_with_weights, b_skipprob):
        skip_dict, variant_probs = build(self.tree, traces_with_weights)
        skip_probs = {self.tree: 1.0, self.a: 0.0, self.b: b_skipprob}
        size_share = (smovetotal(self.b, skip_dict, variant_probs)
                      / smovetotal(self.tree, skip_dict, variant_probs))
        value = voidsalign(self.b, self.tree, skip_dict, variant_probs, skip_probs)
        return value, size_share

    def test_b_always_recorded(self):
        value, size_share = self._voidsalign_b({('a', 'b'): 1.0}, 0.0)
        self.assertAlmostEqual(size_share, 0.5, places=6)
        self.assertAlmostEqual(value, 0.0, places=6)

    def test_b_half_recorded(self):
        value, size_share = self._voidsalign_b(
            {('a', 'b'): 0.5, ('a',): 0.5}, 0.5)
        self.assertAlmostEqual(size_share, 0.5, places=6)
        self.assertAlmostEqual(value, 0.25, places=6)

    def test_b_never_recorded(self):
        value, size_share = self._voidsalign_b({('a',): 1.0}, 1.0)
        self.assertAlmostEqual(size_share, 0.5, places=6)
        # the whole point: full absence scores the FULL size share, not
        # zero - the bug this form fixes.
        self.assertAlmostEqual(value, 0.5, places=6)


class LumpedSizeProportionalityTest(unittest.TestCase):
    '''
    A wholly-absent multi-activity block's SIZE share is proportional
    to its aligncost(<>, .) (min_activity_count), not flattened to the
    same weight as a single missing leaf - the lumping fix working.
    seq(a, seq(x,y), c) on <a,c> lumps the entirely-unwitnessed seq(x,y)
    into one skip move of weight 2; seq(a, z, c) on <a,c> skips the
    single leaf z at weight 1 - same surrounding pattern, sizes differ
    2x, and so must the resulting size shares. Also pins x/y (the
    lumped node's own children) at exactly 0 - they have no execution
    at all under the lumped skip, not an inherited share of it.
    '''

    def test_composite_block_scores_double_a_single_leaf(self):
        a1 = leaf(Activity, 'a', '1')
        x = leaf(Activity, 'x', '2')
        y = leaf(Activity, 'y', '3')
        c1 = leaf(Activity, 'c', '4')
        xy = Sequence(None, [x, y])
        xy.id = '5'
        x.set_parent(xy)
        y.set_parent(xy)
        composite_tree = Sequence(None, [a1, xy, c1])
        composite_tree.id = '6'
        a1.set_parent(composite_tree)
        xy.set_parent(composite_tree)
        c1.set_parent(composite_tree)
        composite_skip_dict, composite_probs = build(composite_tree, {('a', 'c'): 1.0})
        composite_share = (smovetotal(xy, composite_skip_dict, composite_probs)
                            / smovetotal(composite_tree, composite_skip_dict, composite_probs))

        a2 = leaf(Activity, 'a', '1')
        z = leaf(Activity, 'z', '2')
        c2 = leaf(Activity, 'c', '3')
        leaf_tree = Sequence(None, [a2, z, c2])
        leaf_tree.id = '4'
        a2.set_parent(leaf_tree)
        z.set_parent(leaf_tree)
        c2.set_parent(leaf_tree)
        leaf_skip_dict, leaf_probs = build(leaf_tree, {('a', 'c'): 1.0})
        leaf_share = (smovetotal(z, leaf_skip_dict, leaf_probs)
                      / smovetotal(leaf_tree, leaf_skip_dict, leaf_probs))

        # composite: smatch 2 (a,c), smovecount = 2 + min_activity_count(xy) = 4 -> share 2/4
        # leaf:      smatch 2 (a,c), smovecount = 2 + min_activity_count(z)  = 3 -> share 1/3
        self.assertAlmostEqual(composite_share, 2 / 4, places=6)
        self.assertAlmostEqual(leaf_share, 1 / 3, places=6)
        self.assertAlmostEqual(composite_share / leaf_share, (2 / 4) / (1 / 3), places=6)

        # the same proportionality survives through voidsalign when both
        # blocks are given the same skipprob.
        composite_skip_probs = {composite_tree: 1.0, xy: 1.0, x: 1.0, y: 1.0}
        composite_value = voidsalign(xy, composite_tree, composite_skip_dict,
                                      composite_probs, composite_skip_probs)
        leaf_value = voidsalign(z, leaf_tree, leaf_skip_dict, leaf_probs,
                                 {leaf_tree: 1.0, z: 1.0})
        self.assertAlmostEqual(composite_value, composite_share, places=6)
        self.assertAlmostEqual(leaf_value, leaf_share, places=6)
        self.assertGreater(composite_value, leaf_value)

        # x and y themselves have NO execution in this alignment at all
        # (coveragemass.executions gives a node beneath a lumped skip no
        # execution there, per Definition [Executions] - the lumped node
        # alone carries the void) - their smovetotal, and hence
        # voidsalign, is 0 regardless of their own skip_prob. Pinned so a
        # regression reintroducing lump inheritance fails loudly here.
        self.assertEqual(smovetotal(x, composite_skip_dict, composite_probs), 0.0)
        self.assertEqual(smovetotal(y, composite_skip_dict, composite_probs), 0.0)
        self.assertEqual(
            voidsalign(x, composite_tree, composite_skip_dict, composite_probs,
                       composite_skip_probs),
            0.0)
        self.assertEqual(
            voidsalign(y, composite_tree, composite_skip_dict, composite_probs,
                       composite_skip_probs),
            0.0)


class WhollySilentSubprocessTest(unittest.TestCase):
    '''
    A wholly-silent subprocess (min_activity_count 0 - here a bare Tau
    leaf) contributes 0 to its own smovetotal AND to its ancestors' -
    its model move is wrapped as a zero-cost TauPath, not a Skip, so
    smovecount never counts it (see smovecount's own docstring).
    '''

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.tau = leaf(Tau, 'tau', '2', cost=0)
        self.tree = Sequence(None, [self.a, self.tau])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.tau.set_parent(self.tree)
        self.skip_dict, self.variant_probs = build(self.tree, {('a',): 1.0})

    def test_silent_leaf_contributes_nothing_to_its_own_total(self):
        self.assertEqual(smovetotal(self.tau, self.skip_dict, self.variant_probs), 0.0)

    def test_silent_leaf_contributes_nothing_to_the_root_total(self):
        # root total = smatchcount(a) alone = 1, not inflated by tau
        self.assertAlmostEqual(
            smovetotal(self.tree, self.skip_dict, self.variant_probs), 1.0, places=6)

    def test_voidsalign_is_zero_regardless_of_skipprob(self):
        skip_probs = {self.tree: 1.0, self.a: 0.0, self.tau: 1.0}
        self.assertAlmostEqual(
            voidsalign(self.tau, self.tree, self.skip_dict, self.variant_probs, skip_probs),
            0.0, places=6)


class SilentChoiceExcludedFromRootTotalTest(unittest.TestCase):
    '''
    model seq(a, xor(b, tau)), log <a>: sync a, silent (TauPath) move on
    the choice - the choice's own execution has smovecount 0 (no sync,
    no Skip move), so it contributes nothing to the root's smovetotal
    either; the root's total is exactly a's own smatchcount.
    '''

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tau = leaf(Tau, 'tau', '3', cost=0)
        self.choice = Xor(None, [self.b, self.tau])
        self.choice.id = '4'
        self.b.set_parent(self.choice)
        self.tau.set_parent(self.choice)
        self.tree = Sequence(None, [self.a, self.choice])
        self.tree.id = '5'
        self.a.set_parent(self.tree)
        self.choice.set_parent(self.tree)
        self.skip_dict, self.variant_probs = build(self.tree, {('a',): 1.0})

    def test_root_total_ignores_the_silent_move(self):
        self.assertAlmostEqual(
            smovetotal(self.tree, self.skip_dict, self.variant_probs), 1.0, places=6)

    def test_choice_total_is_zero(self):
        self.assertEqual(smovetotal(self.choice, self.skip_dict, self.variant_probs), 0.0)


class MissingVariantAlignmentTest(unittest.TestCase):
    '''
    A variant absent from skip_dict (Delta_sigma = {}, eg a timed-out
    search) contributes zero for its own weight share rather than being
    renormalised across the remaining variants - the same
    1/|Delta_sigma| = 0 convention voidmass_terms already uses.
    '''

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.variant_probs = {('a', 'b'): 0.5, ('a',): 0.5}
        # only the perfectly-fitting variant has any recorded alignment
        self.skip_dict = {variant_key(('a', 'b')): align(self.tree, ('a', 'b'))}

    def test_missing_variant_contributes_zero_not_renormalised(self):
        # perfect-fit variant contributes 0.5 * smovecount(2) = 1.0 to
        # the root total; the missing variant contributes nothing (not
        # another 1.0, which renormalising away the missing entry would
        # give).
        self.assertAlmostEqual(
            smovetotal(self.tree, self.skip_dict, self.variant_probs), 1.0, places=6)


class ZeroDenominatorConventionTest(unittest.TestCase):
    '''voidsalign = 0 where smovetotal(tree) = 0 - no ZeroDivisionError,
    and no assumption that skip_probs[pt] is itself meaningful there.'''

    def test_no_alignments_anywhere_gives_zero(self):
        a = leaf(Activity, 'a', '1')
        tree = Sequence(None, [a])
        tree.id = '2'
        a.set_parent(tree)
        skip_dict = {}
        variant_probs = {('a',): 1.0}
        skip_probs = {tree: 1.0, a: 1.0}
        self.assertEqual(voidsalign(a, tree, skip_dict, variant_probs, skip_probs), 0.0)


class SharedExecutionsCacheTest(unittest.TestCase):
    '''A shared executions_cache (coveragemass.make_executions_cache) must
    give the same result as no cache at all.'''

    def test_cache_and_no_cache_agree(self):
        a1 = leaf(Activity, 'a', '1')
        x = leaf(Activity, 'x', '2')
        y = leaf(Activity, 'y', '3')
        c1 = leaf(Activity, 'c', '4')
        xy = Sequence(None, [x, y])
        xy.id = '5'
        x.set_parent(xy)
        y.set_parent(xy)
        tree = Sequence(None, [a1, xy, c1])
        tree.id = '6'
        a1.set_parent(tree)
        xy.set_parent(tree)
        c1.set_parent(tree)
        skip_dict, variant_probs = build(tree, {('a', 'c'): 1.0})
        skip_probs = {tree: 1.0, xy: 1.0, x: 1.0, y: 1.0}

        cache = make_executions_cache(tree)
        for node in (tree, xy, x, y):
            with self.subTest(node=node):
                self.assertAlmostEqual(
                    voidsalign(node, tree, skip_dict, variant_probs, skip_probs,
                               executions_cache=cache),
                    voidsalign(node, tree, skip_dict, variant_probs, skip_probs),
                    places=6)
