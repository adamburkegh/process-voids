'''
Regression tests for process_voids.coveragemass's voidmass functions
(voidmass-brief.md variants 1 and 2), pinned against the brief's own
worked reference table (voidage.py) on the payment running example -
computed here with the real aligner (not the brief's hand-picked
alignments), following test_coverage_by_alignment.py's RunningExampleTest
fixture pattern exactly, since that already uses the same tree/variants/
weights. Confirmed to match the brief's reference numbers exactly - see
session notes; no aligner-vs-brief divergence found for this example
(contrast test_coverage_by_alignment.py's TiesAcrossOptimalAlignmentsTest,
where the real aligner did diverge from a hand-written expectation).
'''

import unittest

from skipalignments import Activity, Sequence, Aligner

from lab.fixtures import build_running_example_log, build_running_example_tree
from process_voids import pvoid
from process_voids.coveragemass import voidmass_table

ACT_COST = 100000

Aligner.set_level_incentive(0)


def variant_key(trace):
    return ', '.join(trace)


def align(tree, trace):
    states, _ = Aligner(tree).align_normal_form(list(trace), [ACT_COST] * len(trace), True, timeout=100)
    return states


class RunningExampleVoidmassTest(unittest.TestCase):
    '''Pins voidmass-brief.md's reference table exactly (see voidage.py).'''

    def setUp(self):
        self.tree = build_running_example_tree()
        self.o, self.approval, self.sched, self.p = self.tree.children
        self.a, self.e = self.approval.children
        self.s, self.tau = self.sched.children
        self.variant_probs = {
            ('o', 'a', 's', 'p'): 2 / 6,
            ('o', 'a', 'e', 'a', 's', 'p'): 1 / 6,
            ('o', 's', 'p'): 2 / 6,
            ('o', 'a', 'p'): 1 / 6,
        }
        self.skip_dict = {
            variant_key(variant): align(self.tree, variant)
            for variant in self.variant_probs
        }
        self.table = voidmass_table(self.tree, self.skip_dict, self.variant_probs)

    def test_reference_table(self):
        cases = [
            ('N (root)', self.tree, 0.333, 4.167, 0.080, 0.080),
            ('approval', self.approval, 0.333, 1.333, 0.080, 0.250),
            ('a', self.a, 0.333, 1.167, 0.080, 0.286),
            ('e', self.e, 0.0, None, 0.0, 0.0),
            ('o', self.o, 0.0, None, 0.0, 0.0),
            ('sched', self.sched, 0.0, None, 0.0, 0.0),
            ('s', self.s, 0.0, None, 0.0, 0.0),
            ('p', self.p, 0.0, None, 0.0, 0.0),
        ]
        for name, node, deficit_exp, movecount_exp, variant2_exp, variant1_exp in cases:
            with self.subTest(node=name):
                row = self.table[node]
                self.assertAlmostEqual(row['deficit'], deficit_exp, places=3)
                if movecount_exp is not None:
                    self.assertAlmostEqual(row['movecount'], movecount_exp, places=3)
                self.assertAlmostEqual(row['voidmass_process'], variant2_exp, places=3)
                self.assertAlmostEqual(row['voidmass_subprocess'], variant1_exp, places=3)

    def test_variant2_is_additive_over_approval(self):
        # approval = a + e under variant 2 (process-relative)
        self.assertAlmostEqual(
            self.table[self.approval]['voidmass_process'],
            self.table[self.a]['voidmass_process'] + self.table[self.e]['voidmass_process'],
            places=6)

    def test_variant2_root_equals_sum_of_children(self):
        total = sum(self.table[c]['voidmass_process'] for c in self.tree.children)
        self.assertAlmostEqual(self.table[self.tree]['voidmass_process'], total, places=6)

    def test_variant1_is_not_additive(self):
        # parent (approval, 0.250) is LESS than its own child (a, 0.286) -
        # the whole point of variant 1 being scale-free, not size-preserving
        self.assertLess(
            self.table[self.approval]['voidmass_subprocess'],
            self.table[self.a]['voidmass_subprocess'])

    def test_variants_agree_at_root(self):
        # the two divisors coincide at the root
        root = self.table[self.tree]
        self.assertAlmostEqual(root['voidmass_process'], root['voidmass_subprocess'], places=6)


class SizeSensitivityTest(unittest.TestCase):
    '''
    E1 from the brief: two subprocesses of different size (2 and 8
    activities), ablated at the same RATE (50% missing each) - not
    wholly missing. Variant 1 (scale-free) must score them equally;
    variant 2 (size-preserving) must score the larger one roughly 4x
    the smaller.

    Deliberately NOT total ablation (0% observed): when a whole
    subtree goes entirely unwitnessed, the aligner's own optimal-
    alignment normal form collapses it into a single lump Skip move
    (see coveragemass.executions' docstring) regardless of how many
    activities are inside it - movecount/deficit then count alignment
    moves, not model activities, and both subtrees score identically
    (found empirically while building this test - a real limit of
    building on skip-alignments' move-counting, not a bug in
    voidmass_terms). Partial ablation avoids the lump: the aligner
    names each leaf individually once at least one sibling is
    observed, so movecount properly tracks activity count again.
    '''

    def setUp(self):
        small_acts = [Activity(None, f's{i}', ACT_COST) for i in range(2)]
        for i, act in enumerate(small_acts):
            act.id = f'small{i}'
        self.small = Sequence(None, small_acts)
        self.small.id = 'small'
        for act in small_acts:
            act.set_parent(self.small)

        big_acts = [Activity(None, f'b{i}', ACT_COST) for i in range(8)]
        for i, act in enumerate(big_acts):
            act.id = f'big{i}'
        self.big = Sequence(None, big_acts)
        self.big.id = 'big'
        for act in big_acts:
            act.set_parent(self.big)

        self.other = Activity(None, 'other', ACT_COST)
        self.other.id = 'other'

        self.tree = Sequence(None, [self.other, self.small, self.big])
        self.tree.id = 'root'
        self.other.set_parent(self.tree)
        self.small.set_parent(self.tree)
        self.big.set_parent(self.tree)

        # 50% ablation rate for both: 1 of 2 small activities observed,
        # 4 of 8 big activities observed
        trace = ('other', 's0', 'b0', 'b1', 'b2', 'b3')
        self.variant_probs = {trace: 1.0}
        self.skip_dict = {
            variant_key(v): align(self.tree, v) for v in self.variant_probs
        }
        self.table = voidmass_table(self.tree, self.skip_dict, self.variant_probs)

    def test_variant1_scale_free_equal_scores(self):
        self.assertAlmostEqual(
            self.table[self.small]['voidmass_subprocess'],
            self.table[self.big]['voidmass_subprocess'], places=6)
        self.assertAlmostEqual(self.table[self.small]['voidmass_subprocess'], 0.5, places=6)

    def test_variant2_size_preserving_big_is_4x_small(self):
        small_v2 = self.table[self.small]['voidmass_process']
        big_v2 = self.table[self.big]['voidmass_process']
        self.assertGreater(small_v2, 0.0)
        self.assertAlmostEqual(big_v2 / small_v2, 4.0, places=6)

    def test_total_ablation_breaks_size_preservation_via_lumping(self):
        '''
        Documents the failure mode found above: with the SAME two
        subprocesses wholly missing (0% observed instead of 50%), the
        aligner lumps each into one Skip move and variant 2 can no
        longer tell them apart by size. Pinned so this known limit
        doesn't silently change if skip-alignments' lumping behaviour
        ever does.
        '''
        variant_probs = {('other',): 1.0}
        skip_dict = {variant_key(v): align(self.tree, v) for v in variant_probs}
        table = voidmass_table(self.tree, skip_dict, variant_probs)
        small_v2 = table[self.small]['voidmass_process']
        big_v2 = table[self.big]['voidmass_process']
        self.assertAlmostEqual(table[self.small]['movecount'], 1.0, places=6)
        self.assertAlmostEqual(table[self.big]['movecount'], 1.0, places=6)
        self.assertAlmostEqual(big_v2 / small_v2, 1.0, places=6)


class VoidageTest(unittest.TestCase):
    '''
    Variants 3/4 (voidage = skip_prob * voidmass) on the running
    example, using real skip_probs from the full ebi-backed
    DerivationPipeline (pvoid.skipprob) - the same pipeline
    coverage_by_alignment already depends on, so this needs no new
    infrastructure, just the same skip_dict/variant_probs as the other
    tests plus this one extra pipeline run.
    '''

    @classmethod
    def setUpClass(cls):
        # setUpClass, not setUp: this fixture includes a real ebi-backed
        # pvoid.skipprob pipeline run, shared by every test in the class -
        # tests must not mutate it.
        cls.tree = build_running_example_tree()
        cls.o, cls.approval, cls.sched, cls.p = cls.tree.children
        cls.a, cls.e = cls.approval.children
        cls.s, cls.tau = cls.sched.children
        cls.variant_probs = {
            ('o', 'a', 's', 'p'): 2 / 6,
            ('o', 'a', 'e', 'a', 's', 'p'): 1 / 6,
            ('o', 's', 'p'): 2 / 6,
            ('o', 'a', 'p'): 1 / 6,
        }
        cls.skip_dict = {
            variant_key(variant): align(cls.tree, variant)
            for variant in cls.variant_probs
        }

        log = build_running_example_log()
        dv = pvoid.skipprob(log, cls.tree, 'var/lab/test_voidage.slpn')
        cls.skip_probs = dv.skip_probs

        cls.table = voidmass_table(cls.tree, cls.skip_dict, cls.variant_probs,
                                    skip_probs=cls.skip_probs)

    def test_voidage_is_skipprob_times_voidmass(self):
        for node in (self.tree, self.approval, self.a, self.e, self.o,
                     self.sched, self.s, self.p):
            with self.subTest(node=node.id):
                row = self.table[node]
                sp = self.skip_probs[node]
                self.assertAlmostEqual(row['voidage_subprocess'],
                                        sp * row['voidmass_subprocess'], places=9)
                self.assertAlmostEqual(row['voidage_process'],
                                        sp * row['voidmass_process'], places=9)

    def test_voidage_additivity_holds_here_because_only_one_child_contributes(self):
        '''
        NOT a general property - a degenerate case. 'e' has zero deficit
        throughout this example (see reference table), so approval's
        entire voidmass_process comes from 'a' alone, and skip_prob(a)
        happens to equal skip_prob(approval) exactly here (both 1/3 -
        plausibly structural: 'a' is the Loop's mandatory do-child, so
        skipping the Loop at all structurally implies skipping 'a' - see
        _mandatorily_implies). With only one nonzero-deficit child,
        additivity survives the multiply by construction regardless of
        that coincidence. See VoidageAdditivityLossTest for a
        non-degenerate case (two siblings both with nonzero deficit and
        different skip_probs), where the brief's claim does hold.
        '''
        approval_v4 = self.table[self.approval]['voidage_process']
        children_sum = (self.table[self.a]['voidage_process']
                         + self.table[self.e]['voidage_process'])
        self.assertAlmostEqual(approval_v4, children_sum, places=6)

    def test_root_skipprob_near_zero_but_voidmass_nonzero(self):
        '''
        E4 from the brief: root skip_prob should be near zero (whole
        traces mostly conform) while root voidmass is clearly nonzero
        (deficits exist locally) - the gap that gives voidmass a
        headline where skipprob alone reads as "everything's fine".
        '''
        root = self.table[self.tree]
        self.assertLess(self.skip_probs[self.tree], 0.05)
        self.assertGreater(root['voidmass_process'], 0.05)


class VoidageAdditivityLossTest(unittest.TestCase):
    '''
    A non-degenerate case for the brief's additivity-loss claim: two
    sibling leaves (x, y) BOTH with nonzero deficit and DIFFERENT
    skip_probs - x missing from 1/5 of traces, y missing from 2/5,
    so skip_prob(x) != skip_prob(y). Model: seq(x, y).
    '''

    @classmethod
    def setUpClass(cls):
        # setUpClass, not setUp - see VoidageTest.setUpClass for why.
        from process_voids import dtlog

        cls.x = Activity(None, 'x', ACT_COST)
        cls.x.id = '1'
        cls.y = Activity(None, 'y', ACT_COST)
        cls.y.id = '2'
        cls.tree = Sequence(None, [cls.x, cls.y])
        cls.tree.id = '3'
        cls.x.set_parent(cls.tree)
        cls.y.set_parent(cls.tree)

        traces = (['x:0 y:1'] * 2) + (['y:0'] * 1) + (['x:0'] * 2)
        names = [f'c{i}' for i in range(len(traces))]
        log = dtlog.convert_timed(*traces, names=names, time_unit='hours')

        cls.variant_probs = {
            ('x', 'y'): 2 / 5,
            ('y',): 1 / 5,
            ('x',): 2 / 5,
        }
        cls.skip_dict = {
            variant_key(variant): align(cls.tree, variant)
            for variant in cls.variant_probs
        }
        dv = pvoid.skipprob(log, cls.tree, 'var/lab/test_voidage_additivity.slpn')
        cls.skip_probs = dv.skip_probs
        cls.table = voidmass_table(cls.tree, cls.skip_dict, cls.variant_probs,
                                    skip_probs=cls.skip_probs)

    def test_skip_probs_actually_differ(self):
        self.assertNotAlmostEqual(self.skip_probs[self.x], self.skip_probs[self.y], places=6)

    def test_both_children_have_nonzero_deficit(self):
        self.assertGreater(self.table[self.x]['deficit'], 0.0)
        self.assertGreater(self.table[self.y]['deficit'], 0.0)

    def test_voidage_additivity_is_lost(self):
        root_v4 = self.table[self.tree]['voidage_process']
        children_sum = (self.table[self.x]['voidage_process']
                         + self.table[self.y]['voidage_process'])
        self.assertNotAlmostEqual(root_v4, children_sum, places=6)


if __name__ == '__main__':
    unittest.main()
