'''
Hand-derived regression tests for process_voids.voidsalign2 (Definitions
[Coverage by Skip-Weighted Alignment Correspondence] / [Void by
Skip-Weighted Alignment Correspondence]), following
test_voidsalign.py's fixture pattern: the real (ebi-free) skip-aligner
builds each skip_dict, and expected values are hand-derived from the
alignment it is known to find. skip_probs is supplied directly, so what
is under test is the mass arithmetic rather than skipprob itself.

The definition's mass is conditioned on observation at both levels - an
execution with no synchronous move is excluded, and a trace whose
alignments hold no such execution is dropped from the average and from
the W that normalises it - and each observing alignment carries
1/|Upsilon_sigma| of its trace's weight, not 1/|O_sigma|. That is
coveragemass.observed_alignment_mass's arithmetic exactly; what differs
here is the ratio inside it, smatchcount/smovecount, which weights a
lumped skip by the size of the subprocess it stands for.
'''

import unittest
from types import SimpleNamespace

from skipalignments import Activity, Aligner, Sequence, Skip, Tau, Xor

from process_voids.coveragemass import make_executions_cache
from process_voids.voidsalign2 import (
    coversalign2, observed_skip_weighted_mass, voidsalign2,
)

ACT_COST = 100000

Aligner.set_level_incentive(0)


def leaf(cls, name, node_id, cost=ACT_COST):
    node = cls(None, name, cost)
    node.id = node_id
    return node


def variant_key(trace):
    return ', '.join(trace)


def align(tree, trace):
    states, _ = Aligner(tree).align_normal_form(
        list(trace), [ACT_COST] * len(trace), True, timeout=100)
    return states


def build(tree, traces_with_weights):
    skip_dict = {}
    variant_probs = {}
    for trace, weight in traces_with_weights.items():
        skip_dict[variant_key(trace)] = align(tree, trace)
        variant_probs[trace] = weight
    return skip_dict, variant_probs


def seq_abc():
    """seq(a, seq(b, c)) - a two-activity subprocess that lumps."""
    a = leaf(Activity, 'a', '1')
    b = leaf(Activity, 'b', '2')
    c = leaf(Activity, 'c', '3')
    bc = Sequence(None, [b, c])
    bc.id = '4'
    b.set_parent(bc)
    c.set_parent(bc)
    tree = Sequence(None, [a, bc])
    tree.id = '5'
    a.set_parent(tree)
    bc.set_parent(tree)
    return tree, a, bc, b, c


class ObservedMassTest(unittest.TestCase):

    def test_fully_recorded_subprocess_has_mass_one(self):
        tree, a, bc, b, c = seq_abc()
        skip_dict, variant_probs = build(tree, {('a', 'b', 'c'): 1.0})
        self.assertAlmostEqual(
            observed_skip_weighted_mass(bc, skip_dict, variant_probs), 1.0, places=9)

    def test_a_wholly_skipped_subprocess_has_no_observed_mass(self):
        '''
        W = 0: no alignment holds an execution of bc with a synchronous
        move, so the mass is 0 by the definition's own W = 0 case rather
        than by averaging zeros.
        '''
        tree, a, bc, b, c = seq_abc()
        skip_dict, variant_probs = build(tree, {('a',): 1.0})
        self.assertEqual(
            observed_skip_weighted_mass(bc, skip_dict, variant_probs), 0.0)

    def test_traces_with_no_observed_execution_are_dropped_not_scored_zero(self):
        '''
        Half the log skips bc entirely. Those traces leave the average
        AND W, so the mass stays at what the observing half shows (1.0)
        rather than falling to 0.5 - the absence is carried by the
        (1 - skipprob) factor instead, which is what keeps coverage
        linear in absence rather than quadratic.
        '''
        tree, a, bc, b, c = seq_abc()
        skip_dict, variant_probs = build(tree, {('a', 'b', 'c'): 0.5, ('a',): 0.5})
        self.assertAlmostEqual(
            observed_skip_weighted_mass(bc, skip_dict, variant_probs), 1.0, places=9)

    def test_a_partially_recorded_subprocess_is_weighted_by_skip_size(self):
        '''
        <a, b> records b but not c. bc's execution is then one sync move
        (b) plus a skip of c: smatchcount 1, smovecount 1 + 1 = 2, so the
        ratio is 1/2 - the same as an unweighted count here, since the
        skipped part is a single activity.
        '''
        tree, a, bc, b, c = seq_abc()
        skip_dict, variant_probs = build(tree, {('a', 'b'): 1.0})
        self.assertAlmostEqual(
            observed_skip_weighted_mass(bc, skip_dict, variant_probs), 0.5, places=9)

    def test_a_lumped_skip_counts_for_the_size_it_stands_for(self):
        '''
        seq(a, seq(b, c, d)) against <a, b>: the execution of bcd holds
        one sync move and a lumped skip of seq(c, d), which smovecount
        weights 2, not 1. So the ratio is 1/3, where an unweighted
        movecount would read 1/2.
        '''
        a = leaf(Activity, 'a', '1')
        b = leaf(Activity, 'b', '2')
        c = leaf(Activity, 'c', '3')
        d = leaf(Activity, 'd', '4')
        cd = Sequence(None, [c, d])
        cd.id = '5'
        c.set_parent(cd)
        d.set_parent(cd)
        bcd = Sequence(None, [b, cd])
        bcd.id = '6'
        b.set_parent(bcd)
        cd.set_parent(bcd)
        tree = Sequence(None, [a, bcd])
        tree.id = '7'
        a.set_parent(tree)
        bcd.set_parent(tree)

        skip_dict, variant_probs = build(tree, {('a', 'b'): 1.0})
        self.assertAlmostEqual(
            observed_skip_weighted_mass(bcd, skip_dict, variant_probs), 1 / 3, places=9)

    def test_a_silent_skip_adds_no_weight(self):
        '''
        aligncost(empty, msub) of a wholly silent subprocess is 0, so a
        skip over it contributes nothing to smovecount and needs no
        separate silent-move exclusion.
        '''
        a = leaf(Activity, 'a', '1')
        tau = leaf(Tau, 'skip-it', '2', cost=0)
        choice = Xor(None, [a, tau])
        choice.id = '3'
        a.set_parent(choice)
        tau.set_parent(choice)
        b = leaf(Activity, 'b', '4')
        tree = Sequence(None, [choice, b])
        tree.id = '5'
        choice.set_parent(tree)
        b.set_parent(tree)

        skip_dict, variant_probs = build(tree, {('a', 'b'): 1.0})
        self.assertAlmostEqual(
            observed_skip_weighted_mass(tree, skip_dict, variant_probs), 1.0, places=9)

    def test_shared_executions_cache_gives_the_same_answer(self):
        tree, a, bc, b, c = seq_abc()
        skip_dict, variant_probs = build(tree, {('a', 'b'): 1.0})
        cache = make_executions_cache(tree)
        self.assertAlmostEqual(
            observed_skip_weighted_mass(bc, skip_dict, variant_probs,
                                        executions_cache=cache),
            observed_skip_weighted_mass(bc, skip_dict, variant_probs), places=9)


class AveragingTest(unittest.TestCase):
    '''
    The two averages the real aligner's own output does not exercise on
    these small fixtures: 1/|Upsilon_sigma| where a variant has tied
    optimal alignments, and 1/|P_sigma,delta| where one alignment holds
    several executions of the same node. Alignments are hand-built here
    (a path is a list of (log_elem, model_elem) pairs, the shape
    voidmass_pn already constructs for its own skip_dict) so the tie and
    the repeat are controlled rather than hoped for.
    '''

    def setUp(self):
        self.b = leaf(Activity, 'b', '1')
        self.c = leaf(Activity, 'c', '2')
        self.bc = Sequence(None, [self.b, self.c])
        self.bc.id = '3'
        self.b.set_parent(self.bc)
        self.c.set_parent(self.bc)
        self.d = leaf(Activity, 'd', '4')
        self.tree = Sequence(None, [self.bc, self.d])
        self.tree.id = '5'
        self.bc.set_parent(self.tree)
        self.d.set_parent(self.tree)

    def alignment(self, *moves):
        return SimpleNamespace(path=list(moves))

    def sync(self, node):
        return (node.name, node)

    def skip(self, node):
        return ('>>', Skip(node, ACT_COST))

    def test_tied_alignments_carry_their_share_of_the_trace_not_of_the_observers(self):
        '''
        Variant A has two tied optimal alignments and only one observes
        bc; variant B has a single alignment that observes it fully. The
        observing alignment of A carries 1/|Upsilon| = 1/2 of A's weight,
        so the mass is 0.625/0.75 = 5/6. Sharing over |O_sigma| instead
        would give A's observation full weight and read 0.75, letting
        the unobserving alignment vanish rather than dilute.
        '''
        observing = self.alignment(self.sync(self.b), self.skip(self.c),
                                    self.sync(self.d))
        not_observing = self.alignment(self.skip(self.bc), self.sync(self.d))
        skip_dict = {'a-variant': [observing, not_observing],
                     'b-variant': [self.alignment(self.sync(self.b), self.sync(self.c),
                                                   self.sync(self.d))]}
        variant_probs = {('a-variant',): 0.5, ('b-variant',): 0.5}

        self.assertAlmostEqual(
            observed_skip_weighted_mass(self.bc, skip_dict, variant_probs),
            5 / 6, places=9)

    def test_several_executions_in_one_alignment_are_averaged_not_pooled(self):
        '''
        One alignment holds two executions of bc - moves in two runs
        separated by d - with ratios 1/2 and 1. The definition averages
        per execution (3/4), where pooling the counts across both would
        give 2/3.
        '''
        two_runs = self.alignment(
            self.sync(self.b), self.skip(self.c),
            self.sync(self.d),
            self.sync(self.b), self.sync(self.c))
        skip_dict = {'v': [two_runs]}
        variant_probs = {('v',): 1.0}

        self.assertAlmostEqual(
            observed_skip_weighted_mass(self.bc, skip_dict, variant_probs),
            0.75, places=9)

    def test_an_execution_with_no_synchronous_move_leaves_the_average(self):
        '''
        Same shape, but the first run is a bare skip of c with no sync
        move of its own - it is not an observation of bc, so it leaves P
        and the mass is the surviving execution's 1.0 rather than 1/2.
        '''
        two_runs = self.alignment(
            self.skip(self.c),
            self.sync(self.d),
            self.sync(self.b), self.sync(self.c))
        skip_dict = {'v': [two_runs]}
        variant_probs = {('v',): 1.0}

        self.assertAlmostEqual(
            observed_skip_weighted_mass(self.bc, skip_dict, variant_probs),
            1.0, places=9)


class CoverageAndVoidTest(unittest.TestCase):

    def test_coverage_is_one_minus_skipprob_times_mass(self):
        tree, a, bc, b, c = seq_abc()
        skip_dict, variant_probs = build(tree, {('a', 'b'): 1.0})
        skip_probs = {node: 0.25 for node in (tree, a, bc, b, c)}
        mass = observed_skip_weighted_mass(bc, skip_dict, variant_probs)
        self.assertAlmostEqual(
            coversalign2(bc, skip_dict, variant_probs, skip_probs),
            0.75 * mass, places=9)

    def test_void_is_one_minus_coverage(self):
        tree, a, bc, b, c = seq_abc()
        skip_dict, variant_probs = build(tree, {('a', 'b'): 1.0})
        skip_probs = {node: 0.25 for node in (tree, a, bc, b, c)}
        self.assertAlmostEqual(
            voidsalign2(bc, skip_dict, variant_probs, skip_probs),
            1 - coversalign2(bc, skip_dict, variant_probs, skip_probs), places=9)

    def test_never_missing_reads_zero_void(self):
        tree, a, bc, b, c = seq_abc()
        skip_dict, variant_probs = build(tree, {('a', 'b', 'c'): 1.0})
        skip_probs = {node: 0.0 for node in (tree, a, bc, b, c)}
        self.assertAlmostEqual(
            voidsalign2(bc, skip_dict, variant_probs, skip_probs), 0.0, places=9)

    def test_always_missing_reads_one_void(self):
        '''
        The case the superseded form could not report: a wholly missing
        subprocess. W is 0, so coverage is 0 and void is 1, rather than
        the inverted U that returned to zero at total absence.
        '''
        tree, a, bc, b, c = seq_abc()
        skip_dict, variant_probs = build(tree, {('a',): 1.0})
        skip_probs = {tree: 0.0, a: 0.0, bc: 1.0, b: 1.0, c: 1.0}
        self.assertAlmostEqual(
            voidsalign2(bc, skip_dict, variant_probs, skip_probs), 1.0, places=9)

    def test_missing_from_half_the_traces_reads_half_void(self):
        '''
        Recorded completely where recorded at all, absent from half the
        log: mass 1, skipprob 1/2, so coverage is 1/2 and void is 1/2 -
        linear in absence, the property the conditioning buys.
        '''
        tree, a, bc, b, c = seq_abc()
        skip_dict, variant_probs = build(tree, {('a', 'b', 'c'): 0.5, ('a',): 0.5})
        skip_probs = {tree: 0.0, a: 0.0, bc: 0.5, b: 0.5, c: 0.5}
        self.assertAlmostEqual(
            voidsalign2(bc, skip_dict, variant_probs, skip_probs), 0.5, places=9)


class RunningExampleTest(unittest.TestCase):
    '''
    Adam's hand-worked table for the payment running example
    (lab.fixtures), through the full ebi-backed pipeline rather than
    hand-picked alignments - skip probabilities included, so this pins
    the metric end to end rather than the mass arithmetic alone.
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
                                'var/lab/test_voidsalign2.slpn')

    def _row(self, node):
        mass = observed_skip_weighted_mass(node, self.dv.skip_dict_backup, self.dv.pl)
        cover = coversalign2(node, self.dv.skip_dict_backup, self.dv.pl, self.dv.skip_probs)
        void = voidsalign2(node, self.dv.skip_dict_backup, self.dv.pl, self.dv.skip_probs)
        return self.dv.skip_probs[node], mass, cover, void

    def test_reference_table(self):
        cases = [
            ('N (root)', 'tree', 0.0000, 0.9167, 0.9167, 0.0833),
            ('o', 'o', 0.0000, 1.0000, 1.0000, 0.0000),
            ('loop', 'approval', 0.3333, 1.0000, 0.6667, 0.3333),
            ('a', 'a', 0.0000, 1.0000, 1.0000, 0.0000),
            ('e', 'e', 0.0000, 1.0000, 1.0000, 0.0000),
            # The choice's skip probability is 1/6, not 0: sigma6
            # (<o, a, p>) omits s, and the aligner takes the silent
            # branch, which it records as a skip of the Xor itself - so
            # the mass lands on the Xor while s and the tau leaf read 0
            # (see tests/process_voids/test_masked_skip_probs.py).
            ('xor', 'sched', 1 / 6, 1.0000, 5 / 6, 1 / 6),
            ('s', 's', 0.0000, 1.0000, 1.0000, 0.0000),
            ('p', 'p', 0.0000, 1.0000, 1.0000, 0.0000),
        ]
        for name, attr, skipprob, mass, cover, void in cases:
            with self.subTest(node=name):
                node = getattr(self, attr)
                actual = self._row(node)
                for label, expected, got in zip(
                        ('skipprob', 'mass', 'coversalign', 'voidsalign'),
                        (skipprob, mass, cover, void), actual):
                    self.assertAlmostEqual(got, expected, places=4,
                                           msg=f'{name}: {label}')


if __name__ == '__main__':
    unittest.main()
