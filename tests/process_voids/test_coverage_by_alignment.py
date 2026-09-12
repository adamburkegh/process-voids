
import unittest
from types import SimpleNamespace

from skipalignments import Activity, Tau, Sequence, Xor, And, Loop, Aligner, Skip

from lab.fixtures import build_running_example_tree
from process_voids.coveragemass import (
    alignment_mass, make_executions_cache, observed_alignment_mass,
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


class PerfectFitTest(unittest.TestCase):
    """M1: model seq(a,b), log <a,b> x1."""

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.skip_dict, self.variant_probs = build(self.tree, {('a', 'b'): 1.0})

    def test_masses(self):
        cases = [(self.tree, 1.0), (self.a, 1.0), (self.b, 1.0)]
        for node, expected in cases:
            with self.subTest(node=node):
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs),
                    expected, places=6)


class OneModelMoveTest(unittest.TestCase):
    """M2: model seq(a,b), log <a> x1. Alignment: sync a, model move b."""

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.skip_dict, self.variant_probs = build(self.tree, {('a',): 1.0})

    def test_masses(self):
        cases = [(self.tree, 0.5), (self.a, 1.0), (self.b, 0.0)]
        for node, expected in cases:
            with self.subTest(node=node):
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs),
                    expected, places=6)


class RatiosDoNotComposeTest(unittest.TestCase):
    """M3: model seq(a,b,c), log <a> x1. Parent ratio != mean of children."""

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.c = leaf(Activity, 'c', '3')
        self.tree = Sequence(None, [self.a, self.b, self.c])
        self.tree.id = '4'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.c.set_parent(self.tree)
        self.skip_dict, self.variant_probs = build(self.tree, {('a',): 1.0})

    def test_masses(self):
        cases = [(self.tree, 1 / 3), (self.a, 1.0), (self.b, 0.0), (self.c, 0.0)]
        for node, expected in cases:
            with self.subTest(node=node):
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs),
                    expected, places=6)


class SilentMovesExcludedTest(unittest.TestCase):
    """
    M4: model seq(a, xor(b,tau)), log <a> x1. Alignment: sync a, silent
    move. Root movecount must be 1 (not 2) - the primary regression test
    for silent handling.
    """

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

    def test_root_and_a(self):
        cases = [(self.tree, 1.0), (self.a, 1.0)]
        for node, expected in cases:
            with self.subTest(node=node):
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs),
                    expected, places=6)

    def test_choice_convention_split(self):
        self.assertAlmostEqual(
            alignment_mass(self.choice, self.skip_dict, self.variant_probs, 'zero'),
            0.0, places=6)
        self.assertAlmostEqual(
            alignment_mass(self.choice, self.skip_dict, self.variant_probs, 'renormalised'),
            1.0, places=6)


class FilterVsZeroConventionTest(unittest.TestCase):
    """M5: model seq(a, xor(b,tau)), log <a,b> x1, <a> x1."""

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
        self.skip_dict, self.variant_probs = build(
            self.tree, {('a', 'b'): 0.5, ('a',): 0.5})

    def test_choice_and_b_convention_split(self):
        for node in (self.choice, self.b):
            with self.subTest(node=node):
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs, 'zero'),
                    0.5, places=6)
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs, 'renormalised'),
                    1.0, places=6)


class MultipleExecutionsAveragedTest(unittest.TestCase):
    """M6: model seq(o, loop(a,b)), log <o,a,b,a> x1."""

    def setUp(self):
        self.o = leaf(Activity, 'o', '1')
        self.a = leaf(Activity, 'a', '2')
        self.b = leaf(Activity, 'b', '3')
        self.loop = Loop(None, [self.a, self.b])
        self.loop.id = '4'
        self.a.set_parent(self.loop)
        self.b.set_parent(self.loop)
        self.tree = Sequence(None, [self.o, self.loop])
        self.tree.id = '5'
        self.o.set_parent(self.tree)
        self.loop.set_parent(self.tree)
        self.skip_dict, self.variant_probs = build(
            self.tree, {('o', 'a', 'b', 'a'): 1.0})

    def test_masses(self):
        cases = [(self.loop, 1.0), (self.a, 1.0), (self.b, 1.0)]
        for node, expected in cases:
            with self.subTest(node=node):
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs),
                    expected, places=6)


class RarityUnderRenormalisationTest(unittest.TestCase):
    """
    M7: model seq(o, loop(a,b)), log <o,a> x9, <o,a,b,a> x1. b is rare but
    perfectly corroborated when it occurs; renormalisation must not
    conflate rarity with poor recording.
    """

    def setUp(self):
        self.o = leaf(Activity, 'o', '1')
        self.a = leaf(Activity, 'a', '2')
        self.b = leaf(Activity, 'b', '3')
        self.loop = Loop(None, [self.a, self.b])
        self.loop.id = '4'
        self.a.set_parent(self.loop)
        self.b.set_parent(self.loop)
        self.tree = Sequence(None, [self.o, self.loop])
        self.tree.id = '5'
        self.o.set_parent(self.tree)
        self.loop.set_parent(self.tree)
        self.skip_dict, self.variant_probs = build(self.tree, {
            ('o', 'a'): 9 / 10,
            ('o', 'a', 'b', 'a'): 1 / 10,
        })

    def test_b(self):
        self.assertAlmostEqual(
            alignment_mass(self.b, self.skip_dict, self.variant_probs, 'zero'),
            0.1, places=6)
        self.assertAlmostEqual(
            alignment_mass(self.b, self.skip_dict, self.variant_probs, 'renormalised'),
            1.0, places=6)


class TiesAcrossOptimalAlignmentsTest(unittest.TestCase):
    """
    M8: model seq(o, xor(a,b)), log <o> x1. Spec envisioned two tied
    optimal alignments (model move a; model move b). The real aligner
    instead collapses a fully-skipped Xor into ONE composite-level block
    entry (confirmed empirically - see chat), so there is no per-child
    tie to find: n_states == 1, and neither a nor b is named individually
    in the path. choice still gets ratio 0 (one valid, non-vacuous
    execution: the block Skip). a/b individually are fully vacuous
    (never named at all), which zero-convention still reports as 0
    (matching the original spec) but renormalised reports as 1 (the
    "never exercised, default to fully covered" convention - NOT the
    spec's 0, since that assumed a decomposition the aligner doesn't do).
    The skip is the Xor's own, reported by skip_prob(choice): a and b
    have no execution in this alignment, so this trace gives neither
    alignment_mass nor skip_prob anything to say about them.
    """

    def setUp(self):
        self.o = leaf(Activity, 'o', '1')
        self.a = leaf(Activity, 'a', '2')
        self.b = leaf(Activity, 'b', '3')
        self.choice = Xor(None, [self.a, self.b])
        self.choice.id = '4'
        self.a.set_parent(self.choice)
        self.b.set_parent(self.choice)
        self.tree = Sequence(None, [self.o, self.choice])
        self.tree.id = '5'
        self.o.set_parent(self.tree)
        self.choice.set_parent(self.tree)
        self.skip_dict, self.variant_probs = build(self.tree, {('o',): 1.0})

    def test_only_one_alignment_found(self):
        self.assertEqual(len(self.skip_dict[variant_key(('o',))]), 1)

    def test_choice_mass_is_zero(self):
        for convention in ('zero', 'renormalised'):
            with self.subTest(convention=convention):
                self.assertAlmostEqual(
                    alignment_mass(self.choice, self.skip_dict, self.variant_probs,
                                    convention),
                    0.0, places=6)

    def test_children_are_vacuous_not_tied(self):
        for node in (self.a, self.b):
            with self.subTest(node=node):
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs, 'zero'),
                    0.0, places=6)
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs,
                                    'renormalised'),
                    1.0, places=6)


class LogMovesExcludedTest(unittest.TestCase):
    """M9: model seq(a,b), log <a,z,b> x1 with z not in act(model)."""

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.skip_dict, self.variant_probs = build(
            self.tree, {('a', 'z', 'b'): 1.0})

    def test_root_mass_unaffected_by_log_move(self):
        self.assertAlmostEqual(
            alignment_mass(self.tree, self.skip_dict, self.variant_probs),
            1.0, places=6)


class RunningExampleTest(unittest.TestCase):
    """
    Reference values for the payment-approval running example:
    N = seq(o, loop(a,e), xor(s,tau), p), over its four distinct trace
    variants (weights derived from the six-trace fixture: sigma1/5 share
    <o,a,s,p>, sigma3/4 share <o,s,p>).

    When the whole loop goes unwitnessed (<o,s,p>), the normal form lumps
    it into one Skip on the loop - approval's execution, and neither
    child's. So a (the Loop's do-child) has no execution in <o,s,p>, e
    (its redo-child) none outside <o,a,e,a,s,p>, and s (an Xor branch,
    like M8) none in <o,a,p>: each reads differently under the two
    conventions.
    """

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

    def test_reference_values(self):
        cases = [
            ('N', self.tree, 11 / 12, 11 / 12),
            ('o', self.o, 1.0, 1.0),
            ('approval', self.approval, 2 / 3, 2 / 3),
            ('a', self.a, 2 / 3, 1.0),
            ('e', self.e, 1 / 6, 1.0),
            ('sched', self.sched, 5 / 6, 1.0),
            ('s', self.s, 5 / 6, 1.0),
            ('p', self.p, 1.0, 1.0),
        ]
        for name, node, expected_zero, expected_renorm in cases:
            with self.subTest(node=name):
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs, 'zero'),
                    expected_zero, places=4)
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs,
                                    'renormalised'),
                    expected_renorm, places=4)


if __name__ == '__main__':
    unittest.main()


class LumpedLoopDescendantsTest(unittest.TestCase):
    """
    Model seq(o, loop(seq(x,y), z)), log <o> x1: the whole loop goes
    unwitnessed, lumped into one Skip on the loop - the loop's own
    execution. Nothing beneath it has one: not seq(x,y) or x and y,
    which every traversal of the loop traverses, and not the redo-child
    z. All four are vacuous: 0.0 under the zero convention, 1.0
    renormalised.
    """

    def setUp(self):
        self.o = leaf(Activity, 'o', '1')
        self.x = leaf(Activity, 'x', '2')
        self.y = leaf(Activity, 'y', '3')
        self.do = Sequence(None, [self.x, self.y])
        self.do.id = '4'
        self.x.set_parent(self.do)
        self.y.set_parent(self.do)
        self.z = leaf(Activity, 'z', '5')
        self.loop = Loop(None, [self.do, self.z])
        self.loop.id = '6'
        self.do.set_parent(self.loop)
        self.z.set_parent(self.loop)
        self.tree = Sequence(None, [self.o, self.loop])
        self.tree.id = '7'
        self.o.set_parent(self.tree)
        self.loop.set_parent(self.tree)
        self.skip_dict, self.variant_probs = build(self.tree, {('o',): 1.0})

    def test_masses(self):
        cases = [
            ('loop', self.loop, 0.0, 0.0),
            ('seq(x,y)', self.do, 0.0, 1.0),
            ('x', self.x, 0.0, 1.0),
            ('y', self.y, 0.0, 1.0),
            ('z', self.z, 0.0, 1.0),
        ]
        for name, node, expected_zero, expected_renorm in cases:
            with self.subTest(node=name):
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs, 'zero'),
                    expected_zero, places=6)
                self.assertAlmostEqual(
                    alignment_mass(node, self.skip_dict, self.variant_probs,
                                    'renormalised'),
                    expected_renorm, places=6)

    def test_shared_cache_gives_identical_results_to_no_cache(self):
        # lab.exp_disco_degrade scores every node in a report row through
        # one shared make_executions_cache - it must not change any number
        cache = make_executions_cache(self.tree)
        for node in (self.tree, self.o, self.loop, self.do, self.x, self.y, self.z):
            for convention in ('zero', 'renormalised'):
                with self.subTest(node=node.id, convention=convention):
                    self.assertAlmostEqual(
                        alignment_mass(node, self.skip_dict, self.variant_probs, convention,
                                       executions_cache=cache),
                        alignment_mass(node, self.skip_dict, self.variant_probs, convention),
                        places=9)


class TimedOutVariantRatioTest(unittest.TestCase):
    """
    M: model seq(a,b), two trace variants each weight 0.5: <a,b> aligns
    normally (perfect fit, ratio 1.0), <a> is a STAND-IN for a variant
    whose alignment search timed out to zero alignments - skip_dict maps
    it to an empty list explicitly (not simply absent from skip_dict;
    voidmass_pn.voidmass_table_pn always inserts a key for every variant
    it iterates, empty or not - see that module).

    timed_out_ratio, when given, treats such a variant as contributing a
    SYNTHETIC per-variant ratio (not real matchcount/movecount data) to
    the SAME weighted-average machinery every other variant uses -
    1.0 = as if it matched perfectly (deficit=0, the "lower voidage /
    upper coverage" bound), 0.0 = as if it matched nothing (the "upper
    voidage / lower coverage" bound). By default (timed_out_ratio=None)
    such a variant is excluded from the weighted sum entirely, the same
    as a variant genuinely absent from skip_dict.
    """

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)

        fit_states = align(self.tree, ('a', 'b'))
        self.skip_dict = {
            variant_key(('a', 'b')): fit_states,
            variant_key(('a',)): [],  # timed out - zero alignments, not absent
        }
        self.variant_probs = {('a', 'b'): 0.5, ('a',): 0.5}

    def test_default_excludes_the_timed_out_variant_entirely(self):
        # Only the fitting variant's weight (0.5) counts at all, and its
        # own ratio is 1.0.
        self.assertAlmostEqual(
            alignment_mass(self.tree, self.skip_dict, self.variant_probs, 'zero'),
            0.5, places=6)

    def test_ratio_one_gives_the_upper_bound(self):
        # Both variants now contribute their full weight at ratio 1.0.
        self.assertAlmostEqual(
            alignment_mass(self.tree, self.skip_dict, self.variant_probs, 'zero',
                            timed_out_ratio=1.0),
            1.0, places=6)

    def test_ratio_zero_gives_the_lower_bound(self):
        # The timed-out variant now contributes its weight at ratio 0.0.
        self.assertAlmostEqual(
            alignment_mass(self.tree, self.skip_dict, self.variant_probs, 'zero',
                            timed_out_ratio=0.0),
            0.5, places=6)

    def test_renormalised_convention_also_honours_timed_out_ratio(self):
        self.assertAlmostEqual(
            alignment_mass(self.tree, self.skip_dict, self.variant_probs, 'renormalised',
                            timed_out_ratio=0.0),
            0.5, places=6)
        self.assertAlmostEqual(
            alignment_mass(self.tree, self.skip_dict, self.variant_probs, 'renormalised',
                            timed_out_ratio=1.0),
            1.0, places=6)


class ObservedAlignmentMassTest(unittest.TestCase):
    """
    defn:move-coverage's conditioned mass: an execution with no
    synchronous move is excluded, and a trace whose alignments hold no
    observed execution drops out of W rather than contributing a zero to
    it. Model seq(a,b) with b recorded in half the traces reads 1.0,
    where the unconditioned mass (what salign_coverage still uses) reads
    0.5.
    """

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)

    def test_half_recorded_submodel_reads_one(self):
        skip_dict, variant_probs = build(self.tree, {('a', 'b'): 0.5, ('a',): 0.5})
        self.assertAlmostEqual(
            observed_alignment_mass(self.b, skip_dict, variant_probs), 1.0, places=6)
        self.assertAlmostEqual(
            alignment_mass(self.b, skip_dict, variant_probs, 'zero'), 0.5, places=6)

    def test_never_recorded_submodel_reads_zero(self):
        skip_dict, variant_probs = build(self.tree, {('a',): 1.0})
        self.assertAlmostEqual(
            observed_alignment_mass(self.b, skip_dict, variant_probs), 0.0, places=6)

    def test_fully_recorded_submodel_reads_one(self):
        skip_dict, variant_probs = build(self.tree, {('a', 'b'): 1.0})
        self.assertAlmostEqual(
            observed_alignment_mass(self.b, skip_dict, variant_probs), 1.0, places=6)

    def test_partly_recorded_execution_keeps_its_ratio(self):
        # At the root, <a> aligns as sync a then skip b: one execution,
        # matchcount 1 of movecount 2. It was observed, so it counts.
        skip_dict, variant_probs = build(self.tree, {('a',): 1.0})
        self.assertAlmostEqual(
            observed_alignment_mass(self.tree, skip_dict, variant_probs), 0.5, places=6)


class ObservedMassWeightsEveryAlignmentTest(unittest.TestCase):
    """
    W weights each observing alignment by 1/|Gamma_sigma|, not 1/|O_sigma|:
    a variant whose tied alignments disagree about whether the submodel
    was observed contributes only the share of its weight the observing
    ones carry. Hand-built paths, since this pins the arithmetic of the
    outer sum rather than anything the aligner decides.

    Variant <a,b> (weight 0.5) gets two tied alignments: one fits (root
    ratio 1), the other blames the log, skipping both activities with no
    synchronous move at all, so the root is unobserved there. Variant <a>
    (weight 0.5) gets one alignment - sync a, skip b - a root ratio of
    0.5, observed.

        mass = (0.5 * 1/2 * 1 + 0.5 * 1 * 0.5) / (0.5 * 1/2 + 0.5 * 1)
             = 0.5 / 0.75 = 2/3

    Weighting the observing alignments by 1/|O_sigma| would give 0.75
    instead, and the unconditioned mass gives 0.5.
    """

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)

        fitting = [('a', self.a), ('b', self.b)]
        unobserved = [('a', '>>'),
                      ('>>', Skip(self.a, self.a.skip_cost)),
                      ('>>', Skip(self.b, self.b.skip_cost))]
        half = [('a', self.a), ('>>', Skip(self.b, self.b.skip_cost))]
        self.skip_dict = {
            variant_key(('a', 'b')): [SimpleNamespace(path=fitting),
                                       SimpleNamespace(path=unobserved)],
            variant_key(('a',)): [SimpleNamespace(path=half)],
        }
        self.variant_probs = {('a', 'b'): 0.5, ('a',): 0.5}

    def test_an_unobserved_tied_alignment_keeps_its_share_out_of_the_average(self):
        self.assertAlmostEqual(
            observed_alignment_mass(self.tree, self.skip_dict, self.variant_probs),
            2 / 3, places=6)


class ObservedMassTimedOutVariantTest(unittest.TestCase):
    """
    A variant whose alignment search timed out (in skip_dict with zero
    alignments) counts as a single synthetic observed unit at its full
    weight when timed_out_ratio is given, so 0.0 and 1.0 bracket the
    value the cell would otherwise have had. With the default None it is
    dropped entirely, which under this mass renormalises the remaining
    weight instead of scoring the variant as a zero.
    """

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.skip_dict = {
            variant_key(('a', 'b')): align(self.tree, ('a', 'b')),
            variant_key(('a',)): [],  # timed out - zero alignments, not absent
        }
        self.variant_probs = {('a', 'b'): 0.5, ('a',): 0.5}

    def test_ratio_zero_is_the_lower_bound(self):
        self.assertAlmostEqual(
            observed_alignment_mass(self.tree, self.skip_dict, self.variant_probs,
                                     timed_out_ratio=0.0),
            0.5, places=6)

    def test_ratio_one_is_the_upper_bound(self):
        self.assertAlmostEqual(
            observed_alignment_mass(self.tree, self.skip_dict, self.variant_probs,
                                     timed_out_ratio=1.0),
            1.0, places=6)

    def test_default_drops_the_variant_and_renormalises(self):
        self.assertAlmostEqual(
            observed_alignment_mass(self.tree, self.skip_dict, self.variant_probs),
            1.0, places=6)
