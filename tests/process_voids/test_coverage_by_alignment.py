
import unittest

from skipalignments import Activity, Tau, Sequence, Xor, And, Loop, Aligner

from lab.fixtures import build_running_example_tree
from process_voids.coveragemass import alignment_mass

ACT_COST = 100000

Aligner.set_level_incentive(0)


def leaf(cls, name, node_id, cost=ACT_COST):
    node = cls(None, name, cost)
    node.id = node_id
    return node


def variant_key(trace):
    return ', '.join(trace)


def align(tree, trace):
    states, _ = Aligner(tree).align2(list(trace), [ACT_COST] * len(trace), True, timeout=100)
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
    This is why: skip_prob(a)/skip_prob(b) - not alignment_mass - is what
    correctly reports these as skipped once the ancestor-lump-skip
    (1-skip_prob) factor is applied; the two concerns are cleanly split
    by the metric's own (1-skip_prob)*mass factorisation.
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

    a is the Loop's mandatory do-child, itself under a mandatory Sequence
    position - so when the whole loop goes unwitnessed (<o,s,p>), that
    lump skip unambiguously implicates a too, and a's mass should NOT
    change under renormalisation. e (the Loop's optional redo-child) and
    s (an Xor branch, inherently ambiguous like M8) both should change,
    since neither position is unambiguously implicated by an ancestor
    lump skip.
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
            ('a', self.a, 2 / 3, 2 / 3),
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


class HeldOutTwoLevelMandatoryChainTest(unittest.TestCase):
    """
    Held-out prediction, not from the original spec: model
    seq(o, loop(seq(x,y), z)), log <o> x1 (whole loop unwitnessed, one
    lump Skip(Loop)). x and y sit two mandatory levels deep beneath the
    lump (Sequence child, then Loop do-child) and should inherit it
    (0.0 both conventions, not vacuous). z is the loop's redo-child -
    genuinely ambiguous - and should stay vacuous (0.0 zero, 1.0
    renormalised). Predicted by hand from _mandatorily_implies before
    running this test.
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

    def test_predicted_values(self):
        cases = [
            ('loop', self.loop, 0.0, 0.0),
            ('seq(x,y)', self.do, 0.0, 0.0),
            ('x', self.x, 0.0, 0.0),
            ('y', self.y, 0.0, 0.0),
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
