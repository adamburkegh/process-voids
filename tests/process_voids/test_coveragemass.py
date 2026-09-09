
import sys
import unittest

from skipalignments import *
from process_voids.coveragemass import *
from process_voids.slpn import *



sys.stdout.reconfigure(encoding='utf-8')



ACTIVITY_COST = 10000

def activity(label,aid):
    act = Activity(None, label, ACTIVITY_COST)
    act.id = str(aid)
    return act

def set_parent(nl,parent):
    for node in nl:
        node.set_parent(parent)

class CoverageMassTest(unittest.TestCase):


    def test_update_activity_weights(self):
        a = activity('a',1)
        b = activity('b',2)
        c = activity('c',3)
        #
        choice = Xor(None, [b,c])
        choice.id = '4'
        set_parent( [b,c], choice) 
        #
        seq = Sequence( None, [a,choice] )
        set_parent( [a,choice], seq )
        seq.id = '5'
        #
        tree = seq
        #
        slpn = StochasticLabelledPetriNet()
        slpn.addTransition('1', 3)
        slpn.addTransition('2', 2)
        slpn.addTransition('3', 1)
        #
        update_activity_weights(tree,slpn)
        self.assertEqual( a.weight, 3 )
        self.assertEqual( b.weight, 2 )
        self.assertEqual( c.weight, 1 )


    def test_infer_operator_weights(self):
        a = activity('a',1)
        b = activity('b',2)
        c = activity('c',3)
        #
        choice = Xor(None, [b,c])
        set_parent( [b,c], choice) 
        choice.id = '4'
        #
        seq = Sequence( None, [a,choice] )
        seq.id = '5'
        set_parent( [a,choice], seq )
        #
        tree = seq
        a.weight, b.weight, c.weight = 3, 2, 1
        #
        infer_operator_weights(tree)
        self.assertEqual( choice.weight, 3)
        self.assertEqual( seq.weight, 3)

    def test_mass_by_weight(self):
        a = activity('a',1)
        b = activity('b',2)
        c = activity('c',3)
        #
        choice = Xor(None, [b,c])
        choice.id = '4'
        set_parent( [b,c], choice)
        #
        seq = Sequence( None, [a,choice] )
        seq.id = '5'
        set_parent( [a,choice], seq )
        #
        tree = seq
        a.weight, b.weight, c.weight = 3, 2, 1
        #
        infer_operator_weights(tree)
        skip_probs = { a: 0.1, b: 0.9, c: 0, choice: 0.1, seq: 0.2 }
        self.assertEqual( mass_by_weight( a, skip_probs), 0.9 )
        # choice = 0.1* 2/3 + 1.0 * 1/3 ~= 0.399
        self.assertAlmostEqual( mass_by_weight(choice, skip_probs ), 0.399,
                                delta = 0.002 )
        # seq    = (0.9 + 0.399) / 2
        self.assertAlmostEqual( mass_by_weight(tree, skip_probs ), 0.6495,
                                delta = 0.002)

    def test_voidage_by_weight(self):
        a = activity('a',1)
        b = activity('b',2)
        c = activity('c',3)
        #
        choice = Xor(None, [b,c])
        choice.id = '4'
        set_parent( [b,c], choice)
        #
        seq = Sequence( None, [a,choice] )
        seq.id = '5'
        set_parent( [a,choice], seq )
        #
        tree = seq
        a.weight, b.weight, c.weight = 3, 2, 1
        #
        infer_operator_weights(tree)
        skip_probs = { a: 0.1, b: 0.9, c: 0, choice: 0.1, seq: 0.2 }
        # same fixture as test_mass_by_weight - voidage_by_weight is
        # exactly 1 - mass_by_weight at every node, leaf or not
        for node in (a, choice, tree):
            with self.subTest(node=node):
                self.assertAlmostEqual(
                    voidage_by_weight(node, skip_probs),
                    1 - mass_by_weight(node, skip_probs), delta=1e-9)
        self.assertAlmostEqual(voidage_by_weight(a, skip_probs), 0.1, delta=1e-9)
        self.assertAlmostEqual(voidage_by_weight(choice, skip_probs), 0.601, delta=0.002)
        self.assertAlmostEqual(voidage_by_weight(tree, skip_probs), 0.3505, delta=0.002)

    def test_transfer_pt_weights(self):
        a = activity('a',1)
        b = activity('b',2)
        c = activity('c',3)
        #
        choice = Xor(None, [b,c])
        set_parent( [b,c], choice) 
        #
        seq = Sequence( None, [a,choice] )
        set_parent( [a,choice], seq )
        #
        tree = seq
        #
        slpn = StochasticLabelledPetriNet()
        slpn.addTransition('1', 3)
        slpn.addTransition('2', 2)
        slpn.addTransition('3', 1)
        #
        transfer_pt_weights(tree,slpn)
        self.assertEqual( choice.weight, 3)
        self.assertEqual( seq.weight, 3)


class MandatoryNodeCountTest(unittest.TestCase):
    """
    has_silent_alternative / mandatory_node_count / total_node_count -
    the tree-structural diagnostic for how much of a discovered tree
    the void/coverage metrics can actually speak about (see
    coveragemass.py's Mandatory Node Count section docstring).
    """

    def test_leaf_directly_under_xor_tau_has_a_silent_alternative(self):
        a = activity('a', 1)
        tau = Tau(None, 'tau', 0)
        choice = Xor(None, [tau, a])
        set_parent([tau, a], choice)

        self.assertTrue(has_silent_alternative(a))
        self.assertFalse(has_silent_alternative(choice))  # choice itself has no parent

    def test_xor_without_tau_is_a_real_choice_not_a_silent_alternative(self):
        b = activity('b', 1)
        c = activity('c', 2)
        choice = Xor(None, [b, c])
        set_parent([b, c], choice)

        self.assertFalse(has_silent_alternative(b))
        self.assertFalse(has_silent_alternative(c))

    def test_silent_alternative_further_up_the_ancestor_chain_still_counts(self):
        # Xor(Tau, Sequence(a, b)) - a and b's immediate parent is the
        # Sequence, not the Xor, but the whole Sequence (and therefore
        # both a and b) can be silently skipped via the Xor above it.
        a = activity('a', 1)
        b = activity('b', 2)
        seq = Sequence(None, [a, b])
        set_parent([a, b], seq)
        tau = Tau(None, 'tau', 0)
        choice = Xor(None, [tau, seq])
        set_parent([tau, seq], choice)

        self.assertTrue(has_silent_alternative(a))
        self.assertTrue(has_silent_alternative(b))
        self.assertTrue(has_silent_alternative(seq))

    def test_loop_redo_child_has_a_silent_alternative_even_without_tau(self):
        do = activity('do', 1)
        redo = activity('redo', 2)
        loop = Loop(None, [do, redo])
        set_parent([do, redo], loop)

        self.assertFalse(has_silent_alternative(do))
        self.assertTrue(has_silent_alternative(redo))

    def test_mandatory_and_total_counts_on_the_inductive_degenerate_pattern(self):
        # Sequence(Xor(Tau, a), Xor(Tau, b)) - the classic Inductive
        # Miner noise_threshold=0.0 pattern every leaf gets wrapped in.
        a = activity('a', 1)
        tau1 = Tau(None, 'tau1', 0)
        choice1 = Xor(None, [tau1, a])
        set_parent([tau1, a], choice1)

        b = activity('b', 2)
        tau2 = Tau(None, 'tau2', 0)
        choice2 = Xor(None, [tau2, b])
        set_parent([tau2, b], choice2)

        seq = Sequence(None, [choice1, choice2])
        set_parent([choice1, choice2], seq)

        # Non-Tau nodes: seq, choice1, a, choice2, b = 5 total.
        self.assertEqual(total_node_count(seq), 5)
        # Mandatory: seq (no parent), choice1, choice2 (their own
        # parent isn't an Xor) - neither leaf a nor b is mandatory,
        # since each sits directly under an Xor with a Tau sibling.
        self.assertEqual(mandatory_node_count(seq), 3)

    def test_root_is_always_mandatory(self):
        a = activity('a', 1)
        self.assertFalse(has_silent_alternative(a))  # a has no parent at all here
        self.assertEqual(mandatory_node_count(a), 1)
        self.assertEqual(total_node_count(a), 1)

    def test_tau_nodes_excluded_from_both_counts(self):
        a = activity('a', 1)
        tau = Tau(None, 'tau', 0)
        choice = Xor(None, [tau, a])
        set_parent([tau, a], choice)

        # Non-Tau nodes: choice, a = 2. Mandatory: choice only.
        self.assertEqual(total_node_count(choice), 2)
        self.assertEqual(mandatory_node_count(choice), 1)


