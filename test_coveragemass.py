
import sys
import unittest

from processtree import *
from coveragemass import *
from slpn import *



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

    def test_coverage_mass(self):
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
        self.assertEqual( coverage_mass( a, skip_probs), 0.9 )
        # choice = 0.1* 2/3 + 1.0 * 1/3 ~= 0.399
        self.assertAlmostEqual( coverage_mass(choice, skip_probs ), 0.399, 
                                delta = 0.002 )
        # seq    = (0.9 + 0.399) / 2
        self.assertAlmostEqual( coverage_mass(tree, skip_probs ), 0.6495,
                                delta = 0.002)

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


