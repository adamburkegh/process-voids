'''
coveragemass.executions against Definition [Executions]: an execution of
pt is made only of moves whose model element lies in pt's own subtree.

Skip-alignments' normal form lumps an entirely unwitnessed subtree into
one Skip on its coarsest node. A lump on a strict ancestor of pt is not
a move in pt's subtree, so it is no execution of pt: pt was neither
traversed nor skipped in that alignment, and the ancestor alone carries
the void.

The skip probabilities these executions are multiplied with follow the
same rule: EngineWholeTreeLumpTest pins the engine on WholeTreeLumpTest's
model, and test_masked_skip_probs.py checks the two agree node by node
on lumps below the root.
'''

import tempfile
import unittest
from pathlib import Path

from skipalignments import Activity, Sequence, Skip, Aligner

from lab.fixtures import build_running_example_tree
from process_voids import dtlog, pvoid
from process_voids.coveragemass import executions

ACT_COST = 100000

Aligner.set_level_incentive(0)


def leaf(name, node_id):
    node = Activity(None, name, ACT_COST)
    node.id = node_id
    return node


def sequence(node_id, children):
    node = Sequence(None, children)
    node.id = node_id
    for child in children:
        child.set_parent(node)
    return node


def nested_sequence():
    '''seq(a, seq(b, c)), as (root, a, seq(b, c), b, c).'''
    a, b, c = leaf('a', '1'), leaf('b', '2'), leaf('c', '3')
    bc = sequence('4', [b, c])
    root = sequence('5', [a, bc])
    return root, a, bc, b, c


def only_alignment(tree, trace):
    states, _ = Aligner(tree).align_normal_form(list(trace), [ACT_COST] * len(trace), True,
                                                timeout=100)
    assert len(states) == 1, f'expected one optimal alignment, got {len(states)}'
    return states[0].path


def lumps(path):
    return [model_elem for _log_elem, model_elem in path if isinstance(model_elem, Skip)]


class WholeTreeLumpTest(unittest.TestCase):
    '''
    Model seq(a, seq(b, c)), trace <z>: z is a log move and the whole
    tree goes unwitnessed, lumped into one Skip on the root. Only the
    root has an execution - not b or c, even though every traversal of
    the root traverses them.
    '''

    def setUp(self):
        self.root, self.a, self.bc, self.b, self.c = nested_sequence()
        self.path = only_alignment(self.root, ('z',))

    def test_alignment_is_one_log_move_and_one_root_lump(self):
        self.assertEqual([log_elem for log_elem, model_elem in self.path if model_elem == '>>'],
                         ['z'])
        self.assertEqual([lump.node for lump in lumps(self.path)], [self.root])
        self.assertEqual(len(self.path), 2)

    def test_root_has_the_lump_as_its_one_execution(self):
        execs = executions(self.path, self.root)
        self.assertEqual(len(execs), 1)
        self.assertEqual([model_elem for _log_elem, model_elem in execs[0]], lumps(self.path))

    def test_descendants_have_no_execution(self):
        for name, node in (('a', self.a), ('seq(b,c)', self.bc), ('b', self.b), ('c', self.c)):
            with self.subTest(node=name):
                self.assertEqual(executions(self.path, node), [])


class LumpedLoopTest(unittest.TestCase):
    '''
    The running example seq(o, loop(a, e), xor(s, tau), p) on <o, s, p>:
    the approval loop goes unwitnessed, lumped into one Skip on the loop.
    Neither child has an execution - not the do-child a, which every
    traversal of the loop traverses, and not the redo-child e.
    '''

    def setUp(self):
        self.tree = build_running_example_tree()
        _o, self.approval, _sched, _p = self.tree.children
        self.a, self.e = self.approval.children
        self.path = only_alignment(self.tree, ('o', 's', 'p'))

    def test_alignment_lumps_the_loop(self):
        self.assertEqual([lump.node for lump in lumps(self.path)], [self.approval])

    def test_loop_has_the_lump_as_its_one_execution(self):
        execs = executions(self.path, self.approval)
        self.assertEqual(len(execs), 1)
        self.assertEqual([model_elem for _log_elem, model_elem in execs[0]], lumps(self.path))

    def test_children_have_no_execution(self):
        for name, node in (('a', self.a), ('e', self.e)):
            with self.subTest(node=name):
                self.assertEqual(executions(self.path, node), [])


class EngineWholeTreeLumpTest(unittest.TestCase):
    '''
    The engine's skip_probs on WholeTreeLumpTest's model, over the log
    [<a, b, c>, <z>]. Skip probability conditions on traversal: b is
    executed only in <a, b, c>, where it is present, so skip_prob(b) = 0.
    The root is skipped in <z>, so skip_prob(root) = 1/2. An engine that
    propagated the root's skip to its descendants would give 1/2 for
    every node.
    '''

    @classmethod
    def setUpClass(cls):
        cls.root, cls.a, cls.bc, cls.b, cls.c = nested_sequence()
        log = dtlog.convert_timed('a:0 b:1 c:2', 'z:0', names=['c0', 'c1'], time_unit='hours')
        # pvoid.skipprob's SLPN file is only an intermediate
        with tempfile.TemporaryDirectory() as tmp:
            cls.dv = pvoid.skipprob(log, cls.root, str(Path(tmp) / 'whole_tree.slpn'))

    def test_engine_alignment_lumps_the_root(self):
        self.assertEqual([[lump.node for lump in lumps(state.path)]
                          for state in self.dv.skip_dict_backup['z']],
                         [[self.root]])

    def test_root_skipped_in_half_the_traces(self):
        self.assertAlmostEqual(self.dv.skip_probs[self.root], 0.5, places=6)

    def test_descendants_never_skipped_where_executed(self):
        for name, node in (('a', self.a), ('seq(b,c)', self.bc), ('b', self.b), ('c', self.c)):
            with self.subTest(node=name):
                self.assertAlmostEqual(self.dv.skip_probs[node], 0.0, places=6)


if __name__ == '__main__':
    unittest.main()
