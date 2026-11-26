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

from skipalignments import Activity, And, Loop, Sequence, Skip, Aligner

from lab.fixtures import build_running_example_tree
from process_voids import dtlog, pvoid
from process_voids.coveragemass import executions, matchcount, movecount

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


def operator(cls, node_id, children):
    node = cls(None, children)
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


def shapes(path, node):
    '''(matchcount, movecount) per execution of `node`, in path order.'''
    return [(matchcount(ex), movecount(ex)) for ex in executions(path, node)]


def model_moves(path, node):
    '''How many of `path`'s moves lie in `node`'s subtree at all.'''
    ids = set()

    def walk(n):
        ids.add(id(n))
        for child in getattr(n, 'children', []) or []:
            walk(child)
    walk(node)
    found = 0
    for _log_elem, model_elem in path:
        if model_elem == '>>':
            continue
        subject = getattr(model_elem, 'node', model_elem)
        if id(subject) in ids:
            found += 1
    return found


class LooseSubsequenceTest(unittest.TestCase):
    '''
    Definition [Executions] clause (2): the moves of an execution are a
    LOOSE subsequence of the skip alignment, and only moves that are
    themselves part of a traversal of the node may divide it. A log move
    is part of no traversal, and nor is a parallel sibling's move, so
    neither splits an execution however it is positioned.
    '''

    def test_a_log_move_between_two_moves_does_not_divide_the_traversal(self):
        '''seq(a, b) on <a, q, b>: q is in no model position, so it sits
        between the two moves of the one traversal of the root.'''
        a, b = leaf('a', '1'), leaf('b', '2')
        root = sequence('3', [a, b])
        path = only_alignment(root, ('a', 'q', 'b'))

        self.assertEqual([log_elem for log_elem, model_elem in path if model_elem == '>>'], ['q'])
        self.assertEqual(shapes(path, root), [(2, 2)])

    def test_a_parallel_siblings_move_does_not_divide_the_traversal(self):
        '''and(seq(a, b), c) on <a, c, b>: c's move is part of a
        traversal of c, not of seq(a, b), so seq(a, b) is traversed
        once.'''
        a, b, c = leaf('a', '1'), leaf('b', '2'), leaf('c', '3')
        body = sequence('4', [a, b])
        root = operator(And, '5', [body, c])
        path = only_alignment(root, ('a', 'c', 'b'))

        self.assertEqual(shapes(path, body), [(2, 2)])
        self.assertEqual(shapes(path, root), [(3, 3)])

    def test_a_skip_divided_from_its_traversal_by_a_log_move_still_counts(self):
        '''The witness of the normal-form order findings, at the level
        that reads it: seq(x, y, z) missing z, with an unmodelled event
        between. One execution of three moves, two matched - so the
        absence of z is in an execution that has a match, where
        voidsalign's mass can see it.'''
        o, x, y, z, e, p = (leaf(name, str(i + 1)) for i, name
                            in enumerate(('o', 'x', 'y', 'z', 'e', 'p')))
        body = sequence('7', [x, y, z])
        approval = operator(Loop, '8', [body, e])
        root = sequence('9', [o, approval, p])
        path = only_alignment(root, ('o', 'x', 'y', 'q', 'p'))

        self.assertEqual(shapes(path, body), [(2, 3)])

    def test_grouping_preserves_the_move_totals(self):
        '''Whatever the partition, every move in the node's subtree is
        in exactly one execution. This is why the summed definitions
        (voidmass_process) cannot move when the grouping does.'''
        a, b = leaf('a', '1'), leaf('b', '2')
        root = sequence('3', [a, b])
        path = only_alignment(root, ('a', 'q', 'b'))

        self.assertEqual(sum(movecount(ex) for ex in executions(path, root)),
                         model_moves(path, root))


class MaximalTraversalTest(unittest.TestCase):
    '''
    Definition [Executions] clause (1): an execution projects into the
    node's skip language, which is what makes it maximal for a non-loop
    node - two traversals concatenated do not project. So ignoring
    unrelated moves must not run two real traversals together.
    '''

    def setUp(self):
        o, x, y, e, p = (leaf(name, str(i + 1)) for i, name
                         in enumerate(('o', 'x', 'y', 'e', 'p')))
        self.body = sequence('6', [x, y])
        approval = operator(Loop, '7', [self.body, e])
        self.root = sequence('8', [o, approval, p])
        self.approval = approval

    def test_two_iterations_of_a_loop_body_remain_two_executions(self):
        path = only_alignment(self.root, ('o', 'x', 'y', 'e', 'x', 'y', 'p'))
        self.assertEqual(shapes(path, self.body), [(2, 2), (2, 2)])

    def test_two_iterations_split_even_with_a_log_move_between_them(self):
        '''The log move is ignored, but xyxy is still not one traversal
        of seq(x, y).'''
        path = only_alignment(self.root, ('o', 'x', 'y', 'e', 'q', 'x', 'y', 'p'))
        self.assertEqual(shapes(path, self.body), [(2, 2), (2, 2)])

    def test_a_partly_recorded_iteration_keeps_its_own_skip(self):
        '''The second iteration records x but not y, so its skip of y
        belongs to that traversal and not to the complete one.'''
        path = only_alignment(self.root, ('o', 'x', 'y', 'e', 'x', 'p'))
        self.assertEqual(shapes(path, self.body), [(2, 2), (1, 2)])


if __name__ == '__main__':
    unittest.main()
