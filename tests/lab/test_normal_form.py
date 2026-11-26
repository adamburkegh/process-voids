'''
process_voids.normal_form: which adjacent moves of a skip alignment may
be transposed, and reordering one to a chosen precedence.

The independence table (Brief, normal-form order invariance; derived
from the Definition of Skip Alignment, Bar et al. BPM 2025):

    log / skip                  always independent
    sync / skip, skip / skip    independent when both orderings of the
                                model path are in the skip language
    log / sync, log / log       independent when they carry the same
                                activity

Two moves that both consume a log event and carry different activities
are never independent: the trace pins their order.
'''

import unittest

from skipalignments import Activity, And, Loop, Sequence, Skip, Tau, Xor

from lab.normal_form import (
    LOG_FIRST, SKIP_FIRST, independent, move_kind, reorder,
)

ACT_COST = 100000


def leaf(name, node_id, cls=Activity):
    node = cls(None, name, 0 if cls is Tau else ACT_COST)
    node.id = node_id
    return node


def parent(cls, children, node_id):
    node = cls(None, list(children))
    node.id = node_id
    for child in children:
        child.set_parent(node)
    return node


def concurrent_tree():
    '''seq( and(a, b), c ) - a and b concurrent, c after both.'''
    a, b, c = leaf('a', '1'), leaf('b', '2'), leaf('c', '3')
    return parent(Sequence, [parent(And, [a, b], '4'), c], '5'), a, b, c


def sequential_tree():
    '''seq(x, y, z) - nothing concurrent.'''
    x, y, z = leaf('x', '1'), leaf('y', '2'), leaf('z', '3')
    return parent(Sequence, [x, y, z], '4'), x, y, z


def sync(node):
    return (f'{node.name}0', node)


def skip(node):
    return ('>>', Skip(node, ACT_COST))


def log_move(name):
    return (f'{name}0', '>>')


class MoveKindTest(unittest.TestCase):

    def test_the_three_kinds(self):
        tree, a, b, _c = concurrent_tree()
        self.assertEqual(move_kind(sync(a)), 'sync')
        self.assertEqual(move_kind(skip(b)), 'skip')
        self.assertEqual(move_kind(log_move('q')), 'log')


class IndependenceTest(unittest.TestCase):

    def setUp(self):
        self.tree, self.a, self.b, self.c = concurrent_tree()
        self.seq_tree, self.x, self.y, self.z = sequential_tree()

    def test_a_log_move_and_a_skip_are_always_independent(self):
        '''A log move belongs to no execution and consumes an event the
        model does not have, so nothing orders it against a skip.'''
        self.assertTrue(independent(log_move('q'), skip(self.b), self.tree))
        self.assertTrue(independent(skip(self.b), log_move('q'), self.tree))

    def test_a_sync_and_a_skip_in_concurrent_branches_are_independent(self):
        '''The witness case: a and b are children of the same And, so
        both orderings are in the skip language.'''
        self.assertTrue(independent(sync(self.a), skip(self.b), self.tree))
        self.assertTrue(independent(skip(self.b), sync(self.a), self.tree))

    def test_a_sync_and_a_skip_in_sequence_are_not_independent(self):
        '''Sequence pins the order of x and y, so swapping them leaves
        the skip language.'''
        self.assertFalse(independent(sync(self.x), skip(self.y), self.seq_tree))

    def test_two_skips_in_concurrent_branches_are_independent(self):
        self.assertTrue(independent(skip(self.a), skip(self.b), self.tree))

    def test_two_skips_in_sequence_are_not(self):
        self.assertFalse(independent(skip(self.x), skip(self.y), self.seq_tree))

    def test_two_consuming_moves_with_different_activities_are_never_independent(self):
        '''The trace pins their order, whatever the model allows - even
        for a and c, which are not concurrent anyway, and for the
        concurrent pair a and b.'''
        self.assertFalse(independent(sync(self.a), sync(self.b), self.tree))
        self.assertFalse(independent(sync(self.a), log_move('q'), self.tree))

    def test_two_consuming_moves_with_the_same_activity_are_independent(self):
        '''Same label, so the projection onto the trace is unchanged.'''
        self.assertTrue(independent(sync(self.a), log_move('a'), self.tree))
        self.assertTrue(independent(log_move('a'), log_move('a'), self.tree))


class ReorderTest(unittest.TestCase):

    def setUp(self):
        self.tree, self.a, self.b, self.c = concurrent_tree()
        # the aligner's normal form for <a, c>: sync before skip
        self.delta = [sync(self.a), skip(self.b), sync(self.c)]

    def _shape(self, path):
        return [move_kind(m) for m in path]

    def test_log_first_leaves_the_aligners_normal_form_alone(self):
        self.assertEqual(reorder(self.delta, self.tree, LOG_FIRST), self.delta)

    def test_skip_first_moves_the_skip_ahead_of_the_independent_sync(self):
        self.assertEqual(self._shape(reorder(self.delta, self.tree, SKIP_FIRST)),
                         ['skip', 'sync', 'sync'])

    def test_reordering_never_moves_a_skip_past_a_dependent_sync(self):
        '''c is after the And in sequence, so the skip of b cannot
        overtake c/c however the precedence is set.'''
        reordered = reorder(self.delta, self.tree, SKIP_FIRST)
        self.assertEqual([m[1].node.name if move_kind(m) == 'skip' else m[1].name
                          for m in reordered], ['b', 'a', 'c'])

    def test_reordering_is_idempotent(self):
        once = reorder(self.delta, self.tree, SKIP_FIRST)
        self.assertEqual(reorder(once, self.tree, SKIP_FIRST), once)

    def test_a_sequential_alignment_is_unchanged_under_either_precedence(self):
        seq_tree, x, y, z = sequential_tree()
        path = [sync(x), skip(y), sync(z)]
        self.assertEqual(reorder(path, seq_tree, SKIP_FIRST), path)
        self.assertEqual(reorder(path, seq_tree, LOG_FIRST), path)

    def test_reordering_preserves_the_multiset_of_moves(self):
        reordered = reorder(self.delta, self.tree, SKIP_FIRST)
        self.assertEqual(sorted(map(str, reordered)), sorted(map(str, self.delta)))


if __name__ == '__main__':
    unittest.main()
