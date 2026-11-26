'''
Task C of the normal-form order brief: is exec(m, msub, .) - the
partition of alignment positions into per-traversal executions - stable
when independent moves are transposed?

The case the brief could not close on paper is two model-side moves
belonging to DIFFERENT traversals of a loop, where a transposition might
reassign a move from one execution to another.

For that case the implementation says it cannot arise. Two model-side
moves are independent only where both orderings of the model path are in
the skip language, which for a process tree means their lowest common
ancestor is an And and they descend through different branches of it
(process_voids.normal_form). Two moves in different traversals of one
loop have that loop as their lowest common ancestor, and a Loop is not an
And, so no transposition crosses a traversal boundary. On that reading
the hypothesis holds for model-side moves - with the caveat in
normal_form's docstring, that the And criterion is this module's reading
of the skip language.

The case the brief treats as immediate held only once exec followed
Definition [Executions]. It used to group a non-Loop node's moves by
runs of consecutive POSITIONS, so a log move between two of a node's
moves split one execution in two and voidsalign moved with it - the
reading pinned on branch claude/normal-form-order-invariance. An
execution is a loose subsequence, and only moves that are themselves
part of a traversal of the node divide it, so that split was a defect
rather than evidence about the precedence. LogMoveTranspositionTest now
pins the invariance over the same witness trace.
'''

import unittest

from skipalignments import Activity, Aligner, Loop, Sequence

from process_voids.coveragemass import executions, matchcount, movecount
from lab.normal_form import (
    LOG_FIRST, SKIP_FIRST, independent, move_kind, reorder,
)
from process_voids.voidsalign3 import smatchcount, smovecount, voidsalign3
from process_voids.voidmass_pn import voidmass_table_pn

ACT_COST = 100000

Aligner.set_level_incentive(0)


def build_loop_tree():
    '''seq( o, loop(seq(x, y, z), e), p ) - a loop whose body can be
    partly skipped, with no And anywhere.'''
    o, x, y, z, e, p = (Activity(None, name, ACT_COST) for name in ('o', 'x', 'y', 'z', 'e', 'p'))
    for node, node_id in ((o, '1'), (x, '2'), (y, '3'), (z, '4'), (e, '5'), (p, '6')):
        node.id = node_id
    body = Sequence(None, [x, y, z])
    body.id = '7'
    for node in (x, y, z):
        node.set_parent(body)
    loop = Loop(None, [body, e])
    loop.id = '8'
    body.set_parent(loop)
    e.set_parent(loop)
    tree = Sequence(None, [o, loop, p])
    tree.id = '9'
    for node in (o, loop, p):
        node.set_parent(tree)
    return tree, loop, body


def align(tree, trace):
    states, _ = Aligner(tree).align_normal_form(
        list(trace), [ACT_COST] * len(trace), True, timeout=100)
    return states


def pairs_of_differing_kind(path):
    return [(i, move_kind(path[i]), move_kind(path[i + 1]))
            for i in range(len(path) - 1)
            if move_kind(path[i]) != move_kind(path[i + 1])]


def independent_pairs(path):
    return [(i, first, second) for i, first, second in pairs_of_differing_kind(path)
            if independent(path[i], path[i + 1])]


def execution_shapes(path, node):
    '''(matchcount, movecount, smatchcount, smovecount) per execution.'''
    return [(matchcount(ex), movecount(ex), smatchcount(ex), smovecount(ex))
            for ex in executions(path, node)]


class LoopTraversalTest(unittest.TestCase):
    '''
    Two traversals of loop(seq(x, y, z), e), the second missing y - the
    brief's suggested probe. No transposition is available across the
    traversal boundary, so exec cannot shift.
    '''

    @classmethod
    def setUpClass(cls):
        cls.tree, cls.loop, cls.body = build_loop_tree()
        cls.states = align(cls.tree, ('o', 'x', 'y', 'z', 'e', 'x', 'z', 'p'))
        cls.path = cls.states[0].path

    def test_the_trace_traverses_the_loop_twice_and_skips_part_of_one(self):
        self.assertEqual(len(self.states), 1)
        shapes = execution_shapes(self.path, self.body)
        self.assertEqual(len(shapes), 2)
        self.assertEqual(sorted(shapes), [(2, 3, 2, 3), (3, 3, 3, 3)])

    def test_no_model_side_transposition_is_available_at_all(self):
        '''The tree has no And, so every independent pair of differing
        kinds would have to involve a log move - and there are none
        here, since every event is in the model.'''
        self.assertEqual(independent_pairs(self.path), [])

    def test_reordering_is_the_identity_under_either_precedence(self):
        self.assertEqual(reorder(self.path, self.tree, SKIP_FIRST), list(self.path))
        self.assertEqual(reorder(self.path, self.tree, LOG_FIRST), list(self.path))

    def test_exec_and_the_voidsalign_counts_are_therefore_unchanged(self):
        flipped = reorder(self.path, self.tree, SKIP_FIRST)
        for node in (self.loop, self.body):
            with self.subTest(node=node.id):
                self.assertEqual(execution_shapes(flipped, node),
                                 execution_shapes(self.path, node))


class LogMoveTranspositionTest(unittest.TestCase):
    '''
    The transposition that does exist, over the witness that used to
    break C1.

    An unmodelled event beside a skip: log/skip is independent always,
    per the brief's own table. The brief treats this case as immediate,
    because a log move belongs to no execution, and that is now also
    what the code says: coveragemass._group_executions divides a node's
    moves only where the next traversal of that node begins, and a log
    move is part of no traversal of it.

    Either way round, the body seq(x, y, z) reads one execution of three
    moves with two matched, so voidsalign reads 1/3. Under the previous
    consecutive-position grouping the log move sitting between the
    body's moves split that into one whole execution and one with no
    match at all, and since voidsalign's mass counts only executions
    that have a match, the skip was dropped and the void fell to 0.
    '''

    @classmethod
    def setUpClass(cls):
        cls.tree, cls.loop, cls.body = build_loop_tree()
        # q is in no model position; z is never recorded in the one traversal
        cls.states = align(cls.tree, ('o', 'x', 'y', 'q', 'p'))
        cls.path = cls.states[0].path

    def test_an_independent_log_skip_pair_exists(self):
        found = independent_pairs(self.path)
        self.assertTrue(found, f'no independent pair in {[move_kind(m) for m in self.path]}')
        self.assertTrue(all({first, second} == {'log', 'skip'} for _i, first, second in found))

    def test_the_two_precedences_order_it_differently(self):
        log_first = reorder(self.path, self.tree, LOG_FIRST)
        skip_first = reorder(self.path, self.tree, SKIP_FIRST)
        self.assertNotEqual(log_first, skip_first)
        self.assertEqual(sorted(map(str, log_first)), sorted(map(str, skip_first)))

    def test_exec_is_unchanged_when_the_log_move_moves_inside(self):
        '''
        The sequence body, where the split used to show: one traversal
        of seq(x, y, z) recording two of its three moves, wherever the
        unmodelled event sits. [(2, 2, 2, 2), (0, 1, 0, 1)] was the
        defective reading.
        '''
        self.assertEqual(execution_shapes(reorder(self.path, self.tree, SKIP_FIRST),
                                          self.body),
                         [(2, 3, 2, 3)])
        self.assertEqual(execution_shapes(reorder(self.path, self.tree, LOG_FIRST),
                                          self.body),
                         [(2, 3, 2, 3)])

    def test_the_loop_itself_is_unaffected(self):
        log_first = reorder(self.path, self.tree, LOG_FIRST)
        skip_first = reorder(self.path, self.tree, SKIP_FIRST)
        self.assertEqual(execution_shapes(skip_first, self.loop),
                         execution_shapes(log_first, self.loop))

    def test_voidsalign_is_invariant_to_the_precedence(self):
        '''C1, on the witness that falsified it while exec grouped by
        consecutive positions: 1/3 either way, where log_first used to
        read 0. This assertion is the proof of the grouping fix.'''
        values = {}
        for name, precedence in (('log_first', LOG_FIRST), ('skip_first', SKIP_FIRST)):
            values[name] = self._voidsalign(reorder(self.path, self.tree, precedence))
        self.assertAlmostEqual(values['skip_first'], 1 / 3, places=12)
        self.assertAlmostEqual(values['log_first'], 1 / 3, places=12)

    def _voidsalign(self, path):
        from types import SimpleNamespace

        from process_voids.coveragemass import _variant_key

        variant = ('o', 'x', 'y', 'q', 'p')
        skip_probs = {node: 0.0 for node in
                      [self.tree, self.loop, self.body, *self.body.children,
                       *self.tree.children, *self.loop.children]}
        skip_dict = {_variant_key(variant): [SimpleNamespace(path=path)]}
        return voidsalign3(self.body, skip_dict, {variant: 1.0}, skip_probs)


class VoidmassProcessReadsNoNormalFormTest(unittest.TestCase):
    '''
    C2, pinned structurally rather than by a comparison run: the
    classical metric builds its own alignments and takes no
    skip-alignment input, so there is no precedence for it to depend on.
    '''

    def test_voidmass_table_pn_takes_no_skip_alignment_argument(self):
        import inspect

        names = set(inspect.signature(voidmass_table_pn).parameters)
        self.assertNotIn('skip_dict', names)
        self.assertNotIn('paths', names)
        self.assertIn('variant_probs', names)


if __name__ == '__main__':
    unittest.main()
