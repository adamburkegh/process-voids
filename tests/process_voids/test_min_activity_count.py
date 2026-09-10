import unittest

from skipalignments import Activity, Tau, Sequence, Xor, And, Loop

from process_voids.coveragemass import min_activity_count, min_activity_count_by_node


def leaf(cls, name, node_id, cost=100000):
    node = cls(None, name, cost)
    node.id = node_id
    return node


class MinActivityCountTest(unittest.TestCase):
    """
    min_activity_count(node) - the minimum number of labelled activities
    any traversal of node's subtree performs (paper's aligncost(empty,
    node)), a purely structural count independent of any alignment cost
    configuration (MM_COST etc) - Sequence/And sum every child (all must
    run), Xor takes the cheapest (smallest-count) child, Loop counts
    only its do-child (the redo-child can always iterate zero times),
    Tau contributes 0 (not a labelled activity), Activity contributes 1.
    """

    def test_activity_leaf_is_one(self):
        a = leaf(Activity, 'a', '1')
        self.assertEqual(min_activity_count(a), 1)

    def test_tau_leaf_is_zero(self):
        tau = leaf(Tau, 'tau', '1', cost=0)
        self.assertEqual(min_activity_count(tau), 0)

    def test_sequence_sums_every_child(self):
        a, b = leaf(Activity, 'a', '1'), leaf(Activity, 'b', '2')
        tree = Sequence(None, [a, b])
        a.set_parent(tree)
        b.set_parent(tree)
        self.assertEqual(min_activity_count(tree), 2)

    def test_and_sums_every_child(self):
        a, b = leaf(Activity, 'a', '1'), leaf(Activity, 'b', '2')
        tree = And(None, [a, b])
        a.set_parent(tree)
        b.set_parent(tree)
        self.assertEqual(min_activity_count(tree), 2)

    def test_xor_takes_the_cheapest_child(self):
        a = leaf(Activity, 'a', '1')
        bc = Sequence(None, [leaf(Activity, 'b', '2'), leaf(Activity, 'c', '3')])
        for c in bc.children:
            c.set_parent(bc)
        tree = Xor(None, [a, bc])
        a.set_parent(tree)
        bc.set_parent(tree)
        self.assertEqual(min_activity_count(tree), 1)

    def test_xor_with_tau_branch_is_zero(self):
        a = leaf(Activity, 'a', '1')
        tau = leaf(Tau, 'tau', '2', cost=0)
        tree = Xor(None, [a, tau])
        a.set_parent(tree)
        tau.set_parent(tree)
        self.assertEqual(min_activity_count(tree), 0)

    def test_loop_counts_only_the_do_child(self):
        do = leaf(Activity, 'do', '1')
        redo = leaf(Activity, 'redo', '2')
        tree = Loop(None, [do, redo])
        do.set_parent(tree)
        redo.set_parent(tree)
        self.assertEqual(min_activity_count(tree), 1)

    def test_nested_structure(self):
        # seq(xor(a, seq(b,c)), loop(d,e)) -> min(1, 2) + 1 = 2
        a = leaf(Activity, 'a', '1')
        bc = Sequence(None, [leaf(Activity, 'b', '2'), leaf(Activity, 'c', '3')])
        for c in bc.children:
            c.set_parent(bc)
        choice = Xor(None, [a, bc])
        a.set_parent(choice)
        bc.set_parent(choice)
        d, e = leaf(Activity, 'd', '4'), leaf(Activity, 'e', '5')
        loop = Loop(None, [d, e])
        d.set_parent(loop)
        e.set_parent(loop)
        tree = Sequence(None, [choice, loop])
        choice.set_parent(tree)
        loop.set_parent(tree)
        self.assertEqual(min_activity_count(tree), 2)


class MinActivityCountByNodeTest(unittest.TestCase):
    """min_activity_count_by_node(tree) - {node: min_activity_count(node)}
    for every node in tree, computed bottom-up in one pass (used for the
    voidmass_table_pn timed-out-variant ceiling, where every node in the
    tree needs its own count, not just the root)."""

    def setUp(self):
        self.a = leaf(Activity, 'a', '1')
        self.b = leaf(Activity, 'b', '2')
        self.tree = Sequence(None, [self.a, self.b])
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)
        self.counts = min_activity_count_by_node(self.tree)

    def test_every_node_present(self):
        self.assertEqual(set(self.counts), {self.tree, self.a, self.b})

    def test_counts_match_min_activity_count(self):
        for node in (self.tree, self.a, self.b):
            with self.subTest(node=node):
                self.assertEqual(self.counts[node], min_activity_count(node))


if __name__ == '__main__':
    unittest.main()
