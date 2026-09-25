'''
process_voids.pvoid's public surface: one function per void metric, each
returning its value at every node, and a CLI printing one metric as a
tree. Values are checked against the payment_partial values pinned in
test_payment_partial_worked.
'''

import unittest

from lab.fixtures import build_payment_partial_log, build_payment_partial_tree
from process_voids import pvoid


def _nodes(tree):
    root = tree
    o, choice, s, p = root.children
    tau, loop = choice.children
    return {'root': root, 'xor': choice, 'loop': loop, 's': s, 'p': p}


class MetricFunctionsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        log = build_payment_partial_log()
        cls.results = {}
        for name in ('voidsalign', 'voidsat', 'voidmass_process'):
            tree = build_payment_partial_tree()
            cls.results[name] = (_nodes(tree), getattr(pvoid, name)(log, tree))

    def _check(self, name, expected):
        nodes, values = self.results[name]
        for node_name, value in expected.items():
            with self.subTest(metric=name, node=node_name):
                self.assertAlmostEqual(values[nodes[node_name]], value, places=12)

    def test_every_node_has_a_value(self):
        for name, (nodes, values) in self.results.items():
            with self.subTest(metric=name):
                self.assertEqual(len(values), 9)

    def test_voidsalign(self):
        self._check('voidsalign', {'root': 5 / 96, 'xor': 7 / 24, 'loop': 1 / 18,
                                   's': 1 / 8, 'p': 0.0})

    def test_voidsat(self):
        self._check('voidsat', {'xor': 73 / 264, 'loop': 7 / 198, 's': 1 / 8, 'p': 0.0})

    def test_voidmass_process(self):
        self._check('voidmass_process', {'xor': 1 / 34, 'loop': 1 / 34, 's': 1 / 34,
                                         'p': 0.0})

    def test_metrics_names_the_three_functions(self):
        self.assertEqual(pvoid.METRICS, {'voidsalign': pvoid.voidsalign,
                                         'voidsat': pvoid.voidsat,
                                         'voidmass_process': pvoid.voidmass_process})


class ShowTreeTest(unittest.TestCase):

    def setUp(self):
        self.tree = build_payment_partial_tree()
        nodes = _nodes(self.tree)
        all_nodes = [self.tree, *self.tree.children, *nodes['xor'].children,
                     *nodes['loop'].children]
        self.skip_probs = {n: 0.0 for n in all_nodes}
        self.skip_probs[nodes['xor']] = 0.25
        self.values = {n: 0.0 for n in all_nodes}
        self.values[nodes['loop']] = 0.5
        self.lines = pvoid.show_tree(self.tree, self.skip_probs, self.values).splitlines()

    def test_one_line_per_node(self):
        self.assertEqual(len(self.lines), 9)

    def test_each_line_is_skip_probability_then_the_metric(self):
        '''Two numbers per node; the xor's skip probability is bracketed,
        since it can be traversed silently.'''
        xor_line, loop_line = self.lines[2], self.lines[4]
        self.assertEqual(xor_line.split(' : ', 1)[1], '[ 0.25 ], 0.0')
        self.assertEqual(loop_line.split(' : ', 1)[1], '0.0, 0.5')


class CliTest(unittest.TestCase):

    def test_default_metric_is_voidsalign(self):
        self.assertEqual(pvoid.parse_args(['log.xes', 'model.ptml']).metric, 'voidsalign')

    def test_metric_names_carry_no_version(self):
        self.assertEqual(pvoid.parse_args(['l', 'm', '--metric', 'voidsat']).metric, 'voidsat')
        with self.assertRaises(SystemExit):
            pvoid.parse_args(['l', 'm', '--metric', 'voidsalign3'])


if __name__ == '__main__':
    unittest.main()
