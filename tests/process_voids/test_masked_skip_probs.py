'''
Skip probabilities of a node masked inside a lumped skip, and the
metrics that multiply them, on the payment running example
(lab.fixtures), through the full ebi-backed pipeline (pvoid.skipprob).

The optimal alignments of sigma3 and sigma4 (<o, s, p>) skip the whole
approval loop as one lumped move. The loop's do-child a sits inside that
lump: it has no execution of its own in those alignments, so those
traces are on neither side of a's skip-probability ratio, and wherever a
does execute it synchronises. skip_prob(a) is therefore 0, while the
lumped loop's own skip_prob is 2/6.

These are integration pins at the process-voids boundary; skip-alignments
owns the engine's own regression tests. The metrics pinned here are
products skip_prob(node) * mass. A masked node has no execution in the
lumped alignments, for coveragemass.executions as for the skip
probability, so neither factor counts those traces -
SkipProbCountsTheSameTracesAsExecutionsTest checks that for every node.
'''

import tempfile
import unittest
from pathlib import Path

from skipalignments import Activity, Sequence

from lab.fixtures import build_running_example_log, build_running_example_tree
from process_voids import dtlog, pvoid
from process_voids.coveragemass import executions, matchcount, voidmass_table, voidsat

ACT_COST = 100000


class RunningExampleMaskedDoChildTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.log = build_running_example_log()
        cls.tree = build_running_example_tree()
        cls.approval = cls.tree.children[1]
        cls.a = cls.approval.children[0]
        # Temp dir for the intermediate SLPN: nothing here creates var/lab/.
        with tempfile.TemporaryDirectory() as tmp:
            cls.dv = pvoid.skipprob(cls.log, cls.tree,
                                    str(Path(tmp) / 'test_masked_skip_probs.slpn'))

    def test_lumped_loop_keeps_its_skip_probability(self):
        self.assertAlmostEqual(self.dv.skip_probs[self.approval], 2 / 6, places=9)

    def test_masked_do_child_has_zero_skip_probability(self):
        self.assertEqual(self.dv.skip_probs[self.a], 0.0)

    def test_voidage_at_masked_do_child_is_zero(self):
        table = voidmass_table(self.tree, self.dv.skip_dict_backup, self.dv.pl,
                               skip_probs=self.dv.skip_probs)
        self.assertEqual(table[self.a]['voidage_subprocess'], 0.0)
        self.assertEqual(table[self.a]['voidage_process'], 0.0)

    def test_voidsat_at_masked_do_child_is_zero(self):
        self.assertEqual(voidsat(self.a, self.tree, self.dv, self.log), 0.0)


class SkipProbCountsTheSameTracesAsExecutionsTest(unittest.TestCase):
    '''
    Every skip_prob * mass metric pairs skip-alignments' skip_prob with a
    mass built on coveragemass.executions. The two factors are only
    consistent if they count the same traces: for every node n,

        skip_prob(n) == sum_sigma P(sigma) [n's execution in sigma is a skip]
                        / sum_sigma P(sigma) [n has an execution in sigma]

    with executions taken from coveragemass.executions over
    dv.skip_dict_backup, and an execution with no synchronous move
    counted as a skip. The fixture has one optimal alignment per trace
    and at most one execution per node per trace, so the ratio is well
    defined.

    Model: seq(a, seq(b, seq(c, d)), e) over <a,b,c,d,e> x2, <a,b,e>,
    <a,e>. <a,e> lumps seq(b, seq(c, d)); <a,b,e> lumps seq(c, d).
    '''

    @classmethod
    def setUpClass(cls):
        def activity(label, node_id):
            node = Activity(None, label, ACT_COST)
            node.id = node_id
            return node

        def sequence(children, node_id):
            node = Sequence(None, children)
            node.id = node_id
            for child in children:
                child.set_parent(node)
            return node

        a, b, c, d, e = (activity(label, str(i + 1)) for i, label in enumerate('abcde'))
        cd = sequence([c, d], '6')
        bcd = sequence([b, cd], '7')
        cls.tree = sequence([a, bcd, e], '8')
        cls.nodes = {'a': a, 'b': b, 'c': c, 'd': d, 'e': e,
                     'seq(c, d)': cd, 'seq(b, seq(c, d))': bcd, 'root': cls.tree}

        traces = ['a:0 b:1 c:2 d:3 e:4'] * 2 + ['a:0 b:1 e:2', 'a:0 e:1']
        log = dtlog.convert_timed(*traces, names=[f'c{i}' for i in range(len(traces))],
                                  time_unit='hours')
        # Temp dir for the intermediate SLPN - see RunningExampleMaskedDoChildTest.
        with tempfile.TemporaryDirectory() as tmp:
            cls.dv = pvoid.skipprob(log, cls.tree,
                                    str(Path(tmp) / 'test_masked_skip_probs_nested.slpn'))

    def execution_skip_ratio(self, node):
        skipped = executed = 0.0
        for variant, prob in self.dv.pl.items():
            states = self.dv.skip_dict_backup[', '.join(variant)]
            self.assertEqual(len(states), 1, f'{variant}: expected one optimal alignment')
            node_executions = executions(states[0].path, node)
            self.assertLessEqual(len(node_executions), 1,
                                 f'{variant}: expected at most one execution of {node.id}')
            if node_executions:
                executed += prob
                if matchcount(node_executions[0]) == 0:
                    skipped += prob
        return skipped / executed

    def test_skip_prob_matches_execution_skip_ratio(self):
        for name, node in self.nodes.items():
            with self.subTest(node=name):
                self.assertAlmostEqual(self.dv.skip_probs[node],
                                       self.execution_skip_ratio(node), places=9)


if __name__ == '__main__':
    unittest.main()
