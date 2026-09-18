'''
The payment_partial worked example (lab.fixtures), pinned value by value
through the full ebi-backed pipeline: its skip alignments, the SLPN's
split of the one tie, skip probabilities, and voidsalign3, voidsat2 and
voidmass_process at three nodes: the approval loop, the xor that makes
it optional, and s.

The loop is masked under the xor - where approval is absent the xor
took its silent branch - so its skip probability is 0 and its two
matchprob metrics read their mass alone. The xor carries that absence
(skipprob 1/4) with the same mass, so it is where the detection term
acts. s is mandatory and skipped outright in <o a p> (skipprob 1/8).

    N = seq( o, xor( tau, loop(a, e) ), s, p )

    variant        cases     weight
    o a s p        s1, s5    1/4
    o a e a s p    s2        1/8
    o s p          s3, s4    1/4
    o a p          s6        1/8
    o a e s p      s7, s8    1/4

<o a e s p> has two optimal alignments: a partial traversal of the loop
whose closing a is a skip, and a log move discarding e. The SLPN splits
them 3/14 : 11/14, which is what gives skipprob(a) = 1/28. The metrics
do not use that split: each tied alignment carries 1/|Upsilon_sigma| of
its trace, so the two count a half each.

LoopTieIsAMatterOfSizesTest pins why the tie arises: this loop's do- and
redo-parts are one activity each, and a loop whose parts differ in size
does not tie.
'''

import unittest
from datetime import datetime, timedelta

from skipalignments import Activity, Aligner, Loop, Sequence, Skip

from lab.fixtures import build_payment_partial_log, build_payment_partial_tree
from lab.run import CLASSICAL_ALIGNMENT_TIMEOUT
from process_voids.coveragemass import _variant_key, executions
from process_voids.metric_context import CellContext
from process_voids.voidmass_pn import build_id_net
from process_voids.voidsalign3 import smatchcount, smovecount, voidsalign3
from process_voids.voidsat2 import adratio, admass, voidsat2

PARTIAL = ('o', 'a', 'e', 's', 'p')

ACT_COST = 100000

Aligner.set_level_incentive(0)

T0 = datetime(2026, 1, 1, 12, 0, 0)


def timed_trace(*name_hours):
    return [{'concept:name': name, 'time:timestamp': T0 + timedelta(hours=hours)}
            for name, hours in name_hours]


SIGMA7 = timed_trace(('o', 0), ('a', 1), ('e', 5), ('s', 10), ('p', 11))
SIGMA8 = timed_trace(('o', 0), ('a', 1), ('e', 5), ('s', 6), ('p', 7))


def move(log_elem, model_elem):
    '''One move as 'log/model': sk(x) for a skip, tau(Kind) for a silent path.'''
    log_side = '>>' if log_elem == '>>' else log_elem.rstrip('0123456789')
    if isinstance(model_elem, str):
        model_side = model_elem
    elif isinstance(model_elem, Skip):
        model_side = f'sk({model_elem.node.name})'
    elif type(model_elem).__name__ == 'TauPath':
        model_side = f'tau({type(model_elem.node).__name__})'
    else:
        model_side = model_elem.name
    return f'{log_side}/{model_side}'


def shape(path):
    return tuple(move(log_elem, model_elem) for log_elem, model_elem in path)


class PaymentPartialWorkedTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tree = build_payment_partial_tree()
        cls.o, cls.choice, cls.s, cls.p = cls.tree.children
        cls.tau, cls.loop = cls.choice.children
        cls.a, cls.e = cls.loop.children
        cls.log = build_payment_partial_log()
        cls.ctx = CellContext(log=cls.log, tree=cls.tree,
                              slpn_path='var/lab/test_payment_partial_worked.slpn',
                              ppt_weights=None, listeners=[],
                              classical_net=build_id_net(cls.tree),
                              classical_timeout=CLASSICAL_ALIGNMENT_TIMEOUT)
        cls.dv = cls.ctx.stage('dv')
        cls.classical, _ = cls.ctx.stage('classical')
        cls.paths = {k: [s.path for s in v] for k, v in cls.dv.skip_dict_backup.items()}

    def _partial_paths(self):
        '''(partial traversal, log-move explanation) for <o a e s p>.'''
        paths = self.paths[_variant_key(PARTIAL)]
        with_skip = [p for p in paths if any(isinstance(m, Skip) for _, m in p)]
        without = [p for p in paths if not any(isinstance(m, Skip) for _, m in p)]
        return with_skip[0], without[0]

    # --- inputs -----------------------------------------------------------

    def test_variant_probabilities(self):
        self.assertEqual(self.dv.pl, {
            ('o', 'a', 's', 'p'): 1 / 4, ('o', 'a', 'e', 'a', 's', 'p'): 1 / 8,
            ('o', 's', 'p'): 1 / 4, ('o', 'a', 'p'): 1 / 8, PARTIAL: 1 / 4})

    def test_skip_alignments(self):
        '''
        As the worked example assumes, with one correction of notation:
        <o s p> takes the xor's silent branch, a tau path costing
        nothing, not a skip of the xor.
        '''
        expected = {
            ('o', 'a', 's', 'p'): {('o/o', 'a/a', 's/s', 'p/p')},
            ('o', 'a', 'e', 'a', 's', 'p'): {('o/o', 'a/a', 'e/e', 'a/a', 's/s', 'p/p')},
            ('o', 's', 'p'): {('o/o', '>>/tau(Xor)', 's/s', 'p/p')},
            ('o', 'a', 'p'): {('o/o', 'a/a', '>>/sk(s)', 'p/p')},
            PARTIAL: {('o/o', 'a/a', 'e/e', '>>/sk(a)', 's/s', 'p/p'),
                      ('o/o', 'a/a', 'e/>>', 's/s', 'p/p')},
        }
        for variant, shapes in expected.items():
            with self.subTest(variant=variant):
                paths = self.paths[_variant_key(variant)]
                self.assertEqual(len(paths), len(shapes))
                self.assertEqual({shape(p) for p in paths}, shapes)

    def test_the_slpn_splits_the_tie_three_fourteenths_to_eleven(self):
        '''
        cond_prob_'s states are keyed by their indexed log labels. The
        log-move explanation's path leaves the discarded e2 out
        altogether, so without the indices it would read as <o a s p>.
        '''
        by_log_labels = {tuple(log_elem for log_elem, _ in state.path): prob
                         for state, prob in self.dv.cond_prob_[self.tree].items()}
        self.assertAlmostEqual(by_log_labels[('o0', 'a1', 'e2', '>>', 's3', 'p4')], 3 / 14,
                               places=12)
        self.assertAlmostEqual(by_log_labels[('o0', 'a1', 's3', 'p4')], 11 / 14, places=12)

    def test_skip_probabilities(self):
        expected = [(self.tree, 0.0), (self.o, 0.0), (self.choice, 1 / 4), (self.tau, 0.0),
                    (self.loop, 0.0), (self.a, 1 / 28), (self.e, 0.0), (self.s, 1 / 8),
                    (self.p, 0.0)]
        for node, prob in expected:
            with self.subTest(node=node.id):
                self.assertAlmostEqual(self.dv.skip_probs[node], prob, places=12)

    # --- voidsalign3 ------------------------------------------------------

    def test_voidsalign3_ratio_of_the_partial_traversal_is_two_thirds(self):
        partial, log_move = self._partial_paths()
        [execution] = list(executions(partial, self.loop))
        self.assertEqual((smatchcount(execution), smovecount(execution)), (2, 3))
        [execution] = list(executions(log_move, self.loop))
        self.assertEqual((smatchcount(execution), smovecount(execution)), (1, 1))

    def test_voidsalign3_at_the_loop(self):
        '''
        <o s p> never enters the loop, so it is not observed there and
        leaves the average; obscount is 3/4. Every other ratio is 1
        except the partial traversal's 2/3, at half of <o a e s p>'s
        1/4: mass = (1/4 + 1/8 + 1/8 + 1/8 * (2/3 + 1)) / (3/4) = 17/18.
        matchprob(loop) = 1, so void = 1/18.
        '''
        self.assertAlmostEqual(
            voidsalign3(self.loop, self.dv.skip_dict_backup, self.dv.pl, self.dv.skip_probs),
            1 / 18, places=12)

    # --- voidsat2 ---------------------------------------------------------

    def test_voidsat2_ratios_of_sigma7_and_sigma8(self):
        '''
        The partial traversal's skip of a shares the interval before s
        with s/s. sigma7: a/a 1h, e/e 4h, sk(a) 2.5h, so 5/7.5 = 2/3.
        sigma8: sk(a) 0.5h, so 5/5.5 = 10/11. Same variant, same
        alignment, different ratio - what no move-based metric can see.
        The log-move explanation reads 1 for both.
        '''
        partial, log_move = self._partial_paths()
        self.assertAlmostEqual(adratio(self.loop, self.tree, partial, SIGMA7), 2 / 3, places=12)
        self.assertAlmostEqual(adratio(self.loop, self.tree, partial, SIGMA8), 10 / 11, places=12)
        self.assertAlmostEqual(adratio(self.loop, self.tree, log_move, SIGMA7), 1.0, places=12)
        self.assertAlmostEqual(adratio(self.loop, self.tree, log_move, SIGMA8), 1.0, places=12)

    def test_voidsat2_at_the_loop(self):
        '''
        Per trace, not per variant. s3 and s4 never enter the loop, so
        obscount = 6. s1, s2, s5, s6 read 1; s7 reads (2/3 + 1)/2 and
        s8 (10/11 + 1)/2. admass = (4 + 5/6 + 21/22) / 6 = 191/198, and
        matchprob(loop) = 1, so void = 7/198.
        '''
        self.assertAlmostEqual(admass(self.loop, self.tree, self.log, self.paths),
                               191 / 198, places=12)
        self.assertAlmostEqual(
            voidsat2(self.loop, self.tree, self.log, self.paths, self.dv.skip_probs),
            7 / 198, places=12)

    # --- voidmass_process -------------------------------------------------

    def test_voidmass_process_at_the_loop(self):
        '''
        Classical alignments, same tie, same halves. The only deficit in
        the loop is the partial traversal's model move for a: 1/4 * 1/2
        = 1/8. Root moves: 4 * 1/4 + 6 * 1/8 + 3 * 1/4 + 4 * 1/8
        + (6 + 4)/2 * 1/4 = 17/4 (the silent move in <o s p> is not a
        move). voidmass_process = (1/8) / (17/4) = 1/34, and with no
        timeouts the lower and upper bounds agree.
        '''
        self.assertEqual(self.classical.timed_out_count, 0)
        loop_row = self.classical.table[self.loop]
        self.assertAlmostEqual(loop_row['deficit_lower'], 1 / 8, places=12)
        self.assertAlmostEqual(self.classical.table[self.tree]['movecount'], 17 / 4, places=12)
        self.assertAlmostEqual(loop_row['voidmass_process_lower'], 1 / 34, places=12)
        self.assertAlmostEqual(loop_row['voidmass_process_upper'], 1 / 34, places=12)

    # --- the xor: the same mass, with matchprob acting --------------------

    def test_xor(self):
        '''
        <o s p> takes the xor's silent branch: no synchronous move, so
        it leaves both averages and the xor's masses are the loop's.
        matchprob(xor) = 3/4:
          voidsalign3 = 1 - (3/4)(17/18)   = 7/24
          voidsat2    = 1 - (3/4)(191/198) = 73/264
        voidmass_process has no matchprob term, and the xor's only
        deficit is the loop's, so it reads 1/34 at both nodes.
        '''
        self._assert_node(self.choice, voidsalign3_=7 / 24, voidsat2_=73 / 264,
                          voidmass_process=1 / 34)

    # --- s: a mandatory leaf, skipped outright ----------------------------

    def test_s(self):
        '''
        Every execution of s that has a match is s/s alone, ratio 1 by
        moves and by duration; <o a p>'s sk(s) has no match, so it is
        not observed and is carried by matchprob(s) = 7/8 instead:
          voidsalign3 = voidsat2 = 1 - (7/8)(1) = 1/8
        voidmass_process: sk(s) is a model move in <o a p> at 1/8, so
        (1/8) / (17/4) = 1/34.
        '''
        self._assert_node(self.s, voidsalign3_=1 / 8, voidsat2_=1 / 8,
                          voidmass_process=1 / 34)

    def _assert_node(self, node, voidsalign3_, voidsat2_, voidmass_process):
        with self.subTest(metric='voidsalign3'):
            self.assertAlmostEqual(
                voidsalign3(node, self.dv.skip_dict_backup, self.dv.pl, self.dv.skip_probs),
                voidsalign3_, places=12)
        with self.subTest(metric='voidsat2'):
            self.assertAlmostEqual(
                voidsat2(node, self.tree, self.log, self.paths, self.dv.skip_probs),
                voidsat2_, places=12)
        row = self.classical.table[node]
        for bound in ('lower', 'upper'):
            with self.subTest(metric=f'voidmass_process_{bound}'):
                self.assertAlmostEqual(row[f'voidmass_process_{bound}'], voidmass_process,
                                       places=12)


def _node(node, node_id, children=()):
    node.id = node_id
    for child in children:
        child.set_parent(node)
    return node


def _loop_tree(do_names, redo_name):
    '''seq( o, loop(do, redo), p ), do a single activity or a sequence.'''
    ids = iter(str(i) for i in range(1, 20))
    leaf = lambda name: _node(Activity(None, name, ACT_COST), next(ids))
    o = leaf('o')
    do_leaves = [leaf(name) for name in do_names]
    do = do_leaves[0] if len(do_leaves) == 1 else _node(
        Sequence(None, do_leaves), next(ids), do_leaves)
    redo = leaf(redo_name)
    loop = _node(Loop(None, [do, redo]), next(ids), [do, redo])
    p = leaf('p')
    return _node(Sequence(None, [o, loop, p]), next(ids), [o, loop, p])


def _align(tree, trace):
    states, _ = Aligner(tree).align_normal_form(
        list(trace), [ACT_COST] * len(trace), True, timeout=100)
    return {shape(s.path) for s in states}


class LoopTieIsAMatterOfSizesTest(unittest.TestCase):
    '''
    A loop traversal missing its closing do-part ties with discarding
    the redo only when the two cost the same: one model move per missing
    do-activity against one log move per redo and recorded do event.
    '''

    def test_equal_parts_tie(self):
        '''loop(a, e) on <o a e p>: one model move for a, or one log move for e.'''
        self.assertEqual(_align(_loop_tree(['a'], 'e'), ('o', 'a', 'e', 'p')), {
            ('o/o', 'a/a', 'e/e', '>>/sk(a)', 'p/p'),
            ('o/o', 'a/a', 'e/>>', 'p/p'),
        })

    def test_unequal_parts_do_not_tie(self):
        '''
        loop(seq(x, y, z), e) on <o x y z e x z p>: completing the second
        traversal costs one model move for y; discarding it costs log
        moves for e, x and z. The partial traversal is the only optimal
        alignment.
        '''
        self.assertEqual(
            _align(_loop_tree(['x', 'y', 'z'], 'e'), ('o', 'x', 'y', 'z', 'e', 'x', 'z', 'p')),
            {('o/o', 'x/x', 'y/y', 'z/z', 'e/e', 'x/x', '>>/sk(y)', 'z/z', 'p/p')})


if __name__ == '__main__':
    unittest.main()
