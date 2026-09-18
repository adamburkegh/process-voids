"""
Lemma [Additivity] for Definition [Void by Process-Relative Alignment
Moves] - voidmass_process, from classical alignments
(process_voids.voidmass_pn.voidmass_table_pn).

For an antichain C of nodes whose leaves cover every labelled leaf of m,
the values at C sum to the value at the root. The denominator is the
root's move total at every node, so this holds exactly when the
numerators - deficits - partition: every move lies in the executions of
exactly one node of C. Classical alignments put every move on a leaf, so
they do.

Property-based: hypothesis generates the tree, the log and the cut, so
the lemma is checked over tree shapes and cuts no hand-written case
would reach - choices, loops, parallel blocks, silent leaves left out of
the cut, and log moves for activities the model does not have.

Only the no-timeout case. A timed-out variant adds the same upper-bound
allowance at every node, so _upper is not additive by construction; see
voidmass_pn.voidmass_table_pn.
"""

import itertools
import unittest

from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from skipalignments import Activity, And, Loop, Sequence, Tau, Xor

from process_voids.voidmass_pn import build_id_net, voidmass_table_pn

COST = 100000
TOLERANCE = 1e-9
OUT_OF_ALPHABET = 'z'

_OPERATORS = {'seq': Sequence, 'xor': Xor, 'and': And, 'loop': Loop}

# A tree shape: ('act', ()) or ('tau', ()) at a leaf, else (operator,
# children). Loops take exactly a do-child and a redo-child.
_LEAF = st.sampled_from([('act', ()), ('act', ()), ('act', ()), ('tau', ())])


def _operator(children):
    # Sequence weighted up: its children are mandatory, which is what
    # makes a missing activity a model move - deficit - rather than a
    # branch the aligner simply does not take.
    return st.one_of(
        st.tuples(st.sampled_from(['seq', 'seq', 'xor', 'and']),
                  st.lists(children, min_size=2, max_size=3).map(tuple)),
        st.tuples(st.just('loop'), st.lists(children, min_size=2, max_size=2).map(tuple)),
    )


# Rooted at an operator: a lone leaf admits only the trivial cut.
SHAPES = st.recursive(_LEAF, _operator, max_leaves=5).filter(lambda s: s[0] in _OPERATORS)


def _count_activities(shape):
    kind, children = shape
    if kind == 'act':
        return 1
    return sum(_count_activities(child) for child in children)


def build_tree(shape, labels):
    """A skip-alignments process tree for `shape`, naming activity leaves
    from `labels` in order and giving every node a distinct id."""
    ids = (str(i) for i in itertools.count())
    taus = (f'tau{i}' for i in itertools.count())

    def build(shape):
        kind, children = shape
        if kind == 'act':
            node = Activity(None, next(labels), COST)
        elif kind == 'tau':
            node = Tau(None, next(taus), 0)
        else:
            built = [build(child) for child in children]
            node = _OPERATORS[kind](None, built)
            for child in built:
                child.set_parent(node)
        node.id = next(ids)
        return node

    return build(shape)


def children_of(node):
    return getattr(node, 'children', None) or []


def is_silent(node):
    return isinstance(node, Tau)


@st.composite
def cuts(draw, tree):
    """
    An antichain covering every labelled leaf: at each node below the
    root, either take it whole or descend into its children. The root is
    always descended from, since [root] alone sums to the root
    trivially. A silent leaf reached by descending may be left out
    altogether, since the cover need only reach the labelled leaves.
    """
    def cut(node, whole_allowed):
        kids = children_of(node)
        if not kids:
            if is_silent(node) and draw(st.booleans()):
                return []
            return [node]
        if whole_allowed and draw(st.booleans()):
            return [node]
        return [member for child in kids for member in cut(child, True)]
    return cut(tree, whole_allowed=False)


def _traces(alphabet):
    """
    Traces of two kinds. An in-order subset of the alphabet - the model's
    activities in leaf order, some left out - makes a left-out mandatory
    activity a model move, so deficits are common. A noisy draw adds
    reordering, repeats and an out-of-alphabet activity, whose log moves
    count at no node.
    """
    in_order = st.lists(st.booleans(), min_size=len(alphabet), max_size=len(alphabet)).map(
        lambda keep: [a for a, k in zip(alphabet, keep) if k]).filter(bool)
    noisy = st.lists(st.sampled_from(alphabet + [OUT_OF_ALPHABET]), min_size=1, max_size=4)
    return st.lists(st.one_of(in_order, noisy), min_size=1, max_size=3)


def variant_probs_of(traces):
    counts = {}
    for trace in traces:
        counts[tuple(trace)] = counts.get(tuple(trace), 0) + 1
    return {variant: n / len(traces) for variant, n in counts.items()}


def voidmass_process_table(tree, traces):
    net, im, fm, activity_to_id, tau_ids, id_loop_list = build_id_net(tree)
    result = voidmass_table_pn(tree, variant_probs_of(traces), net, im, fm,
                               activity_to_id, tau_ids, id_loop_list=id_loop_list,
                               timeout=30)
    assume(result.timed_out_count == 0)
    return result.table


@st.composite
def models_logs_and_cuts(draw):
    shape = draw(SHAPES)
    n_activities = _count_activities(shape)
    assume(n_activities >= 1)
    alphabet = [chr(ord('a') + i) for i in range(n_activities)]
    tree = build_tree(shape, iter(alphabet))
    traces = draw(_traces(alphabet))
    return tree, traces, draw(cuts(tree))


class RepeatedLabelsKnownFailureTest(unittest.TestCase):
    """
    The lemma holds for the definition, and for this implementation only
    where leaf labels are distinct - which is why the property above
    draws them so.

    The definition attributes each move to exactly one execution.
    voidmass_pn.terms_by_node instead maps a move back to its activity
    LABEL and credits every node whose leaves carry that label, so two
    leaves sharing a label both claim every move on it. Hypothesis shrank
    the first counterexample to this: seq(a, a) against <a>. One a
    synchronises and the other is a model move, so the root reads deficit
    1 over 2 moves, 0.5; by the definition the model move belongs to one
    leaf and the two leaves sum to 0.5, but each is credited with it and
    they sum to 1.0.

    Pinned at what the code reads today, so the failure stays visible and
    the test turns red if attribution changes.
    """

    def test_seq_a_a_against_a_counts_the_missing_a_at_both_leaves(self):
        first = build_tree(('act', ()), iter(['a']))
        second = build_tree(('act', ()), iter(['a']))
        tree = Sequence(None, [first, second])
        first.id, second.id, tree.id = 'first', 'second', 'root'
        first.set_parent(tree)
        second.set_parent(tree)

        table = voidmass_process_table(tree, [['a']])

        self.assertAlmostEqual(table[tree]['voidmass_process_lower'], 0.5)
        self.assertAlmostEqual(table[first]['voidmass_process_lower'], 0.5)
        self.assertAlmostEqual(table[second]['voidmass_process_lower'], 0.5)


class AdditivityTest(unittest.TestCase):

    @settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    @given(models_logs_and_cuts())
    def test_values_over_a_covering_antichain_sum_to_the_roots(self, case):
        tree, traces, cut = case
        table = voidmass_process_table(tree, traces)
        total = sum(table[node]['voidmass_process_lower'] for node in cut)
        self.assertAlmostEqual(total, table[tree]['voidmass_process_lower'],
                               delta=TOLERANCE,
                               msg=f'tree {tree}, traces {traces}, cut {[n.id for n in cut]}')

    @settings(max_examples=30, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    @given(models_logs_and_cuts(), st.data())
    def test_a_cut_missing_a_labelled_leaf_with_a_deficit_falls_short(self, case, data):
        """
        The cover condition is necessary, and the test above has teeth:
        leave out a labelled leaf that carries deficit and the sum is
        short of the root. Were the generator producing only perfectly
        aligned logs, hypothesis could find no such leaf and would fail
        this test as unsatisfiable rather than let the property pass on
        zeros.
        """
        tree, traces, cut = case
        table = voidmass_process_table(tree, traces)
        leaves = [node for node in cut if not children_of(node) and not is_silent(node)]
        assume(leaves)
        dropped = data.draw(st.sampled_from(leaves))
        assume(table[dropped]['voidmass_process_lower'] > 0)
        total = sum(table[node]['voidmass_process_lower'] for node in cut if node is not dropped)
        self.assertLess(total, table[tree]['voidmass_process_lower'] - TOLERANCE)

    @settings(max_examples=30, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    @given(models_logs_and_cuts())
    def test_without_a_timeout_the_bounds_coincide(self, case):
        """So _lower above is the value itself, not merely a bound."""
        tree, traces, _cut = case
        table = voidmass_process_table(tree, traces)
        for node, row in table.items():
            self.assertAlmostEqual(row['voidmass_process_lower'],
                                   row['voidmass_process_upper'], delta=TOLERANCE)


if __name__ == '__main__':
    unittest.main()
