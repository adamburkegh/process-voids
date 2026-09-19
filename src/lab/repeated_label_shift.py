"""
How far classical per-node voidmass_process moves on a tree whose leaves
repeat activity labels, between the label-based computation and the
leaf-based one.

Label-based is what process_voids.voidmass_pn did before the
repeated-label fix: a move was credited to every node whose leaves carry
its label, and tied optimal alignments were deduplicated by a signature
of (kind, label), which merges alignments differing only in which
same-labelled leaf fired. Leaf-based credits a move to the one leaf
whose transition fired, and deduplicates by (kind, leaf). Where no label
repeats, the two agree at every node.

Both are computed here from one pass of alignments, independently of
voidmass_table_pn, so the comparison holds before and after the fix
lands. Only the no-timeout value is compared: a variant whose alignment
search times out is counted and left out of both.

    python -m lab.repeated_label_shift <log> <tree-cache.pkl>
    python -m lab.repeated_label_shift <log> <tree-cache.pkl> \
        --degradations activity_gradual trace --levels 0.3 0.5

The tree is read from a tree-cache pair file, so a result can be
checked against the exact tree it was scored against. A discovered tree
typically fits its own undegraded log, leaving no model moves to
misattribute, so the shift shows under degradation - the same
lab.degradation dimensions and levels the sweeps use.
"""

import argparse
import os
import pickle
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import pm4py_config as pm4py

from lab.degradation import DEGRADATIONS
from process_voids.metric_context import _variant_probs as _pipeline_variant_probs
from process_voids.voidmass_pn import align_variant_all, build_id_net, leaf_by_transition

LABEL = 'label'
LEAF = 'leaf'


def _is_silent(node):
    return type(node).__name__ == 'Tau'


def counted_moves(alignment, leaves):
    """(kind, leaf) for every synchronous or model move on an activity
    leaf - the moves deficit and movecount count."""
    moves = []
    for t in alignment:
        trace_side, model_side = t.label
        leaf = leaves.get(t.name[1])
        if model_side == '>>' or leaf is None or _is_silent(leaf):
            continue
        moves.append(('model' if trace_side == '>>' else 'sync', leaf))
    return moves


def signature(alignment, leaves, by):
    """A tied alignment's causal story, by LABEL or by LEAF."""
    if by == LABEL:
        return tuple((kind, leaf.name) for kind, leaf in counted_moves(alignment, leaves))
    return tuple((kind, leaf.id) for kind, leaf in counted_moves(alignment, leaves))


def _nodes(tree):
    yield tree
    for child in tree.children:
        yield from _nodes(child)


@dataclass
class Shift:
    """{node: (deficit, movecount)} under each computation, plus the
    variants left out because their search timed out."""
    label: dict = field(default_factory=dict)
    leaf: dict = field(default_factory=dict)
    timed_out: int = 0

    def voidmass_process(self, by, node, tree):
        table = self.label if by == LABEL else self.leaf
        root_movecount = table[tree][1]
        return table[node][0] / root_movecount if root_movecount else 0.0


def shift_table(tree, variant_probs, timeout=30):
    """Shift for `tree` over `variant_probs` ({activity tuple: weight})."""
    net, im, fm, activity_to_id, tau_ids, loops = build_id_net(tree)
    leaves = leaf_by_transition(net, tree)
    nodes = list(_nodes(tree))
    subtree = {node: {id(n) for n in _nodes(node)} for node in nodes}
    labels = {node: {n.name for n in _nodes(node) if not n.children} for node in nodes}
    result = Shift(label={n: (0.0, 0.0) for n in nodes}, leaf={n: (0.0, 0.0) for n in nodes})

    for variant, weight in variant_probs.items():
        raw = align_variant_all(list(variant), net, im, fm, activity_to_id, tau_ids,
                                id_loop_list=loops, timeout=timeout)
        if not raw:
            result.timed_out += 1
            continue
        for by, table in ((LABEL, result.label), (LEAF, result.leaf)):
            stories = {}
            for alignment in raw:
                stories.setdefault(signature(alignment, leaves, by), alignment)
            share = weight / len(stories)
            for alignment in stories.values():
                moves = counted_moves(alignment, leaves)
                for node in nodes:
                    if by == LABEL:
                        mine = [kind for kind, leaf in moves if leaf.name in labels[node]]
                    else:
                        mine = [kind for kind, leaf in moves if id(leaf) in subtree[node]]
                    d, m = table[node]
                    table[node] = (d + share * mine.count('model'), m + share * len(mine))
    return result


def variant_probs_of_log(log):
    """
    {activity tuple: probability} for a log DataFrame - the pipeline's own
    construction, not a reimplementation. Where events in a case share a
    timestamp, their order depends on how the sort is done, so a
    reimplementation can build different variants from the same log; this
    has to match the values it is compared against exactly.
    """
    return _pipeline_variant_probs(log)


def shift_rows(tree, shift):
    """One row per node: its value under each computation and the change."""
    rows = []
    for node in _nodes(tree):
        old = shift.voidmass_process(LABEL, node, tree)
        new = shift.voidmass_process(LEAF, node, tree)
        rows.append({'node_id': node.id, 'node_type': type(node).__name__,
                     'is_root': node is tree, 'label_based': old, 'leaf_based': new,
                     'change': new - old})
    return pd.DataFrame(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('log', help='XES log path')
    parser.add_argument('tree_cache', help='tree-cache pair file (.pkl)')
    parser.add_argument('--degradations', nargs='+', default=[],
                        help=f'lab.degradation dimensions (choices: {sorted(DEGRADATIONS)}); '
                             'none compares the undegraded log only')
    parser.add_argument('--levels', type=float, nargs='+', default=[0.3, 0.5, 0.7])
    parser.add_argument('--timeout', type=int, default=30,
                        help='per-variant alignment search budget, seconds')
    parser.add_argument('--out', help='CSV to write the per-node rows to')
    args = parser.parse_args(argv)

    with open(args.tree_cache, 'rb') as f:
        tree, _ppt_weights = pickle.load(f)
    base = pm4py.read_xes(os.fspath(args.log))
    cells = [(None, 0.0, base)] + [
        (dim, level, DEGRADATIONS[dim](base, level)[0])
        for dim in args.degradations for level in args.levels]

    print(f'{Path(args.tree_cache).name}')
    frames = []
    for dim, level, log in cells:
        shift = shift_table(tree, variant_probs_of_log(log), timeout=args.timeout)
        rows = shift_rows(tree, shift).assign(degradation_dim=dim, degradation_level=level)
        frames.append(rows)
        root = rows[rows['is_root']].iloc[0]
        moved = rows[rows['change'].abs() > 1e-9]
        print(f'  {dim or "undegraded"} {level:.1f}: root {root.label_based:.4f} -> '
              f'{root.leaf_based:.4f} ({root.change:+.4f}); {len(moved)} of {len(rows)} '
              f'nodes move, largest {moved["change"].abs().max() if len(moved) else 0.0:.4f}; '
              f'{shift.timed_out} variants timed out')
    if args.out:
        pd.concat(frames).to_csv(args.out, index=False)
        print(f'  wrote {args.out}')


if __name__ == '__main__':
    main()
