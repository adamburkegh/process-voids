"""
How far the void metrics move when the skip-alignment normal form's move
precedence is reversed.

skip-alignments returns one representative of each class of transposable
alignments, ordered log, then synchronous, then skip. Every optimal
alignment reachable by transposing INDEPENDENT adjacent moves
(process_voids.normal_form) is equally valid, so the precedence is a
choice the implementation makes rather than something the evidence
decides. This reorders each returned alignment to the reverse precedence
and rescores, which answers whether a metric depends on that choice.

    python -m lab.normal_form_compare
    python -m lab.normal_form_compare --fixtures payment_partial --out var/nf.csv

voidmass_process is not compared: it is computed from classical
alignments (process_voids.voidmass_pn builds its own), and never reads a
skip-alignment normal form, so there is no second value to compare
against. test_normal_form_invariance.py pins that it takes no
skip-alignment input.
"""

import argparse
from types import SimpleNamespace

import pandas as pd

from lab.fixtures import (
    build_partial_sequence_log, build_partial_sequence_tree,
    build_payment_partial_log, build_payment_partial_tree,
    build_running_example_log, build_running_example_tree,
)
from process_voids import pvoid
from lab.normal_form import LOG_FIRST, SKIP_FIRST, independent, move_kind, reorder
from process_voids.voidsalign3 import voidsalign3
from process_voids.voidsat2 import voidsat2

FIXTURES = {
    'payment_approval': (build_running_example_log, build_running_example_tree),
    'payment_partial': (build_payment_partial_log, build_payment_partial_tree),
    'partial_sequence': (build_partial_sequence_log, build_partial_sequence_tree),
}


def _nodes(tree):
    yield tree
    for child in getattr(tree, 'children', []) or []:
        yield from _nodes(child)


def independent_pairs(paths_by_variant):
    """
    Every (variant key, position, kinds) where two adjacent moves are
    independent and of different kinds - the transpositions that could
    move a value at all. An empty list means a fixture cannot show a
    difference, whatever the metrics do.
    """
    found = []
    for key, paths in paths_by_variant.items():
        for path in paths:
            for i in range(len(path) - 1):
                first, second = path[i], path[i + 1]
                if move_kind(first) == move_kind(second):
                    continue
                if independent(first, second):
                    found.append((key, i, (move_kind(first), move_kind(second))))
    return found


def reordered(paths_by_variant, tree, precedence):
    return {key: [reorder(path, tree, precedence) for path in paths]
            for key, paths in paths_by_variant.items()}


def _as_states(paths_by_variant):
    """voidsalign3 reads a skip_dict of objects carrying .path."""
    return {key: [SimpleNamespace(path=path) for path in paths]
            for key, paths in paths_by_variant.items()}


def comparison_rows(tree, log, paths_by_variant, variant_probs, skip_probs):
    """
    One row per (node, metric): its value under each precedence and the
    change. `paths_by_variant` is the aligner's own output, keyed by
    coveragemass._variant_key.
    """
    under = {}
    for name, precedence in (('log_first', LOG_FIRST), ('skip_first', SKIP_FIRST)):
        paths = reordered(paths_by_variant, tree, precedence)
        states = _as_states(paths)
        under[name] = {
            'voidsalign3': {node: voidsalign3(node, states, variant_probs, skip_probs)
                            for node in _nodes(tree)},
            'voidsat2': {node: voidsat2(node, tree, log, paths, skip_probs)
                         for node in _nodes(tree)},
        }

    rows = []
    for node in _nodes(tree):
        for metric in ('voidsalign3', 'voidsat2'):
            first = under['log_first'][metric][node]
            second = under['skip_first'][metric][node]
            rows.append({'node_id': node.id, 'node_type': type(node).__name__,
                         'is_root': node is tree, 'metric': metric,
                         'log_first': first, 'skip_first': second,
                         'change': second - first})
    return pd.DataFrame(rows)


def compare_fixture(name, slpn_dir='var/lab'):
    """(rows, independent pairs) for one fixture, through the real pipeline."""
    build_log, build_tree = FIXTURES[name]
    log, tree = build_log(), build_tree()
    dv = pvoid.skipprob(log, tree, f'{slpn_dir}/normal_form_compare_{name}.slpn')
    paths = {key: [state.path for state in states]
             for key, states in dv.skip_dict_backup.items()}
    rows = comparison_rows(tree, log, paths, dv.pl, dv.skip_probs)
    return rows, independent_pairs(paths)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--fixtures', nargs='+', default=sorted(FIXTURES),
                        choices=sorted(FIXTURES))
    parser.add_argument('--out', help='CSV to write the per-node rows to')
    args = parser.parse_args(argv)

    frames = []
    for name in args.fixtures:
        rows, pairs = compare_fixture(name)
        frames.append(rows.assign(fixture=name))
        moved = rows[rows['change'].abs() > 1e-12]
        print(f'{name}: {len(pairs)} independent adjacent pair(s) of differing kind; '
              f'{len(moved)} of {len(rows)} (node, metric) values move')
        for _, row in moved.iterrows():
            print(f'    {row.metric} at node {row.node_id} ({row.node_type}): '
                  f'{row.log_first:.6f} -> {row.skip_first:.6f} ({row.change:+.6f})')
    if args.out:
        pd.concat(frames).to_csv(args.out, index=False)
        print(f'wrote {args.out}')


if __name__ == '__main__':
    main()
