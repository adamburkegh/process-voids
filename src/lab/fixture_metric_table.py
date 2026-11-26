"""
Every metric lab.run computes, at every node of each hand-built fixture.

    python -m lab.fixture_metric_table --out var/lab/fixture_metrics.csv
    python -m lab.fixture_metric_table --fixtures payment_partial
    python -m lab.fixture_metric_table --diff before.csv after.csv

The fixture-scale counterpart to a run: the same ProcessMetric roster and
the same CellContext staging as lab.run, over the small trees in
lab.fixtures rather than a discovered one. Discovery is not reproducible
across processes with PYTHONHASHSEED unset, so these trees are what two
runs of different code can be compared on.

That is what it is for - checking whether a change to the shared metric
machinery moves a value anywhere, including at the nodes no test pins.
Run it before and after the change and diff the two CSVs.
"""

import argparse

import pandas as pd

from lab.fixtures import (
    build_partial_sequence_log, build_partial_sequence_tree,
    build_payment_partial_log, build_payment_partial_tree,
    build_running_example_log, build_running_example_tree,
)
from lab.run import ALL_METRICS
from process_voids.metric_context import CellContext, METRIC_ERROR, score_all
from process_voids.voidmass_pn import build_id_net

FIXTURES = {
    'payment_approval': (build_running_example_log, build_running_example_tree),
    'payment_partial': (build_payment_partial_log, build_payment_partial_tree),
    'partial_sequence': (build_partial_sequence_log, build_partial_sequence_tree),
}


def _nodes(tree):
    yield tree
    for child in getattr(tree, 'children', []) or []:
        yield from _nodes(child)


def fixture_rows(name, metrics=ALL_METRICS, slpn_dir='var/lab'):
    """One row per node of fixture `name`, through the real pipeline."""
    build_log, build_tree = FIXTURES[name]
    log, tree = build_log(), build_tree()
    ctx = CellContext(log=log, tree=tree,
                      slpn_path=f'{slpn_dir}/fixture_metric_table_{name}.slpn',
                      classical_net=build_id_net(tree))
    for stage_id in {stage for metric in metrics for stage in metric.needs}:
        ctx.stage(stage_id)
    rows = []
    for node in _nodes(tree):
        values = score_all(ctx, metrics, node=node)
        rows.append({'fixture': name, 'node_id': node.id,
                     'node_type': type(node).__name__,
                     **{metric_id: (None if value is METRIC_ERROR else value)
                        for metric_id, value in values.items()}})
    return rows


KEY = ['fixture', 'node_id']


def diff_rows(before, after, tolerance=1e-12):
    """
    One row per (fixture, node, metric) whose value differs between two
    tables of this shape - the answer a before-and-after pair is run
    for. A metric missing from either table is reported as such rather
    than silently skipped.
    """
    metrics = [column for column in after.columns
               if column not in KEY and column != 'node_type']
    first, second = before.set_index(KEY), after.set_index(KEY)
    rows = []
    for metric in metrics:
        if metric not in first.columns:
            rows.append({'metric': metric, 'note': 'absent before'})
            continue
        for key in second.index:
            was, now = first.at[key, metric], second.at[key, metric]
            if pd.isna(was) and pd.isna(now):
                continue
            if pd.notna(was) and pd.notna(now) and abs(now - was) <= tolerance:
                continue
            fixture, node_id = key
            rows.append({'fixture': fixture, 'node_id': node_id,
                         'node_type': second.at[key, 'node_type'],
                         'metric': metric, 'before': was, 'after': now,
                         'change': now - was, 'note': ''})
    return pd.DataFrame(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--fixtures', nargs='+', default=sorted(FIXTURES),
                        choices=sorted(FIXTURES))
    parser.add_argument('--out', default='var/lab/fixture_metrics.csv',
                        help='CSV to write the per-node rows to')
    parser.add_argument('--diff', nargs=2, metavar=('BEFORE', 'AFTER'),
                        help='compare two tables already written, and compute nothing')
    args = parser.parse_args(argv)

    if args.diff:
        moved = diff_rows(pd.read_csv(args.diff[0]), pd.read_csv(args.diff[1]))
        print('no (node, metric) value differs' if moved.empty
              else moved.to_string(index=False))
        return

    rows = [row for name in args.fixtures for row in fixture_rows(name)]
    frame = pd.DataFrame(rows)
    frame.to_csv(args.out, index=False)
    print(frame.to_string(index=False))
    print(f'wrote {args.out} ({len(rows)} rows)')


if __name__ == '__main__':
    main()
