"""
Experiment: interval surprise - a coverage-void metric needing no
alignments and no tunable parameters (process_voids.surprise).

A tree is discovered once per (log, combo) from the original, undegraded
log - the same fixed-reference-model design as exp_disco_degrade.py.
Each degradation dimension x level then degrades the log and recomputes
surprise (both attribution schemes) against that same fixed tree,
reporting per-node totals alongside root-level headline figures.

Usage:
    python -m lab.exp_surprise <log_path> [<log_path> ...] \\
        [--combos inductive toothpaste] [--degradations activity trace] \\
        [--levels 0.0 0.1 0.2]
"""

import argparse
import logging
import pickle
import time
from pathlib import Path

import pandas as pd
import pm4py_config as pm4py

from lab.degradation import DEGRADATIONS
from lab.discovery import COMBOS
from lab.logconfig import configure
from process_voids.coveragemass import log_to_traces
from process_voids.surprise import (
    event_surprise, observed_intervals, surprise_totals,
    compute_predecessors, predecessor_totals,
)

logger = logging.getLogger(__name__)


def _format_dropped(dropped, limit=50):
    """(dropped_str, dropped_count) - see exp_disco_degrade._format_dropped."""
    items = sorted(str(d) for d in dropped)
    if len(items) > limit:
        return f'{len(items)} items', len(items)
    return ', '.join(items), len(items)


CELL_COLS = ['log', 'combo', 'degradation_dim', 'degradation_level']

# The metric ids this experiment emits, per (distribution in {self,
# baseline}) - see lab.metric_registry, whose drift test imports these
# directly rather than re-deriving them from the CSV output.
NODE_METRIC_KEYS = ('containment_bits', 'predecessor_bits')
SUMMARY_METRIC_KEYS = ('headline_bits', 'bits_per_event')


def _merge_write(df, path, cell_cols=CELL_COLS):
    """
    Upsert `df` into the CSV at `path` by `cell_cols`: existing rows on
    disk whose cell (log, combo, dim, level) matches a row in `df` are
    dropped entirely and replaced by df's rows for that cell, everything
    else on disk is kept, and the combined result is written back.

    Purging by cell rather than by each row's own full key matters
    because a cell's rows are produced atomically (both distribution
    variants together, or a single error row with no distribution set)
    - keying by the full row would leave a stale error row behind after
    a rerun succeeds (its distribution is NaN, so it never matches the
    new self/baseline rows' keys), or leave orphaned node rows behind
    for node_ids that no longer exist if the tree changed between runs.

    This is the only way results ever reach disk here - callers never
    need to juggle separate output paths or merge runs by hand to add a
    combo/level without re-running (and possibly clobbering or
    duplicating) what was already there.

    An empty df (nothing computed this run - e.g. every cell errored)
    has no columns to key by, so there is nothing to purge or add -
    leave whatever's on disk untouched rather than crash on cell_cols.
    """
    path = Path(path)
    if df.empty:
        return pd.read_csv(path) if path.exists() else df
    cells = set(df[cell_cols].fillna('').apply(tuple, axis=1))
    if path.exists():
        existing = pd.read_csv(path)
        existing_cells = existing[cell_cols].fillna('').apply(tuple, axis=1)
        existing = existing[~existing_cells.isin(cells)]
        combined = pd.concat([existing, df], ignore_index=True)
    else:
        combined = df
    path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(path, index=False)
    return combined


TREE_CACHE_DIR = Path('var/lab/tree_cache')


def _discover_cached(log_name, combo_name, combo, base_log):
    """
    The discovered tree depends only on (log, combo), never on
    degradation dim/level - cache it to disk so a separate run against
    the same log doesn't pay discovery's cost again (toothpaste's
    external-subprocess discovery in particular can be slow). Keyed by
    filename only, not log content - delete the cache file (or the
    whole var/lab/tree_cache/ dir) if the underlying log changes.
    """
    cache_path = TREE_CACHE_DIR / f'{log_name}__{combo_name}.pkl'
    if cache_path.exists():
        with open(cache_path, 'rb') as f:
            return pickle.load(f)
    tree = combo.discover(base_log).tree
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, 'wb') as f:
        pickle.dump(tree, f)
    return tree


def _node_rows(log_name, combo_name, dim, level, distribution, containment, pred_totals):
    rows = []
    for node in set(containment) | set(pred_totals):
        metric_values = (containment.get(node, 0.0), pred_totals.get(node, 0.0))
        rows.append({
            'log': log_name,
            'combo': combo_name,
            'degradation_dim': dim,
            'degradation_level': level,
            'distribution': distribution,
            'node_id': node.id,
            'node_type': type(node).__name__,
            'alphabet': ','.join(sorted(set(node.get_leaf_labels()))),
            **dict(zip(NODE_METRIC_KEYS, metric_values)),
        })
    return rows


def _compute_variant(tree, predecessors, traces, obs):
    rows = event_surprise(traces, obs=obs)
    containment, out_of_alphabet_count = surprise_totals(tree, rows)
    pred_totals, ambiguous_count, unattributable_count = \
        predecessor_totals(tree, rows, predecessors)
    n_events = len(rows)
    headline = containment[tree]
    summary_values = (headline, headline / n_events if n_events else None)
    return (containment, pred_totals), {
        'n_events': n_events,
        **dict(zip(SUMMARY_METRIC_KEYS, summary_values)),
        'ambiguous_event_count': ambiguous_count,
        'unattributable_event_count': unattributable_count,
        'out_of_alphabet_event_count': out_of_alphabet_count,
    }


def _compute_cell(tree, predecessors, log, base_obs=None):
    """
    {distribution -> (node_totals, summary_fields)} for one (tree, log)
    pairing. 'self' always estimates the tail distribution from `log`
    itself (the metric's normal, self-contained mode - see surprise.py's
    module docstring on why this is self-limiting under heavy loss).
    When base_obs is given (the distribution estimated from the
    *undegraded* log), an additional 'baseline' variant reruns the same
    degraded events against that fixed, uncontaminated distribution -
    isolating whether the metric's response to degradation is real
    signal or an artifact of the estimator degrading along with the log.
    """
    traces = log_to_traces(log)
    variants = {'self': _compute_variant(tree, predecessors, traces, None)}
    if base_obs is not None:
        variants['baseline'] = _compute_variant(tree, predecessors, traces, base_obs)
    return variants


def run_surprise(log_paths, combos=COMBOS, degradations=DEGRADATIONS, levels=(0.0,),
                  out_csv='var/lab/results/exp_surprise.csv',
                  summary_csv='var/lab/results/exp_surprise_summary.csv'):
    logger.info('Experiment: surprise | logs=%s | combos=%s | degradations=%s | levels=%s',
                [Path(p).stem for p in log_paths], list(combos), list(degradations), levels)

    node_rows = []
    summary_rows = []
    for log_path in log_paths:
        log_name = Path(log_path).stem
        base_log = pm4py.read_xes(log_path)
        base_obs = observed_intervals(log_to_traces(base_log))

        for combo_name, combo in combos.items():
            cell = f'{log_name} / {combo_name}'
            started = time.monotonic()
            try:
                tree = _discover_cached(log_name, combo_name, combo, base_log)
            except NotImplementedError:
                for dim in degradations:
                    for level in levels:
                        summary_rows.append({
                            'log': log_name, 'combo': combo_name, 'degradation_dim': dim,
                            'degradation_level': level, 'status': 'not_implemented'})
                logger.info('%s - skipped (not_implemented)', cell)
                continue
            except Exception as e:
                for dim in degradations:
                    for level in levels:
                        summary_rows.append({
                            'log': log_name, 'combo': combo_name, 'degradation_dim': dim,
                            'degradation_level': level, 'status': f'discovery error: {e}'})
                logger.warning('%s - discovery error: %s', cell, e)
                continue
            logger.info('%s - discovered in %.1fs', cell, time.monotonic() - started)

            predecessors = compute_predecessors(tree)

            # level 0.0 drops nothing regardless of dimension, so it's the
            # same (tree, log) computation under every dim - compute it
            # once here and reuse across dims, mirroring exp_disco_degrade.
            zero_level_variants = None
            if 0.0 in levels:
                zero_level_variants = _compute_cell(tree, predecessors, base_log,
                                                     base_obs=base_obs)

            for dim, degrade in degradations.items():
                for level in levels:
                    sub_cell = f'{cell} / {dim} / {level}'
                    started = time.monotonic()
                    try:
                        if level == 0.0 and zero_level_variants is not None:
                            variants = zero_level_variants
                            dropped_str, dropped_count = '', 0
                        else:
                            degraded_log, dropped = degrade(base_log, level)
                            dropped_str, dropped_count = _format_dropped(dropped)
                            variants = _compute_cell(tree, predecessors, degraded_log,
                                                      base_obs=base_obs)

                        for distribution, ((containment, pred_totals), fields) in variants.items():
                            node_rows.extend(_node_rows(log_name, combo_name, dim, level,
                                                         distribution, containment, pred_totals))
                            summary_rows.append({
                                'log': log_name, 'combo': combo_name, 'degradation_dim': dim,
                                'degradation_level': level, 'distribution': distribution,
                                'dropped_count': dropped_count, 'status': 'ok', **fields,
                            })
                        self_fields = variants['self'][1]
                        logger.info('%s - done in %.1fs (%d events, %.2f bits)',
                                    sub_cell, time.monotonic() - started,
                                    self_fields['n_events'], self_fields['headline_bits'])
                    except Exception as e:
                        summary_rows.append({
                            'log': log_name, 'combo': combo_name, 'degradation_dim': dim,
                            'degradation_level': level, 'status': f'error: {e}'})
                        logger.warning('%s - error: %s', sub_cell, e)

    node_df = _merge_write(pd.DataFrame(node_rows), out_csv)
    summary_df = _merge_write(pd.DataFrame(summary_rows), summary_csv)
    return node_df, summary_df


def top_nodes(node_df, column, n=15):
    """The n rows with the highest `column` value, one row per (log, combo, dim, level, distribution)."""
    if node_df.empty:
        return node_df
    return (node_df.sort_values(column, ascending=False)
                    .groupby(['log', 'combo', 'degradation_dim', 'degradation_level',
                              'distribution'], sort=False)
                    .head(n))


def format_node_table(df, column):
    cols = ['log', 'combo', 'degradation_dim', 'degradation_level', 'distribution',
            'node_id', 'node_type', 'alphabet', column]
    return df[cols].to_string(index=False) if not df.empty else '  (no nodes)'


def print_report(node_df, summary_df, top_n=15):
    print('SUMMARY')
    print(summary_df.to_string(index=False))

    print(f'\nTOP {top_n} NODES BY CONTAINMENT')
    print(format_node_table(top_nodes(node_df, 'containment_bits', top_n), 'containment_bits'))

    print(f'\nTOP {top_n} NODES BY PREDECESSOR ATTRIBUTION')
    print(format_node_table(top_nodes(node_df, 'predecessor_bits', top_n), 'predecessor_bits'))


def main():
    configure()
    logger.info('Starting exp_surprise')
    parser = argparse.ArgumentParser(
        description='Interval surprise experiment: per-node containment and '
                     'predecessor attribution for a discovered tree against its (possibly '
                     'degraded) log.')
    parser.add_argument('logs', nargs='+', help='XES log path(s)')
    parser.add_argument('--combos', nargs='+', default=list(COMBOS),
                         help=f'Discovery combos to use (default: all of {list(COMBOS)})')
    parser.add_argument('--degradations', nargs='+', default=list(DEGRADATIONS),
                         help=f'Degradation dimensions to use (default: all of {list(DEGRADATIONS)})')
    parser.add_argument('--levels', type=float, nargs='+', default=[0.0],
                         help='Degradation levels in [0,1] (default: [0.0], ie undegraded)')
    parser.add_argument('--out', default='var/lab/results/exp_surprise.csv')
    parser.add_argument('--summary-out', default='var/lab/results/exp_surprise_summary.csv')
    parser.add_argument('--top', type=int, default=15,
                         help='How many top nodes to print per (log, combo, dim, level) (default: 15)')
    args = parser.parse_args()

    combos = {name: COMBOS[name] for name in args.combos}
    degradations = {name: DEGRADATIONS[name] for name in args.degradations}
    node_df, summary_df = run_surprise(args.logs, combos=combos, degradations=degradations,
                                        levels=args.levels, out_csv=args.out,
                                        summary_csv=args.summary_out)

    # node_df/summary_df are the full accumulated files (merge_write returns
    # everything on disk, not just what this call computed) - scope the
    # printed report to what was actually requested this run, so it doesn't
    # balloon with every prior run's results as the files accumulate.
    log_names = [Path(p).stem for p in args.logs]
    def _requested(df):
        return df[df['log'].isin(log_names) & df['combo'].isin(combos) &
                   df['degradation_dim'].isin(degradations) & df['degradation_level'].isin(args.levels)]

    print_report(_requested(node_df), _requested(summary_df), top_n=args.top)
    logger.info('Wrote %s and %s', args.out, args.summary_out)


if __name__ == '__main__':
    main()
