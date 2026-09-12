"""
Experiment: interval surprise - a coverage-void metric needing no
alignments and no tunable parameters (process_voids.surprise).

A tree is discovered once per (log, combo) from the original, undegraded
log - the same fixed-reference-model design as exp_disco_degrade.py.
Each degradation dimension x level then degrades the log and recomputes
surprise (both attribution schemes) against that same fixed tree,
reporting per-node totals alongside root-level headline figures. Each
metric is reported under two separate ids/columns, not a shared id plus
a distribution column: unsuffixed ('containment_bits', 'predecessor_bits',
'headline_bits', 'bits_per_event') estimates the tail distribution from
the log being scored ('self', the metric's normal, deployable mode);
the '_baseline' suffix reruns the same events against the *undegraded*
log's distribution instead - a benchmark-only oracle comparison not
available under real missing-activity scenarios, since no undegraded
reference log exists there. See NODE_METRIC_KEYS/NODE_METRIC_BASELINE_KEYS/
SUMMARY_METRIC_KEYS/SUMMARY_METRIC_BASELINE_KEYS below and
_compute_cell's docstring.

Usage:
    python -m lab.exp_surprise <log_path> [<log_path> ...] \\
        [--combos inductive toothpaste] [--degradations activity trace] \\
        [--levels 0.0 0.1 0.2]
"""

import argparse
import logging
import time
from pathlib import Path

import pandas as pd
import pm4py_config as pm4py

from lab.degradation import DEGRADATIONS
from lab.discovery import COMBOS, discover_cached
from lab.logconfig import configure
from lab.timing import Timer
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

# The metric ids this experiment emits - see lab.metric_registry, whose
# drift test imports these directly rather than re-deriving them from
# the CSV output. 'self' (unsuffixed) and 'baseline' (_baseline
# suffix) are separate ids, not a shared id plus a distribution column
# - they're different quantities (self is the deployable metric,
# baseline is a benchmark-only oracle comparison against the undegraded
# log's distribution - see _compute_cell), not two readings of the same
# one.
NODE_METRIC_KEYS = ('containment_bits', 'predecessor_bits')
NODE_METRIC_BASELINE_KEYS = ('containment_bits_baseline', 'predecessor_bits_baseline')
SUMMARY_METRIC_KEYS = ('headline_bits', 'bits_per_event')
SUMMARY_METRIC_BASELINE_KEYS = ('headline_bits_baseline', 'bits_per_event_baseline')


def _merge_write(df, path, cell_cols=CELL_COLS):
    """
    Upsert `df` into the CSV at `path` by `cell_cols`: existing rows on
    disk whose cell (log, combo, dim, level) matches a row in `df` are
    dropped entirely and replaced by df's rows for that cell, everything
    else on disk is kept, and the combined result is written back.

    Purging by cell rather than by each row's own full key matters
    because a cell's rows are produced atomically (self and baseline
    columns together per node, or a single error row with everything
    else NaN) - keying by the full row would leave a stale error row
    behind after a rerun succeeds (its metric columns are NaN, so it
    never matches the new row's key), or leave orphaned node rows
    behind for node_ids that no longer exist if the tree changed
    between runs.

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


def _node_rows(log_name, combo_name, dim, level, self_totals, baseline_totals):
    """One row per tree node, self and baseline bits as separate columns
    (NODE_METRIC_KEYS / NODE_METRIC_BASELINE_KEYS) rather than a shared
    id plus a distribution column - see the module-level constants'
    docstring."""
    self_containment, self_pred = self_totals
    baseline_containment, baseline_pred = baseline_totals
    rows = []
    all_nodes = (set(self_containment) | set(self_pred)
                 | set(baseline_containment) | set(baseline_pred))
    for node in all_nodes:
        self_values = (self_containment.get(node, 0.0), self_pred.get(node, 0.0))
        baseline_values = (baseline_containment.get(node, 0.0), baseline_pred.get(node, 0.0))
        rows.append({
            'log': log_name,
            'combo': combo_name,
            'degradation_dim': dim,
            'degradation_level': level,
            'node_id': node.id,
            'node_type': type(node).__name__,
            'alphabet': ','.join(sorted(set(node.get_leaf_labels()))),
            **dict(zip(NODE_METRIC_KEYS, self_values)),
            **dict(zip(NODE_METRIC_BASELINE_KEYS, baseline_values)),
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


def _compute_cell(tree, predecessors, log, base_obs):
    """
    (self_totals, self_fields, baseline_totals, baseline_fields) for one
    (tree, log) pairing. 'self' estimates the tail distribution from
    `log` itself (the metric's normal, self-contained mode - see
    surprise.py's module docstring on why this is self-limiting under
    heavy loss). 'baseline' reruns the same events against base_obs (the
    distribution estimated from the *undegraded* log) instead - a
    benchmark-only oracle comparison isolating whether the metric's
    response to degradation is real signal or an artifact of the
    self-estimator degrading along with the log; not available under
    real missing-activity scenarios, where no undegraded reference log
    exists.
    """
    traces = log_to_traces(log)
    self_totals, self_fields = _compute_variant(tree, predecessors, traces, None)
    baseline_totals, baseline_fields = _compute_variant(tree, predecessors, traces, base_obs)
    return self_totals, self_fields, baseline_totals, baseline_fields


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
                tree, _ = discover_cached(log_name, combo_name, combo, base_log)
            except NotImplementedError:
                for dim in degradations:
                    for level in levels:
                        summary_rows.append({
                            'log': log_name, 'combo': combo_name, 'degradation_dim': dim,
                            'degradation_level': level, 'status': 'not_implemented',
                            'elapsed_s': None})
                logger.info('%s - skipped (not_implemented)', cell)
                continue
            except Exception as e:
                logger.exception('%s - discovery exception', cell)
                for dim in degradations:
                    for level in levels:
                        summary_rows.append({
                            'log': log_name, 'combo': combo_name, 'degradation_dim': dim,
                            'degradation_level': level,
                            'status': f'discovery error: {type(e).__name__}: {e}',
                            'elapsed_s': None})
                continue
            logger.info('%s - discovered in %.1fs', cell, time.monotonic() - started)

            predecessors = compute_predecessors(tree)

            # level 0.0 drops nothing regardless of dimension, so it's the
            # same (tree, log) computation under every dim - compute it
            # once here and reuse across dims, mirroring exp_disco_degrade.
            zero_level_cell = None
            zero_level_elapsed_s = None
            if 0.0 in levels:
                with Timer() as t:
                    zero_level_cell = _compute_cell(tree, predecessors, base_log, base_obs)
                zero_level_elapsed_s = t.elapsed_s

            for dim, degrade in degradations.items():
                for level in levels:
                    sub_cell = f'{cell} / {dim} / {level}'
                    with Timer() as t:
                        try:
                            if level == 0.0 and zero_level_cell is not None:
                                self_totals, self_fields, baseline_totals, baseline_fields = zero_level_cell
                                dropped_str, dropped_count = '', 0
                            else:
                                degraded_log, dropped = degrade(base_log, level)
                                dropped_str, dropped_count = _format_dropped(dropped)
                                self_totals, self_fields, baseline_totals, baseline_fields = \
                                    _compute_cell(tree, predecessors, degraded_log, base_obs)

                            node_rows.extend(_node_rows(log_name, combo_name, dim, level,
                                                         self_totals, baseline_totals))

                            baseline_summary_values = (baseline_fields['headline_bits'],
                                                        baseline_fields['bits_per_event'])
                            summary_row = {
                                'log': log_name, 'combo': combo_name, 'degradation_dim': dim,
                                'degradation_level': level, 'dropped_count': dropped_count,
                                'status': 'ok', **self_fields,
                                **dict(zip(SUMMARY_METRIC_BASELINE_KEYS, baseline_summary_values)),
                            }
                        except Exception as e:
                            logger.exception('%s - exception', sub_cell)
                            summary_row = {
                                'log': log_name, 'combo': combo_name, 'degradation_dim': dim,
                                'degradation_level': level,
                                'status': f'error: {type(e).__name__}: {e}',
                            }
                    # level 0.0 reuses the shared computation above (same
                    # rationale as zero_level_metrics elsewhere in this
                    # project) - report that computation's own elapsed_s
                    # rather than the near-zero time this iteration itself
                    # took to just look the cached result up.
                    elapsed_s = (zero_level_elapsed_s if level == 0.0 and zero_level_cell is not None
                                 else t.elapsed_s)
                    summary_row['elapsed_s'] = elapsed_s
                    summary_rows.append(summary_row)
                    if summary_row['status'] == 'ok':
                        logger.info('%s - done in %.1fs (%d events, %.2f bits)',
                                    sub_cell, elapsed_s,
                                    self_fields['n_events'], self_fields['headline_bits'])
                    else:
                        logger.warning('%s - done in %.1fs (status=%s)',
                                        sub_cell, elapsed_s, summary_row['status'])

    node_df = _merge_write(pd.DataFrame(node_rows), out_csv)
    summary_df = _merge_write(pd.DataFrame(summary_rows), summary_csv)
    return node_df, summary_df


def top_nodes(node_df, column, n=15):
    """The n rows with the highest `column` value, one row per (log, combo, dim, level)."""
    if node_df.empty:
        return node_df
    return (node_df.sort_values(column, ascending=False)
                    .groupby(['log', 'combo', 'degradation_dim', 'degradation_level'], sort=False)
                    .head(n))


def format_node_table(df, column):
    cols = ['log', 'combo', 'degradation_dim', 'degradation_level',
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
