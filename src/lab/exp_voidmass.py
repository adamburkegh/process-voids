"""
Experiment: voidmass/voidage dose-response to ablating a target subprocess.

A tree is discovered once per (log, combo) from the original,
undegraded log - the same fixed-reference-model design used elsewhere
in this package. A target subprocess's events are then removed from a
growing, EXPLICIT number of cases (not a fraction - see
lab.degradation.degrade_target_subprocess for why), and the full
ebi-backed alignment pipeline (pvoid.skipprob) is recomputed at each
level, reporting all four voidmass/voidage variants per node alongside
the target's rank among all nodes.

This needs the expensive alignment pipeline (unlike lab.exp_surprise,
which was built specifically to avoid it) - expect this to be far
slower per cell, especially on large logs.

Usage:
    python -m lab.exp_voidmass <log_path> --target a [--combos inductive] \\
        [--n-drops 0 1 2 3] [--out ...] [--summary-out ...]
"""

import argparse
import logging
import time
from pathlib import Path

import pandas as pd
import pm4py_config as pm4py

from lab.degradation import degrade_target_subprocess
from lab.discovery import COMBOS, discover_cached
from lab.logconfig import configure
from process_voids import pvoid
from process_voids.coveragemass import voidmass_table

logger = logging.getLogger(__name__)

CELL_COLS = ['log', 'combo', 'target', 'n_drop_cases']

# Metric ids this script's per-node rows carry, in the order
# run_voidmass_doseresponse actually writes them - see lab.metric_registry,
# which this tuple's membership is checked against.
NODE_METRIC_KEYS = ('skip_prob', 'deficit', 'movecount', 'voidmass_subprocess',
                    'voidmass_process', 'voidage_subprocess', 'voidage_process')

# Metric ids this script's summary rows carry - see NODE_METRIC_KEYS.
SUMMARY_METRIC_KEYS = ('target_voidmass_subprocess', 'target_voidmass_process',
                       'target_voidage_subprocess', 'target_voidage_process',
                       'target_rank_voidmass_process', 'target_rank_voidage_process',
                       'n_optimal_alignments')


def _merge_write(df, path, cell_cols=CELL_COLS):
    """
    Upsert by cell - see lab.exp_surprise._merge_write for the full
    rationale. An empty df (e.g. target not found anywhere) has no
    columns to key by - leave whatever's on disk untouched.
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


def _target_nodes(tree, target_activities):
    """Leaves of tree whose name is in target_activities."""
    return [leaf for leaf in tree.get_leafs()
            if getattr(leaf, 'name', None) in target_activities]


def _rank_descending(table, node, column):
    """1-indexed rank of node among all of table's nodes, sorted by column descending."""
    value = table[node][column]
    higher = sum(1 for row in table.values() if row[column] > value)
    return higher + 1


def _n_optimal_alignments(skip_dict):
    """Total optimal alignments found across all variants (the sum of
    |Gamma_sigma|)."""
    return sum(len(states) for states in skip_dict.values())


def run_voidmass_doseresponse(log_paths, target_activities, combos=COMBOS, n_drops=None,
                               out_csv='var/lab/results/exp_voidmass.csv',
                               summary_csv='var/lab/results/exp_voidmass_summary.csv'):
    target_activities = set(target_activities)
    target_name = ','.join(sorted(target_activities))
    logger.info('Experiment: voidmass dose-response | logs=%s | target=%s | combos=%s',
                [Path(p).stem for p in log_paths], target_name, list(combos))

    node_rows = []
    summary_rows = []
    for log_path in log_paths:
        log_name = Path(log_path).stem
        base_log = pm4py.read_xes(log_path)

        for combo_name, combo in combos.items():
            cell = f'{log_name} / {combo_name}'
            started = time.monotonic()
            try:
                tree, _ = discover_cached(log_name, combo_name, combo, base_log)
            except Exception as e:
                logger.warning('%s - discovery error: %s', cell, e)
                continue
            logger.info('%s - discovered in %.1fs', cell, time.monotonic() - started)

            target_nodes = _target_nodes(tree, target_activities)
            if not target_nodes:
                logger.warning('%s - target %s not found in discovered tree, skipping',
                                cell, target_name)
                continue

            max_drops = len(base_log.loc[
                base_log['concept:name'].isin(target_activities), 'case:concept:name'
            ].unique())
            levels = n_drops if n_drops is not None else list(range(max_drops + 1))

            for n in levels:
                sub_cell = f'{cell} / n_drop_cases={n}'
                started = time.monotonic()
                try:
                    if n == 0:
                        degraded_log, dropped = base_log, set()
                    else:
                        degraded_log, dropped = degrade_target_subprocess(
                            base_log, target_activities, n)

                    slpn_path = f'var/lab/voidmass_{log_name}_{combo_name}_n{n}.slpn'
                    dv = pvoid.skipprob(degraded_log, tree, slpn_path)
                    table = voidmass_table(tree, dv.skip_dict_backup, dv.pl,
                                            skip_probs=dv.skip_probs)

                    for node, row in table.items():
                        node_rows.append({
                            'log': log_name, 'combo': combo_name, 'target': target_name,
                            'n_drop_cases': n, 'node_id': node.id,
                            'node_type': type(node).__name__,
                            'alphabet': ','.join(sorted(set(node.get_leaf_labels()))),
                            'skip_prob': dv.skip_probs[node],
                            **row,
                        })

                    target_node = target_nodes[0]
                    summary_rows.append({
                        'log': log_name, 'combo': combo_name, 'target': target_name,
                        'n_drop_cases': n, 'dropped_case_count': len(dropped),
                        'status': 'ok',
                        'target_voidmass_subprocess': table[target_node]['voidmass_subprocess'],
                        'target_voidmass_process': table[target_node]['voidmass_process'],
                        'target_voidage_subprocess': table[target_node]['voidage_subprocess'],
                        'target_voidage_process': table[target_node]['voidage_process'],
                        'target_rank_voidmass_process': _rank_descending(
                            table, target_node, 'voidmass_process'),
                        'target_rank_voidage_process': _rank_descending(
                            table, target_node, 'voidage_process'),
                        'n_nodes': len(table),
                        'n_optimal_alignments': _n_optimal_alignments(dv.skip_dict_backup),
                        'elapsed_s': time.monotonic() - started,
                    })
                    logger.info('%s - done in %.1fs', sub_cell, time.monotonic() - started)
                except Exception as e:
                    summary_rows.append({
                        'log': log_name, 'combo': combo_name, 'target': target_name,
                        'n_drop_cases': n, 'status': f'error: {e}'})
                    logger.warning('%s - error: %s', sub_cell, e)

    node_df = _merge_write(pd.DataFrame(node_rows), out_csv)
    summary_df = _merge_write(pd.DataFrame(summary_rows), summary_csv)
    return node_df, summary_df


def main():
    configure()
    logger.info('Starting exp_voidmass')
    parser = argparse.ArgumentParser(
        description='Voidmass/voidage dose-response: ablate a target subprocess from a '
                     'growing number of cases and track all four variants per node.')
    parser.add_argument('logs', nargs='+', help='XES log path(s)')
    parser.add_argument('--target', nargs='+', required=True,
                         help='Target subprocess activity label(s)')
    parser.add_argument('--combos', nargs='+', default=list(COMBOS),
                         help=f'Discovery combos to use (default: all of {list(COMBOS)})')
    parser.add_argument('--n-drops', type=int, nargs='+', default=None,
                         help='Explicit case-drop counts (default: every level from 0 to '
                              'the number of cases containing the target)')
    parser.add_argument('--out', default='var/lab/results/exp_voidmass.csv')
    parser.add_argument('--summary-out', default='var/lab/results/exp_voidmass_summary.csv')
    args = parser.parse_args()

    combos = {name: COMBOS[name] for name in args.combos}
    node_df, summary_df = run_voidmass_doseresponse(
        args.logs, args.target, combos=combos, n_drops=args.n_drops,
        out_csv=args.out, summary_csv=args.summary_out)

    log_names = [Path(p).stem for p in args.logs]
    requested = summary_df[summary_df['log'].isin(log_names) & summary_df['combo'].isin(combos)]
    print(requested.to_string(index=False))
    logger.info('Wrote %s and %s', args.out, args.summary_out)


if __name__ == '__main__':
    main()
