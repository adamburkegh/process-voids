"""
Experiment: discovered model vs (degraded) log.

For each log, a model is discovered once per combo from the original,
undegraded log - that discovered model is the fixed reference for that
(log, combo). Each degradation dimension x level then degrades the log
and checks it against that same fixed model, producing a dose-response
curve in the coverage metrics as degradation increases.

Computes the full metric roster at the ROOT of the (log, combo)'s
discovered tree, plus a full per-node breakdown (see below) covering
every node including ablation targets like the claims fixture's
appeal_seq/loop_block/assess - lab.claims_fixture's CLAIMS_COMBOS/
CLAIMS_DEGRADATIONS registers that fixture's known tree and its named
ablation targets as combo/degradation-dimension names usable here, so
there's no separate per-target concept or script needed: compute_metrics'
skip-
alignments-based weight_coverage/weight_voidage/skipprob/salign_coverage,
plus process_voids.voidmass_pn's classical-alignment voidmass_deficit/
voidmass_movecount/voidmass_subprocess/voidmass_process/
alignment_coverage_pn. duration_coverage stays off the roster (dropped
- see lab.metrics).

Writes TWO CSVs per run (node_out_csv defaults to inserting '_nodes'
before out_csv's extension): the root-level one above, one row per
(log, combo, dim, level), plus a per-node one - every metric above
again, but evaluated at EVERY node in the discovered tree (Activity,
Tau, and composite Sequence/Xor/And/Loop alike), one row per (log,
combo, dim, level, node_id) - see _node_rows/PER_NODE_METRIC_KEYS. The
root-level row was always a lossy collapse of data already computed
per-node internally (voidmass_table_pn builds a full table; only
vm_table[tree] was ever kept) - this project is fundamentally about
subprocess-level voids, so throwing that away by default was the wrong
call. node_skip_prob is a new id there, not a reuse of the root CSV's
'skipprob' (which is mean_skipprob's average over every Activity leaf
in the whole tree, not any one node's own value).

NOTE on the classical-alignment metrics specifically: they're computed
via align_pn_all with id_loop_list=[] (skip-alignments' own
insert_cycle_checks gap - see session notes), so a discovered tree with
a tau-skippable loop could in principle make that search cycle; it's
still bounded by the per-variant timeout (default 100s, see
voidmass_table_pn), just wastefully so, not hung outright. Worth
watching for surprisingly slow cells on a real, unfamiliar discovered
tree rather than assuming it's always cheap the way it was on the
small synthetic fixtures this was validated against.

run_disco_degrade() takes explicit logs/combos/degradations/levels and
runs their Cartesian product - the full parameter catalog lives in
lab/params.py, and specific reproducible subset selections live in
lab/runs.py (RUNS), runnable by name:

    python -m lab.exp_disco_degrade --run smoke
    python -m lab.exp_disco_degrade logs/rtfm_fine_appeal.xes.gz --levels 0.0 0.5
"""

import argparse
import logging
import time
from dataclasses import replace
from datetime import datetime
from pathlib import Path

import pandas as pd
import pm4py_config as pm4py

from lab.claims_fixture import CLAIMS_COMBOS, CLAIMS_DEGRADATIONS
from lab.logconfig import configure, enable_skipalignments_debug
from lab.metrics import compute_metrics
from lab.params import ALL_COMBOS, ALL_DEGRADATIONS, ALL_LEVELS
from lab.runs import Experiment, RUNS
from lab.timing import Timer
from process_voids.coveragemass import (
    TREE_METRIC_KEYS, mandatory_node_count, total_node_count,
    mass_by_weight, voidage_by_weight, coverage_by_alignment,
)
from process_voids.voidmass_pn import build_id_net, voidmass_table_pn, coverage_by_alignment_pn

logger = logging.getLogger(__name__)

CLASSICAL_ALIGNMENT_TIMEOUT = 100


def _format_dropped(dropped, limit=50):
    """
    (dropped_str, dropped_count). Dropped can be huge for trace-wise
    degradation at high levels (one entry per dropped case) - listing
    all of them bloats the CSV by orders of magnitude for no benefit,
    so beyond `limit` items we just report the count.
    """
    items = sorted(str(d) for d in dropped)
    if len(items) > limit:
        return f'{len(items)} items', len(items)
    return ', '.join(items), len(items)


def _log_stats(log):
    """(n_cases, n_variants), or (None, None) if log isn't a DataFrame."""
    if not hasattr(log, 'groupby'):
        return None, None
    n_cases = log['case:concept:name'].nunique()
    n_variants = log.groupby('case:concept:name')['concept:name'] \
                     .apply(tuple).nunique()
    return n_cases, n_variants


def _variant_probs(log):
    """{trace variant (activity tuple): probability} from a log's case
    frequencies - duplicated from lab.exp_claims_degrade rather than
    imported, per this package's small-helper convention."""
    n_cases = log['case:concept:name'].nunique()
    variants = {}
    for _case, group in log.groupby('case:concept:name', sort=False):
        variant = tuple(group.sort_values('time:timestamp')['concept:name'])
        variants[variant] = variants.get(variant, 0) + 1
    return {v: c / n_cases for v, c in variants.items()}


CLASSICAL_METRIC_KEYS = ('voidmass_deficit', 'voidmass_movecount',
                          'voidmass_subprocess', 'voidmass_process',
                          'alignment_coverage_pn')

# weight_coverage/weight_voidage/salign_coverage are the SAME quantities
# (same functions, same registry entries) as lab.metrics.METRIC_KEYS -
# just evaluated at an arbitrary node instead of only the tree root,
# same as CLASSICAL_METRIC_KEYS/TREE_METRIC_KEYS already are per-node in
# voidmass_table_pn / mandatory_node_count. node_skip_prob is new: the
# root-level 'skipprob' column is mean_skipprob's mean over every
# Activity leaf in the WHOLE tree regardless of subtree (see that
# function's docstring - it ignores its own tree argument for scoping),
# not this specific node's own probability, so it can't be reused here
# without silently changing what the id means - see lab.metric_registry
# for why every id must mean exactly one thing.
PER_NODE_METRIC_KEYS = ('weight_coverage', 'weight_voidage', 'node_skip_prob', 'salign_coverage')


def _classical_metrics(tree, log, net, im, fm, activity_to_id, tau_ids, id_loop_list, dv,
                        timeout=CLASSICAL_ALIGNMENT_TIMEOUT):
    """
    (root-level classical-alignment metrics dict, full per-node vm_table)
    - process_voids.voidmass_pn, reusing dv.skip_probs (already computed
    by compute_metrics) rather than deriving a separate skip-probability
    estimate. vm_table is returned alongside the root-only dict (not
    just discarded) so callers can also build a full per-node breakdown
    - see _node_rows - without paying for a second, redundant
    voidmass_table_pn call."""
    variant_probs = _variant_probs(log)
    vm_table = voidmass_table_pn(tree, variant_probs, net, im, fm, activity_to_id,
                                  tau_ids, id_loop_list=id_loop_list, timeout=timeout)
    root_row = vm_table[tree]
    values = (
        root_row['deficit'],
        root_row['movecount'],
        root_row['voidmass_subprocess'],
        root_row['voidmass_process'],
        coverage_by_alignment_pn(tree, dv.skip_probs[tree], vm_table),
    )
    return dict(zip(CLASSICAL_METRIC_KEYS, values)), vm_table


def _node_rows(log_name, combo_name, dim, level, dv, vm_table):
    """
    One row per node in vm_table (every node in the tree - Activity,
    Tau, and composite Sequence/Xor/And/Loop nodes alike), the full
    per-node breakdown the root-level row in `rows` collapses away.
    Every metric here is the SAME quantity/id as the root-level row,
    just evaluated at that specific node - see PER_NODE_METRIC_KEYS.
    """
    node_rows = []
    for node, classical_row in vm_table.items():
        per_node_values = (
            mass_by_weight(node, dv.skip_probs),
            voidage_by_weight(node, dv.skip_probs),
            dv.skip_probs[node],
            coverage_by_alignment(node, dv),
        )
        classical_values = (
            classical_row['deficit'],
            classical_row['movecount'],
            classical_row['voidmass_subprocess'],
            classical_row['voidmass_process'],
            coverage_by_alignment_pn(node, dv.skip_probs[node], vm_table),
        )
        tree_values = (mandatory_node_count(node), total_node_count(node))
        node_rows.append({
            'log': log_name, 'combo': combo_name,
            'degradation_dim': dim, 'degradation_level': level,
            'node_id': node.id, 'node_type': type(node).__name__,
            'alphabet': ','.join(sorted(set(node.get_leaf_labels()))),
            **dict(zip(PER_NODE_METRIC_KEYS, per_node_values)),
            **dict(zip(CLASSICAL_METRIC_KEYS, classical_values)),
            **dict(zip(TREE_METRIC_KEYS, tree_values)),
        })
    return node_rows


def run_disco_degrade(log_paths, combos=ALL_COMBOS, degradations=ALL_DEGRADATIONS,
                     levels=ALL_LEVELS, out_csv='var/lab/results/exp_disco_degrade.csv',
                     node_out_csv=None):
    logger.info('Experiment: disco_degrade | logs=%s | combos=%s | '
                'degradations=%s | levels=%s',
                [Path(p).stem for p in log_paths], list(combos),
                list(degradations), levels)

    rows = []
    node_rows = []
    for log_path in log_paths:
        log_name = Path(log_path).stem
        base_log = pm4py.read_xes(log_path)
        n_cases, n_variants = _log_stats(base_log)
        logger.info('Log %s: %s cases, %s variants', log_name, n_cases, n_variants)

        for combo_name, combo in combos.items():
            started_discover = time.monotonic()
            try:
                result = combo.discover(base_log)
                tree, ppt_weights = result.tree, result.ppt_weights
                discover_status = 'ok'
            except NotImplementedError:
                tree, ppt_weights = None, None
                discover_status = 'not_implemented'
            except Exception as e:
                logger.exception('Discovery: log=%s combo=%s - exception', log_name, combo_name)
                tree, ppt_weights = None, None
                discover_status = f'discovery error: {type(e).__name__}: {e}'
            logger.info('Discovery: log=%s combo=%s -> %s (%.1fs)',
                        log_name, combo_name, discover_status,
                        time.monotonic() - started_discover)

            # build_id_net is a cheap pm4py conversion (no alignment
            # computation) - safe to do once per (log, combo) up front
            # and reuse across every degradation level below, same as
            # the tree/ppt_weights fixed-reference-model design.
            classical_net = build_id_net(tree) if tree is not None else None

            # mandatory_node_count/total_node_count depend only on the
            # discovered tree's structure, not the (possibly degraded)
            # log - same fixed-per-(log, combo) computation as above, and
            # unaffected by whether a particular cell's metric
            # computation below succeeds or errors, so merged into every
            # row for this (log, combo) unconditionally once known.
            tree_metrics = (dict(zip(TREE_METRIC_KEYS,
                                      (mandatory_node_count(tree), total_node_count(tree))))
                            if tree is not None else {k: None for k in TREE_METRIC_KEYS})

            # Level 0.0 applies zero drops regardless of dimension, so it's
            # the same (log, tree) computation under every dim - compute it
            # once here and reuse the result across dims below, rather than
            # redoing an identical (and potentially very expensive) run per
            # dimension.
            zero_level_status = None
            zero_level_metrics = None
            zero_level_elapsed_s = None
            zero_level_dv = None
            zero_level_vm_table = None
            if tree is not None and 0.0 in levels:
                cell = f'{log_name} / {combo_name} / (all dims) / 0.0'
                logger.debug('%s - starting', cell)
                slpn_path = f'var/lab/disco_degrade_{log_name}_{combo_name}_level0.slpn'
                with Timer() as t:
                    try:
                        zero_level_metrics, dv = compute_metrics(
                            base_log, tree, slpn_path, ppt_weights=ppt_weights, return_dv=True)
                        net, im, fm, activity_to_id, tau_ids, id_loop_list = classical_net
                        classical_metrics, vm_table = _classical_metrics(
                            tree, base_log, net, im, fm, activity_to_id, tau_ids, id_loop_list, dv)
                        zero_level_metrics.update(classical_metrics)
                        zero_level_dv, zero_level_vm_table = dv, vm_table
                        zero_level_status = 'ok'
                    except Exception as e:
                        logger.exception('%s - exception', cell)
                        zero_level_status = f'error: {type(e).__name__}: {e}'
                zero_level_elapsed_s = t.elapsed_s
                if zero_level_status == 'ok':
                    logger.info('%s - done in %.1fs', cell, zero_level_elapsed_s)
                else:
                    logger.warning('%s - done in %.1fs (status=%s)',
                                    cell, zero_level_elapsed_s, zero_level_status)

            for dim, degrade in degradations.items():
                for level in levels:
                    cell = f'{log_name} / {combo_name} / {dim} / {level}'

                    if tree is None:
                        row = {
                            'log': log_name, 'combo': combo_name,
                            'degradation_dim': dim, 'degradation_level': level,
                            'dropped': '', 'dropped_count': 0, 'status': discover_status,
                            'elapsed_s': None,
                            'weight_coverage': None, 'weight_voidage': None, 'skipprob': None,
                            'salign_coverage': None,
                            **{k: None for k in CLASSICAL_METRIC_KEYS},
                            **tree_metrics,
                        }
                        rows.append(row)
                        logger.info('%s - skipped (%s)', cell, discover_status)
                        continue

                    if level == 0.0:
                        row = {
                            'log': log_name, 'combo': combo_name,
                            'degradation_dim': dim, 'degradation_level': level,
                            'dropped': '', 'dropped_count': 0, 'status': zero_level_status,
                            'elapsed_s': zero_level_elapsed_s,
                            **tree_metrics,
                        }
                        if zero_level_status == 'ok':
                            row.update(zero_level_metrics)
                            node_rows.extend(_node_rows(log_name, combo_name, dim, level,
                                                         zero_level_dv, zero_level_vm_table))
                        else:
                            row.update(weight_coverage=None, weight_voidage=None,
                                       skipprob=None, salign_coverage=None,
                                       **{k: None for k in CLASSICAL_METRIC_KEYS})
                        rows.append(row)
                        logger.debug('%s - reused level-0.0 result', cell)
                        continue

                    logger.debug('%s - starting', cell)

                    degraded_log, dropped = degrade(base_log, level)
                    dropped_str, dropped_count = _format_dropped(dropped)
                    row = {
                        'log': log_name,
                        'combo': combo_name,
                        'degradation_dim': dim,
                        'degradation_level': level,
                        'dropped': dropped_str,
                        'dropped_count': dropped_count,
                        **tree_metrics,
                    }
                    slpn_path = (f'var/lab/disco_degrade_{log_name}_{combo_name}_'
                                 f'{dim}_{level}.slpn')
                    with Timer() as t:
                        try:
                            metrics, dv = compute_metrics(
                                degraded_log, tree, slpn_path,
                                ppt_weights=ppt_weights, return_dv=True)
                            net, im, fm, activity_to_id, tau_ids, id_loop_list = classical_net
                            classical_metrics, vm_table = _classical_metrics(
                                tree, degraded_log, net, im, fm, activity_to_id, tau_ids,
                                id_loop_list, dv)
                            metrics.update(classical_metrics)
                            row['status'] = 'ok'
                            row.update(metrics)
                            node_rows.extend(_node_rows(log_name, combo_name, dim, level,
                                                         dv, vm_table))
                        except Exception as e:
                            logger.exception('%s - exception', cell)
                            row.update(status=f'error: {type(e).__name__}: {e}',
                                        weight_coverage=None, weight_voidage=None, skipprob=None,
                                        salign_coverage=None,
                                        **{k: None for k in CLASSICAL_METRIC_KEYS})
                    row['elapsed_s'] = t.elapsed_s
                    rows.append(row)

                    if row['status'] == 'ok':
                        logger.info('%s - done in %.1fs', cell, t.elapsed_s)
                    else:
                        logger.warning('%s - done in %.1fs (status=%s)',
                                        cell, t.elapsed_s, row['status'])

    df = pd.DataFrame(rows)
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    if node_out_csv is None:
        out_path = Path(out_csv)
        node_out_csv = str(out_path.with_name(f'{out_path.stem}_nodes{out_path.suffix}'))
    node_df = pd.DataFrame(node_rows)
    Path(node_out_csv).parent.mkdir(parents=True, exist_ok=True)
    node_df.to_csv(node_out_csv, index=False)
    logger.info('Wrote %d node rows to %s', len(node_df), node_out_csv)

    return df, node_df


def _timestamped(path):
    """
    Insert a YYYYMMDD-HHMMSS stamp before the extension, eg
    'var/lab/results/rtfm.csv' -> 'var/lab/results/rtfm_20260908-121549.csv'.

    Applied unconditionally in main() below, whether --out was given or
    not: run_disco_degrade (unlike exp_claims_degrade/exp_voidmass/
    exp_surprise) has no _merge_write - it's a plain df.to_csv overwrite
    every call, and every cell here is expensive (a real log, a real
    discovered tree, the full classical-alignment pass) - a same-name
    rerun silently clobbering the last one is exactly the mistake this
    exists to make structurally impossible rather than something to
    remember to avoid by hand each time.
    """
    path = Path(path)
    stamp = datetime.now().strftime('%Y%m%d-%H%M%S')
    return str(path.with_name(f'{path.stem}_{stamp}{path.suffix}'))


def _resolve_experiment(args, parser):
    """
    The Experiment this invocation resolves to, from either --run or ad
    hoc logs/--combos/--levels - the single place both --dry-run and
    the real run derive their configuration from, so the two can never
    show/do different things.
    """
    if args.run:
        if args.run not in RUNS:
            parser.error(f'Unknown run {args.run!r}; choices: {sorted(RUNS)}')
        out_csv = args.out or f'var/lab/results/{args.run}.csv'
        return replace(RUNS[args.run], out_csv=out_csv)

    if not args.logs:
        parser.error('logs are required unless --run is given')

    # NAMED_COMBOS/NAMED_DEGRADATIONS (real-discovery + claims-fixture
    # registrations merged) are for ad hoc --combos/--degradations
    # LOOKUP only - ALL_COMBOS/ALL_DEGRADATIONS (real-discovery only)
    # stay the actual default when neither flag is given, so a bare
    # `--run full`/no-flags invocation never silently picks up
    # claims_known (which ignores its log argument entirely and always
    # returns the same known tree, nonsensical against a real log).
    named_combos = {**ALL_COMBOS, **CLAIMS_COMBOS}
    named_degradations = {**ALL_DEGRADATIONS, **CLAIMS_DEGRADATIONS}

    if args.combos:
        for name in args.combos:
            if name not in named_combos:
                parser.error(f'Unknown combo {name!r}; choices: {sorted(named_combos)}')
        combos = {name: named_combos[name] for name in args.combos}
    else:
        combos = ALL_COMBOS

    if args.degradations:
        for name in args.degradations:
            if name not in named_degradations:
                parser.error(f'Unknown degradation {name!r}; '
                              f'choices: {sorted(named_degradations)}')
        degradations = {name: named_degradations[name] for name in args.degradations}
    else:
        degradations = ALL_DEGRADATIONS

    levels = args.levels or ALL_LEVELS
    out_csv = args.out or 'var/lab/results/exp_disco_degrade.csv'
    return Experiment(name='ad hoc', log_paths=args.logs, combos=combos,
                       degradations=degradations, levels=levels, out_csv=out_csv)


def main():
    configure()
    logger.info('Starting exp_disco_degrade')
    parser = argparse.ArgumentParser(
        description='Dose-response coverage experiment: fixed discovered '
                     'model vs degraded log.')
    parser.add_argument('logs', nargs='*', help='XES log path(s) (ignored if --run given)')
    parser.add_argument('--run', help='Named run from lab.runs.RUNS to execute '
                                       '(overrides logs/--levels/--combos/--degradations)')
    parser.add_argument('--combos', nargs='+',
                         help='Discovery combos to use (default: all of '
                              f'{list(ALL_COMBOS)}) - also accepts claims-fixture '
                              f'combos ({list(CLAIMS_COMBOS)}), validated in '
                              '_resolve_experiment rather than here so both '
                              'registries can be checked together.')
    parser.add_argument('--degradations', nargs='+',
                         help='Degradation dimensions to use (default: all of '
                              f'{list(ALL_DEGRADATIONS)}) - also accepts claims-fixture '
                              f'ablation targets ({list(CLAIMS_DEGRADATIONS)}).')
    parser.add_argument('--levels', type=float, nargs='+',
                         help=f'Degradation levels in [0,1] (default: {ALL_LEVELS})')
    parser.add_argument('--out', help='Output CSV path base (default depends on --run, '
                                       'else var/lab/results/exp_disco_degrade.csv) - a '
                                       'timestamp is always inserted before the extension, '
                                       'so this is never overwritten by a later run.')
    parser.add_argument('--verbose', action='store_true',
                         help='Enable skip-alignments debug logging '
                              '(waste ratios, per-variant timing)')
    parser.add_argument('--dry-run', action='store_true',
                         help='Print the resolved experiment configuration '
                              '(logs/combos/degradations/levels/cell count) and '
                              'exit without computing anything.')
    args = parser.parse_args()

    if args.verbose:
        enable_skipalignments_debug()

    experiment = _resolve_experiment(args, parser)

    if args.dry_run:
        print(experiment.describe())
        return

    out_csv = _timestamped(experiment.out_csv)
    df, node_df = run_disco_degrade(log_paths=experiment.log_paths, combos=experiment.combos,
                                     degradations=experiment.degradations,
                                     levels=experiment.levels, out_csv=out_csv)
    print(df)
    logger.info("Wrote %s (%d rows) and %d node rows", out_csv, len(df), len(node_df))


if __name__ == '__main__':
    main()
