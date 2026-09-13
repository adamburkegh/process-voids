"""
Consolidated experiment entry point: cells = log x combo x degradation
dim x level. combos/degradations both come from lab.params' catalogs
(ALL_COMBOS/ALL_DEGRADATIONS), or lab.claims_fixture's registrations for
the claims fixture, or a named --run (lab.runs.RUNS). degradations={}
(or not given) runs one cell per (log, combo) with no dim/level at all -
the "fixed model, no degradation" case.

Metric ids and their meaning live in lab.metric_registry, not here -
ALL_METRICS below is which of those ids this module currently computes
and how, via process_voids.metric_context.ProcessMetric/CellContext -
see _compute_cell. --metrics restricts a run to a subset, by bare id or
by GROUPS name; only the stages a selected metric's `needs` actually
requires get triggered per cell.

Writes THREE CSVs per run (node_out_csv/timings_out_csv default to
inserting '_nodes'/'_timings' before out_csv's extension): the root-level
one (one row per cell), a per-node one (every selected metric again, one
row per (cell, node_id) - see NODE_ROW_COLUMNS), and a long-form timing
one (one row per stage/metric actually computed in a cell - see
lab.timing.TimingListener).

    python -m lab.run --run smoke
    python -m lab.run logs/rtfm_fine_appeal.xes.gz --levels 0.0 0.5
    python -m lab.run logs/rtfm_fine_appeal.xes.gz --metrics voidsat classical
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
from lab.discovery import discover_cached
from lab.logconfig import configure, enable_skipalignments_debug
from lab.metrics import mean_leaf_skipprob
from lab.params import ALL_COMBOS, ALL_DEGRADATIONS, ALL_LEVELS
from lab.runs import Experiment, RUNS
from lab.timing import Timer, TimingListener
from process_voids.coveragemass import (
    TREE_METRIC_KEYS, mandatory_node_count, total_node_count,
    mass_by_weight, voidage_by_weight, coverage_by_alignment, voidsat,
)
from process_voids.metric_context import CellContext, ProcessMetric, score_all, METRIC_ERROR
from process_voids.voidmass_pn import build_id_net, coverage_by_alignment_pn
from process_voids.voidsalign import voidsalign

logger = logging.getLogger(__name__)

CLASSICAL_ALIGNMENT_TIMEOUT = 100


def _format_dropped(dropped, limit=50):
    """(dropped_str, dropped_count) - beyond `limit` items, just the count."""
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


CLASSICAL_METRIC_KEYS = ('voidmass_deficit_lower', 'voidmass_deficit_upper',
                          'voidmass_movecount', 'voidmass_movecount_bound',
                          'voidmass_subprocess_lower', 'voidmass_subprocess_upper',
                          'voidmass_process_lower', 'voidmass_process_upper',
                          'alignment_coverage_pn2_lower', 'alignment_coverage_pn2_upper')

# Per-cell bookkeeping, root CSV only (not metrics, so not in
# lab.metric_registry).
TIMEOUT_DIAGNOSTIC_KEYS = ('timed_out_count', 'timed_out_weight')

PER_NODE_METRIC_KEYS = ('weight_coverage', 'weight_voidage', 'skipprob', 'salign_coverage',
                         'voidsalign')

ALIGNED_DURATION_METRIC_KEYS = ('voidsat',)

# --metrics group names -> the ids each expands to.
GROUPS = {
    'skip_alignment': PER_NODE_METRIC_KEYS,
    'classical': CLASSICAL_METRIC_KEYS,
    'aligned_duration': ALIGNED_DURATION_METRIC_KEYS,
}

NODE_ROW_COLUMNS = (
    ['log', 'combo', 'degradation_dim', 'degradation_level',
     'node_id', 'node_type', 'alphabet']
    + list(PER_NODE_METRIC_KEYS) + list(CLASSICAL_METRIC_KEYS)
    + list(ALIGNED_DURATION_METRIC_KEYS) + list(TREE_METRIC_KEYS)
)


def _classical_field(field):
    def compute(ctx, node):
        result, _variant_probs = ctx.stage('classical')
        return result.table[node][field]
    return compute


def _alignment_coverage_pn(timed_out_ratio):
    def compute(ctx, node):
        result, variant_probs = ctx.stage('classical')
        return coverage_by_alignment_pn(node, ctx.stage('dv').skip_probs[node], result.skip_dict,
                                        variant_probs, executions_cache=ctx.stage('executions_cache'),
                                        timed_out_ratio=timed_out_ratio)
    return compute


# Every metric this runner can score, both at the tree root and at every
# node - see _compute_cell. mean_leaf_skipprob and TIMEOUT_DIAGNOSTIC_KEYS
# are root-only bookkeeping, not in this list - see _compute_cell.
ALL_METRICS = [
    ProcessMetric(id='weight_coverage', scope='node', needs=('dv',),
                  compute=lambda ctx, node: mass_by_weight(node, ctx.stage('dv').skip_probs)),
    ProcessMetric(id='weight_voidage', scope='node', needs=('dv',),
                  compute=lambda ctx, node: voidage_by_weight(node, ctx.stage('dv').skip_probs)),
    ProcessMetric(id='skipprob', scope='node', needs=('dv',),
                  compute=lambda ctx, node: ctx.stage('dv').skip_probs[node]),
    ProcessMetric(id='salign_coverage', scope='node', needs=('dv', 'executions_cache'),
                  compute=lambda ctx, node: coverage_by_alignment(
                      node, ctx.stage('dv'), executions_cache=ctx.stage('executions_cache'))),
    ProcessMetric(id='voidsalign', scope='node', needs=('dv', 'executions_cache'),
                  compute=lambda ctx, node: voidsalign(
                      node, ctx.tree, ctx.stage('dv').skip_dict_backup, ctx.stage('dv').pl,
                      ctx.stage('dv').skip_probs, executions_cache=ctx.stage('executions_cache'))),
    ProcessMetric(id='voidmass_deficit_lower', scope='node', needs=('classical',),
                  compute=_classical_field('deficit_lower')),
    ProcessMetric(id='voidmass_deficit_upper', scope='node', needs=('classical',),
                  compute=_classical_field('deficit_upper')),
    ProcessMetric(id='voidmass_movecount', scope='node', needs=('classical',),
                  compute=_classical_field('movecount')),
    ProcessMetric(id='voidmass_movecount_bound', scope='node', needs=('classical',),
                  compute=_classical_field('movecount_bound')),
    ProcessMetric(id='voidmass_subprocess_lower', scope='node', needs=('classical',),
                  compute=_classical_field('voidmass_subprocess_lower')),
    ProcessMetric(id='voidmass_subprocess_upper', scope='node', needs=('classical',),
                  compute=_classical_field('voidmass_subprocess_upper')),
    ProcessMetric(id='voidmass_process_lower', scope='node', needs=('classical',),
                  compute=_classical_field('voidmass_process_lower')),
    ProcessMetric(id='voidmass_process_upper', scope='node', needs=('classical',),
                  compute=_classical_field('voidmass_process_upper')),
    ProcessMetric(id='alignment_coverage_pn2_lower', scope='node',
                  needs=('classical', 'dv', 'executions_cache'),
                  compute=_alignment_coverage_pn(0.0)),
    ProcessMetric(id='alignment_coverage_pn2_upper', scope='node',
                  needs=('classical', 'dv', 'executions_cache'),
                  compute=_alignment_coverage_pn(1.0)),
    ProcessMetric(id='voidsat', scope='node', needs=('dv', 'aligned_duration_cache'),
                  compute=lambda ctx, node: voidsat(
                      node, ctx.tree, ctx.stage('dv'), ctx.log,
                      cache=ctx.stage('aligned_duration_cache'))),
]

_METRICS_BY_ID = {m.id: m for m in ALL_METRICS}

MEAN_LEAF_SKIPPROB_METRIC = ProcessMetric(
    id='mean_leaf_skipprob', scope='root', needs=('dv',),
    compute=lambda ctx, node: mean_leaf_skipprob(ctx.tree, ctx.stage('dv').skip_probs))

def _null_metric_values(metrics):
    """
    {metric.id: None for every metric in `metrics`} plus the root-only
    extras - one shared fallback for a cell that never got as far as
    computing anything (a stage failure propagating out of _compute_cell,
    or a discovery failure).

    Takes `metrics` rather than closing over ALL_METRICS so a caller
    that restricts the roster (run()'s own `metrics` parameter) gets a
    null fallback matching what it actually scores - an error row
    carrying a key no successful row in the same run has (or vice versa)
    is exactly the column-inconsistency a restricted roster must avoid.
    """
    values = {m.id: None for m in metrics}
    values.update({k: None for k in TIMEOUT_DIAGNOSTIC_KEYS})
    values['mean_leaf_skipprob'] = None
    return values


# The default fallback, matching every registered metric.
NULL_METRIC_VALUES = _null_metric_values(ALL_METRICS)


def _resolve_metrics(metrics):
    """
    ALL_METRICS filtered to `metrics` - or ALL_METRICS unchanged when
    `metrics` is None/empty. Accepts two shapes: a list of ProcessMetric
    objects already (passed through as-is, order preserved - the shape
    a caller that built its own inclusive/exclusive roster already has,
    eg lab.exp_disco_degrade's --exclude-metric complement), or a list of
    bare ids and/or GROUPS names (resolved by lookup, duplicates
    collapsed, ALL_METRICS' own order preserved - the CLI --metrics
    shape).
    """
    if not metrics:
        return ALL_METRICS
    if all(isinstance(item, ProcessMetric) for item in metrics):
        return list(metrics)
    selected_ids = set()
    for name in metrics:
        if name in GROUPS:
            selected_ids.update(GROUPS[name])
        elif name in _METRICS_BY_ID:
            selected_ids.add(name)
        else:
            raise ValueError(f'Unknown metric or group {name!r}; choices: '
                              f'{sorted(_METRICS_BY_ID)} or a group in {sorted(GROUPS)}')
    return [m for m in ALL_METRICS if m.id in selected_ids]


def _all_nodes(tree):
    """
    Every node in tree (Activity/Tau leaves and composite Sequence/Xor/
    And/Loop nodes alike), root first - independent of any stage, so the
    per-node loop below doesn't force the classical stage (an expensive
    real alignment search) to run just to get the node list, when no
    selected metric actually needs it. Same node set/order as
    voidmass_pn.voidmass_table_pn's own walk.
    """
    yield tree
    for child in getattr(tree, 'children', []):
        yield from _all_nodes(child)


def _compute_cell(log_name, combo_name, dim, level, log, tree, ppt_weights, classical_net,
                  slpn_path, timing_rows, metrics=ALL_METRICS):
    """
    Runs one cell's metric computation via a fresh CellContext, returning
    (root_metrics, node_rows). Every metric in `metrics` is scored once
    per node in the tree (root included), so the root dict is that same
    scoring for node=tree, not a separate computation.

    Every stage a selected metric needs is triggered directly here, not
    left to be reached lazily from inside a metric's own compute() - a
    stage reached only from within a ProcessMetric closure has its
    failure isolated to that one metric (per node) by CellContext.score,
    rather than failing the whole cell. dv/classical failing (the
    ebi/alignment-search calls - the realistic failure mode) should
    still fail the whole cell uniformly, so both propagate out of this
    function when needed; the caller's own try/except records that as a
    cell-wide error row. A stage no selected metric needs is never
    triggered at all.

    metrics: the ProcessMetric subset to score (see _resolve_metrics).

    timing_rows: extended in place with this cell's TimingListener rows,
    stamped with this cell's (log, combo, degradation_dim,
    degradation_level) cell key.
    """
    listener = TimingListener()
    ctx = CellContext(log=log, tree=tree, slpn_path=slpn_path, ppt_weights=ppt_weights,
                      listeners=[listener], classical_net=classical_net,
                      classical_timeout=CLASSICAL_ALIGNMENT_TIMEOUT)
    try:
        needed_stages = {stage for m in metrics for stage in m.needs}
        for stage_id in needed_stages:
            ctx.stage(stage_id)

        selected_ids = {m.id for m in metrics}
        node_rows = []
        for node in _all_nodes(tree):
            values = score_all(ctx, metrics, node=node)
            node_rows.append({
                'log': log_name, 'combo': combo_name,
                'degradation_dim': dim, 'degradation_level': level,
                'node_id': node.id, 'node_type': type(node).__name__,
                'alphabet': ','.join(sorted(set(node.get_leaf_labels()))),
                **{k: None for k in NULL_METRIC_VALUES if k not in selected_ids},
                **{k: (None if v is METRIC_ERROR else v) for k, v in values.items()},
                'mandatory_node_count': mandatory_node_count(node),
                'total_node_count': total_node_count(node),
            })

        # _all_nodes yields the root first, so node_rows[0] is the tree
        # root's row - the root dict is that same scoring, not a second
        # call, plus the root-only extras below.
        root_metrics = {m.id: node_rows[0][m.id] for m in metrics}
        if 'dv' in needed_stages:
            mean_leaf_value = ctx.score(MEAN_LEAF_SKIPPROB_METRIC)
            root_metrics['mean_leaf_skipprob'] = (None if mean_leaf_value is METRIC_ERROR
                                                  else mean_leaf_value)
        else:
            root_metrics['mean_leaf_skipprob'] = None
        if 'classical' in needed_stages:
            result, _variant_probs = ctx.stage('classical')  # already computed, memoised
            root_metrics['timed_out_count'] = result.timed_out_count
            root_metrics['timed_out_weight'] = result.timed_out_weight
        else:
            root_metrics['timed_out_count'] = None
            root_metrics['timed_out_weight'] = None
        return root_metrics, node_rows
    finally:
        timing_rows.extend({'log': log_name, 'combo': combo_name, 'degradation_dim': dim,
                            'degradation_level': level, **row}
                           for row in listener.rows)


def _compute_no_degradation_cell(log_name, combo_name, log, tree, ppt_weights, classical_net,
                                 slpn_path, timing_rows, metrics):
    """
    One cell for a (log, combo) pair with no degradation dimension at
    all - the "fixed model, as-is log" case (dim=level=None in the
    written row). Delegates to _compute_cell with dim/level both None,
    which is already a valid cell key everywhere else in this module -
    no separate row shape needed.
    """
    return _compute_cell(log_name, combo_name, None, None, log, tree, ppt_weights,
                         classical_net, slpn_path, timing_rows, metrics=metrics)


def run(log_paths, combos=ALL_COMBOS, degradations=ALL_DEGRADATIONS, levels=ALL_LEVELS,
       metrics=None, out_csv='var/lab/results/run.csv', node_out_csv=None, timings_out_csv=None):
    """
    Runs every (log, combo, degradation_dim, degradation_level) cell -
    or, if `degradations` is empty, one (log, combo) cell with no
    degradation at all - and writes the root/node/timings CSVs. Returns
    (root_df, node_df, timings_df).
    """
    selected_metrics = _resolve_metrics(metrics)
    null_metric_values = _null_metric_values(selected_metrics)
    logger.info('Experiment: run | logs=%s | combos=%s | degradations=%s | levels=%s | metrics=%s',
                [Path(p).stem for p in log_paths], list(combos), list(degradations), levels,
                [m.id for m in selected_metrics])

    rows = []
    node_rows = []
    timing_rows = []
    for log_path in log_paths:
        log_name = Path(log_path).stem
        base_log = pm4py.read_xes(log_path)
        n_cases, n_variants = _log_stats(base_log)
        logger.info('Log %s: %s cases, %s variants', log_name, n_cases, n_variants)

        for combo_name, combo in combos.items():
            started_discover = time.monotonic()
            try:
                tree, ppt_weights = discover_cached(
                    log_name, combo_name, combo, base_log)
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

            classical_net = build_id_net(tree) if tree is not None else None
            tree_metrics = (dict(zip(TREE_METRIC_KEYS,
                                      (mandatory_node_count(tree), total_node_count(tree))))
                            if tree is not None else {k: None for k in TREE_METRIC_KEYS})

            if not degradations:
                cell = f'{log_name} / {combo_name} / (no degradation)'
                if tree is None:
                    rows.append({
                        'log': log_name, 'combo': combo_name,
                        'degradation_dim': None, 'degradation_level': None,
                        'dropped': '', 'dropped_count': 0, 'status': discover_status,
                        'elapsed_s': None, **null_metric_values, **tree_metrics,
                    })
                    logger.info('%s - skipped (%s)', cell, discover_status)
                    continue
                slpn_path = f'var/lab/run_{log_name}_{combo_name}.slpn'
                row = {'log': log_name, 'combo': combo_name,
                       'degradation_dim': None, 'degradation_level': None,
                       'dropped': '', 'dropped_count': 0, **tree_metrics}
                with Timer() as t:
                    try:
                        cell_metrics, cell_node_rows = _compute_no_degradation_cell(
                            log_name, combo_name, base_log, tree, ppt_weights, classical_net,
                            slpn_path, timing_rows, selected_metrics)
                        row['status'] = 'ok'
                        row.update(cell_metrics)
                        node_rows.extend(cell_node_rows)
                    except Exception as e:
                        logger.exception('%s - exception', cell)
                        row.update(status=f'error: {type(e).__name__}: {e}', **null_metric_values)
                row['elapsed_s'] = t.elapsed_s
                rows.append(row)
                if row['status'] == 'ok':
                    logger.info('%s - done in %.1fs', cell, t.elapsed_s)
                else:
                    logger.warning('%s - done in %.1fs (status=%s)', cell, t.elapsed_s, row['status'])
                continue

            zero_level_status = None
            zero_level_metrics = None
            zero_level_elapsed_s = None
            zero_level_node_rows = None
            if tree is not None and 0.0 in levels:
                cell = f'{log_name} / {combo_name} / (all dims) / 0.0'
                logger.debug('%s - starting', cell)
                slpn_path = f'var/lab/run_{log_name}_{combo_name}_level0.slpn'
                with Timer() as t:
                    try:
                        # Computed HERE, once, rather than once per dim
                        # below - mass_by_weight/voidage_by_weight read
                        # tree.weight/child.weight directly off the
                        # shared, mutable tree object, which the 'dv'
                        # stage's transfer_pt_weights overwrites on every
                        # OTHER cell's call too. A second dim reusing
                        # this "shared" level-0.0 result would otherwise
                        # recompute weight_coverage/weight_voidage LIVE,
                        # after an unrelated nonzero-level cell of the
                        # first dim has already reset those attributes.
                        # dim is stamped in per use below, since these
                        # rows are otherwise identical regardless of
                        # which dim reuses them.
                        zero_level_metrics, zero_level_node_rows = _compute_cell(
                            log_name, combo_name, None, 0.0, base_log, tree, ppt_weights,
                            classical_net, slpn_path, timing_rows, metrics=selected_metrics)
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
                        rows.append({
                            'log': log_name, 'combo': combo_name,
                            'degradation_dim': dim, 'degradation_level': level,
                            'dropped': '', 'dropped_count': 0, 'status': discover_status,
                            'elapsed_s': None, **null_metric_values, **tree_metrics,
                        })
                        logger.info('%s - skipped (%s)', cell, discover_status)
                        continue

                    if level == 0.0:
                        row = {
                            'log': log_name, 'combo': combo_name,
                            'degradation_dim': dim, 'degradation_level': level,
                            'dropped': '', 'dropped_count': 0, 'status': zero_level_status,
                            'elapsed_s': zero_level_elapsed_s, **tree_metrics,
                        }
                        if zero_level_status == 'ok':
                            row.update(zero_level_metrics)
                            node_rows.extend({**r, 'degradation_dim': dim}
                                             for r in zero_level_node_rows)
                        else:
                            row.update(null_metric_values)
                        rows.append(row)
                        logger.debug('%s - reused level-0.0 result', cell)
                        continue

                    logger.debug('%s - starting', cell)
                    degraded_log, dropped = degrade(base_log, level)
                    dropped_str, dropped_count = _format_dropped(dropped)
                    row = {
                        'log': log_name, 'combo': combo_name,
                        'degradation_dim': dim, 'degradation_level': level,
                        'dropped': dropped_str, 'dropped_count': dropped_count, **tree_metrics,
                    }
                    slpn_path = (f'var/lab/run_{log_name}_{combo_name}_{dim}_{level}.slpn')
                    with Timer() as t:
                        try:
                            cell_metrics, cell_node_rows = _compute_cell(
                                log_name, combo_name, dim, level, degraded_log, tree, ppt_weights,
                                classical_net, slpn_path, timing_rows, metrics=selected_metrics)
                            row['status'] = 'ok'
                            row.update(cell_metrics)
                            node_rows.extend(cell_node_rows)
                        except Exception as e:
                            logger.exception('%s - exception', cell)
                            row.update(status=f'error: {type(e).__name__}: {e}', **null_metric_values)
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
    node_df = pd.DataFrame(node_rows, columns=NODE_ROW_COLUMNS)
    Path(node_out_csv).parent.mkdir(parents=True, exist_ok=True)
    node_df.to_csv(node_out_csv, index=False)
    logger.info('Wrote %d node rows to %s', len(node_df), node_out_csv)

    if timings_out_csv is None:
        out_path = Path(out_csv)
        timings_out_csv = str(out_path.with_name(f'{out_path.stem}_timings{out_path.suffix}'))
    timings_df = pd.DataFrame(timing_rows)
    Path(timings_out_csv).parent.mkdir(parents=True, exist_ok=True)
    timings_df.to_csv(timings_out_csv, index=False)
    logger.info('Wrote %d timing rows to %s', len(timings_df), timings_out_csv)

    return df, node_df, timings_df


def _timestamped(path):
    """Insert a YYYYMMDD-HHMMSS stamp before the extension."""
    path = Path(path)
    stamp = datetime.now().strftime('%Y%m%d-%H%M%S')
    return str(path.with_name(f'{path.stem}_{stamp}{path.suffix}'))


def _resolve_experiment(args, parser):
    """The Experiment this invocation resolves to, from either --run or
    ad hoc logs/--combos/--levels - the single place both --dry-run and
    the real run derive their configuration from."""
    if args.run:
        if args.run not in RUNS:
            parser.error(f'Unknown run {args.run!r}; choices: {sorted(RUNS)}')
        out_csv = args.out or f'var/lab/results/{args.run}.csv'
        return replace(RUNS[args.run], out_csv=out_csv)

    if not args.logs:
        parser.error('logs are required unless --run is given')

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
    elif args.no_degradation:
        degradations = {}
    else:
        degradations = ALL_DEGRADATIONS

    levels = args.levels or ALL_LEVELS
    out_csv = args.out or 'var/lab/results/run.csv'
    return Experiment(name='ad hoc', log_paths=args.logs, combos=combos,
                       degradations=degradations, levels=levels, out_csv=out_csv)


def main():
    configure()
    logger.info('Starting lab.run')
    parser = argparse.ArgumentParser(
        description='Consolidated experiment runner: cells = log x combo x '
                     'degradation dim x level.')
    parser.add_argument('logs', nargs='*', help='XES log path(s) (ignored if --run given)')
    parser.add_argument('--run', help='Named run from lab.runs.RUNS to execute '
                                       '(overrides logs/--levels/--combos/--degradations)')
    parser.add_argument('--combos', nargs='+',
                         help='Discovery combos to use (default: all of '
                              f'{list(ALL_COMBOS)}) - also accepts claims-fixture '
                              f'combos ({list(CLAIMS_COMBOS)}).')
    parser.add_argument('--degradations', nargs='+',
                         help='Degradation dimensions to use (default: all of '
                              f'{list(ALL_DEGRADATIONS)}) - also accepts claims-fixture '
                              f'ablation targets ({list(CLAIMS_DEGRADATIONS)}).')
    parser.add_argument('--no-degradation', action='store_true',
                         help='Run each log/combo as-is, with no degradation dimension at '
                              'all (one row per log/combo). Ignored if --degradations is given.')
    parser.add_argument('--levels', type=float, nargs='+',
                         help=f'Degradation levels in [0,1] (default: {ALL_LEVELS})')
    parser.add_argument('--metrics', nargs='+',
                         help='Metric ids or group names to compute (default: everything - '
                              f'groups: {sorted(GROUPS)}).')
    parser.add_argument('--out', default=None,
                         help='Output CSV path base (default depends on --run, else '
                              'var/lab/results/run.csv) - a timestamp is always inserted '
                              'before the extension.')
    parser.add_argument('--verbose', action='store_true',
                         help='Enable skip-alignments debug logging '
                              '(waste ratios, per-variant timing)')
    parser.add_argument('--dry-run', action='store_true',
                         help='Print the resolved experiment configuration and exit '
                              'without computing anything.')
    args = parser.parse_args()

    if args.verbose:
        enable_skipalignments_debug()

    experiment = _resolve_experiment(args, parser)

    if args.dry_run:
        print(experiment.describe())
        return

    out_csv = _timestamped(experiment.out_csv)
    df, node_df, timings_df = run(
        log_paths=experiment.log_paths, combos=experiment.combos,
        degradations=experiment.degradations, levels=experiment.levels,
        metrics=args.metrics, out_csv=out_csv)
    print(df)
    logger.info("Wrote %s (%d rows), %d node rows and %d timing rows",
                out_csv, len(df), len(node_df), len(timings_df))


if __name__ == '__main__':
    main()
