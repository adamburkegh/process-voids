"""
Experiment: discovered model vs (degraded) log.

For each log, a model is discovered once per combo from the original,
undegraded log - that discovered model is the fixed reference for that
(log, combo). Each degradation dimension x level then degrades the log
and checks it against that same fixed model, producing a dose-response
curve in the coverage metrics as degradation increases.

Every metric in ALL_METRICS (process_voids.metric_context.ProcessMetric
declarations - skip-alignments-based weight_coverage/weight_voidage/
skipprob/salign_coverage/voidsalign, process_voids.voidmass_pn's classical-
alignment CLASSICAL_METRIC_KEYS (voidmass deficit/movecount/subprocess/
process and alignment_coverage_pn2 - see that constant), and
process_voids.coveragemass's skip-alignment-based, real-elapsed-time
voidsat) is scored once per node
in the discovered tree (Activity, Tau, and composite Sequence/Xor/And/Loop
alike) via a per-cell CellContext - see _compute_cell. The root-level row
is that same scoring at the tree root, not a separate computation; the
per-node CSV keeps every node's scores, since this project is about
subprocess-level voids. duration_coverage is not on the roster
(product-only - see lab.metric_registry). skipprob here is
dv.skip_probs[node] directly (skip-alignments' own published definition)
- not lab.metrics.mean_leaf_skipprob (a different, root-only statistic,
scored separately - see MEAN_LEAF_SKIPPROB_METRIC).

lab.claims_fixture's CLAIMS_COMBOS/CLAIMS_DEGRADATIONS registers that
fixture's known tree and its named ablation targets (appeal_seq/
loop_block/assess) as combo/degradation-dimension names usable here, so
there's no separate per-target concept or script needed.

Writes THREE CSVs per run (node_out_csv/timings_out_csv default to
inserting '_nodes'/'_timings' before out_csv's extension): the root-level
one above (one row per (log, combo, dim, level)), a per-node one (every
metric above again, one row per (log, combo, dim, level, node_id) - see
NODE_ROW_COLUMNS), and a long-form timing one (one row per stage/metric
actually computed in a cell - see lab.timing.TimingListener).

NOTE on the classical-alignment metrics specifically: build_id_net's
id_loop_list is passed through to align_pn_all, so its cycle guard can
engage on a tau-skippable loop, and each variant's search is bounded
by CLASSICAL_ALIGNMENT_TIMEOUT. A variant that times out is bounded
rather than estimated (the _lower/_upper columns), and one that runs
close to the timeout is logged as a warning - worth watching for on a
real, unfamiliar discovered tree rather than assuming it's always cheap
the way it is on small synthetic fixtures.

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


# _lower/_upper: align_variant_all can return zero alignments for a
# variant (a per-variant timeout, seen in real runs). Its real
# contribution is unknown, so each cell reports provable bounds on the
# value it would have had without the timeout - see voidmass_table_pn
# and coveragemass.alignment_mass for the derivations.
# voidmass_movecount is the OBSERVED total from completed variants;
# voidmass_movecount_bound is the denominator the bounds divide by.
CLASSICAL_METRIC_KEYS = ('voidmass_deficit_lower', 'voidmass_deficit_upper',
                          'voidmass_movecount', 'voidmass_movecount_bound',
                          'voidmass_subprocess_lower', 'voidmass_subprocess_upper',
                          'voidmass_process_lower', 'voidmass_process_upper',
                          'alignment_coverage_pn2_lower', 'alignment_coverage_pn2_upper')

# Per-cell bookkeeping, root CSV only (not metrics, so not in
# lab.metric_registry): how many variants' alignment search timed out,
# and their summed probability. The weight is what matters - 0.03% of a
# log timing out is immaterial, 20% is not - which the count can't show.
TIMEOUT_DIAGNOSTIC_KEYS = ('timed_out_count', 'timed_out_weight')

# weight_coverage/weight_voidage/skipprob/salign_coverage/voidsalign are
# the SAME quantities (same functions/lookups, same registry entries) as
# lab.metrics.METRIC_KEYS - just evaluated at an arbitrary node instead
# of only the tree root, same as CLASSICAL_METRIC_KEYS/TREE_METRIC_KEYS
# already are per-node in voidmass_table_pn / mandatory_node_count.
# skipprob = dv.skip_probs[node] directly (skip-alignments' own
# definition) - NOT lab.metrics.mean_leaf_skipprob, a different,
# unrelated statistic (see lab.metrics' module docstring). voidsalign
# shares dv.skip_dict_backup/dv.pl with salign_coverage's own
# alignment_mass call - same lumped skip-alignment normal form, just a
# different move-weighting (see process_voids.voidsalign's own module
# docstring).
PER_NODE_METRIC_KEYS = ('weight_coverage', 'weight_voidage', 'skipprob', 'salign_coverage',
                         'voidsalign')

# voidsat is skip-alignments-based like PER_NODE_METRIC_KEYS (reuses
# dv.skip_dict_backup/dv.skip_probs, not the classical-alignment path),
# but needs the REAL log (real per-trace timestamps - duration is a
# per-instance quantity, see coveragemass.admass) rather than just dv,
# so it gets its own small key group/helper instead of folding into
# either PER_NODE_METRIC_KEYS's or CLASSICAL_METRIC_KEYS's existing
# computation, which don't carry the log through to where it's needed.
ALIGNED_DURATION_METRIC_KEYS = ('voidsat',)


# Every key _compute_cell's node row dicts carry, in the same order -
# used to give the per-node CSV a real header even when node_rows is empty
# (every cell in a run errored) rather than pandas' default empty-columns
# DataFrame, which writes a headerless file that raises EmptyDataError in
# any downstream reader expecting an empty-but-columned frame instead.
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


# Every metric this experiment scores, both at the tree root and at
# every node - see _compute_cell. Each is the SAME quantity/id whether
# evaluated at the root or a subprocess node - the root row is just this
# list scored at node=tree, not a separately-derived computation (see
# module docstring). mean_leaf_skipprob and TIMEOUT_DIAGNOSTIC_KEYS are
# root-only bookkeeping, not in this list - see _compute_cell.
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

MEAN_LEAF_SKIPPROB_METRIC = ProcessMetric(
    id='mean_leaf_skipprob', scope='root', needs=('dv',),
    compute=lambda ctx, node: mean_leaf_skipprob(ctx.tree, ctx.stage('dv').skip_probs))

def _null_metric_values(metrics):
    """
    {metric.id: None for every metric in `metrics`} plus the root-only
    extras - one shared fallback for a cell that never got as far as
    computing anything (a stage failure propagating out of _compute_cell,
    or a discovery failure).

    Takes `metrics` rather than closing over ALL_METRICS so a caller that
    excludes a metric (run_disco_degrade's own `metrics` parameter) gets
    a null fallback matching what it actually scores - an error row
    carrying a key no successful row in the same run has (or vice versa)
    is exactly the column-inconsistency excluding a metric must avoid.
    """
    values = {m.id: None for m in metrics}
    values.update({k: None for k in TIMEOUT_DIAGNOSTIC_KEYS})
    values['mean_leaf_skipprob'] = None
    return values


# The default fallback, matching every registered metric.
NULL_METRIC_VALUES = _null_metric_values(ALL_METRICS)


def _compute_cell(log_name, combo_name, dim, level, log, tree, ppt_weights, classical_net,
                  slpn_path, timing_rows, metrics=ALL_METRICS):
    """
    Runs one cell's full metric computation via a fresh CellContext,
    returning (root_metrics, node_rows). Every metric in `metrics`
    (default ALL_METRICS - see run_disco_degrade's own `metrics`
    parameter for why a caller might narrow this) is scored once per
    node in the classical stage's table (root included), so the root
    dict is that same scoring for node=tree, not a separate computation.

    A stage failure (the dv/classical alignment pipelines - the realistic
    failure mode) propagates out of this function; the caller's own
    try/except records that as a cell-wide error row. A single
    ProcessMetric's own compute() failing
    is isolated per (metric, node) by CellContext.score and surfaces as
    METRIC_ERROR -> None in the row it belongs to, not a cell-wide
    failure - see process_voids.metric_context.

    timing_rows: extended in place with this cell's TimingListener rows,
    stamped with this cell's (log, combo, degradation_dim,
    degradation_level) cell key - see lab.timing.TimingListener.
    """
    listener = TimingListener()
    ctx = CellContext(log=log, tree=tree, slpn_path=slpn_path, ppt_weights=ppt_weights,
                      listeners=[listener], classical_net=classical_net,
                      classical_timeout=CLASSICAL_ALIGNMENT_TIMEOUT)
    try:
        # Both stages triggered directly, outside any ProcessMetric's
        # compute() - a stage accessed only from within a metric closure
        # (via score_all/ctx.score below) has its failure caught and
        # isolated to just that metric, per node, rather than failing the
        # whole cell. dv/classical failing (the ebi/alignment-search
        # calls - the realistic failure mode) should still fail the
        # whole cell uniformly.
        result, _variant_probs = ctx.stage('classical')
        ctx.stage('dv')

        node_rows = []
        for node in result.table:
            values = score_all(ctx, metrics, node=node)
            node_rows.append({
                'log': log_name, 'combo': combo_name,
                'degradation_dim': dim, 'degradation_level': level,
                'node_id': node.id, 'node_type': type(node).__name__,
                'alphabet': ','.join(sorted(set(node.get_leaf_labels()))),
                **{k: (None if v is METRIC_ERROR else v) for k, v in values.items()},
                'mandatory_node_count': mandatory_node_count(node),
                'total_node_count': total_node_count(node),
            })

            # result.table's own _walk inserts the root first (see
            # voidmass_pn), so node_rows[0] is the tree root's row - the
            # root dict is that same scoring, not a second call, plus the
            # two root-only extras.
        root_metrics = ({m.id: node_rows[0][m.id] for m in metrics} if node_rows
                        else _null_metric_values(metrics))
        mean_leaf_value = ctx.score(MEAN_LEAF_SKIPPROB_METRIC)
        root_metrics['mean_leaf_skipprob'] = (None if mean_leaf_value is METRIC_ERROR
                                              else mean_leaf_value)
        root_metrics['timed_out_count'] = result.timed_out_count
        root_metrics['timed_out_weight'] = result.timed_out_weight
        return root_metrics, node_rows
    finally:
        # Flushed even when a stage failure propagates out of this
        # function below - the listener has already recorded that
        # stage's own stage_failed row by the time the exception reaches
        # here, and it's the one useful diagnostic a caller gets for a
        # cell that otherwise contributes only an error-status row.
        timing_rows.extend({'log': log_name, 'combo': combo_name, 'degradation_dim': dim,
                            'degradation_level': level, **row}
                           for row in listener.rows)


def run_disco_degrade(log_paths, combos=ALL_COMBOS, degradations=ALL_DEGRADATIONS,
                     levels=ALL_LEVELS, metrics=ALL_METRICS,
                     out_csv='var/lab/results/exp_disco_degrade.csv',
                     node_out_csv=None, timings_out_csv=None):
    """
    metrics: the ProcessMetrics to score (default ALL_METRICS) - narrow
    this to skip a metric entirely for this run, eg one that's known-
    expensive at a log's scale (voidsat's per-trace cost) or known-wrong
    pending an upstream fix, where re-collecting it later is cheap. An
    excluded metric gets no timing row and its CSV column is absent from
    every scored row (NaN once written to the fixed NODE_ROW_COLUMNS
    schema) - never a KeyError, and never present on some rows but not
    others, since the null-value fallback for a failed/skipped cell is
    derived from this SAME `metrics` argument (see _null_metric_values).
    """
    null_metric_values = _null_metric_values(metrics)
    logger.info('Experiment: disco_degrade | logs=%s | combos=%s | '
                'degradations=%s | levels=%s | metrics=%s',
                [Path(p).stem for p in log_paths], list(combos),
                list(degradations), levels, [m.id for m in metrics])

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
            zero_level_node_rows = None
            if tree is not None and 0.0 in levels:
                cell = f'{log_name} / {combo_name} / (all dims) / 0.0'
                logger.debug('%s - starting', cell)
                slpn_path = f'var/lab/disco_degrade_{log_name}_{combo_name}_level0.slpn'
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
                        # first dim has already reset those attributes to
                        # weights estimated from a degraded log - silently
                        # corrupting the per-node CSV's weight-derived
                        # columns for every dim after the first. dim is
                        # stamped in per use below, since these rows are
                        # otherwise identical regardless of which dim
                        # reuses them.
                        zero_level_metrics, zero_level_node_rows = _compute_cell(
                            log_name, combo_name, None, 0.0, base_log, tree, ppt_weights,
                            classical_net, slpn_path, timing_rows, metrics=metrics)
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
                            **null_metric_values,
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
                            # Reuse the SAME precomputed rows for every
                            # dim, just stamped with this dim's name -
                            # not a fresh _compute_cell call (see where
                            # zero_level_node_rows is built, above, for
                            # why recomputing here would be wrong, not
                            # just wasteful).
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
                            cell_root_metrics, cell_node_rows = _compute_cell(
                                log_name, combo_name, dim, level, degraded_log, tree, ppt_weights,
                                classical_net, slpn_path, timing_rows, metrics=metrics)
                            row['status'] = 'ok'
                            row.update(cell_root_metrics)
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
    """
    Insert a YYYYMMDD-HHMMSS stamp before the extension, eg
    'var/lab/results/rtfm.csv' -> 'var/lab/results/rtfm_20260908-121549.csv'.

    Applied unconditionally in main() below, whether --out was given or
    not: run_disco_degrade (unlike exp_voidmass/exp_surprise) has no
    _merge_write - it's a plain df.to_csv overwrite
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
    parser.add_argument('--exclude-metric', nargs='+', default=[], metavar='METRIC_ID',
                         help='Skip scoring the given metric id(s) entirely for this run '
                              f'(choices: {[m.id for m in ALL_METRICS]}) - for a metric '
                              "that's known-expensive at a log's scale (eg voidsat's "
                              'per-trace cost) or known-wrong pending an upstream fix, '
                              'where re-collecting it later is cheap. An excluded metric '
                              'gets no timing row and its CSV column is entirely empty '
                              'for this run.')
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

    all_metric_ids = {m.id for m in ALL_METRICS}
    unknown = [m for m in args.exclude_metric if m not in all_metric_ids]
    if unknown:
        parser.error(f'Unknown --exclude-metric {unknown}; choices: {sorted(all_metric_ids)}')
    metrics = [m for m in ALL_METRICS if m.id not in args.exclude_metric]

    experiment = _resolve_experiment(args, parser)

    if args.dry_run:
        print(experiment.describe())
        if args.exclude_metric:
            print(f'  excluding metrics: {args.exclude_metric}')
        return

    out_csv = _timestamped(experiment.out_csv)
    df, node_df, timings_df = run_disco_degrade(
        log_paths=experiment.log_paths, combos=experiment.combos,
        degradations=experiment.degradations, levels=experiment.levels,
        metrics=metrics, out_csv=out_csv)
    print(df)
    logger.info("Wrote %s (%d rows), %d node rows and %d timing rows",
                out_csv, len(df), len(node_df), len(timings_df))


if __name__ == '__main__':
    main()
