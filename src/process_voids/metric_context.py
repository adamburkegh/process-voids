"""
A single (log, tree) cell's shared, lazily-computed stages, and a
ProcessMetric declaration that scores against them. A metric declares the
stage ids it needs (ProcessMetric.needs, read through CellContext.stage)
rather than a caller threading the right values through by hand at every
call site.

A CellContext is built fresh per cell and never reused across cells - see
STAGES' 'dv' entry, which mutates the tree's own .weight attributes in
place (process_voids.coveragemass.transfer_pt_weights); a stage-dependent
metric is only correct if scored before a LATER cell's own 'dv' stage
runs again, which holds as long as one context is fully scored before the
next is built.

Reference-scoped values that don't depend on the cell's own (possibly
degraded) log - a discovered tree, a classical net, an undegraded-log
surprise distribution - are NOT stages here: they're computed once by
whichever runner is iterating cells and passed into CellContext.__init__
as refs, so this module stays scoped to one cell's own lazy stages.
"""

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from process_voids import pvoid, slpn_importer
from process_voids.coveragemass import (
    transfer_pt_weights, make_executions_cache, make_aligned_duration_cache, log_to_traces,
)
from process_voids.surprise import event_surprise
from process_voids.voidmass_pn import voidmass_table_pn


METRIC_ERROR = object()  # sentinel: distinct from a genuine None metric value


@dataclass(frozen=True)
class ProcessMetric:
    id: str
    scope: str          # 'node' | 'root'
    needs: tuple         # stage ids this metric's compute() reads via ctx.stage(...)
    compute: Callable    # (ctx, node) -> value; node is ctx.tree when scope='root'


def _dv_stage(ctx):
    Path(ctx.slpn_path).parent.mkdir(parents=True, exist_ok=True)
    dv = pvoid.skipprob(ctx.log, ctx.tree, ctx.slpn_path, ppt_weights=ctx.ppt_weights)
    slpn = slpn_importer.read_slpn(ctx.slpn_path)
    transfer_pt_weights(ctx.tree, slpn)
    return dv


def _executions_cache_stage(ctx):
    return make_executions_cache(ctx.tree)


def _traces_stage(ctx):
    return log_to_traces(ctx.log)


def _aligned_duration_cache_stage(ctx):
    return make_aligned_duration_cache(ctx.tree, ctx.log, traces=ctx.stage('traces'))


def _surprise_self_stage(ctx):
    return event_surprise(ctx.stage('traces'), obs=None)


def _variant_probs(log):
    """{trace variant (activity tuple): probability} from a log's case
    frequencies - duplicated from lab.exp_disco_degrade rather than
    imported, per this package's small-helper convention."""
    n_cases = log['case:concept:name'].nunique()
    variants = {}
    for _case, group in log.groupby('case:concept:name', sort=False):
        variant = tuple(group.sort_values('time:timestamp')['concept:name'])
        variants[variant] = variants.get(variant, 0) + 1
    return {v: c / n_cases for v, c in variants.items()}


def _classical_stage(ctx):
    """
    (VoidmassPnResult, variant_probs) for this cell's (tree, log) - the
    classical-alignment pipeline. ctx.refs['classical_net'] is the
    (net, im, fm, activity_to_id, tau_ids, id_loop_list) tuple from
    voidmass_pn.build_id_net(tree) - reference-scoped (depends only on
    tree structure, not the cell's own log), so the runner computes it
    once per (log, combo) and passes it in rather than this stage
    rebuilding it every cell. ctx.refs['classical_timeout'] (optional)
    overrides voidmass_table_pn's default per-variant timeout.
    """
    net, im, fm, activity_to_id, tau_ids, id_loop_list = ctx.refs['classical_net']
    variant_probs = _variant_probs(ctx.log)
    kwargs = {'id_loop_list': id_loop_list}
    if 'classical_timeout' in ctx.refs:
        kwargs['timeout'] = ctx.refs['classical_timeout']
    result = voidmass_table_pn(ctx.tree, variant_probs, net, im, fm, activity_to_id,
                               tau_ids, **kwargs)
    return result, variant_probs


STAGES = {
    'dv': _dv_stage,
    'executions_cache': _executions_cache_stage,
    'traces': _traces_stage,
    'aligned_duration_cache': _aligned_duration_cache_stage,
    'surprise_self': _surprise_self_stage,
    'classical': _classical_stage,
}


class CellContext:
    """
    One cell's (log, tree) worth of lazily-computed, memoised stages, plus
    ProcessMetric scoring with lifecycle events and per-metric error
    isolation.

    refs: reference-scoped values the runner already computed for this
    (log, combo) - e.g. a discovered tree's classical net, or the
    undegraded-log surprise distribution - available to a ProcessMetric's
    compute() as ctx.refs[...], never recomputed here.

    listeners: callables invoked as listener(event, ctx, id_, node, **extra)
    for 'stage_started'/'stage_finished'/'stage_failed'/'metric_started'/
    'metric_finished'/'metric_failed'. 'stage_finished'/'metric_finished'
    carry elapsed_s; 'stage_failed'/'metric_failed' carry elapsed_s and
    exception. Fired around the computing call only - a memoised stage's
    later accesses (success OR failure) fire nothing further.
    """

    STAGES = STAGES

    def __init__(self, log, tree, slpn_path=None, ppt_weights=None, listeners=(), **refs):
        self.log = log
        self.tree = tree
        self.slpn_path = slpn_path
        self.ppt_weights = ppt_weights
        self.listeners = tuple(listeners)
        self.refs = refs
        self._stages = {}

    def _emit(self, event, id_, node, **extra):
        for listener in self.listeners:
            listener(event, self, id_, node, **extra)

    def stage(self, stage_id):
        """
        Computes and memoises STAGES[stage_id](self) on first access. A
        stage that raises is memoised as a failure too - its exception is
        re-raised (not recomputed) on every later access within this
        context's lifetime, so an expensive, failing stage (eg the
        ebi-backed 'dv' stage) runs at most once per cell even though
        every metric needing it fails independently via score()'s own
        exception handling.
        """
        if stage_id not in self._stages:
            self._emit('stage_started', stage_id, None)
            started = time.monotonic()
            try:
                value = self.STAGES[stage_id](self)
            except Exception as e:
                self._emit('stage_failed', stage_id, None,
                           elapsed_s=time.monotonic() - started, exception=e)
                self._stages[stage_id] = ('error', e)
                raise
            self._emit('stage_finished', stage_id, None, elapsed_s=time.monotonic() - started)
            self._stages[stage_id] = ('ok', value)
        status, payload = self._stages[stage_id]
        if status == 'error':
            raise payload
        return payload

    def score(self, metric, node=None):
        node = self.tree if node is None else node
        self._emit('metric_started', metric.id, node)
        started = time.monotonic()
        try:
            value = metric.compute(self, node)
        except Exception as e:
            self._emit('metric_failed', metric.id, node,
                       elapsed_s=time.monotonic() - started, exception=e)
            return METRIC_ERROR
        self._emit('metric_finished', metric.id, node, elapsed_s=time.monotonic() - started)
        return value


def score_all(ctx, metrics, node=None):
    """
    {metric.id: ctx.score(metric, node) for metric in metrics} - the
    row-assembly counterpart to CellContext.score. A ProcessMetric stays
    a single declared quantity; this is the one place "loop several
    metrics over one node" lives, so a caller building a CSV row doesn't
    hand-roll that loop (and its own error handling) at every call site.
    """
    return {metric.id: ctx.score(metric, node) for metric in metrics}
