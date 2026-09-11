"""
Shared timing helpers for experiment entrypoints.

Timer: one measurement per row/cell, not per individual metric function -
each lab/exp_*.py wraps its own per-cell computation block in one Timer
and reports that single elapsed_s. Still used directly by scripts that
haven't moved onto process_voids.metric_context's CellContext.

TimingListener: per-metric/per-stage timing for a CellContext, driven by
its lifecycle events rather than a hand-placed Timer block.
"""

import time


class Timer:
    """
    Context manager measuring wall-clock seconds for its block.
    .elapsed_s is set on __exit__ regardless of whether the block
    raised, so a caller whose try/except lives inside the `with` still
    gets a timing on the error path, not just the success path.

    Usage:
        with Timer() as t:
            ...compute a cell's metrics...
        row['elapsed_s'] = t.elapsed_s
    """
    def __init__(self):
        self.elapsed_s = None

    def __enter__(self):
        self._started = time.monotonic()
        return self

    def __exit__(self, *exc_info):
        self.elapsed_s = time.monotonic() - self._started
        return False


class TimingListener:
    """
    process_voids.metric_context.CellContext lifecycle listener collecting
    one long-form row per stage/metric actually computed this cell
    (metric_or_stage, seconds, status) - a memoised stage's later accesses
    add nothing, since CellContext only fires *_finished/*_failed around
    the computing call. The runner adds its own cell-identifying columns
    (log, combo, degradation_dim, degradation_level, ...) when it flushes
    .rows into a _timings CSV; construct a fresh listener per cell rather
    than clearing .rows between cells, so a listener's rows are always one
    cell's worth.

    Usage:
        listener = TimingListener()
        ctx = CellContext(..., listeners=[listener])
        ...score every metric for this cell...
        timings_rows.extend({**cell_key, **row} for row in listener.rows)
    """
    def __init__(self):
        self.rows = []

    def __call__(self, event, ctx, id_, node, **extra):
        if event in ('stage_finished', 'metric_finished', 'metric_failed'):
            self.rows.append({
                'metric_or_stage': id_,
                'seconds': extra['elapsed_s'],
                'status': 'error' if event == 'metric_failed' else 'ok',
            })
