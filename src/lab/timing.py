"""
Shared timing helper for experiment entrypoints - one measurement per
row/cell, not per individual metric function. compute_metrics,
voidmass_table_pn, mandatory_node_count etc. stay untimed themselves;
each lab/exp_*.py wraps its own per-cell computation block in one
Timer and reports that single elapsed_s, so timing is consistent
across scripts without every metric having to instrument itself.
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
