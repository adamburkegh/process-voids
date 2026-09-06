"""
Log degradation for dose-response experiments.

Independent dimensions of degradation, each parameterised by a level
in [0, 1]:
  - activity-wise: drop a fraction of distinct activity labels, removing
    all their events from the log (an activity becomes wholly unobserved)
  - activity-wise (gradual): same fixed drop order and same one-at-a-time
    commitment, but the "currently being dropped" activity's own events
    are removed as a growing fraction rather than all at once - see
    degrade_activity_wise_gradual for why this exists alongside the
    original (kept separately, not replaced, for traceability against
    stats already gathered with the step version).
  - trace-wise: drop a fraction of whole cases from the log

Selection is a seeded shuffle, so a given (log, level, seed) is
reproducible and lower levels are subsets of what higher levels drop.
"""

import random
from typing import Set, Tuple

import pandas as pd


def degrade_activity_wise(log: pd.DataFrame, level: float,
                           seed: int = 42) -> Tuple[pd.DataFrame, Set[str]]:
    activities = sorted(log['concept:name'].unique())
    rng = random.Random(seed)
    rng.shuffle(activities)
    n_drop = round(level * len(activities))
    dropped = set(activities[:n_drop])
    return log[~log['concept:name'].isin(dropped)].copy(), dropped


def degrade_activity_wise_gradual(log: pd.DataFrame, level: float,
                                   seed: int = 42) -> Tuple[pd.DataFrame, Set[str]]:
    """
    Like degrade_activity_wise, but ramps continuously instead of
    stepping. Activities are still fully dropped one at a time in a
    fixed shuffled order - so lower levels are still subsets of what
    higher levels drop, and dropped (fully-removed activities only,
    same meaning as the step version) is still a faithful record of
    which activities have been committed - but the "currently being
    dropped" activity's own events are removed as a growing fraction
    (a seeded shuffle of just that activity's own occurrences) rather
    than disappearing all at once, giving a piecewise-linear
    dose-response curve instead of a k-step staircase (k = number of
    distinct activities). Useful precisely when k is small, where the
    step version is too coarse to show any response between levels.
    """
    activities = sorted(log['concept:name'].unique())
    rng = random.Random(seed)
    rng.shuffle(activities)
    k = len(activities)
    progress = level * k
    n_full = min(int(progress), k)
    frac = progress - n_full
    dropped = set(activities[:n_full])

    mask = log['concept:name'].isin(dropped)
    if n_full < k and frac > 0:
        current = activities[n_full]
        current_idx = list(log.index[log['concept:name'] == current])
        random.Random(f'{seed}:{current}').shuffle(current_idx)
        n_drop_current = round(frac * len(current_idx))
        mask = mask | log.index.isin(current_idx[:n_drop_current])
    return log[~mask].copy(), dropped


def degrade_trace_wise(log: pd.DataFrame, level: float,
                        seed: int = 42) -> Tuple[pd.DataFrame, Set[str]]:
    cases = sorted(log['case:concept:name'].unique())
    rng = random.Random(seed)
    rng.shuffle(cases)
    n_drop = round(level * len(cases))
    dropped = set(cases[:n_drop])
    return log[~log['case:concept:name'].isin(dropped)].copy(), dropped


DEGRADATIONS = {
    'activity': degrade_activity_wise,
    'activity_gradual': degrade_activity_wise_gradual,
    'trace': degrade_trace_wise,
}
