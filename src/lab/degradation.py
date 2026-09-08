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
    # 'activity' (the step version) dropped from the default roster -
    # activity_gradual supersedes it for sweep purposes (same fixed drop
    # order, but a continuous ramp instead of a k-step staircase, so it
    # never wastes/collapses levels the way the step version can on a
    # small alphabet - see degrade_activity_wise_gradual's own
    # docstring). degrade_activity_wise itself is untouched and still
    # directly importable for anything that specifically wants the step
    # behaviour.
    'activity_gradual': degrade_activity_wise_gradual,
    'trace': degrade_trace_wise,
}


def degrade_target_subprocess(log: pd.DataFrame, target_activities: Set[str],
                               n_drop_cases: int, seed: int = 42,
                               exclude_cases: Set[str] = None) -> Tuple[pd.DataFrame, Set[str]]:
    '''
    Remove every event whose activity is in target_activities from
    exactly n_drop_cases cases - an explicit count, not a fraction (see
    voidmass-brief.md's E2: fractions rounding to the same integer count
    at low levels produced duplicate, uninformative rows). Only cases
    that actually contain at least one target-activity event are
    eligible - dropping from a case that never had the target would
    silently waste a dose-response level. Selection is a seeded shuffle
    of eligible cases, so lower counts are always a subset of what
    higher counts drop. Only the target's own events are removed from a
    dropped case - everything else in that case is untouched.

    exclude_cases removes cases from eligibility entirely, before the
    shuffle - for a case that would confound the experiment if it were
    ever ablated (eg it's the log's one deliberately-deviated trace, so
    dropping it would remove that deviation as an accidental side effect
    of the ablation rather than the ablation itself producing the
    result - see exp_claims_degrade.py). Not the same as "eligible but
    never happens to be picked": excluded cases are never candidates at
    any n_drop_cases, deliberately, not by chance of the shuffle.

    Not part of DEGRADATIONS: this parameterises by (target, count), not
    a single [0,1] level, so it doesn't fit that registry's shape - used
    directly by lab/exp_voidmass.py's dose-response sweep instead.
    '''
    eligible_mask = log['concept:name'].isin(target_activities)
    eligible_cases = sorted(log.loc[eligible_mask, 'case:concept:name'].unique())
    if exclude_cases:
        eligible_cases = [c for c in eligible_cases if c not in exclude_cases]
    if n_drop_cases > len(eligible_cases):
        raise ValueError(f'n_drop_cases={n_drop_cases} exceeds '
                          f'{len(eligible_cases)} eligible cases')
    rng = random.Random(seed)
    rng.shuffle(eligible_cases)
    dropped_cases = set(eligible_cases[:n_drop_cases])
    drop_mask = log['case:concept:name'].isin(dropped_cases) & eligible_mask
    return log[~drop_mask].copy(), dropped_cases
