"""
Log degradation for dose-response experiments.

Independent dimensions of degradation, each parameterised by a level
in [0, 1]:
  - activity-wise: drop a fraction of distinct activity labels, removing
    all their events from the log (an activity becomes wholly unobserved)
  - activity-wise (gradual): same fixed drop order and same one-at-a-time
    commitment, but the "currently being dropped" activity's own events
    are removed as a growing fraction rather than all at once - see
    degrade_activity_wise_gradual for when this matters.
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


def degrade_activity_by_frequency(log: pd.DataFrame, level: float) -> Tuple[pd.DataFrame, Set[str]]:
    '''
    Activities dropped rarest-first by event count, one at a time, until
    roughly `level` of the log's total EVENT volume is gone - unlike
    degrade_activity_wise(_gradual), whose dose is a count of distinct
    activity labels regardless of how much event volume each one
    represents. level=0.5 here means "half the log's events are gone",
    not "half the distinct labels are gone", so the dose corresponds to
    an actual, meaningful quantity - how much of the log's behaviour has
    disappeared - rather than a count that a single very common or very
    rare activity can dominate arbitrarily.

    No seed, unlike the other activity-wise degradations: frequency
    order is a real, fixed property of the log, not an arbitrary choice
    that needs averaging over multiple random orderings to be
    informative (a seeded shuffle's specific draw is itself a confound -
    see the mechanism this replaces). Ties in frequency are broken by
    activity label for reproducibility. The one remaining randomised
    step - which of the boundary activity's own events are the ones
    dropped to hit the target fraction exactly - uses a fixed,
    activity-keyed seed, matching degrade_activity_wise_gradual's own
    convention for the same problem.

    dropped: activities fully removed (same convention as
    degrade_activity_wise_gradual - the boundary activity being
    partially dropped to hit the target exactly is not included).
    '''
    counts = log['concept:name'].value_counts()
    activities_rarest_first = sorted(counts.index, key=lambda a: (counts[a], a))
    target_events = level * len(log)

    dropped = set()
    events_dropped = 0
    mask = pd.Series(False, index=log.index)
    for activity in activities_rarest_first:
        activity_count = counts[activity]
        if events_dropped + activity_count <= target_events:
            dropped.add(activity)
            events_dropped += activity_count
            mask = mask | (log['concept:name'] == activity)
        else:
            remaining = target_events - events_dropped
            if remaining > 0:
                current_idx = list(log.index[log['concept:name'] == activity])
                random.Random(f'freq_gradual:{activity}').shuffle(current_idx)
                n_drop_current = round(remaining)
                mask = mask | log.index.isin(current_idx[:n_drop_current])
            break
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
    # The step version (degrade_activity_wise) is not on the default
    # roster - activity_gradual has the same fixed drop order, but a
    # continuous ramp instead of a k-step staircase, so it never wastes/
    # collapses levels the way the step version can on a small alphabet
    # (see degrade_activity_wise_gradual's own docstring).
    # degrade_activity_wise stays directly importable for anything that
    # specifically wants the step behaviour.
    'activity_gradual': degrade_activity_wise_gradual,
    'activity_frequency_gradual': degrade_activity_by_frequency,
    'trace': degrade_trace_wise,
}


def degrade_target_subprocess(log: pd.DataFrame, target_activities: Set[str],
                               n_drop_cases: int, seed: int = 42,
                               exclude_cases: Set[str] = None) -> Tuple[pd.DataFrame, Set[str]]:
    '''
    Remove every event whose activity is in target_activities from
    exactly n_drop_cases cases - an explicit count, not a fraction: on a
    small eligible-case pool, neighbouring low fractions round to the
    same integer count and give duplicate, uninformative dose-response
    rows. Only cases that actually contain at least one target-activity
    event are
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
    result - see lab.claims_fixture's CLAIMS_EXCLUDE_CASES). Not the
    same as "eligible but never happens to be picked": excluded cases
    are never candidates at any n_drop_cases, deliberately, not by
    chance of the shuffle.

    Not part of DEGRADATIONS: this parameterises by (target, count), not
    a single [0,1] level, so it doesn't fit that registry's shape - used
    directly by lab.exp_voidmass's dose-response sweep, and wrapped as
    [0,1]-level dimensions by lab.claims_fixture's CLAIMS_DEGRADATIONS.
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
