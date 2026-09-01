"""
Log degradation for dose-response experiments.

Two independent dimensions of degradation, each parameterised by a level
in [0, 1]:
  - activity-wise: drop a fraction of distinct activity labels, removing
    all their events from the log (an activity becomes wholly unobserved)
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
    'trace': degrade_trace_wise,
}
