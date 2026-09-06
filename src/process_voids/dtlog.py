"""
Parses delimited trace strings into a pm4py-format event log (a pandas
DataFrame with case:concept:name / concept:name / time:timestamp
columns).

Forked from the `convert("a b c", "a b", "a")` convention in pmkoalas'
dtlog module (https://github.com/PMKoalas/pmkoalas, Adam Peter Banham,
MIT licence) and rewritten against pm4py/pandas directly, since nothing
else in this codebase uses koalas' own Trace/EventLog types.
"""

import datetime
from pathlib import Path
from typing import List, Optional

import pandas as pd

DEFAULT_DELIM = ' '
DEFAULT_EVENT_DELIM = ':'
DEFAULT_BASE_TIME = datetime.datetime(2026, 1, 1)

TIME_UNITS = {
    'seconds': datetime.timedelta(seconds=1),
    'minutes': datetime.timedelta(minutes=1),
    'hours': datetime.timedelta(hours=1),
    'days': datetime.timedelta(days=1),
}

_COLUMNS = ['case:concept:name', 'concept:name', 'time:timestamp']


def _case_names(n: int, names: Optional[List[str]]) -> List[str]:
    if names is not None:
        if len(names) != n:
            raise ValueError(f'names has {len(names)} entries for {n} traces')
        return list(names)
    return [f'case_{i + 1}' for i in range(n)]


def convert(*traces: str, delimiter: str = DEFAULT_DELIM,
            names: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Each trace string is a delimiter-separated sequence of activity
    labels, eg convert("o a s p", "o s p"). Events get a synthetic
    one-tick-per-event ordering; use convert_timed for real timestamps.
    """
    case_names = _case_names(len(traces), names)
    rows = []
    tick = 0
    for case, trace in zip(case_names, traces):
        for activity in (trace.split(delimiter) if trace else []):
            rows.append({
                'case:concept:name': case,
                'concept:name': activity,
                'time:timestamp': DEFAULT_BASE_TIME + datetime.timedelta(seconds=tick),
            })
            tick += 1
    return pd.DataFrame(rows, columns=_COLUMNS)


def convert_timed(*traces: str, delimiter: str = DEFAULT_DELIM,
                   event_delim: str = DEFAULT_EVENT_DELIM,
                   names: Optional[List[str]] = None,
                   base_time: datetime.datetime = DEFAULT_BASE_TIME,
                   time_unit: str = 'hours') -> pd.DataFrame:
    """
    Each trace string is a delimiter-separated sequence of
    "activity<event_delim>offset" pairs, eg
    convert_timed("o:0 a:1 s:2 p:3"), matching the a:t notation used in
    the paper's running example.
    """
    if time_unit not in TIME_UNITS:
        raise ValueError(f'time_unit must be one of {sorted(TIME_UNITS)}')
    unit = TIME_UNITS[time_unit]

    case_names = _case_names(len(traces), names)
    rows = []
    for case, trace in zip(case_names, traces):
        for event in (trace.split(delimiter) if trace else []):
            activity, offset = event.split(event_delim)
            rows.append({
                'case:concept:name': case,
                'concept:name': activity,
                'time:timestamp': base_time + float(offset) * unit,
            })
    log = pd.DataFrame(rows, columns=_COLUMNS)
    return log.sort_values(['case:concept:name', 'time:timestamp']).reset_index(drop=True)


def write_xes(log: pd.DataFrame, path: str) -> None:
    """
    Write a pm4py-format log DataFrame to an XES file.

    pm4py's writer assumes nanosecond-resolution timestamps; pandas now
    defaults to microsecond resolution, which silently corrupts the
    written times (a 1000x factor - eg 2026-01-01 round-trips as
    1970-01-21) unless corrected here first.

    astype('datetime64[ns]') also refuses to drop timezone info outright
    (raises rather than silently converting), so a tz-aware column is
    normalised to naive UTC first - matches what pm4py.read_xes hands
    back for logs with an explicit UTC offset.
    """
    import pm4py_config as pm4py
    log = log.copy()
    timestamps = log['time:timestamp']
    if timestamps.dt.tz is not None:
        timestamps = timestamps.dt.tz_convert('UTC').dt.tz_localize(None)
    log['time:timestamp'] = timestamps.astype('datetime64[ns]')
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    pm4py.write_xes(log, path)
