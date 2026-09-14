"""
One row per run, appended to a persistent CSV.

Every run-level fact - which combos/degradations/levels were used, which
metrics were scored, the process-voids and skip-alignments versions - is
otherwise recorded only as text, in that run's own
var/lab/logs/<script>_<timestamp>.log. That log stays; this file answers
a different question. A result CSV has to be readable a year later
without finding, and correctly matching by timestamp proximity, a
separate text log, so the facts that change how its columns should be
read live in a file that can be joined to it - on out_csv.

Appended, never upserted. lab.exp_surprise._merge_write purges the rows
matching a key before writing, which is right for results (a rerun
supersedes the cell it recomputed) and wrong here: a rerun writing to the
same out_csv is a second run, and erasing the first one's row would
discard exactly the history this file exists to keep.

These are administrative columns, not metrics, so they have no
lab.metric_registry entries - the registry is the dictionary for result
columns.
"""

import logging
import os
import time
from pathlib import Path

import pandas as pd

from lab.logconfig import (
    _version_line, dependency_version_line, git_version, project_version)

logger = logging.getLogger(__name__)

DEFAULT_RUN_HISTORY_CSV = 'var/lab/results/run_history.csv'

# Recorded when PYTHONHASHSEED is not set. Not '' - an empty cell reads
# back from CSV as NaN, indistinguishable from a column that was never
# written, which is the confusion this file exists to remove.
SEED_UNSET = 'unset'

RUN_HISTORY_COLUMNS = (
    'run_timestamp',
    'out_csv',            # join key back to the per-cell result rows
    'run_name',
    'logs',
    'combos',
    'degradations',
    'levels',
    'metrics',
    'excluded_metrics',
    'pythonhashseed',
    'process_voids_version',
    'skipalignments_version',
)


def _joined(values):
    return ';'.join(str(v) for v in values)


def run_history_row(out_csv, log_paths, combos, degradations, levels, metric_ids,
                    run_name=''):
    """
    One run's row. `metric_ids` is what was actually scored, so
    excluded_metrics is its complement against the full roster - which
    makes the row truthful however the selection was expressed
    (lab.run's inclusive --metrics, or exp_disco_degrade's
    --exclude-metric).

    Logs are recorded by Path().stem, the same spelling the result rows
    use for 'log', so history joins against them without massaging.
    out_csv is recorded with forward slashes rather than this machine's
    separator, so the file stays legible - and joinable against a path
    typed by hand - away from Windows.

    The seed is read from the environment rather than from whatever was
    requested: PYTHONHASHSEED is consumed by the interpreter at startup,
    so what a flag asked for and what this process actually ran under
    can differ, and only the latter explains the result.
    """
    from lab.run import ALL_METRICS

    scored = list(metric_ids)
    excluded = [m.id for m in ALL_METRICS if m.id not in set(scored)]
    return {
        'run_timestamp': time.strftime('%Y-%m-%dT%H:%M:%S'),
        'out_csv': Path(out_csv).as_posix(),
        'run_name': run_name,
        'logs': _joined(Path(p).stem for p in log_paths),
        'combos': _joined(combos),
        'degradations': _joined(degradations),
        'levels': _joined(levels),
        'metrics': _joined(scored),
        'excluded_metrics': _joined(excluded),
        'pythonhashseed': os.environ.get('PYTHONHASHSEED') or SEED_UNSET,
        'process_voids_version': _version_line(project_version(), *git_version()),
        'skipalignments_version': dependency_version_line('skipalignments'),
    }


def append_run_history(path, row):
    """
    Append `row` to the CSV at `path`, creating it with
    RUN_HISTORY_COLUMNS if it doesn't exist yet.

    Never raises. A completed sweep is expensive - losing its results
    because the bookkeeping file couldn't be written would be much the
    worse failure - so a write failure is logged as a warning and the
    caller carries on.
    """
    try:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame([row], columns=list(RUN_HISTORY_COLUMNS))
        header = not path.exists()
        df.to_csv(path, mode='a', header=header, index=False)
    except Exception:
        logger.warning('Could not write run history to %s - the run itself is '
                       'unaffected', path, exc_info=True)


def record_run(path, **row_kwargs):
    """
    Build this run's row and append it - the single entry point a runner
    calls.

    Building the row is inside the same guard as writing it: it shells
    out to git for both packages' commit state, so it can fail for
    reasons that have nothing to do with the run that just succeeded.
    """
    try:
        row = run_history_row(**row_kwargs)
    except Exception:
        logger.warning('Could not build the run history row - the run itself is '
                       'unaffected', exc_info=True)
        return
    append_run_history(path, row)
