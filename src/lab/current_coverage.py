'''
"Are we current" for a log: which of today's LIVE metrics does its
MOST RECENT sweep actually have.

lab.collection_report answers a different question on purpose - "has
this metric ever been collected for this log, across every result CSV
ever written", including retired ids and long-stale runs (its own
docstring: "not 'what's its latest value' or 'is it still live'"). That
view is right for an archive; it is the wrong view for deciding what to
rerun, since a log can show Done there from a sweep run before a metric
even existed in its current form (eg voidsalign2 vs voidsalign3). This
module looks only at each log's newest file and only at metrics
presently live in the registry - the same "did we lose track of a gap"
question this session kept re-answering by hand, reading a CSV's header
and comparing it against what the registry currently expects.

"Most recent" is by file modification time, not by parsing a
_timestamped() filename - the two normally agree, but mtime is what
actually reflects when the data was produced regardless of naming.
'''

import argparse
from pathlib import Path

import pandas as pd

from lab.collection_report import _root_level_result_csvs, to_markdown_table
from lab.metric_registry import METRICS
from lab.plots import _RUNNER

RESULTS_DIR = 'var/lab/results'


def live_metric_ids():
    '''Metric ids presently live in the registry AND emitted by
    exp_disco_degrade - the runner this module's result CSVs come from.
    Excludes exp_surprise/exp_voidmass's own metrics, which live in a
    separate pipeline with a separate output schema and would otherwise
    show as a permanent, uninformative wall of blanks here.'''
    return [metric_id for metric_id, metric in METRICS.items()
            if metric.status == 'live' and _RUNNER in metric.scripts]


def latest_csv_per_log(results_dir=RESULTS_DIR):
    '''{log_name: Path} - the most recently modified root-level result
    CSV containing that log's rows, one entry per log name seen in any
    scanned CSV's 'log' column. A CSV with no 'log' column contributes
    nothing, same as lab.collection_report.'''
    latest = {}  # log_name -> (mtime, path)
    for csv_path in _root_level_result_csvs(results_dir):
        df = pd.read_csv(csv_path, usecols=lambda c: True)
        if 'log' not in df.columns:
            continue
        mtime = csv_path.stat().st_mtime
        for log_name in df['log'].unique():
            current = latest.get(log_name)
            if current is None or mtime > current[0]:
                latest[log_name] = (mtime, csv_path)
    return {log_name: path for log_name, (_mtime, path) in latest.items()}


def current_coverage(results_dir=RESULTS_DIR):
    '''
    DataFrame indexed by live metric id, one boolean column per log:
    True if that log's newest result CSV has at least one non-null
    value for that metric. A metric column absent from that file reads
    False, even if an older file for the same log once had it.
    '''
    metric_ids = live_metric_ids()
    latest = latest_csv_per_log(results_dir)
    columns = {}
    for log_name, csv_path in sorted(latest.items()):
        df = pd.read_csv(csv_path)
        group = df[df['log'] == log_name]
        columns[log_name] = [
            metric_id in group.columns and group[metric_id].notna().any()
            for metric_id in metric_ids
        ]
    return pd.DataFrame(columns, index=metric_ids)


def source_files(results_dir=RESULTS_DIR):
    '''{log_name: Path} - which file current_coverage's columns are
    reporting on, for a caller that wants to show recency alongside
    the table (eg the CLI's own printout).'''
    return latest_csv_per_log(results_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Which live metrics each log's newest result CSV has.")
    parser.add_argument('--results-dir', default=RESULTS_DIR,
                         help='directory of root-level result CSVs (default: %(default)s); '
                              "the published results/ folder works too")
    args = parser.parse_args()

    df = current_coverage(args.results_dir)
    print(to_markdown_table(df))
    print()
    print('as of:')
    for log_name, path in sorted(source_files(args.results_dir).items()):
        print(f'  {log_name}: {path}')


if __name__ == '__main__':
    main()
