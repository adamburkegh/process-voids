'''
Reports which lab.metric_registry ids have real (non-null) values
recorded for which logs, across every root-level result CSV ever
written - "has this been collected at all", not "what's its latest
value" or "is it still live" (retired ids can still show Done from
before they were retired).

A result CSV's 'log' column (not its filename) says which log a row
belongs to - ad hoc/probe/smoke runs count exactly like a named
--run sweep. *_nodes.csv/*_timings.csv companions are not scanned
separately (see _root_level_result_csvs); a CSV with no 'log' column,
or none of whose columns match a registry id (eg lab.mass_term_probe's
own bespoke output), contributes nothing - not an error.
'''

from pathlib import Path

import pandas as pd

from lab.metric_registry import METRICS

RESULTS_DIR = 'var/lab/results'
METRIC_IDS = list(METRICS)


def _root_level_result_csvs(results_dir=RESULTS_DIR):
    '''Every result CSV except the *_nodes.csv/*_timings.csv companions
    a root-level run also writes - those carry the same metric columns
    at finer granularity or none at all, so scanning them too would only
    double-count.'''
    return sorted(
        p for p in Path(results_dir).glob('*.csv')
        if not p.stem.endswith('_nodes') and not p.stem.endswith('_timings')
    )


def collected_metrics(results_dir=RESULTS_DIR):
    '''
    Returns a DataFrame indexed by metric id, one boolean column per log
    name seen in any scanned CSV's 'log' column: True if some CSV has at
    least one non-null value for that (metric, log) pair.
    '''
    collected = {}  # log -> set of metric ids seen with real data

    for csv_path in _root_level_result_csvs(results_dir):
        df = pd.read_csv(csv_path)
        if 'log' not in df.columns:
            continue
        present_ids = [m for m in METRIC_IDS if m in df.columns]
        if not present_ids:
            continue
        for log_name, group in df.groupby('log'):
            seen = collected.setdefault(log_name, set())
            for metric_id in present_ids:
                if group[metric_id].notna().any():
                    seen.add(metric_id)

    logs = sorted(collected)
    return pd.DataFrame(
        {log: [metric_id in collected[log] for metric_id in METRIC_IDS] for log in logs},
        index=METRIC_IDS,
    )


def to_markdown_table(df):
    '''Hand-rolled Markdown table (no tabulate dependency, not currently
    installed) - True/False rendered as 'Done'/'' per this first cut's
    own spec. LaTeX output can follow later if wanted.'''
    headers = ['metric'] + list(df.columns)
    lines = ['| ' + ' | '.join(headers) + ' |',
             '| ' + ' | '.join(['---'] * len(headers)) + ' |']
    for metric_id, row in df.iterrows():
        cells = [metric_id] + ['Done' if v else '' for v in row]
        lines.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def main():
    df = collected_metrics()
    print(to_markdown_table(df))


if __name__ == '__main__':
    main()
