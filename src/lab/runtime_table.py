'''
metric x (log, combo) -> total wall-clock seconds actually spent
computing it, from *_timings.csv's per-stage rows.

A metric's OWN row is usually near-zero - the real cost sits in shared
stages (classical, dv, ...) several metrics depend on (see lab.run's
ProcessMetric.needs). This sums a metric's own row plus its declared
stage dependencies, so the reported number reads as "what this metric
alone would have cost" - not additive across metrics without double-
counting a shared stage two of them both need.

voidsalign accepts either voidsalign3 or voidsalign2's own row, since
different logs' sweeps can be pinned to different registry versions -
same convention lab.cross_log_plots uses for its columns.

Usage:
    python -m lab.runtime_table \\
        --log bpi2013_closed_problems=var/lab/results/bpi2013_closed_problems_20260917-005007_timings.csv \\
        --log rtfm=var/lab/results/rtfm_20260917-054813_timings.csv \\
        --log bpic2020_rfp=var/lab/results/bpic2020_rfp_20260917-091018_timings.csv \\
        --combos inductive_noise20 toothpaste_noise10
'''

import argparse

import pandas as pd

# metric id -> (display label, own row ids, declared stage dependencies).
# Own row ids from two registry versions (voidsalign3/voidsalign2,
# voidmass_process2/voidmass_process) are harmless to sum together - a
# single file only ever has one of the two.
METRIC_TIMING_SPECS = {
    'voidsalign3': ('Void by Skip Alignment', ('voidsalign3', 'voidsalign2'),
                     ('dv', 'executions_cache')),
    'voidsat2': ('Void by Aligned Durations', ('voidsat2',), ('dv', 'aligned_duration_cache')),
    'voidmass_process': ('Void by Process Relative Moves',
                          ('voidmass_process2_lower', 'voidmass_process2_upper',
                           'voidmass_process_lower', 'voidmass_process_upper'), ('classical',)),
}

DEFAULT_METRICS = ('voidsalign3', 'voidsat2', 'voidmass_process')


def load_timings(log_csvs):
    '''{log_label: DataFrame} - one read per *_timings.csv.'''
    return {label: pd.read_csv(path) for label, path in log_csvs.items()}


def _seconds_for(df, combo, row_ids):
    matched = df[(df['combo'] == combo) & (df['metric_or_stage'].isin(row_ids))]
    return matched['seconds'].sum()


def total_seconds(df, combo, metric_id):
    '''Total seconds spent on `metric_id` for `combo`, across every
    degradation dimension/level in `df` - own row(s) plus declared
    stage dependencies (METRIC_TIMING_SPECS).'''
    _label, own_rows, stage_rows = METRIC_TIMING_SPECS[metric_id]
    return _seconds_for(df, combo, own_rows) + _seconds_for(df, combo, stage_rows)


def format_duration(seconds):
    '''0s / 45s / 2m 5s / 1h 2m - whichever units are non-zero, coarsest
    two.'''
    seconds = int(round(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f'{hours}h {minutes}m'
    if minutes:
        return f'{minutes}m {secs}s'
    return f'{secs}s'


def runtime_table(log_timings, combos, metric_ids=DEFAULT_METRICS):
    '''
    DataFrame indexed by (metric label, combo), one column per log in
    `log_timings`, values in seconds (float - format_duration at
    display time, not here, so the numbers stay usable programmatically).
    '''
    index = [(METRIC_TIMING_SPECS[metric_id][0], combo)
             for metric_id in metric_ids for combo in combos]
    data = {}
    for log_label, df in log_timings.items():
        data[log_label] = [total_seconds(df, combo, metric_id)
                            for metric_id in metric_ids for combo in combos]
    return pd.DataFrame(data, index=pd.MultiIndex.from_tuples(index, names=['metric', 'combo']))


def to_markdown_table(df):
    headers = ['metric', 'combo'] + list(df.columns)
    lines = ['| ' + ' | '.join(headers) + ' |',
             '| ' + ' | '.join(['---'] * len(headers)) + ' |']
    for (metric, combo), row in df.iterrows():
        cells = [metric, combo] + [format_duration(v) for v in row]
        lines.append('| ' + ' | '.join(cells) + ' |')
    return '\n'.join(lines)


def to_latex_table(df):
    ncols = 2 + len(df.columns)
    lines = [
        r'\begin{tabular}{' + 'l' * ncols + '}',
        r'\toprule',
        ' & '.join(['Metric', 'Combo'] + list(df.columns)) + r' \\',
        r'\midrule',
    ]
    for (metric, combo), row in df.iterrows():
        cells = [metric, combo] + [format_duration(v) for v in row]
        lines.append(' & '.join(cells) + r' \\')
    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular}')
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='metric x (log, combo) -> total wall-clock seconds, from *_timings.csv files.')
    parser.add_argument('--log', action='append', required=True, dest='logs',
                         metavar='LABEL=PATH',
                         help='a log to include, given as label=timings_csv_path; repeat per log')
    parser.add_argument('--combos', nargs='+', required=True)
    parser.add_argument('--metrics', nargs='+', default=list(DEFAULT_METRICS),
                         choices=list(METRIC_TIMING_SPECS))
    args = parser.parse_args()

    log_csvs = dict(item.split('=', 1) for item in args.logs)
    log_timings = load_timings(log_csvs)
    df = runtime_table(log_timings, args.combos, args.metrics)

    print('# Markdown\n')
    print(to_markdown_table(df))
    print('\n# LaTeX\n')
    print(to_latex_table(df))


if __name__ == '__main__':
    main()
