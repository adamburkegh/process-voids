'''
Standard deviation across degradation levels, for the trace-degradation
dimension specifically.

The trace-degradation dose-response curves are informative to look at
but too flat/uneventful to earn a figure in a paper. This collapses
each curve (one per metric x log x combo, restricted to
degradation_dim='trace') to its standard deviation across levels - the
same "how much did this line actually move" question a reader would
ask of the graph, answered as one number instead of a plot, so a
minimum and maximum across every line can be quoted in place of showing
them all.

Reuses lab.cross_log_plots' METRIC_SPECS/_series_for/load_logs rather
than redefining them - same metrics, same log-loading, same missing-
column-is-skipped convention, just a different reduction over the
result.

Usage:
    python -m lab.trace_variability \\
        --log bpi2013_closed_problems=var/lab/results/bpi2013_closed_problems_20260917-005007.csv \\
        --log rtfm=var/lab/results/rtfm_20260917-054813.csv \\
        --log bpic2020_rfp=var/lab/results/bpic2020_rfp_20260917-091018.csv \\
        --combos inductive_noise20 toothpaste_noise10
'''

import argparse

import pandas as pd

from lab.cross_log_plots import DEFAULT_METRICS, METRIC_SPECS, _series_for, load_logs


def std_dev_for(df, combo, metric_id):
    '''Standard deviation of `metric_id`'s values across every
    degradation_dim='trace' row for `combo` in `df` - the midpoint
    series for a banded (lower, upper) metric, same as what
    lab.cross_log_plots actually draws as the line. None if the column
    isn't in this log's CSV or there's no matching row, same skip-not-
    error convention as lab.cross_log_plots.'''
    _label, cols = METRIC_SPECS[metric_id]
    series = _series_for(df, combo, 'trace', cols)
    if series is None:
        return None
    if len(cols) == 1:
        values = series[cols[0]]
    else:
        lower = pd.to_numeric(series[cols[0]])
        upper = pd.to_numeric(series[cols[1]])
        values = (lower + upper) / 2
    return values.std()


def variability_table(log_dfs, combos, metric_ids=DEFAULT_METRICS):
    '''DataFrame indexed by (metric label, combo), one column per log in
    `log_dfs`, values are std_dev_for's result - NaN where a log has no
    data for that metric/combo, not an error.'''
    index = [(METRIC_SPECS[metric_id][0], combo)
             for metric_id in metric_ids for combo in combos]
    data = {}
    for log_label, df in log_dfs.items():
        data[log_label] = [std_dev_for(df, combo, metric_id)
                            for metric_id in metric_ids for combo in combos]
    return pd.DataFrame(data, index=pd.MultiIndex.from_tuples(index, names=['metric', 'combo']))


def summarize(df):
    '''(min_value, (metric, combo, log) at the min, max_value, (metric,
    combo, log) at the max) - NaN cells (a log missing that metric)
    excluded from both. stack() flattens the (metric, combo) row index
    and the log column index into one 3-level index, so idxmin/idxmax
    already return the (metric, combo, log) triple directly.'''
    stacked = df.stack()
    min_key = stacked.idxmin()
    max_key = stacked.idxmax()
    return stacked[min_key], min_key, stacked[max_key], max_key


def main():
    parser = argparse.ArgumentParser(
        description="Standard deviation of each metric's values across trace-degradation "
                     'levels, one line per (metric, log, combo) - a numeric stand-in for the '
                     'trace-dimension dose-response plots.')
    parser.add_argument('--log', action='append', required=True, dest='logs',
                         metavar='LABEL=PATH',
                         help='a log to include, given as label=csv_path; repeat per log')
    parser.add_argument('--combos', nargs='+', required=True)
    parser.add_argument('--metrics', nargs='+', default=list(DEFAULT_METRICS),
                         choices=list(METRIC_SPECS))
    args = parser.parse_args()

    log_csvs = dict(item.split('=', 1) for item in args.logs)
    log_dfs = load_logs(log_csvs)
    df = variability_table(log_dfs, args.combos, args.metrics)

    print(df.to_string())
    min_val, min_where, max_val, max_where = summarize(df)
    print(f'\nmin std dev: {min_val:.4f} at {min_where}')
    print(f'max std dev: {max_val:.4f} at {max_where}')


if __name__ == '__main__':
    main()
