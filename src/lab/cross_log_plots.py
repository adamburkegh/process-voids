'''
Cross-log dose-response plots.

lab.plots compares combos or degradation dimensions within ONE log's
own result CSV. This compares the SAME combo and degradation dimension
across DIFFERENT logs' result CSVs, side by side - one line per log,
since each log lives in its own result file and, for some metrics, its
own version of that file (a sweep run before a metric landed has no
column for it at all).

One figure per (combo, degradation_dim) pair requested, one panel per
metric in METRIC_SPECS, one line per log within a panel. A log missing
a metric's column, or with no rows for that combo/dim, is skipped for
that line - not an error - same convention as lab.check_run's missing-
column handling. A <base>_lower/<base>_upper pair is plotted as a
midpoint line with a shaded band, matching lab.plots' banded-panel
convention.

Usage:
    python -m lab.cross_log_plots \\
        --log bpi2013_closed_problems=var/lab/results/bpi2013_closed_problems_20260917-005007.csv \\
        --log rtfm=var/lab/results/rtfm_20260917-054813.csv \\
        --log bpic2020_rfp=var/lab/results/bpic2020_rfp_20260917-091018.csv \\
        --combos inductive_noise20 toothpaste_noise10 \\
        --degradation-dims trace activity_frequency_gradual
'''

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from lab.plots import _exclude_degenerate

# label -> (panel title, columns). A 1-tuple is a plain metric; a
# (lower, upper) pair is one banded panel, plotted as its midpoint with
# a shaded band - see plot_cross_log.
METRIC_SPECS = {
    'voidsalign3': ('voidsalign3', ('voidsalign3',)),
    'voidsat2': ('voidsat2', ('voidsat2',)),
    'voidmass_process': ('voidmass_process', ('voidmass_process_lower', 'voidmass_process_upper')),
}

DEFAULT_METRICS = ('voidsalign3', 'voidsat2', 'voidmass_process')


def load_logs(log_csvs):
    '''{log_label: DataFrame} - one read per CSV, degenerate rows
    (status != 'ok', or degradation_level == 1.0) dropped up front, the
    same filter lab.plots applies before plotting.'''
    return {label: _exclude_degenerate(pd.read_csv(path)) for label, path in log_csvs.items()}


def _series_for(df, combo, degradation_dim, columns):
    '''Rows for one (combo, degradation_dim), sorted by degradation_level
    - or None if any of `columns` isn't in this log's CSV (a sweep
    predating that metric) or no row matches, so the caller can skip
    the line rather than error.'''
    if any(col not in df.columns for col in columns):
        return None
    rows = df[(df['combo'] == combo) & (df['degradation_dim'] == degradation_dim)]
    if rows.empty:
        return None
    return rows.sort_values('degradation_level')


def _draw_series(ax, series, cols, label):
    '''Draws one line for `label` on `ax` from `series` (from
    _series_for) - a plain line for a 1-column metric, a midpoint line
    with a shaded band for a (lower, upper) pair. Shared by both
    panel_by modes below, which differ only in what's held fixed per
    panel and what varies per line.'''
    x = series['degradation_level']
    if len(cols) == 1:
        ax.plot(x, series[cols[0]], marker='o', label=label)
    else:
        lower_col, upper_col = cols
        lower = pd.to_numeric(series[lower_col])
        upper = pd.to_numeric(series[upper_col])
        midpoint = (lower + upper) / 2
        line, = ax.plot(x, midpoint, marker='o', label=label)
        ax.fill_between(x, lower, upper, alpha=0.2, color=line.get_color())


def plot_cross_log(log_dfs, combos, degradation_dims, metric_ids=DEFAULT_METRICS,
                    out_dir='var/lab/plots/cross_log', fmt='png', ylim=None,
                    panel_by='metric'):
    '''
    Writes one figure per (combo, degradation_dim) pair. Returns the
    list of paths written.

    panel_by='metric' (default): figures are one per (combo,
    degradation_dim); one panel per metric in `metric_ids`, one line
    per log in `log_dfs` within it - "how does each log compare on
    this metric, for this ref model".
    panel_by='log': same figures; one panel per log in `log_dfs`, one
    line per metric in `metric_ids` within it - "how do this log's
    metrics compare to each other".
    panel_by='log_metric': figures are one per degradation_dim only
    (combo is no longer a figure axis); one panel per (log, metric)
    pair, one line per combo in `combos` within it - "how do the ref
    models compare, for this log and metric".

    Whichever mode, a log missing a metric's column, or with no rows
    for the combo/dim in question, is skipped for that line - not an
    error - same convention as lab.check_run's missing-column handling.

    ylim: (min, max) applied to every panel, or None (default) to let
    each panel autoscale to its own data - lab.plots' own default,
    since a shared range can flatten a small-range metric like
    voidmass_process (0-0.05ish on some fixtures) to a near-flat line
    against a coverage-style metric near 1. Pass ylim=(0, 1) for
    consistent, directly-comparable axes across panels/figures instead,
    accepting that tradeoff.
    '''
    if panel_by not in ('metric', 'log', 'log_metric'):
        raise ValueError(f"panel_by must be 'metric', 'log' or 'log_metric', got {panel_by!r}")
    specs = [METRIC_SPECS[metric_id] for metric_id in metric_ids]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    def _finish(fig, axes, title, out_path):
        for ax in axes:
            ax.set_xlabel('degradation level')
            if ylim is not None:
                ax.set_ylim(*ylim)
            ax.legend()
        fig.suptitle(title)
        fig.tight_layout()
        fig.savefig(out_path)
        plt.close(fig)
        written.append(out_path)

    if panel_by == 'log_metric':
        panels = [(log_label, df, label, cols)
                  for log_label, df in log_dfs.items() for label, cols in specs]
        for dim in degradation_dims:
            fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 4))
            if len(panels) == 1:
                axes = [axes]
            for ax, (log_label, df, label, cols) in zip(axes, panels):
                for combo in combos:
                    series = _series_for(df, combo, dim, cols)
                    if series is not None:
                        _draw_series(ax, series, cols, combo)
                ax.set_title(f'{log_label} - {label}')
                ax.set_ylabel(label)
            _finish(fig, axes, dim, out_dir / f'{dim}_by_log_metric.{fmt}')
        return written

    panels = specs if panel_by == 'metric' else list(log_dfs.items())

    for combo in combos:
        for dim in degradation_dims:
            fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 4))
            if len(panels) == 1:
                axes = [axes]
            for ax, panel in zip(axes, panels):
                if panel_by == 'metric':
                    label, cols = panel
                    for log_label, df in log_dfs.items():
                        series = _series_for(df, combo, dim, cols)
                        if series is not None:
                            _draw_series(ax, series, cols, log_label)
                    ax.set_ylabel(label)
                else:
                    log_label, df = panel
                    for label, cols in specs:
                        series = _series_for(df, combo, dim, cols)
                        if series is not None:
                            _draw_series(ax, series, cols, label)
                    ax.set_title(log_label)
                    ax.set_ylabel('value')
            suffix = '' if panel_by == 'metric' else '_by_log'
            _finish(fig, axes, f'{combo} - {dim}', out_dir / f'{combo}_{dim}{suffix}.{fmt}')

    return written


def main():
    parser = argparse.ArgumentParser(
        description='Cross-log dose-response plots: same combo and degradation '
                     'dimension, one line per log, each read from its own result CSV.')
    parser.add_argument('--log', action='append', required=True, dest='logs',
                         metavar='LABEL=PATH',
                         help='a log to plot, given as label=csv_path; repeat for each log')
    parser.add_argument('--combos', nargs='+', required=True)
    parser.add_argument('--degradation-dims', nargs='+', required=True)
    parser.add_argument('--metrics', nargs='+', default=list(DEFAULT_METRICS),
                         choices=list(METRIC_SPECS))
    parser.add_argument('--out-dir', default='var/lab/plots/cross_log')
    parser.add_argument('--format', default='png', choices=['pdf', 'png', 'svg'])
    parser.add_argument('--ylim', type=float, nargs=2, default=None, metavar=('MIN', 'MAX'),
                         help='fixed y-axis range for every panel, eg --ylim 0 1; '
                              'default autoscales each panel to its own data')
    parser.add_argument('--panel-by', default='metric', choices=['metric', 'log', 'log_metric'],
                         help="'metric' (default): panel per metric, line per log. "
                              "'log': panel per log, line per metric. 'log_metric': panel "
                              "per (log, metric) pair, line per combo - figures then split "
                              "by degradation_dim only, not combo.")
    args = parser.parse_args()

    log_csvs = dict(item.split('=', 1) for item in args.logs)
    log_dfs = load_logs(log_csvs)
    written = plot_cross_log(log_dfs, args.combos, args.degradation_dims,
                              metric_ids=args.metrics, out_dir=args.out_dir, fmt=args.format,
                              ylim=tuple(args.ylim) if args.ylim else None, panel_by=args.panel_by)
    for path in written:
        print(f'Wrote {path}')


if __name__ == '__main__':
    main()
