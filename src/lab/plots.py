"""
Dose-response plots for degradation-sweep experiment output.

plot_dose_response: exp_disco_degrade's shape (log, degradation_dim,
combo, degradation_level columns) - weight_coverage / skipprob /
salign_coverage / alignment_coverage_pn / voidmass_subprocess /
voidmass_process vs degradation_level, one line per discovery+estimator
combo, faceted by (log, degradation dimension). weight_voidage isn't
plotted separately - it's exactly 1 - weight_coverage, so its own panel
would just be a mirror image with no new signal.

plot_claims_degrade: exp_claims_degrade's shape (target, n_drop_cases
columns) - the same three metrics plus voidmass_subprocess/
voidmass_process vs n_drop_cases, one line per ablation target. Not
sharing a y-axis range with the coverage metrics: voidmass values sit
in a much smaller range (0-0.05ish here) than the coverage-style
metrics (near 1), so a shared 0-1 range would flatten them to nothing.

main() dispatches on the CSV's columns - pass either experiment's
output and it picks the right plot.

Usage:
    python -m lab.plots var/lab/results/exp_disco_degrade.csv
    python -m lab.plots var/lab/results/exp_claims_degrade.csv
    python -m lab.plots var/lab/results/exp_disco_degrade.csv --out-dir var/lab/results/plots
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

METRICS = ['weight_coverage', 'skipprob', 'salign_coverage', 'alignment_coverage_pn',
           'voidmass_subprocess', 'voidmass_process']


def _exclude_degenerate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rows fit to plot: status == 'ok' and degradation_level != 1.0.

    Level 1.0 is excluded on ANY dimension (activity, activity_gradual,
    or trace) because the degraded log is then literally empty (every
    activity or every case dropped), and every metric falls through its
    own "no data" default there rather than measuring anything. Those
    defaults don't even agree with each other (some land at their
    "perfect" identity value, salign_coverage at its "worst"), so the
    level=1.0 point is actively misleading on a dose-response curve, not
    just an uninteresting edge case.

    lab.params.ALL_LEVELS no longer includes 1.0 at all (a wasted,
    uninformative compute point, not just an unplottable one), so this
    filter is now defensive rather than load-bearing for a default run -
    it still matters for older result CSVs on disk, or a call that
    passes an explicit --levels including 1.0.
    """
    return df[(df['status'] == 'ok') & (df['degradation_level'] != 1.0)]


def plot_dose_response(df: pd.DataFrame, out_dir: str = 'var/lab/results/plots',
                        fmt: str = 'png'):
    """
    Writes one figure per (log, degradation_dim), each with a subplot
    per metric in METRICS, one line per combo. Pass fmt='pdf' for
    vector output that drops straight into a LaTeX build. Returns the
    list of paths written. See _exclude_degenerate for what's dropped
    before plotting.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    ok = _exclude_degenerate(df)

    for (log, dim), group in df.groupby(['log', 'degradation_dim']):
        fig, axes = plt.subplots(1, len(METRICS), figsize=(5 * len(METRICS), 4))
        if len(METRICS) == 1:
            axes = [axes]

        group_ok = ok[(ok['log'] == log) & (ok['degradation_dim'] == dim)]
        for ax, metric in zip(axes, METRICS):
            # skipprob is a skip probability (higher = worse); every other
            # metric here is a coverage proxy (higher = better) - plot
            # 1-skipprob so all four panels read the same direction.
            label = '1 - skipprob' if metric == 'skipprob' else metric
            for combo, combo_group in group_ok.groupby('combo'):
                combo_group = combo_group.sort_values('degradation_level')
                if combo_group.empty:
                    continue
                y = (1 - combo_group[metric] if metric == 'skipprob'
                     else combo_group[metric])
                ax.plot(combo_group['degradation_level'], y,
                        marker='o', label=combo)
            ax.set_xlabel(f'{dim}-wise degradation level')
            ax.set_ylabel(label)
            ax.set_ylim(-0.05, 1.05)
            ax.legend()

        fig.suptitle(f'{log} - {dim}-wise degradation')
        fig.tight_layout()

        out_path = out_dir / f'{log}_{dim}.{fmt}'
        fig.savefig(out_path)
        plt.close(fig)
        written.append(out_path)

    return written


CLAIMS_METRICS = ['weight_coverage', 'skipprob', 'salign_coverage', 'alignment_coverage_pn',
                   'voidmass_subprocess', 'voidmass_process']


def plot_claims_degrade(df: pd.DataFrame, out_dir: str = 'var/lab/results/plots',
                         fmt: str = 'png'):
    """
    One figure for exp_claims_degrade output: a subplot per metric in
    CLAIMS_METRICS, one line per ablation target, x-axis n_drop_cases.
    Returns the list of paths written (always one, for consistency with
    plot_dose_response's return shape).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ok = df[df['status'] == 'ok']

    fig, axes = plt.subplots(1, len(CLAIMS_METRICS), figsize=(5 * len(CLAIMS_METRICS), 4))
    if len(CLAIMS_METRICS) == 1:
        axes = [axes]

    for ax, metric in zip(axes, CLAIMS_METRICS):
        # skipprob is a skip probability (higher = worse); every other
        # metric here is a coverage/void proxy in its own natural
        # direction - plot 1-skipprob so its panel reads the same way
        # as weight_coverage/salign_coverage (higher = better);
        # voidmass panels are deliberately left as-is (higher = worse,
        # the opposite direction) since they're not coverage proxies.
        label = '1 - skipprob' if metric == 'skipprob' else metric
        for target, group in ok.groupby('target'):
            group = group.sort_values('n_drop_cases')
            if group.empty:
                continue
            y = (1 - group[metric] if metric == 'skipprob' else group[metric])
            ax.plot(group['n_drop_cases'], y, marker='o', label=target)
        ax.set_xlabel('n_drop_cases')
        ax.set_ylabel(label)
        ax.legend()

    fig.suptitle('Claims fixture - degradation vs the known generating tree')
    fig.tight_layout()

    out_path = out_dir / f'claims_degrade.{fmt}'
    fig.savefig(out_path)
    plt.close(fig)
    return [out_path]


def main():
    parser = argparse.ArgumentParser(
        description='Plot dose-response curves from exp_disco_degrade or '
                     'exp_claims_degrade output (auto-detected from the CSV columns).')
    parser.add_argument('csv', help='experiment result CSV path')
    parser.add_argument('--out-dir', default='var/lab/results/plots')
    parser.add_argument('--format', default='png', choices=['pdf', 'png', 'svg'],
                         help='output format (default: png; pdf for LaTeX inclusion)')
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    if 'target' in df.columns and 'n_drop_cases' in df.columns:
        written = plot_claims_degrade(df, out_dir=args.out_dir, fmt=args.format)
    else:
        written = plot_dose_response(df, out_dir=args.out_dir, fmt=args.format)
    for path in written:
        print(f'Wrote {path}')


if __name__ == '__main__':
    main()
