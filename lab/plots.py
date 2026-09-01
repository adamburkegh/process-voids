"""
Dose-response plots for exp_disco_degrade output: weight_coverage /
skipprob vs degradation_level, one line per discovery+estimator combo,
faceted by (log, degradation dimension).

Usage:
    python -m lab.plots var/lab/results/exp_disco_degrade.csv
    python -m lab.plots var/lab/results/exp_disco_degrade.csv --out-dir var/lab/results/plots
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

METRICS = ['weight_coverage', 'skipprob']


def plot_dose_response(df: pd.DataFrame, out_dir: str = 'var/lab/results/plots',
                        fmt: str = 'png'):
    """
    Writes one figure per (log, degradation_dim), each with a subplot
    per metric in METRICS, one line per combo. Pass fmt='pdf' for
    vector output that drops straight into a LaTeX build. Returns the
    list of paths written.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    ok = df[df['status'] == 'ok']

    for (log, dim), group in df.groupby(['log', 'degradation_dim']):
        fig, axes = plt.subplots(1, len(METRICS), figsize=(5 * len(METRICS), 4))
        if len(METRICS) == 1:
            axes = [axes]

        group_ok = ok[(ok['log'] == log) & (ok['degradation_dim'] == dim)]
        for ax, metric in zip(axes, METRICS):
            for combo, combo_group in group_ok.groupby('combo'):
                combo_group = combo_group.sort_values('degradation_level')
                if combo_group.empty:
                    continue
                ax.plot(combo_group['degradation_level'], combo_group[metric],
                        marker='o', label=combo)
            ax.set_xlabel(f'{dim}-wise degradation level')
            ax.set_ylabel(metric)
            ax.set_ylim(-0.05, 1.05)
            ax.legend()

        fig.suptitle(f'{log} - {dim}-wise degradation')
        fig.tight_layout()

        out_path = out_dir / f'{log}_{dim}.{fmt}'
        fig.savefig(out_path)
        plt.close(fig)
        written.append(out_path)

    return written


def main():
    parser = argparse.ArgumentParser(
        description='Plot dose-response coverage curves from exp_disco_degrade output.')
    parser.add_argument('csv', help='exp_disco_degrade CSV path')
    parser.add_argument('--out-dir', default='var/lab/results/plots')
    parser.add_argument('--format', default='png', choices=['pdf', 'png', 'svg'],
                         help='output format (default: png; pdf for LaTeX inclusion)')
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    written = plot_dose_response(df, out_dir=args.out_dir, fmt=args.format)
    for path in written:
        print(f'Wrote {path}')


if __name__ == '__main__':
    main()
