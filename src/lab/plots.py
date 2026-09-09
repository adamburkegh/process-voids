"""
Dose-response plots for exp_disco_degrade output (log, combo,
degradation_dim, degradation_level columns - the one schema every
degradation sweep in this project produces now, including the claims
fixture via lab.claims_fixture's CLAIMS_COMBOS/CLAIMS_DEGRADATIONS
registration, so there's a single plot function rather than one per
experiment script).

plot_dose_response: weight_coverage / skipprob / salign_coverage /
alignment_coverage_pn / voidmass_subprocess / voidmass_process vs
degradation_level, one figure per (log, facet value), one line per
`line_by` value within it. weight_voidage isn't plotted separately -
it's exactly 1 - weight_coverage, so its own panel would just be a
mirror image with no new signal.

line_by='combo' (default): compare discovery/estimator combos under
one fixed degradation dimension - faceted by (log, degradation_dim).
The natural question for a real discovered-model sweep: "which combo
holds up best as this dimension degrades?"

line_by='degradation_dim': compare degradation dimensions (eg ablation
targets) for one fixed combo - faceted by (log, combo). The natural
question when combo is a single fixed value, eg claims_known: "does
this ablation target behave differently from that one?" - this is what
exp_claims_degrade.py's own dedicated plot function used to do, before
its output converged onto this same schema and made a separate
function unnecessary.

Usage:
    python -m lab.plots var/lab/results/exp_disco_degrade.csv
    python -m lab.plots var/lab/results/claims_degrade.csv --line-by degradation_dim
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
                        fmt: str = 'png', line_by: str = 'combo'):
    """
    Writes one figure per (log, facet value), each with a subplot per
    metric in METRICS, one line per `line_by` value within it - see
    this module's docstring for line_by='combo' vs 'degradation_dim'.
    Pass fmt='pdf' for vector output that drops straight into a LaTeX
    build. Returns the list of paths written. See _exclude_degenerate
    for what's dropped before plotting.

    No fixed y-axis range: voidmass_subprocess/voidmass_process sit in
    a much smaller range (0-0.05ish on the claims fixture) than the
    coverage-style metrics (near 1), so a shared 0-1 range would
    flatten them to nothing - each panel auto-scales to its own data.
    """
    if line_by not in ('combo', 'degradation_dim'):
        raise ValueError(f"line_by must be 'combo' or 'degradation_dim', got {line_by!r}")
    facet_by = 'degradation_dim' if line_by == 'combo' else 'combo'

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    ok = _exclude_degenerate(df)

    for (log, facet_val), group in df.groupby(['log', facet_by]):
        fig, axes = plt.subplots(1, len(METRICS), figsize=(5 * len(METRICS), 4))
        if len(METRICS) == 1:
            axes = [axes]

        group_ok = ok[(ok['log'] == log) & (ok[facet_by] == facet_val)]
        for ax, metric in zip(axes, METRICS):
            # skipprob is a skip probability (higher = worse); every other
            # metric here is a coverage proxy (higher = better) - plot
            # 1-skipprob so all four panels read the same direction.
            label = '1 - skipprob' if metric == 'skipprob' else metric
            for line_val, line_group in group_ok.groupby(line_by):
                line_group = line_group.sort_values('degradation_level')
                if line_group.empty:
                    continue
                y = (1 - line_group[metric] if metric == 'skipprob'
                     else line_group[metric])
                ax.plot(line_group['degradation_level'], y,
                        marker='o', label=line_val)
            ax.set_xlabel(f'{facet_val}-wise degradation level' if facet_by == 'degradation_dim'
                           else 'degradation level')
            ax.set_ylabel(label)
            ax.legend()

        fig.suptitle(f'{log} - {facet_val}')
        fig.tight_layout()

        out_path = out_dir / f'{log}_{facet_val}.{fmt}'
        fig.savefig(out_path)
        plt.close(fig)
        written.append(out_path)

    return written


def main():
    parser = argparse.ArgumentParser(
        description='Plot dose-response curves from exp_disco_degrade output '
                     '(including the claims fixture, run via lab.claims_fixture\'s '
                     'CLAIMS_COMBOS/CLAIMS_DEGRADATIONS registration).')
    parser.add_argument('csv', help='experiment result CSV path')
    parser.add_argument('--out-dir', default='var/lab/results/plots')
    parser.add_argument('--format', default='png', choices=['pdf', 'png', 'svg'],
                         help='output format (default: png; pdf for LaTeX inclusion)')
    parser.add_argument('--line-by', default='combo', choices=['combo', 'degradation_dim'],
                         help="'combo' (default): compare combos under one degradation "
                              "dimension. 'degradation_dim': compare degradation "
                              "dimensions for one combo - use this for claims-fixture "
                              "output, where combo is always the single fixed "
                              "'claims_known' value.")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    written = plot_dose_response(df, out_dir=args.out_dir, fmt=args.format, line_by=args.line_by)
    for path in written:
        print(f'Wrote {path}')


if __name__ == '__main__':
    main()
