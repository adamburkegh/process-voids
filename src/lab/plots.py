"""
Dose-response plots for exp_disco_degrade output (log, combo,
degradation_dim, degradation_level columns - the one schema every
degradation sweep in this project produces, including the claims
fixture via lab.claims_fixture's CLAIMS_COMBOS/CLAIMS_DEGRADATIONS
registration, so there's a single plot function rather than one per
experiment script).

plot_dose_response: weight_coverage / skipprob / salign_coverage /
alignment_coverage_pn / voidmass_subprocess / voidmass_process vs
degradation_level, one figure per (log, facet value), one line per
`line_by` value within it. weight_voidage isn't plotted separately -
it's exactly 1 - weight_coverage, so its own panel would just be a
mirror image with no new signal.

alignment_coverage_pn/voidmass_subprocess/voidmass_process are each a
_lower/_upper bound pair, not a single column (the timed-out-variant
bound - see lab.exp_disco_degrade.CLASSICAL_METRIC_KEYS): plotted as the
midpoint line with a shaded band between the two. A cell with no
timed-out variants has lower == upper, so the band collapses to a plain
line with no special-casing needed.

line_by='combo' (default): compare discovery/estimator combos under
one fixed degradation dimension - faceted by (log, degradation_dim).
The natural question for a real discovered-model sweep: "which combo
holds up best as this dimension degrades?"

line_by='degradation_dim': compare degradation dimensions (eg ablation
targets) for one fixed combo - faceted by (log, combo). The natural
question when combo is a single fixed value, eg claims_known: "does
this ablation target behave differently from that one?"

average_over_nodes + --average-nodes: a THIRD view, alongside (not
instead of) the root-level one above - collapses the per-node CSV
(exp_disco_degrade's *_nodes.csv) into one smoothed curve per cell by
averaging each metric across every non-Tau node. Root-only is one
number with no subprocess information; plotting every node as its own
line is unreadable past a handful of nodes; this is the middle ground.

Usage:
    python -m lab.plots var/lab/results/exp_disco_degrade.csv
    python -m lab.plots var/lab/results/claims_degrade.csv --line-by degradation_dim
    python -m lab.plots var/lab/results/exp_disco_degrade.csv --out-dir var/lab/results/plots
    python -m lab.plots var/lab/results/exp_disco_degrade_nodes.csv --average-nodes
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# (label, columns): columns is a 1-tuple for a plain value or a 2-tuple
# (lower, upper) for a bound pair - see module docstring. label is the
# axis text (and the lookup key for the skipprob-inversion special case
# below), independent of the column name(s) it reads.
METRICS = [
    ('weight_coverage', ('weight_coverage',)),
    ('skipprob', ('skipprob',)),
    ('salign_coverage', ('salign_coverage',)),
    ('alignment_coverage_pn', ('alignment_coverage_pn_lower', 'alignment_coverage_pn_upper')),
    ('voidmass_subprocess', ('voidmass_subprocess_lower', 'voidmass_subprocess_upper')),
    ('voidmass_process', ('voidmass_process_lower', 'voidmass_process_upper')),
]


def _metric_columns():
    """Every column METRICS reads, flattened - the one place the flat
    list lives, so average_over_nodes can't drift from what
    plot_dose_response itself actually plots."""
    return [col for _label, cols in METRICS for col in cols]


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

    lab.params.ALL_LEVELS excludes 1.0 (a wasted, uninformative compute
    point, not just an unplottable one), so for a default run this
    filter is defensive rather than load-bearing - it matters for a
    result CSV from an explicit --levels including 1.0.
    """
    return df[(df['status'] == 'ok') & (df['degradation_level'] != 1.0)]


def average_over_nodes(node_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapses exp_disco_degrade's per-node CSV (one row per (log, combo,
    degradation_dim, degradation_level, node_id)) into one row per
    (log, combo, degradation_dim, degradation_level) by averaging each
    metric across every non-Tau node - a smoothed, subprocess-aware
    middle ground between the root-only view (one number, no subprocess
    information at all) and plotting every node as its own line (too
    noisy to read with more than a handful of nodes). The root-level
    plots remain the primary view - this is an additional one, not a
    replacement.

    Tau nodes are excluded, same convention as
    process_voids.coveragemass's mandatory_node_count/total_node_count:
    a Tau leaf represents "do nothing", not a thing whose coverage/void
    reading should pull the average toward its own degenerate values
    (eg skipprob=1.0, weight_coverage=0.0 on every Tau node, regardless
    of how the rest of the tree is actually behaving).

    Feed the result straight into plot_dose_response - a 'status'='ok'
    column is added so the output matches the root-level CSV's shape
    exactly (the per-node CSV has no status column of its own; a node
    CSV only ever holds successful cells, see run_disco_degrade).
    degradation_level=1.0 rows are dropped before averaging (rather than
    left to _exclude_degenerate downstream) since a single degenerate
    node's reading would otherwise contaminate that cell's average even
    when other nodes in it look fine.
    """
    metric_cols = _metric_columns()
    group_cols = ['log', 'combo', 'degradation_dim', 'degradation_level']

    filtered = node_df[(node_df['node_type'] != 'Tau') & (node_df['degradation_level'] != 1.0)]
    averaged = filtered.groupby(group_cols)[metric_cols].mean().reset_index()
    averaged['status'] = 'ok'
    return averaged


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
        for ax, (metric_label, cols) in zip(axes, METRICS):
            # skipprob is a skip probability (higher = worse); every other
            # metric here is a coverage proxy (higher = better) - plot
            # 1-skipprob so all four panels read the same direction.
            label = '1 - skipprob' if metric_label == 'skipprob' else metric_label
            for line_val, line_group in group_ok.groupby(line_by):
                line_group = line_group.sort_values('degradation_level')
                if line_group.empty:
                    continue
                x = line_group['degradation_level']
                if len(cols) == 1:
                    y = (1 - line_group[cols[0]] if metric_label == 'skipprob'
                         else line_group[cols[0]])
                    ax.plot(x, y, marker='o', label=line_val)
                else:
                    lower_col, upper_col = cols
                    # to_numeric: a column that HELD a non-numeric value
                    # anywhere (eg an excluded level=1.0 row elsewhere in
                    # the same CSV) can stay object-dtype even after
                    # _exclude_degenerate drops that row - fill_between's
                    # own isfinite check can't handle object dtype even
                    # when every remaining value is a real float, unlike
                    # ax.plot below, which tolerates it.
                    lower = pd.to_numeric(line_group[lower_col])
                    upper = pd.to_numeric(line_group[upper_col])
                    midpoint = (lower + upper) / 2
                    line, = ax.plot(x, midpoint, marker='o', label=line_val)
                    ax.fill_between(x, lower, upper, alpha=0.2, color=line.get_color())
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
    parser.add_argument('--average-nodes', action='store_true',
                         help='Treat csv as exp_disco_degrade\'s per-node output '
                              '(the *_nodes.csv file, not the root-level one) and '
                              'average each metric across every non-Tau node per cell '
                              'before plotting - a smoothed, subprocess-aware curve. '
                              'Does not replace plotting the root-level CSV directly, '
                              'just an additional view.')
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    if args.average_nodes:
        df = average_over_nodes(df)
    written = plot_dose_response(df, out_dir=args.out_dir, fmt=args.format, line_by=args.line_by)
    for path in written:
        print(f'Wrote {path}')


if __name__ == '__main__':
    main()
