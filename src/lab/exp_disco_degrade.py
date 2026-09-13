"""
Experiment: discovered model vs (degraded) log.

A thin wrapper over lab.run - the metric roster (ALL_METRICS), cell loop
and CellContext wiring all live there now (see lab.run's own module
docstring). This module keeps its own stable CLI/signature (no
--metrics/--no-degradation - the full roster, dim x level sweep every
time, except --exclude-metric below) and its original default output
filenames, for anyone still invoking it directly:

    python -m lab.exp_disco_degrade --run smoke
    python -m lab.exp_disco_degrade logs/rtfm_fine_appeal.xes.gz --levels 0.0 0.5
    python -m lab.exp_disco_degrade --run rtfm --exclude-metric voidsat

--exclude-metric skips scoring the given metric id(s) entirely for a run
(eg voidsat, whose per-trace cost is prohibitive on a high case-count
log) - the inclusive complement of ALL_METRICS is what actually gets
passed to lab.run.run's own `metrics=`.

lab.claims_fixture's CLAIMS_COMBOS/CLAIMS_DEGRADATIONS registers that
fixture's known tree and its named ablation targets (appeal_seq/
loop_block/assess) as combo/degradation-dimension names usable here.
"""

import argparse
import logging
from dataclasses import replace
from datetime import datetime
from pathlib import Path

from lab import run as lab_run
from lab.claims_fixture import CLAIMS_COMBOS, CLAIMS_DEGRADATIONS
from lab.logconfig import configure, enable_skipalignments_debug
from lab.params import ALL_COMBOS, ALL_DEGRADATIONS, ALL_LEVELS
from lab.run import (
    ALL_METRICS, CLASSICAL_ALIGNMENT_TIMEOUT, CLASSICAL_METRIC_KEYS,
    TIMEOUT_DIAGNOSTIC_KEYS, PER_NODE_METRIC_KEYS,
    ALIGNED_DURATION_METRIC_KEYS, NODE_ROW_COLUMNS,
)
from lab.runs import Experiment, RUNS

logger = logging.getLogger(__name__)


def run_disco_degrade(log_paths, combos=ALL_COMBOS, degradations=ALL_DEGRADATIONS,
                     levels=ALL_LEVELS, metrics=ALL_METRICS,
                     out_csv='var/lab/results/exp_disco_degrade.csv',
                     node_out_csv=None, timings_out_csv=None):
    """
    Every (log, combo, dim, level) cell, scoring `metrics` (default
    ALL_METRICS) - see lab.run.run, which this delegates to unchanged.
    `metrics` is the inclusive roster to score; main()'s own
    --exclude-metric computes this as ALL_METRICS minus the excluded ids.
    """
    return lab_run.run(log_paths, combos=combos, degradations=degradations, levels=levels,
                       metrics=metrics, out_csv=out_csv, node_out_csv=node_out_csv,
                       timings_out_csv=timings_out_csv)


def _timestamped(path):
    """
    Insert a YYYYMMDD-HHMMSS stamp before the extension, eg
    'var/lab/results/rtfm.csv' -> 'var/lab/results/rtfm_20260908-121549.csv'.

    Applied unconditionally in main() below, whether --out was given or
    not: run_disco_degrade has no _merge_write - it's a plain df.to_csv
    overwrite every call, and every cell here is expensive (a real log,
    a real discovered tree, the full classical-alignment pass) - a
    same-name rerun silently clobbering the last one is exactly the
    mistake this exists to make structurally impossible.
    """
    path = Path(path)
    stamp = datetime.now().strftime('%Y%m%d-%H%M%S')
    return str(path.with_name(f'{path.stem}_{stamp}{path.suffix}'))


def _resolve_experiment(args, parser):
    """
    The Experiment this invocation resolves to, from either --run or ad
    hoc logs/--combos/--levels - the single place both --dry-run and
    the real run derive their configuration from, so the two can never
    show/do different things.
    """
    if args.run:
        if args.run not in RUNS:
            parser.error(f'Unknown run {args.run!r}; choices: {sorted(RUNS)}')
        out_csv = args.out or f'var/lab/results/{args.run}.csv'
        return replace(RUNS[args.run], out_csv=out_csv)

    if not args.logs:
        parser.error('logs are required unless --run is given')

    # NAMED_COMBOS/NAMED_DEGRADATIONS (real-discovery + claims-fixture
    # registrations merged) are for ad hoc --combos/--degradations
    # LOOKUP only - ALL_COMBOS/ALL_DEGRADATIONS (real-discovery only)
    # stay the actual default when neither flag is given, so a bare
    # `--run full`/no-flags invocation never silently picks up
    # claims_known (which ignores its log argument entirely and always
    # returns the same known tree, nonsensical against a real log).
    named_combos = {**ALL_COMBOS, **CLAIMS_COMBOS}
    named_degradations = {**ALL_DEGRADATIONS, **CLAIMS_DEGRADATIONS}

    if args.combos:
        for name in args.combos:
            if name not in named_combos:
                parser.error(f'Unknown combo {name!r}; choices: {sorted(named_combos)}')
        combos = {name: named_combos[name] for name in args.combos}
    else:
        combos = ALL_COMBOS

    if args.degradations:
        for name in args.degradations:
            if name not in named_degradations:
                parser.error(f'Unknown degradation {name!r}; '
                              f'choices: {sorted(named_degradations)}')
        degradations = {name: named_degradations[name] for name in args.degradations}
    else:
        degradations = ALL_DEGRADATIONS

    levels = args.levels or ALL_LEVELS
    out_csv = args.out or 'var/lab/results/exp_disco_degrade.csv'
    return Experiment(name='ad hoc', log_paths=args.logs, combos=combos,
                       degradations=degradations, levels=levels, out_csv=out_csv)


def main():
    configure()
    logger.info('Starting exp_disco_degrade')
    parser = argparse.ArgumentParser(
        description='Dose-response coverage experiment: fixed discovered '
                     'model vs degraded log.')
    parser.add_argument('logs', nargs='*', help='XES log path(s) (ignored if --run given)')
    parser.add_argument('--run', help='Named run from lab.runs.RUNS to execute '
                                       '(overrides logs/--levels/--combos/--degradations)')
    parser.add_argument('--combos', nargs='+',
                         help='Discovery combos to use (default: all of '
                              f'{list(ALL_COMBOS)}) - also accepts claims-fixture '
                              f'combos ({list(CLAIMS_COMBOS)}), validated in '
                              '_resolve_experiment rather than here so both '
                              'registries can be checked together.')
    parser.add_argument('--degradations', nargs='+',
                         help='Degradation dimensions to use (default: all of '
                              f'{list(ALL_DEGRADATIONS)}) - also accepts claims-fixture '
                              f'ablation targets ({list(CLAIMS_DEGRADATIONS)}).')
    parser.add_argument('--levels', type=float, nargs='+',
                         help=f'Degradation levels in [0,1] (default: {ALL_LEVELS})')
    parser.add_argument('--exclude-metric', nargs='+', default=[], metavar='METRIC_ID',
                         help='Skip scoring the given metric id(s) entirely for this run '
                              f'(choices: {[m.id for m in ALL_METRICS]}) - for a metric '
                              "that's known-expensive at a log's scale (eg voidsat's "
                              'per-trace cost) or known-wrong pending an upstream fix, '
                              'where re-collecting it later is cheap. An excluded metric '
                              'gets no timing row and its CSV column is entirely empty '
                              'for this run.')
    parser.add_argument('--out', help='Output CSV path base (default depends on --run, '
                                       'else var/lab/results/exp_disco_degrade.csv) - a '
                                       'timestamp is always inserted before the extension, '
                                       'so this is never overwritten by a later run.')
    parser.add_argument('--verbose', action='store_true',
                         help='Enable skip-alignments debug logging '
                              '(waste ratios, per-variant timing)')
    parser.add_argument('--dry-run', action='store_true',
                         help='Print the resolved experiment configuration '
                              '(logs/combos/degradations/levels/cell count) and '
                              'exit without computing anything.')
    args = parser.parse_args()

    if args.verbose:
        enable_skipalignments_debug()

    all_metric_ids = {m.id for m in ALL_METRICS}
    unknown = [m for m in args.exclude_metric if m not in all_metric_ids]
    if unknown:
        parser.error(f'Unknown --exclude-metric {unknown}; choices: {sorted(all_metric_ids)}')
    metrics = [m for m in ALL_METRICS if m.id not in args.exclude_metric]

    experiment = _resolve_experiment(args, parser)

    if args.dry_run:
        print(experiment.describe())
        if args.exclude_metric:
            print(f'  excluding metrics: {args.exclude_metric}')
        return

    out_csv = _timestamped(experiment.out_csv)
    df, node_df, timings_df = run_disco_degrade(
        log_paths=experiment.log_paths, combos=experiment.combos,
        degradations=experiment.degradations, levels=experiment.levels,
        metrics=metrics, out_csv=out_csv)
    print(df)
    logger.info("Wrote %s (%d rows), %d node rows and %d timing rows",
                out_csv, len(df), len(node_df), len(timings_df))


if __name__ == '__main__':
    main()
