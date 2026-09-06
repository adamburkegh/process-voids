"""
Experiment: provided model vs (varying) log.

No degradation and no discovery sweep - a fixed process tree is checked
against one or more logs and coverage is captured for each. Run with no
arguments to try it on the built-in running example (lab/fixtures.py).

Usage:
    python -m lab.exp_premodel
    python -m lab.exp_premodel --model models/rtfm_extra.ptml logs/rtfm_fine_appeal.xes.gz
"""

import argparse
import logging
import time
from pathlib import Path

import pandas as pd
import pm4py_config as pm4py

from lab.fixtures import build_running_example_log, build_running_example_tree
from lab.logconfig import configure, enable_skipalignments_debug
from lab.metrics import compute_metrics
from process_voids.tree import from_pm4py

logger = logging.getLogger(__name__)


def load_tree(model_path):
    pt_pm4py = pm4py.read_ptml(model_path)
    return from_pm4py(pt_pm4py)


def run_premodel(cases, out_csv='var/lab/results/exp_premodel.csv'):
    """cases: list of (label, log, tree) tuples."""
    logger.info('Experiment: premodel | cases=%s', [label for label, _, _ in cases])

    rows = []
    for label, log, tree in cases:
        started = time.monotonic()
        logger.debug('%s - starting', label)

        row = {'log': label}
        slpn_path = f'var/lab/premodel_{label}.slpn'
        try:
            metrics = compute_metrics(log, tree, slpn_path)
            row['status'] = 'ok'
            row.update(metrics)
        except Exception as e:
            row.update(status=f'error: {e}', weight_coverage=None,
                        skipprob=None, duration_coverage=None, alignment_coverage=None)
        rows.append(row)

        elapsed = time.monotonic() - started
        if row['status'] == 'ok':
            logger.info('%s - done in %.1fs', label, elapsed)
        else:
            logger.warning('%s - done in %.1fs (status=%s)', label, elapsed, row['status'])

    df = pd.DataFrame(rows)
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df


def main():
    configure()
    logger.info('Starting exp_premodel')
    parser = argparse.ArgumentParser(
        description='Coverage experiment: provided model vs one or more logs.')
    parser.add_argument('logs', nargs='*', help='XES log path(s)')
    parser.add_argument('--model', help='PTML model to test the logs against '
                                         '(required if any logs are given)')
    parser.add_argument('--out', default='var/lab/results/exp_premodel.csv')
    parser.add_argument('--verbose', action='store_true',
                         help='Enable skip-alignments debug logging '
                              '(waste ratios, per-variant timing)')
    args = parser.parse_args()

    if args.verbose:
        enable_skipalignments_debug()

    if args.logs:
        if not args.model:
            parser.error('--model is required when log paths are given')
        tree = load_tree(args.model)
        cases = [(Path(p).stem, pm4py.read_xes(p), tree) for p in args.logs]
    else:
        cases = [('running_example', build_running_example_log(),
                   build_running_example_tree())]

    df = run_premodel(cases, out_csv=args.out)
    print(df)
    logger.info('Wrote %s', args.out)


if __name__ == '__main__':
    main()
