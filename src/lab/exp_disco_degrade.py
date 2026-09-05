"""
Experiment: discovered model vs (degraded) log.

For each log, a model is discovered once per combo from the original,
undegraded log - that discovered model is the fixed reference for that
(log, combo). Each degradation dimension x level then degrades the log
and checks it against that same fixed model, producing a dose-response
curve in the coverage metrics as degradation increases.

run_disco_degrade() takes explicit logs/combos/degradations/levels and
runs their Cartesian product - the full parameter catalog lives in
lab/params.py, and specific reproducible subset selections live in
lab/runs.py (RUNS), runnable by name:

    python -m lab.exp_disco_degrade --run smoke
    python -m lab.exp_disco_degrade logs/rtfm_fine_appeal.xes.gz --levels 0.0 0.5
"""

import argparse
import logging
import time
from pathlib import Path

import pandas as pd
import pm4py_config as pm4py

from lab.logconfig import configure, enable_skipalignments_debug
from lab.metrics import compute_metrics
from lab.params import ALL_COMBOS, ALL_DEGRADATIONS, ALL_LEVELS

logger = logging.getLogger(__name__)


def _format_dropped(dropped, limit=50):
    """
    (dropped_str, dropped_count). Dropped can be huge for trace-wise
    degradation at high levels (one entry per dropped case) - listing
    all of them bloats the CSV by orders of magnitude for no benefit,
    so beyond `limit` items we just report the count.
    """
    items = sorted(str(d) for d in dropped)
    if len(items) > limit:
        return f'{len(items)} items', len(items)
    return ', '.join(items), len(items)


def _log_stats(log):
    """(n_cases, n_variants), or (None, None) if log isn't a DataFrame."""
    if not hasattr(log, 'groupby'):
        return None, None
    n_cases = log['case:concept:name'].nunique()
    n_variants = log.groupby('case:concept:name')['concept:name'] \
                     .apply(tuple).nunique()
    return n_cases, n_variants


def run_disco_degrade(log_paths, combos=ALL_COMBOS, degradations=ALL_DEGRADATIONS,
                     levels=ALL_LEVELS, out_csv='var/lab/results/exp_disco_degrade.csv'):
    logger.info('Experiment: disco_degrade | logs=%s | combos=%s | '
                'degradations=%s | levels=%s',
                [Path(p).stem for p in log_paths], list(combos),
                list(degradations), levels)

    rows = []
    for log_path in log_paths:
        log_name = Path(log_path).stem
        base_log = pm4py.read_xes(log_path)
        n_cases, n_variants = _log_stats(base_log)
        logger.info('Log %s: %s cases, %s variants', log_name, n_cases, n_variants)

        for combo_name, combo in combos.items():
            started_discover = time.monotonic()
            try:
                result = combo.discover(base_log)
                tree, ppt_weights = result.tree, result.ppt_weights
                discover_status = 'ok'
            except NotImplementedError:
                tree, ppt_weights = None, None
                discover_status = 'not_implemented'
            except Exception as e:
                tree, ppt_weights = None, None
                discover_status = f'discovery error: {e}'
            logger.info('Discovery: log=%s combo=%s -> %s (%.1fs)',
                        log_name, combo_name, discover_status,
                        time.monotonic() - started_discover)

            # Level 0.0 applies zero drops regardless of dimension, so it's
            # the same (log, tree) computation under every dim - compute it
            # once here and reuse the result across dims below, rather than
            # redoing an identical (and potentially very expensive) run per
            # dimension.
            zero_level_status = None
            zero_level_metrics = None
            if tree is not None and 0.0 in levels:
                cell = f'{log_name} / {combo_name} / (all dims) / 0.0'
                started = time.monotonic()
                logger.debug('%s - starting', cell)
                slpn_path = f'var/lab/disco_degrade_{log_name}_{combo_name}_level0.slpn'
                try:
                    zero_level_metrics = compute_metrics(
                        base_log, tree, slpn_path, ppt_weights=ppt_weights)
                    zero_level_status = 'ok'
                except Exception as e:
                    zero_level_status = f'error: {e}'
                elapsed = time.monotonic() - started
                if zero_level_status == 'ok':
                    logger.info('%s - done in %.1fs', cell, elapsed)
                else:
                    logger.warning('%s - done in %.1fs (status=%s)',
                                    cell, elapsed, zero_level_status)

            for dim, degrade in degradations.items():
                for level in levels:
                    cell = f'{log_name} / {combo_name} / {dim} / {level}'

                    if tree is None:
                        row = {
                            'log': log_name, 'combo': combo_name,
                            'degradation_dim': dim, 'degradation_level': level,
                            'dropped': '', 'dropped_count': 0, 'status': discover_status,
                            'weight_coverage': None, 'skipprob': None,
                            'duration_coverage': None, 'alignment_coverage': None,
                        }
                        rows.append(row)
                        logger.info('%s - skipped (%s)', cell, discover_status)
                        continue

                    if level == 0.0:
                        row = {
                            'log': log_name, 'combo': combo_name,
                            'degradation_dim': dim, 'degradation_level': level,
                            'dropped': '', 'dropped_count': 0, 'status': zero_level_status,
                        }
                        if zero_level_status == 'ok':
                            row.update(zero_level_metrics)
                        else:
                            row.update(weight_coverage=None, skipprob=None,
                                       duration_coverage=None, alignment_coverage=None)
                        rows.append(row)
                        logger.debug('%s - reused level-0.0 result', cell)
                        continue

                    started = time.monotonic()
                    logger.debug('%s - starting', cell)

                    degraded_log, dropped = degrade(base_log, level)
                    dropped_str, dropped_count = _format_dropped(dropped)
                    row = {
                        'log': log_name,
                        'combo': combo_name,
                        'degradation_dim': dim,
                        'degradation_level': level,
                        'dropped': dropped_str,
                        'dropped_count': dropped_count,
                    }
                    slpn_path = (f'var/lab/disco_degrade_{log_name}_{combo_name}_'
                                 f'{dim}_{level}.slpn')
                    try:
                        metrics = compute_metrics(
                            degraded_log, tree, slpn_path,
                            ppt_weights=ppt_weights)
                        row['status'] = 'ok'
                        row.update(metrics)
                    except Exception as e:
                        row.update(status=f'error: {e}', weight_coverage=None,
                                    skipprob=None, duration_coverage=None, alignment_coverage=None)
                    rows.append(row)

                    elapsed = time.monotonic() - started
                    if row['status'] == 'ok':
                        logger.info('%s - done in %.1fs', cell, elapsed)
                    else:
                        logger.warning('%s - done in %.1fs (status=%s)',
                                        cell, elapsed, row['status'])

    df = pd.DataFrame(rows)
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df


def main():
    configure()
    logger.info('Starting exp_disco_degrade')
    parser = argparse.ArgumentParser(
        description='Dose-response coverage experiment: fixed discovered '
                     'model vs degraded log.')
    parser.add_argument('logs', nargs='*', help='XES log path(s) (ignored if --run given)')
    parser.add_argument('--run', help='Named run from lab.runs.RUNS to execute '
                                       '(overrides logs/--levels)')
    parser.add_argument('--levels', type=float, nargs='+',
                         help=f'Degradation levels in [0,1] (default: {ALL_LEVELS})')
    parser.add_argument('--out', help='Output CSV path (default depends on --run, '
                                       'else var/lab/results/exp_disco_degrade.csv)')
    parser.add_argument('--verbose', action='store_true',
                         help='Enable skip-alignments debug logging '
                              '(waste ratios, per-variant timing)')
    args = parser.parse_args()

    if args.verbose:
        enable_skipalignments_debug()

    if args.run:
        from lab.runs import RUNS
        if args.run not in RUNS:
            parser.error(f'Unknown run {args.run!r}; choices: {sorted(RUNS)}')
        kwargs = dict(RUNS[args.run])
        kwargs['out_csv'] = args.out or f'var/lab/results/{args.run}.csv'
    else:
        if not args.logs:
            parser.error('logs are required unless --run is given')
        kwargs = dict(log_paths=args.logs)
        if args.levels:
            kwargs['levels'] = args.levels
        if args.out:
            kwargs['out_csv'] = args.out

    df = run_disco_degrade(**kwargs)
    print(df)
    logger.info("Wrote %s", kwargs.get('out_csv', 'var/lab/results/exp_disco_degrade.csv'))


if __name__ == '__main__':
    main()
