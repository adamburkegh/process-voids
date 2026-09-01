"""
Experiment: discovered model vs (degraded) log.

For each log, a model is discovered once per combo from the original,
undegraded log - that discovered model is the fixed reference for that
(log, combo). Each degradation dimension x level then degrades the log
and checks it against that same fixed model, producing a dose-response
curve in the coverage metrics as degradation increases.

run_experiment1() takes explicit logs/combos/degradations/levels and
runs their Cartesian product - the full parameter catalog lives in
lab/params.py, and specific reproducible subset selections live in
lab/runs.py (RUNS), runnable by name:

    python -m lab.exp_disco_degrade --run smoke
    python -m lab.exp_disco_degrade logs/rtfm_fine_appeal.xes.gz --levels 0.0 0.5
"""

import argparse
from pathlib import Path

import pandas as pd
import pm4py

from lab.metrics import compute_metrics
from lab.params import ALL_COMBOS, ALL_DEGRADATIONS, ALL_LEVELS


def run_experiment1(log_paths, combos=ALL_COMBOS, degradations=ALL_DEGRADATIONS,
                     levels=ALL_LEVELS, out_csv='var/lab/results/exp_disco_degrade.csv'):
    rows = []
    for log_path in log_paths:
        log_name = Path(log_path).stem
        base_log = pm4py.read_xes(log_path)

        for combo_name, combo in combos.items():
            try:
                tree = combo.discover(base_log)
                discover_status = 'ok'
            except NotImplementedError:
                tree = None
                discover_status = 'not_implemented'
            except Exception as e:
                tree = None
                discover_status = f'discovery error: {e}'

            for dim, degrade in degradations.items():
                for level in levels:
                    degraded_log, dropped = degrade(base_log, level)
                    row = {
                        'log': log_name,
                        'combo': combo_name,
                        'degradation_dim': dim,
                        'degradation_level': level,
                        'dropped': ', '.join(sorted(str(d) for d in dropped)),
                    }
                    if tree is None:
                        row.update(status=discover_status, weight_coverage=None,
                                   skipprob=None, duration_coverage=None, alignment_coverage=None)
                        rows.append(row)
                        continue
                    slpn_path = (f'var/lab/e1_{log_name}_{combo_name}_'
                                 f'{dim}_{level}.slpn')
                    try:
                        metrics = compute_metrics(
                            degraded_log, tree, slpn_path,
                            has_estimator=combo.has_estimator)
                        row['status'] = 'ok'
                        row.update(metrics)
                    except Exception as e:
                        row.update(status=f'error: {e}', weight_coverage=None,
                                    skipprob=None, duration_coverage=None, alignment_coverage=None)
                    rows.append(row)

    df = pd.DataFrame(rows)
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df


def main():
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
    args = parser.parse_args()

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

    df = run_experiment1(**kwargs)
    print(df)
    print(f"Wrote {kwargs.get('out_csv', 'var/lab/results/exp_disco_degrade.csv')}")


if __name__ == '__main__':
    main()
