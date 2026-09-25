'''
Builds a compact, publishable folder from full-size experiment result
CSVs, so anyone who reruns an experiment can regenerate the same folder.

A publish is a snapshot, not a rolling update: point --out-dir at a dated
subfolder (results/<YYYYMMDD>/csv, alongside a plots/ folder for
lab.cross_log_plots' output and a README.md - see results/20260919 for
the pattern) so a later publish adds a new folder rather than silently
overwriting an earlier one a paper or report already cites.

Per log id it writes:

  <log_id>.csv          the root-level result CSV, verbatim (including the
                        `log` column, which holds the XES file stem as the
                        run wrote it) except for columns that name files
                        outside source control (DROP_COLUMNS). Rows are
                        sorted by cell, so the file does not depend on how
                        many runs it was assembled from.
  <log_id>_timings.csv  the run's <stem>_timings.csv collapsed to one row
                        per (combo, degradation_dim, degradation_level,
                        stage): `seconds` summed, `rows` counting what was
                        summed, `status` 'ok' unless any summed row failed.
                        lab.runtime_table reads it unchanged.

and one runs.csv for the lot: which source file each log came from, with
the versions, seed and configuration lab.run_history recorded for it.

Naming a log id more than once adds a later source. A later source
replaces the earlier ones cell by cell - root row and timings both - which
is how a retry run that recomputed a few failed cells is folded into the
run it retries without editing either file.

The per-node CSVs are not published: the cross-log plots, tables and
coverage scans read only the root-level files.

Usage:
    python -m lab.publish_results \\
        --source rtfm=var/lab/results/rtfm_20260918-102305.csv \\
        --source rtfm=var/lab/results/rtfm_retry_20260918-124136.csv \\
        --source bpic2020_rfp=var/lab/results/bpic2020_rfp_20260917-091018.csv \\
        --out-dir results/20260919/csv
'''

import argparse
from pathlib import Path

import pandas as pd

RESULTS_DIR = 'results'

# Names a gitignored cache file - a reference to something outside source
# control, in a folder meant to be public.
DROP_COLUMNS = ('tree_cache_file',)

CELL_SORT = ['combo', 'degradation_dim', 'degradation_level']
TIMING_GROUP = ['log', 'combo', 'degradation_dim', 'degradation_level', 'metric_or_stage']

# From lab.run_history's row for the run; the join columns are not repeated.
PROVENANCE_COLUMNS = ('run_timestamp', 'combos', 'degradations', 'levels', 'metrics',
                      'pythonhashseed', 'process_voids_version', 'skipalignments_version')


def timings_path(root_csv):
    '''The *_timings.csv companion lab.run writes beside a root-level CSV.'''
    path = Path(root_csv)
    return path.with_name(f'{path.stem}_timings{path.suffix}')


def combine_status(statuses):
    ''''ok' if every status is, else the distinct failures, sorted and
    ';'-joined - a failed stage stays visible after aggregation.'''
    failures = sorted({status for status in statuses if status != 'ok'})
    return ';'.join(failures) if failures else 'ok'


def aggregate_timings(timings):
    '''One row per (log, combo, degradation_dim, degradation_level,
    metric_or_stage). dropna=False keeps the level-0.0 rows: that level is
    computed once and shared across dimensions, so their degradation_dim is
    blank, and grouping on a NaN key would silently discard them.'''
    grouped = timings.groupby(TIMING_GROUP, dropna=False)
    return grouped.agg(seconds=('seconds', 'sum'),
                       rows=('seconds', 'size'),
                       status=('status', combine_status)).reset_index()


def _cell_keys(df):
    '''One string per row naming its (combo, degradation_dim, level) cell.
    A blank dimension keys as '' - NaN never compares equal to NaN.'''
    dim = df['degradation_dim'].fillna('').astype(str)
    return df['combo'].astype(str) + '|' + dim + '|' + df['degradation_level'].map(repr)


def _later_wins(frames):
    '''Concatenate, each later frame replacing the cells it holds.'''
    merged = frames[0]
    for later in frames[1:]:
        replaced = set(_cell_keys(later))
        merged = pd.concat([merged[~_cell_keys(merged).isin(replaced)], later],
                           ignore_index=True)
    return merged


def _read_history(path):
    return pd.read_csv(path) if path is not None and Path(path).exists() else None


def _provenance(log_id, source, history):
    row = {'log_id': log_id, 'source_csv': Path(source).name}
    row.update({column: None for column in PROVENANCE_COLUMNS})
    if history is not None:
        # An out_csv rerun is appended to run_history, not upserted: the
        # last row is the run that wrote the file now on disk.
        names = history['out_csv'].map(lambda path: Path(str(path)).name)
        matches = history[names == row['source_csv']]
        if not matches.empty:
            last = matches.iloc[-1]
            row.update({column: last.get(column) for column in PROVENANCE_COLUMNS})
    return row


def publish(sources, out_dir=RESULTS_DIR, run_history_csv=None):
    '''
    Write out_dir's files for {log_id: [root_csv, ...]} and return every
    path written. A missing *_timings.csv companion raises: the durations
    are part of what is published.

    run_history_csv defaults to run_history.csv beside the first source.
    '''
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if run_history_csv is None:
        first = next(iter(sources.values()))[0]
        run_history_csv = Path(first).parent / 'run_history.csv'
    history = _read_history(run_history_csv)

    written = []
    provenance = []
    for log_id, root_csvs in sources.items():
        roots = _later_wins([pd.read_csv(path) for path in root_csvs])
        roots = roots.drop(columns=[c for c in DROP_COLUMNS if c in roots.columns])
        roots = roots.sort_values(CELL_SORT, kind='stable')
        root_out = out_dir / f'{log_id}.csv'
        roots.to_csv(root_out, index=False)

        timings = _later_wins([pd.read_csv(timings_path(path)) for path in root_csvs])
        timings_out = out_dir / f'{log_id}_timings.csv'
        aggregate_timings(timings).to_csv(timings_out, index=False)

        written += [root_out, timings_out]
        provenance += [_provenance(log_id, path, history) for path in root_csvs]

    runs_out = out_dir / 'runs.csv'
    pd.DataFrame(provenance).to_csv(runs_out, index=False)
    return written + [runs_out]


def main():
    parser = argparse.ArgumentParser(
        description='Build the public results/ folder from full-size result CSVs.')
    parser.add_argument('--source', action='append', required=True, dest='sources',
                         metavar='LOG_ID=PATH',
                         help='a root-level result CSV for a log id; repeat the id to add a '
                              'later source that replaces earlier ones cell by cell')
    parser.add_argument('--out-dir', default=RESULTS_DIR)
    parser.add_argument('--run-history', default=None,
                         help="lab.run_history's CSV (default: run_history.csv beside the "
                              'first source)')
    args = parser.parse_args()

    sources = {}
    for item in args.sources:
        log_id, path = item.split('=', 1)
        sources.setdefault(log_id, []).append(path)
    for path in publish(sources, args.out_dir, args.run_history):
        print(f'Wrote {path}')


if __name__ == '__main__':
    main()
