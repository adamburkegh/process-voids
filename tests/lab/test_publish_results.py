'''
Tests for lab.publish_results.

The published results/ folder must be reproducible from full-size result
CSVs by script, so these tests pin what the script does to them: the root
CSV passes through verbatim apart from the one column that names a
gitignored path, the timings collapse to one row per (combo, dimension,
level, stage) without losing the level-0.0 rows whose dimension is blank,
and a later source overrides an earlier one cell by cell.
'''

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from lab.publish_results import (
    aggregate_timings, combine_status, publish, timings_path)
from lab.runtime_table import total_seconds


def _root_row(combo, dim, level, status='ok', **extra):
    return {'log': 'BPIC2020_rfp', 'combo': combo, 'degradation_dim': dim,
            'degradation_level': level, 'status': status,
            'tree_cache_file': 'var/lab/tree_cache/x.pkl', 'voidsat2': 0.5, **extra}


def _timing_row(combo, dim, level, stage, seconds, status='ok'):
    return {'log': 'BPIC2020_rfp', 'combo': combo, 'degradation_dim': dim,
            'degradation_level': level, 'metric_or_stage': stage,
            'seconds': seconds, 'status': status}


def _write(dir_path, name, root_rows, timing_rows):
    root = Path(dir_path) / name
    pd.DataFrame(root_rows).to_csv(root, index=False)
    pd.DataFrame(timing_rows).to_csv(timings_path(root), index=False)
    return root


class TimingsPathTest(unittest.TestCase):
    def test_inserts_timings_before_extension(self):
        self.assertEqual(timings_path('a/b/run_2026.csv'), Path('a/b/run_2026_timings.csv'))


class CombineStatusTest(unittest.TestCase):
    def test_all_ok_is_ok(self):
        self.assertEqual(combine_status(pd.Series(['ok', 'ok'])), 'ok')

    def test_any_failure_is_reported_not_hidden(self):
        self.assertEqual(combine_status(pd.Series(['ok', 'error: b', 'error: a'])),
                         'error: a;error: b')


class AggregateTimingsTest(unittest.TestCase):
    def test_sums_seconds_per_cell_and_stage(self):
        df = pd.DataFrame([
            _timing_row('c', 'trace', 0.1, 'skipprob', 1.0),
            _timing_row('c', 'trace', 0.1, 'skipprob', 2.5),
            _timing_row('c', 'trace', 0.1, 'classical', 10.0),
        ])
        out = aggregate_timings(df).set_index('metric_or_stage')
        self.assertAlmostEqual(out.loc['skipprob', 'seconds'], 3.5)
        self.assertEqual(out.loc['skipprob', 'rows'], 2)
        self.assertAlmostEqual(out.loc['classical', 'seconds'], 10.0)

    def test_keeps_level_zero_rows_whose_dimension_is_blank(self):
        '''Level 0.0 is computed once and shared across dimensions, so its
        timing rows carry no degradation_dim. A groupby that drops NaN
        keys would silently delete exactly those rows.'''
        df = pd.DataFrame([
            _timing_row('c', float('nan'), 0.0, 'classical', 100.0),
            _timing_row('c', 'trace', 0.1, 'classical', 10.0),
        ])
        out = aggregate_timings(df)
        self.assertEqual(len(out), 2)
        self.assertAlmostEqual(out['seconds'].sum(), 110.0)

    def test_failed_stage_status_survives_aggregation(self):
        df = pd.DataFrame([
            _timing_row('c', 'trace', 0.1, 'dv', 1.0),
            _timing_row('c', 'trace', 0.1, 'dv', 1.0, status='error: boom'),
        ])
        self.assertEqual(aggregate_timings(df)['status'].iloc[0], 'error: boom')


class PublishTest(unittest.TestCase):
    def test_root_csv_verbatim_apart_from_dropped_column(self):
        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as out:
            root = _write(src, 'run_1.csv',
                          [_root_row('inductive_noise20', 'trace', 0.0),
                           _root_row('inductive_noise20', 'trace', 0.1)],
                          [_timing_row('inductive_noise20', 'trace', 0.1, 'dv', 1.0)])
            publish({'bpic2020_rfp': [root]}, out)
            published = pd.read_csv(Path(out) / 'bpic2020_rfp.csv')
            self.assertNotIn('tree_cache_file', published.columns)
            # the log column stays as the run wrote it, not the registry id
            self.assertEqual(set(published['log']), {'BPIC2020_rfp'})
            self.assertEqual(list(published['degradation_level']), [0.0, 0.1])
            self.assertIn('voidsat2', published.columns)

    def test_published_timings_feed_runtime_table_unchanged(self):
        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as out:
            raw = [_timing_row('c', 'trace', 0.1, 'classical', 4.0),
                   _timing_row('c', 'trace', 0.1, 'voidmass_process_lower', 0.5),
                   _timing_row('c', 'trace', 0.2, 'classical', 6.0),
                   _timing_row('c', 'trace', 0.2, 'voidmass_process_lower', 0.5)]
            root = _write(src, 'run_1.csv', [_root_row('c', 'trace', 0.1)], raw)
            publish({'rtfm': [root]}, out)
            published = pd.read_csv(Path(out) / 'rtfm_timings.csv')
            self.assertAlmostEqual(
                total_seconds(published, 'c', 'voidmass_process'),
                total_seconds(pd.DataFrame(raw), 'c', 'voidmass_process'))

    def test_later_source_overrides_earlier_cell_by_cell(self):
        '''A retry run replaces just the cells it recomputed - root row and
        timings both - and leaves every other cell of the first run alone.'''
        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as out:
            first = _write(
                src, 'first.csv',
                [_root_row('c', 'trace', 0.1), _root_row('c', 'trace', 0.2, status='error: boom')],
                [_timing_row('c', 'trace', 0.1, 'dv', 1.0),
                 _timing_row('c', 'trace', 0.2, 'dv', 9.0, status='error: boom')])
            retry = _write(
                src, 'retry.csv',
                [_root_row('c', 'trace', 0.2)],
                [_timing_row('c', 'trace', 0.2, 'dv', 2.0)])
            publish({'rtfm': [first, retry]}, out)
            root = pd.read_csv(Path(out) / 'rtfm.csv').sort_values('degradation_level')
            self.assertEqual(list(root['status']), ['ok', 'ok'])
            timings = pd.read_csv(Path(out) / 'rtfm_timings.csv').set_index('degradation_level')
            self.assertAlmostEqual(timings.loc[0.1, 'seconds'], 1.0)
            self.assertAlmostEqual(timings.loc[0.2, 'seconds'], 2.0)
            self.assertEqual(timings.loc[0.2, 'status'], 'ok')

    def test_runs_csv_joins_run_history_on_file_name(self):
        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as out:
            root = _write(src, 'run_1.csv', [_root_row('c', 'trace', 0.1)],
                          [_timing_row('c', 'trace', 0.1, 'dv', 1.0)])
            history = Path(src) / 'run_history.csv'
            pd.DataFrame([
                {'run_timestamp': '2026-09-01T00:00:00', 'out_csv': 'var/lab/results/run_1.csv',
                 'combos': 'c', 'degradations': 'trace', 'levels': '0.1', 'metrics': 'voidsat2',
                 'pythonhashseed': 'unset', 'process_voids_version': 'old',
                 'skipalignments_version': '0.3.0'},
                # a rerun writing the same out_csv is appended, not upserted
                {'run_timestamp': '2026-09-02T00:00:00', 'out_csv': 'var/lab/results/run_1.csv',
                 'combos': 'c', 'degradations': 'trace', 'levels': '0.1', 'metrics': 'voidsat2',
                 'pythonhashseed': 'unset', 'process_voids_version': 'new',
                 'skipalignments_version': '0.3.0'},
            ]).to_csv(history, index=False)
            publish({'rtfm': [root]}, out, run_history_csv=history)
            runs = pd.read_csv(Path(out) / 'runs.csv')
            self.assertEqual(list(runs['log_id']), ['rtfm'])
            self.assertEqual(list(runs['source_csv']), ['run_1.csv'])
            self.assertEqual(list(runs['process_voids_version']), ['new'])

    def test_runs_csv_tolerates_a_run_older_than_run_history(self):
        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as out:
            root = _write(src, 'run_1.csv', [_root_row('c', 'trace', 0.1)],
                          [_timing_row('c', 'trace', 0.1, 'dv', 1.0)])
            publish({'rtfm': [root]}, out, run_history_csv=Path(src) / 'absent.csv')
            runs = pd.read_csv(Path(out) / 'runs.csv')
            self.assertEqual(list(runs['source_csv']), ['run_1.csv'])
            self.assertTrue(pd.isna(runs['process_voids_version'].iloc[0]))

    def test_runs_csv_has_no_log_column_so_coverage_scans_skip_it(self):
        '''lab.current_coverage reads every root-level CSV that has a 'log'
        column; runs.csv is provenance, not results.'''
        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as out:
            root = _write(src, 'run_1.csv', [_root_row('c', 'trace', 0.1)],
                          [_timing_row('c', 'trace', 0.1, 'dv', 1.0)])
            publish({'rtfm': [root]}, out)
            self.assertNotIn('log', pd.read_csv(Path(out) / 'runs.csv').columns)

    def test_returns_every_path_written(self):
        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as out:
            root = _write(src, 'run_1.csv', [_root_row('c', 'trace', 0.1)],
                          [_timing_row('c', 'trace', 0.1, 'dv', 1.0)])
            written = publish({'rtfm': [root]}, out)
            self.assertEqual(sorted(p.name for p in written),
                             ['rtfm.csv', 'rtfm_timings.csv', 'runs.csv'])

    def test_missing_timings_companion_is_an_error(self):
        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as out:
            root = Path(src) / 'run_1.csv'
            pd.DataFrame([_root_row('c', 'trace', 0.1)]).to_csv(root, index=False)
            with self.assertRaises(FileNotFoundError):
                publish({'rtfm': [root]}, out)


if __name__ == '__main__':
    unittest.main()
