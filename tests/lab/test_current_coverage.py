'''
Tests for lab.current_coverage.

lab.collection_report answers "has this metric ever been collected for
this log, across every result CSV ever written" - deliberately including
retired ids and stale runs. This module answers a different, narrower
question: "for this log's MOST RECENT sweep, which of today's LIVE
metrics does it actually have" - the "are we current" question that
kept getting reconstructed by hand this session (checking a CSV's
header against what the registry currently expects).
'''

import os
import tempfile
import time
import unittest
from pathlib import Path

import pandas as pd

from lab.current_coverage import current_coverage, latest_csv_per_log


def _write_csv(dir_path, name, log, columns, mtime_offset=0):
    path = Path(dir_path) / name
    row = {'log': log, **{c: 0.5 for c in columns}}
    pd.DataFrame([row]).to_csv(path, index=False)
    if mtime_offset:
        stat = path.stat()
        os.utime(path, (stat.st_atime, stat.st_mtime + mtime_offset))
    return path


class LatestCsvPerLogTest(unittest.TestCase):
    def test_picks_the_most_recently_modified_file(self):
        with tempfile.TemporaryDirectory() as d:
            older = _write_csv(d, 'a.csv', 'rtfm', ['skipprob'])
            time.sleep(0.01)
            newer = _write_csv(d, 'b.csv', 'rtfm', ['voidsat2'])
            latest = latest_csv_per_log(d)
            self.assertEqual(latest['rtfm'], newer)

    def test_separate_logs_tracked_independently(self):
        with tempfile.TemporaryDirectory() as d:
            a = _write_csv(d, 'a.csv', 'rtfm', ['skipprob'])
            b = _write_csv(d, 'b.csv', 'bpic2020_rfp', ['skipprob'])
            latest = latest_csv_per_log(d)
            self.assertEqual(latest['rtfm'], a)
            self.assertEqual(latest['bpic2020_rfp'], b)

    def test_nodes_and_timings_companions_ignored(self):
        with tempfile.TemporaryDirectory() as d:
            _write_csv(d, 'a_nodes.csv', 'rtfm', ['skipprob'])
            _write_csv(d, 'a_timings.csv', 'rtfm', ['skipprob'])
            latest = latest_csv_per_log(d)
            self.assertNotIn('rtfm', latest)


class CurrentCoverageTest(unittest.TestCase):
    def test_only_live_metrics_are_rows(self):
        with tempfile.TemporaryDirectory() as d:
            _write_csv(d, 'a.csv', 'rtfm', ['skipprob', 'weight_coverage'])
            df = current_coverage(d)
            # weight_coverage is retired (superseded by weight_voidage) -
            # present in the CSV, but must not appear as a row here.
            self.assertIn('skipprob', df.index)
            self.assertNotIn('weight_coverage', df.index)

    def test_missing_column_reads_false(self):
        with tempfile.TemporaryDirectory() as d:
            _write_csv(d, 'a.csv', 'rtfm', ['skipprob'])
            df = current_coverage(d)
            self.assertFalse(df.loc['voidsat2', 'rtfm'])
            self.assertTrue(df.loc['skipprob', 'rtfm'])

    def test_stale_file_does_not_leak_into_current_view(self):
        '''An older file's metric that the newer file doesn't have must
        read False - unlike lab.collection_report, which would OR them
        together across every file ever written.'''
        with tempfile.TemporaryDirectory() as d:
            _write_csv(d, 'old.csv', 'rtfm', ['voidsalign3'])
            time.sleep(0.01)
            _write_csv(d, 'new.csv', 'rtfm', ['voidsat2'])
            df = current_coverage(d)
            self.assertFalse(df.loc['voidsalign3', 'rtfm'])
            self.assertTrue(df.loc['voidsat2', 'rtfm'])


class MainTest(unittest.TestCase):
    def test_results_dir_flag_reads_that_directory(self):
        import contextlib
        import io
        from unittest.mock import patch
        from lab.current_coverage import main
        with tempfile.TemporaryDirectory() as d:
            _write_csv(d, 'rtfm.csv', 'rtfm', ['skipprob'])
            # provenance file beside the results: no 'log' column, must be skipped
            pd.DataFrame([{'log_id': 'rtfm', 'source_csv': 'x.csv'}]).to_csv(
                Path(d) / 'runs.csv', index=False)
            out = io.StringIO()
            with patch('sys.argv', ['current_coverage', '--results-dir', d]), \
                    contextlib.redirect_stdout(out):
                main()
        self.assertIn('| skipprob | Done |', out.getvalue())
        self.assertIn('as of:', out.getvalue())


if __name__ == '__main__':
    unittest.main()
