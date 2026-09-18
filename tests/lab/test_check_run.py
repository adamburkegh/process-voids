'''
lab.check_run answers the "did this run go the way I expected" question
this session kept re-answering with a fresh python -c snippet each
time: row count, status breakdown, and which logs/combos/degradation
dims a result CSV actually covers.
'''

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from lab.check_run import summarize, format_summary


def _write_csv(rows):
    path = Path(tempfile.mkdtemp()) / 'result.csv'
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


class SummarizeTest(unittest.TestCase):
    def test_counts_rows_and_status(self):
        csv_path = _write_csv([
            {'log': 'a', 'combo': 'x', 'degradation_dim': 'd', 'status': 'ok'},
            {'log': 'a', 'combo': 'x', 'degradation_dim': 'd', 'status': 'ok'},
            {'log': 'a', 'combo': 'y', 'degradation_dim': 'd', 'status': 'error: boom'},
        ])
        summary = summarize(csv_path)
        self.assertEqual(summary['rows'], 3)
        self.assertEqual(summary['status'], {'ok': 2, 'error: boom': 1})

    def test_lists_distinct_logs_combos_and_dims(self):
        csv_path = _write_csv([
            {'log': 'a', 'combo': 'x', 'degradation_dim': 'd1', 'status': 'ok'},
            {'log': 'b', 'combo': 'y', 'degradation_dim': 'd2', 'status': 'ok'},
        ])
        summary = summarize(csv_path)
        self.assertEqual(summary['logs'], ['a', 'b'])
        self.assertEqual(summary['combos'], ['x', 'y'])
        self.assertEqual(summary['degradation_dims'], ['d1', 'd2'])

    def test_missing_columns_are_skipped_not_errors(self):
        # a timings-shaped CSV has no 'combo'/'degradation_dim' columns
        csv_path = _write_csv([{'log': 'a', 'seconds': 1.2, 'status': 'ok'}])
        summary = summarize(csv_path)  # must not raise
        self.assertIsNone(summary['combos'])
        self.assertIsNone(summary['degradation_dims'])
        self.assertEqual(summary['logs'], ['a'])


class FormatSummaryTest(unittest.TestCase):
    def test_omits_none_fields(self):
        text = format_summary({
            'rows': 1, 'status': {'ok': 1}, 'logs': ['a'],
            'combos': None, 'degradation_dims': None,
        })
        self.assertIn('rows: 1', text)
        self.assertIn("status: {'ok': 1}", text)
        self.assertIn("logs: ['a']", text)
        self.assertNotIn('combos', text)
        self.assertNotIn('degradation_dims', text)


if __name__ == '__main__':
    unittest.main()
