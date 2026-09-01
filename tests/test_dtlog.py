
import datetime
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import pm4py

from process_voids import dtlog


def rows(df):
    return list(zip(df['case:concept:name'], df['concept:name']))


class ConvertTest(unittest.TestCase):
    # shaped after pmkoalas' own dtlog tests (test_convert_empty,
    # test_convert_single_trace, ...), adapted for DataFrame output.

    def test_convert_empty(self):
        df = dtlog.convert()
        self.assertEqual(len(df), 0)

    def test_convert_single_trace(self):
        df = dtlog.convert('a')
        self.assertEqual(rows(df), [('case_1', 'a')])

    def test_convert_single_trace_multi_event(self):
        df = dtlog.convert('a b')
        self.assertEqual(rows(df), [('case_1', 'a'), ('case_1', 'b')])

    def test_convert_custom_delimiter(self):
        df = dtlog.convert('a-b', delimiter='-')
        self.assertEqual(rows(df), [('case_1', 'a'), ('case_1', 'b')])

    def test_convert_multiple_traces(self):
        df = dtlog.convert('a b', 'a', 'a b')
        self.assertEqual(rows(df), [
            ('case_1', 'a'), ('case_1', 'b'),
            ('case_2', 'a'),
            ('case_3', 'a'), ('case_3', 'b'),
        ])

    def test_convert_default_case_names(self):
        df = dtlog.convert('a', 'b', 'c')
        self.assertEqual(sorted(df['case:concept:name'].unique()),
                          ['case_1', 'case_2', 'case_3'])

    def test_convert_explicit_names(self):
        df = dtlog.convert('a b', 'a', names=['sigma1', 'sigma4'])
        self.assertEqual(rows(df), [
            ('sigma1', 'a'), ('sigma1', 'b'),
            ('sigma4', 'a'),
        ])

    def test_convert_names_length_mismatch(self):
        with self.assertRaises(ValueError):
            dtlog.convert('a', 'b', names=['sigma1'])

    def test_convert_events_strictly_increasing(self):
        df = dtlog.convert('a b', 'c d')
        self.assertTrue(df['time:timestamp'].is_monotonic_increasing)
        self.assertEqual(len(df['time:timestamp'].unique()), 4)


class ConvertTimedTest(unittest.TestCase):

    def test_basic_offsets(self):
        df = dtlog.convert_timed('o:0 a:1 s:2 p:3', names=['sigma1'])
        base = df.loc[df['concept:name'] == 'o', 'time:timestamp'].iloc[0]
        deltas = {
            row['concept:name']: row['time:timestamp'] - base
            for _, row in df.iterrows()
        }
        self.assertEqual(deltas['a'], datetime.timedelta(hours=1))
        self.assertEqual(deltas['s'], datetime.timedelta(hours=2))
        self.assertEqual(deltas['p'], datetime.timedelta(hours=3))

    def test_multiple_cases_paper_example(self):
        df = dtlog.convert_timed(
            'o:0 a:1 s:2 p:3',
            'o:0 a:1 e:5 a:9 s:10 p:11',
            'o:0 s:10 p:11',
            'o:0 s:1 p:2',
            'o:0 a:1 s:10 p:11',
            names=['sigma1', 'sigma2', 'sigma3', 'sigma4', 'sigma5'],
        )
        self.assertEqual(len(df), 4 + 6 + 3 + 3 + 4)
        self.assertEqual(
            list(df[df['case:concept:name'] == 'sigma3']['concept:name']),
            ['o', 's', 'p'])

    def test_time_units(self):
        hours = dtlog.convert_timed('a:0 b:1', time_unit='hours')
        minutes = dtlog.convert_timed('a:0 b:1', time_unit='minutes')
        seconds = dtlog.convert_timed('a:0 b:1', time_unit='seconds')
        for df, unit in [(hours, datetime.timedelta(hours=1)),
                          (minutes, datetime.timedelta(minutes=1)),
                          (seconds, datetime.timedelta(seconds=1))]:
            delta = df['time:timestamp'].iloc[1] - df['time:timestamp'].iloc[0]
            self.assertEqual(delta, unit)

    def test_invalid_time_unit(self):
        with self.assertRaises(ValueError):
            dtlog.convert_timed('a:0', time_unit='fortnights')

    def test_names_length_mismatch(self):
        with self.assertRaises(ValueError):
            dtlog.convert_timed('a:0', 'b:0', names=['only_one'])

    def test_sorted_by_case_then_time(self):
        df = dtlog.convert_timed('a:5 b:6', 'c:0 d:1', names=['later', 'earlier'])
        self.assertEqual(list(df['case:concept:name']), sorted(df['case:concept:name']))


class WriteXesTest(unittest.TestCase):
    # regression test for the pandas datetime64[us] vs pm4py's assumed
    # datetime64[ns] write bug (2026 dates round-tripped as 1970-01-21)

    def test_round_trip_preserves_timestamps(self):
        df = dtlog.convert_timed('o:0 a:1 s:2 p:3', names=['sigma1'])
        with tempfile.TemporaryDirectory() as tmp:
            path = str(Path(tmp) / 'nested' / 'log.xes')
            dtlog.write_xes(df, path)
            self.assertTrue(Path(path).exists())
            back = pm4py.read_xes(path)

        original = df.set_index('concept:name')['time:timestamp']
        round_tripped = back.set_index('concept:name')['time:timestamp']
        for activity in ['o', 'a', 's', 'p']:
            expected = original[activity]
            actual = pd.Timestamp(round_tripped[activity]).tz_localize(None)
            self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
