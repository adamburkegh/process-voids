import time
import unittest

from lab.timing import Timer, TimingListener
from process_voids.metric_context import CellContext, Metric


class TimerTest(unittest.TestCase):
    def test_elapsed_s_is_none_before_the_block_runs(self):
        t = Timer()
        self.assertIsNone(t.elapsed_s)

    def test_elapsed_s_is_set_and_nonnegative_after_a_successful_block(self):
        with Timer() as t:
            pass
        self.assertIsNotNone(t.elapsed_s)
        self.assertGreaterEqual(t.elapsed_s, 0)

    def test_elapsed_s_reflects_roughly_how_long_the_block_took(self):
        with Timer() as t:
            time.sleep(0.05)
        self.assertGreaterEqual(t.elapsed_s, 0.05)

    def test_elapsed_s_is_still_set_when_the_block_raises(self):
        t = Timer()
        with self.assertRaises(ValueError):
            with t:
                raise ValueError('boom')
        self.assertIsNotNone(t.elapsed_s)

    def test_exception_from_the_block_still_propagates(self):
        with self.assertRaises(ValueError):
            with Timer():
                raise ValueError('boom')


class TimingListenerTest(unittest.TestCase):
    def _ctx(self, listener):
        return CellContext(log='fake_log', tree='fake_tree', listeners=[listener])

    def test_records_one_row_per_finished_metric_with_seconds(self):
        listener = TimingListener()
        ctx = self._ctx(listener)
        metric = Metric(id='m', scope='node', needs=(), compute=lambda c, node: 'ok')
        ctx.score(metric, node='n')
        self.assertEqual(len(listener.rows), 1)
        row = listener.rows[0]
        self.assertEqual(row['metric_or_stage'], 'm')
        self.assertEqual(row['status'], 'ok')
        self.assertGreaterEqual(row['seconds'], 0.0)

    def test_a_failed_metric_still_gets_a_row_marked_error(self):
        listener = TimingListener()
        ctx = self._ctx(listener)
        metric = Metric(id='m', scope='node', needs=(), compute=lambda c, node: 1 / 0)
        ctx.score(metric, node='n')
        self.assertEqual(len(listener.rows), 1)
        self.assertEqual(listener.rows[0]['status'], 'error')

    def test_records_a_row_for_a_computed_stage_too(self):
        listener = TimingListener()
        ctx = self._ctx(listener)
        ctx.STAGES = dict(ctx.STAGES, thing=lambda c: 'value')
        ctx.stage('thing')
        self.assertEqual([r['metric_or_stage'] for r in listener.rows], ['thing'])

    def test_a_memoised_stage_access_does_not_add_a_second_row(self):
        listener = TimingListener()
        ctx = self._ctx(listener)
        ctx.STAGES = dict(ctx.STAGES, thing=lambda c: 'value')
        ctx.stage('thing')
        ctx.stage('thing')
        self.assertEqual(len(listener.rows), 1)

    def test_started_events_do_not_add_rows(self):
        listener = TimingListener()
        ctx = self._ctx(listener)
        metric = Metric(id='m', scope='node', needs=(), compute=lambda c, node: 'ok')
        ctx.score(metric, node='n')
        # one row total: the finished event only, not a separate started row
        self.assertEqual(len(listener.rows), 1)


if __name__ == '__main__':
    unittest.main()
