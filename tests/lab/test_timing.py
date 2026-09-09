import time
import unittest

from lab.timing import Timer


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


if __name__ == '__main__':
    unittest.main()
