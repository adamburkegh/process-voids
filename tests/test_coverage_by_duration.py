
import unittest

from lab.fixtures import build_running_example_log, build_running_example_tree
from process_voids.coveragemass import dur, log_to_traces, sdur, coverage_by_duration


class CoverageByDurationTest(unittest.TestCase):
    """
    Values cross-checked by hand against the running example N = seq(o,
    loop(a,e), xor(s,tau), p) over its six traces:

        sigma1 = o:0 a:1 s:2 p:3          (duration 3)
        sigma2 = o:0 a:1 e:5 a:9 s:10 p:11 (duration 11)
        sigma3 = o:0 s:10 p:11            (duration 11)
        sigma4 = o:0 s:1 p:2              (duration 2)
        sigma5 = o:0 a:1 s:10 p:11        (duration 11)
        sigma6 = o:0 a:1 p:2              (duration 2, the only trace
                                            without s - an urgent payment
                                            that skips scheduling)

    total duration = 40 hours. "mass" below is sdur-total / total-duration,
    ie coverage_by_duration with skip_probs = 0 everywhere, isolating
    the duration-attribution calculation from skip probability. dur()/
    sdur() work in seconds (datetime.total_seconds()), so raw sdur
    totals here are scaled by HOUR; mass values are ratios and are
    unaffected by the unit.
    """

    HOUR = 3600

    def setUp(self):
        self.tree = build_running_example_tree()
        self.o, self.approval, self.sched, self.p = self.tree.children
        self.a, self.e = self.approval.children
        self.s, self.tau = self.sched.children
        self.traces = log_to_traces(build_running_example_log())
        self.total_dur = sum(dur(sigma) for sigma in self.traces)
        self.zero_skip_probs = {node: 0 for node in
                                 [self.tree, self.o, self.approval, self.a,
                                  self.e, self.sched, self.s, self.tau, self.p]}

    def test_trace_durations(self):
        durations = sorted(dur(sigma) / self.HOUR for sigma in self.traces)
        self.assertEqual(durations, [2, 2, 3, 11, 11, 11])
        self.assertEqual(self.total_dur / self.HOUR, 40)

    def test_sdur_and_mass(self):
        # (node, alphabet, has_submodels, expected sdur total, expected mass)
        cases = [
            (self.tree,     {'o', 'a', 'e', 's', 'p'}, True,  40,   1.00),
            (self.o,        {'o'},                     False, 0,    0.00),
            (self.approval, {'a', 'e'},                 True,  12,   0.30),
            (self.a,        {'a'},                     False, 5.5,  0.1375),
            (self.e,        {'e'},                     False, 4,    0.10),
            (self.sched,    {'s'},                      True,  22,   0.55),
            (self.s,        {'s'},                     False, 22,   0.55),
            (self.p,        {'p'},                     False, 6,    0.15),
        ]
        for node, activities, has_submodels, expected_sdur, expected_mass in cases:
            with self.subTest(node=node):
                total_sdur = sum(sdur(activities, sigma, has_submodels)
                                  for sigma in self.traces)
                self.assertAlmostEqual(total_sdur / self.HOUR, expected_sdur,
                                        places=6)
                mass = coverage_by_duration(node, self.traces,
                                             self.zero_skip_probs, self.total_dur)
                self.assertAlmostEqual(mass, expected_mass, places=4)

    def test_skip_probability_discounts_mass(self):
        skip_probs = dict(self.zero_skip_probs)
        skip_probs[self.approval] = 0.4
        mass = coverage_by_duration(self.approval, self.traces, skip_probs,
                                     self.total_dur)
        self.assertAlmostEqual(mass, 0.6 * 0.30, places=4)


if __name__ == '__main__':
    unittest.main()
