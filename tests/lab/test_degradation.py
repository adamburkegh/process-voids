'''
Unit tests for lab.degradation, focused on degrade_activity_wise_gradual
(the piecewise-linear ramp) since the plain step version is simple
enough to be covered indirectly by test_exp_disco_degrade.py's stubs.
'''

import unittest

import pandas as pd

from lab.degradation import (
    degrade_activity_wise, degrade_activity_wise_gradual, degrade_target_subprocess,
)


def make_log(activities, n_per_activity=10):
    '''A log with n_per_activity events per activity, one case per event
    (case identity doesn't matter for activity-wise degradation).'''
    rows = []
    i = 0
    for a in activities:
        for _ in range(n_per_activity):
            rows.append({'case:concept:name': f'c{i}', 'concept:name': a})
            i += 1
    return pd.DataFrame(rows)


class GradualActivityDegradationTest(unittest.TestCase):
    def setUp(self):
        self.activities = ['a', 'b', 'c', 'd']
        self.log = make_log(self.activities, n_per_activity=10)

    def _counts(self, level):
        degraded, dropped = degrade_activity_wise_gradual(self.log, level)
        counts = degraded['concept:name'].value_counts().to_dict()
        return {a: counts.get(a, 0) for a in self.activities}, dropped

    def test_zero_level_drops_nothing(self):
        counts, dropped = self._counts(0.0)
        self.assertEqual(dropped, set())
        self.assertEqual(counts, {a: 10 for a in self.activities})

    def test_full_level_drops_everything(self):
        counts, dropped = self._counts(1.0)
        self.assertEqual(dropped, set(self.activities))
        self.assertEqual(counts, {a: 0 for a in self.activities})

    def test_only_one_activity_is_partial_at_a_time(self):
        # k=4 activities, so level 0.375 = 1.5 activity-units of work:
        # one fully dropped, one at ~50%, two untouched.
        counts, dropped = self._counts(0.375)
        n_zero = sum(1 for c in counts.values() if c == 0)
        n_full = sum(1 for c in counts.values() if c == 10)
        n_partial = sum(1 for c in counts.values() if 0 < c < 10)
        self.assertEqual(n_zero, 1)
        self.assertEqual(n_full, 2)
        self.assertEqual(n_partial, 1)
        self.assertEqual(len(dropped), 1)

    def test_monotonically_nested_across_levels(self):
        # whatever is missing at a lower level must still be missing at
        # every higher level (per-event, not just per-activity)
        _, dropped_low = degrade_activity_wise_gradual(self.log, 0.2)
        degraded_low, _ = degrade_activity_wise_gradual(self.log, 0.2)
        degraded_high, _ = degrade_activity_wise_gradual(self.log, 0.6)
        remaining_low = set(degraded_low.index)
        remaining_high = set(degraded_high.index)
        self.assertTrue(remaining_high.issubset(remaining_low))

    def test_ramp_is_continuous_not_a_step(self):
        # unlike the step version, small level increases within a single
        # activity's ramp should change the event count
        counts_a, _ = self._counts(0.1)
        counts_b, _ = self._counts(0.2)
        self.assertNotEqual(counts_a, counts_b)

    def test_step_version_is_unaffected_by_this_change(self):
        # sanity: the original function still behaves exactly as before
        degraded, dropped = degrade_activity_wise(self.log, 0.5)
        self.assertEqual(len(dropped), 2)
        counts = degraded['concept:name'].value_counts().to_dict()
        self.assertTrue(all(c in (0, 10) for c in counts.values()))


def make_case_log(case_activities):
    '''case_activities: {case_id: [activities...]}.'''
    rows = []
    for case, activities in case_activities.items():
        for a in activities:
            rows.append({'case:concept:name': case, 'concept:name': a})
    return pd.DataFrame(rows)


class TargetSubprocessDegradationTest(unittest.TestCase):
    def setUp(self):
        # 3 cases contain 'x' (the target), 2 don't - eligible = 3
        self.log = make_case_log({
            'c1': ['o', 'x', 'p'],
            'c2': ['o', 'x', 'p'],
            'c3': ['o', 'x', 'p'],
            'c4': ['o', 'p'],
            'c5': ['o', 'p'],
        })

    def test_zero_drops_leaves_log_unchanged(self):
        degraded, dropped = degrade_target_subprocess(self.log, {'x'}, 0)
        self.assertEqual(dropped, set())
        self.assertEqual(len(degraded), len(self.log))
        self.assertEqual((degraded['concept:name'] == 'x').sum(), 3)

    def test_partial_drop_removes_only_target_from_dropped_cases(self):
        degraded, dropped = degrade_target_subprocess(self.log, {'x'}, 2)
        self.assertEqual(len(dropped), 2)
        self.assertTrue(dropped.issubset({'c1', 'c2', 'c3'}))
        # exactly 1 'x' event remains (3 eligible cases - 2 dropped)
        self.assertEqual((degraded['concept:name'] == 'x').sum(), 1)
        # non-target activities in dropped cases are untouched
        for case in dropped:
            case_rows = degraded[degraded['case:concept:name'] == case]
            self.assertEqual(set(case_rows['concept:name']), {'o', 'p'})
        # untouched cases (c4, c5, and whichever of c1-c3 wasn't dropped)
        # still have all their original events
        self.assertEqual(len(degraded), len(self.log) - 2)

    def test_full_drop_removes_target_from_whole_log(self):
        degraded, dropped = degrade_target_subprocess(self.log, {'x'}, 3)
        self.assertEqual(dropped, {'c1', 'c2', 'c3'})
        self.assertEqual((degraded['concept:name'] == 'x').sum(), 0)

    def test_exceeding_eligible_cases_raises(self):
        with self.assertRaises(ValueError):
            degrade_target_subprocess(self.log, {'x'}, 4)

    def test_monotonically_nested_across_counts(self):
        _, dropped_1 = degrade_target_subprocess(self.log, {'x'}, 1)
        _, dropped_2 = degrade_target_subprocess(self.log, {'x'}, 2)
        self.assertTrue(dropped_1.issubset(dropped_2))

    def test_ineligible_cases_never_touched(self):
        for n in (0, 1, 2, 3):
            with self.subTest(n=n):
                _, dropped = degrade_target_subprocess(self.log, {'x'}, n)
                self.assertFalse(dropped & {'c4', 'c5'})


if __name__ == '__main__':
    unittest.main()
