'''
lab.mass_term_probe's actual timing finding is a real-world measurement,
not a checkable invariant - nothing here asserts a specific duration.
The one property worth protecting is the thing that made the finding
trustworthy in the first place: alternating which method runs first each
repetition, ruling out warm-cache/first-call bias as an explanation for
a difference. A future edit that silently stops alternating would still
produce a plausible-looking DataFrame, with nothing else signalling that
the exact confound this script exists to control for had come back.
'''

import unittest
from unittest.mock import patch, MagicMock

import pandas as pd
from skipalignments import Activity

from lab.mass_term_probe import run_probe


def _fake_log():
    '''A minimal real DataFrame, not a mock - _classical_stage's own
    _variant_probs and run_probe's _distinct_variant_lists both do a
    real groupby/sort over whatever log they're given.'''
    return pd.DataFrame([
        {'case:concept:name': 'c1', 'concept:name': 'a',
         'time:timestamp': pd.Timestamp('2026-01-01')},
    ])


class RunProbeAlternatesOrderTest(unittest.TestCase):
    def test_call_order_alternates_by_repetition(self):
        call_order = []

        def fake_align_sk_all(variant_strings, tree, timeout):
            call_order.append('mass_term')
            future = MagicMock()
            future.result.return_value = None
            return [future]

        def fake_voidmass_table_pn(*args, **kwargs):
            call_order.append('classical')
            return None

        tree = Activity(None, 'a', 100000)
        tree.id = '1'
        classical_net = ('NET', 'IM', 'FM', {}, set(), [])

        with patch('lab.mass_term_probe._cell',
                   return_value=(_fake_log(), tree, classical_net)), \
             patch('lab.mass_term_probe.align_sk_all', side_effect=fake_align_sk_all), \
             patch('process_voids.metric_context.voidmass_table_pn',
                   side_effect=fake_voidmass_table_pn):
            # _cell is mocked and ignores its arguments, but run_probe
            # itself resolves ALL_LOGS[log_key] before calling it - a
            # real key, any real key, is needed to get that far.
            run_probe(pairs=[('payment_approval', 'fake_combo')], reps=4)

        expected = []
        for rep in range(4):
            expected.extend(['classical', 'mass_term'] if rep % 2 == 0
                             else ['mass_term', 'classical'])
        self.assertEqual(call_order, expected)


if __name__ == '__main__':
    unittest.main()
