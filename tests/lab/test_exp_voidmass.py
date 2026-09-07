'''
Orchestration tests for lab.exp_voidmass.run_voidmass_doseresponse -
sweep control flow (target-not-found skip, error handling, merge-write
upsert with the new 'target' cell dimension), not the voidmass
computation itself (see tests/process_voids/test_voidmass.py for that).

Mirrors tests/lab/test_exp_disco_degrade.py's mocking style: fake
combos, and here also a fake pvoid.skipprob/voidmass_table so no real
alignment computation runs. _discover_cached is also patched to bypass
its disk cache (var/lab/tree_cache/) - every test class here reuses the
same fake log path/combo name, so without this they'd collide on the
same cache file and pollute each other with a stale tree from whichever
test happened to run first.
'''

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from skipalignments import Activity, Sequence

from lab.discovery import DiscoveryCombo, DiscoveryResult
from lab.exp_voidmass import run_voidmass_doseresponse

FAKE_ROW = {'deficit': 1.0, 'movecount': 2.0, 'voidmass_subprocess': 0.5,
            'voidmass_process': 0.25, 'voidage_subprocess': 0.1, 'voidage_process': 0.05}

FAILING_LOG = object()


def make_tree():
    a = Activity(None, 'a', 100000)
    a.id = '1'
    b = Activity(None, 'b', 100000)
    b.id = '2'
    tree = Sequence(None, [a, b])
    tree.id = '3'
    a.set_parent(tree)
    b.set_parent(tree)
    return tree, a, b


class FakeDv:
    def __init__(self, tree, a, b):
        self.skip_dict_backup = {'a, b': ['STATE1']}
        self.pl = {('a', 'b'): 1.0}
        self.skip_probs = {tree: 0.0, a: 0.0, b: 0.0}


def fake_log():
    return pd.DataFrame([
        {'case:concept:name': 'c1', 'concept:name': 'a'},
        {'case:concept:name': 'c1', 'concept:name': 'b'},
        {'case:concept:name': 'c2', 'concept:name': 'a'},
        {'case:concept:name': 'c2', 'concept:name': 'b'},
    ])


def degrade_stub(log, target_activities, n_drop_cases):
    return log, {f'dropped_case_{i}' for i in range(n_drop_cases)}


def _uncached_discover(log_name, combo_name, combo, base_log):
    return combo.discover(base_log).tree


class BasicSweepTest(unittest.TestCase):
    def setUp(self):
        self.tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        self.tmp_summary = Path(tempfile.mkdtemp()) / 'summary.csv'
        self.tree, self.a, self.b = make_tree()
        self.combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(self.tree))}
        self.fake_table = {self.tree: dict(FAKE_ROW), self.a: dict(FAKE_ROW), self.b: dict(FAKE_ROW)}

    def _run(self, **kwargs):
        with patch('lab.exp_voidmass._discover_cached', side_effect=_uncached_discover), \
             patch('lab.exp_voidmass.pm4py.read_xes', return_value=fake_log()), \
             patch('lab.exp_voidmass.degrade_target_subprocess', side_effect=degrade_stub), \
             patch('lab.exp_voidmass.pvoid.skipprob',
                   return_value=FakeDv(self.tree, self.a, self.b)), \
             patch('lab.exp_voidmass.voidmass_table', return_value=self.fake_table):
            return run_voidmass_doseresponse(
                ['fake_log.xes'], ['a'], combos=self.combos, n_drops=[0, 1],
                out_csv=str(self.tmp_out), summary_csv=str(self.tmp_summary), **kwargs)

    def test_produces_one_summary_row_per_level(self):
        node_df, summary_df = self._run()
        self.assertEqual(len(summary_df), 2)
        self.assertEqual(set(summary_df['n_drop_cases']), {0, 1})
        self.assertTrue((summary_df['status'] == 'ok').all())

    def test_node_rows_cover_every_tree_node_per_level(self):
        node_df, _ = self._run()
        # 3 nodes (tree, a, b) x 2 levels = 6 node rows
        self.assertEqual(len(node_df), 6)

    def test_summary_reports_target_values_from_the_table(self):
        _, summary_df = self._run()
        row = summary_df.iloc[0]
        self.assertAlmostEqual(row['target_voidmass_subprocess'], FAKE_ROW['voidmass_subprocess'])
        self.assertAlmostEqual(row['target_voidage_process'], FAKE_ROW['voidage_process'])

    def test_dropped_case_count_reflects_the_degrade_stub(self):
        _, summary_df = self._run()
        counts = dict(zip(summary_df['n_drop_cases'], summary_df['dropped_case_count']))
        self.assertEqual(counts[0], 0)
        self.assertEqual(counts[1], 1)


class TargetNotFoundTest(unittest.TestCase):
    # A target activity absent from the discovered tree should be
    # skipped entirely - no rows, no crash - since there is nothing to
    # report a dose-response curve for.

    def test_missing_target_produces_no_rows(self):
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        tmp_summary = Path(tempfile.mkdtemp()) / 'summary.csv'
        tree, a, b = make_tree()
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(tree))}

        with patch('lab.exp_voidmass._discover_cached', side_effect=_uncached_discover), \
             patch('lab.exp_voidmass.pm4py.read_xes', return_value=fake_log()):
            node_df, summary_df = run_voidmass_doseresponse(
                ['fake_log.xes'], ['nonexistent_activity'], combos=combos, n_drops=[0],
                out_csv=str(tmp_out), summary_csv=str(tmp_summary))

        self.assertEqual(len(node_df), 0)
        self.assertEqual(len(summary_df), 0)


class ComputeErrorTest(unittest.TestCase):
    # A pvoid.skipprob (or voidmass_table) failure at one level should
    # be recorded as an error row, not crash the whole sweep.

    def test_error_at_one_level_does_not_block_others(self):
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        tmp_summary = Path(tempfile.mkdtemp()) / 'summary.csv'
        tree, a, b = make_tree()
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(tree))}
        fake_table = {tree: dict(FAKE_ROW), a: dict(FAKE_ROW), b: dict(FAKE_ROW)}

        def flaky_skipprob(log, tree_arg, slpn_path):
            if log is FAILING_LOG:
                raise RuntimeError('boom')
            return FakeDv(tree, a, b)

        def degrade(log, target_activities, n_drop_cases):
            return (FAILING_LOG if n_drop_cases == 1 else log), set()

        with patch('lab.exp_voidmass._discover_cached', side_effect=_uncached_discover), \
             patch('lab.exp_voidmass.pm4py.read_xes', return_value=fake_log()), \
             patch('lab.exp_voidmass.degrade_target_subprocess', side_effect=degrade), \
             patch('lab.exp_voidmass.pvoid.skipprob', side_effect=flaky_skipprob), \
             patch('lab.exp_voidmass.voidmass_table', return_value=fake_table):
            node_df, summary_df = run_voidmass_doseresponse(
                ['fake_log.xes'], ['a'], combos=combos, n_drops=[0, 1],
                out_csv=str(tmp_out), summary_csv=str(tmp_summary))

        self.assertEqual(len(summary_df), 2)
        ok_row = summary_df[summary_df['n_drop_cases'] == 0].iloc[0]
        err_row = summary_df[summary_df['n_drop_cases'] == 1].iloc[0]
        self.assertEqual(ok_row['status'], 'ok')
        self.assertIn('error: boom', err_row['status'])
        # the failing level contributes no node rows
        self.assertEqual(len(node_df), 3)


class MergeWriteTargetDimensionTest(unittest.TestCase):
    # Rerunning with a DIFFERENT target must not clobber a previous
    # target's rows - 'target' is part of the cell key, same principle
    # as lab.exp_surprise's merge-write.

    def test_different_targets_coexist(self):
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        tmp_summary = Path(tempfile.mkdtemp()) / 'summary.csv'
        tree, a, b = make_tree()
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(tree))}
        fake_table = {tree: dict(FAKE_ROW), a: dict(FAKE_ROW), b: dict(FAKE_ROW)}

        with patch('lab.exp_voidmass._discover_cached', side_effect=_uncached_discover), \
             patch('lab.exp_voidmass.pm4py.read_xes', return_value=fake_log()), \
             patch('lab.exp_voidmass.degrade_target_subprocess', side_effect=degrade_stub), \
             patch('lab.exp_voidmass.pvoid.skipprob',
                   return_value=FakeDv(tree, a, b)), \
             patch('lab.exp_voidmass.voidmass_table', return_value=fake_table):
            run_voidmass_doseresponse(['fake_log.xes'], ['a'], combos=combos, n_drops=[0],
                                       out_csv=str(tmp_out), summary_csv=str(tmp_summary))
            _, summary_df = run_voidmass_doseresponse(
                ['fake_log.xes'], ['b'], combos=combos, n_drops=[0],
                out_csv=str(tmp_out), summary_csv=str(tmp_summary))

        self.assertEqual(set(summary_df['target']), {'a', 'b'})
        self.assertEqual(len(summary_df), 2)


if __name__ == '__main__':
    unittest.main()
