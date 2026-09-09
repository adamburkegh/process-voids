'''
Orchestration tests for lab.exp_claims_degrade.run_claims_degradation -
sweep control flow (target-node resolution, error handling, merge-write
upsert), not the metrics themselves (covered elsewhere: lab.metrics for
compute_metrics, tests/process_voids/test_voidmass_pn_prototype.py for
voidmass_table_pn). compute_metrics and voidmass_table_pn are mocked
(both run the expensive ebi-backed pipeline); build_id_net is left real
since EbiOccurance.build_petri_net is a cheap pm4py conversion with no
alignment computation.
'''

import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from skipalignments import Activity, Sequence

from lab.exp_claims_degrade import run_claims_degradation

FAKE_METRICS = {'weight_coverage': 0.9, 'weight_voidage': 0.1, 'skipprob': 0.1,
                 'salign_coverage': 0.98}
FAKE_ROW = {'deficit': 1.0, 'movecount': 2.0,
            'voidmass_subprocess': 0.5, 'voidmass_process': 0.25,
            'alignment_mass_pooled': 0.75}

FAILING_LOG = object()


class FakeDv:
    '''Stands in for compute_metrics(..., return_dv=True)'s second
    return value - only .skip_probs is used downstream
    (coverage_by_alignment_pn), so a fixed value for every node is
    enough; defaultdict means any node in whichever fake tree a test
    builds resolves without needing to be listed explicitly.'''
    def __init__(self, value=0.2):
        self.skip_probs = defaultdict(lambda: value)


def make_tree():
    assess = Activity(None, 'assess', 100000)
    assess.id = '1'
    request_docs = Activity(None, 'request_docs', 100000)
    request_docs.id = '2'
    receive_docs = Activity(None, 'receive_docs', 100000)
    receive_docs.id = '3'
    lodge_appeal = Activity(None, 'lodge_appeal', 100000)
    lodge_appeal.id = '4'
    decide_appeal = Activity(None, 'decide_appeal', 100000)
    decide_appeal.id = '5'

    loop_block = Sequence(None, [request_docs, receive_docs])
    loop_block.id = '6'
    request_docs.set_parent(loop_block)
    receive_docs.set_parent(loop_block)

    appeal_seq = Sequence(None, [lodge_appeal, decide_appeal])
    appeal_seq.id = '7'
    lodge_appeal.set_parent(appeal_seq)
    decide_appeal.set_parent(appeal_seq)

    tree = Sequence(None, [assess, loop_block, appeal_seq])
    tree.id = '8'
    assess.set_parent(tree)
    loop_block.set_parent(tree)
    appeal_seq.set_parent(tree)
    return tree, assess, loop_block, appeal_seq


def fake_log():
    return pd.DataFrame([
        {'case:concept:name': 'c1', 'concept:name': 'assess',
         'time:timestamp': pd.Timestamp('2026-01-01')},
        {'case:concept:name': 'c1', 'concept:name': 'request_docs',
         'time:timestamp': pd.Timestamp('2026-01-01 01:00')},
    ])


def fake_ground_truth_csv(tmp_dir, deviated_cases=()):
    path = Path(tmp_dir) / 'ground_truth.csv'
    pd.DataFrame({
        'case': ['c1', 'c2', 'c3'],
        'deviation': ['reordered' if c in deviated_cases else '' for c in ('c1', 'c2', 'c3')],
    }).to_csv(path, index=False)
    return str(path)


def degrade_stub(log, target_activities, n_drop_cases, exclude_cases=None):
    return log, {f'dropped_case_{i}' for i in range(n_drop_cases)}


class BasicSweepTest(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.tmp_out = Path(self.tmp_dir) / 'out.csv'
        self.ground_truth_csv = fake_ground_truth_csv(self.tmp_dir)
        self.tree, self.assess, self.loop_block, self.appeal_seq = make_tree()
        self.fake_table = {self.tree: dict(FAKE_ROW), self.assess: dict(FAKE_ROW),
                            self.loop_block: dict(FAKE_ROW), self.appeal_seq: dict(FAKE_ROW)}

    def _run(self, **kwargs):
        with patch('lab.exp_claims_degrade.build_claims_tree', return_value=self.tree), \
             patch('lab.exp_claims_degrade.pm4py.read_xes', return_value=fake_log()), \
             patch('lab.exp_claims_degrade.degrade_target_subprocess', side_effect=degrade_stub), \
             patch('lab.exp_claims_degrade.compute_metrics', return_value=(dict(FAKE_METRICS), FakeDv())), \
             patch('lab.exp_claims_degrade.voidmass_table_pn', return_value=self.fake_table):
            return run_claims_degradation(
                target_names=['assess', 'loop_block'], n_drops=[0, 1],
                out_csv=str(self.tmp_out), ground_truth_csv=self.ground_truth_csv, **kwargs)

    def test_one_row_per_target_per_level(self):
        df = self._run()
        self.assertEqual(len(df), 4)  # 2 targets x 2 levels
        self.assertEqual(set(df['target']), {'assess', 'loop_block'})
        self.assertTrue((df['status'] == 'ok').all())

    def test_metrics_and_voidmass_columns_present(self):
        df = self._run()
        row = df.iloc[0]
        for col in ('weight_coverage', 'skipprob', 'salign_coverage',
                    'voidmass_deficit', 'voidmass_subprocess', 'voidmass_process',
                    'alignment_coverage_pn'):
            self.assertIn(col, df.columns)
        self.assertAlmostEqual(row['weight_coverage'], FAKE_METRICS['weight_coverage'])
        self.assertAlmostEqual(row['voidmass_deficit'], FAKE_ROW['deficit'])
        # FakeDv's skip_prob is 0.2, FAKE_ROW's alignment_mass_pooled is
        # 0.75 - (1 - 0.2) * 0.75
        self.assertAlmostEqual(row['alignment_coverage_pn'], 0.8 * 0.75)

    def test_dropped_case_count_reflects_the_degrade_stub(self):
        df = self._run()
        counts = dict(zip(zip(df['target'], df['n_drop_cases']), df['dropped_case_count']))
        self.assertEqual(counts[('assess', 0)], 0)
        self.assertEqual(counts[('assess', 1)], 1)


class TargetNodeResolutionTest(unittest.TestCase):
    '''appeal_seq's leaf set {lodge_appeal, decide_appeal} must resolve
    to the Sequence node itself, not some ancestor or a single leaf.'''

    def test_appeal_seq_resolves_to_the_right_node(self):
        tree, assess, loop_block, appeal_seq = make_tree()
        fake_table = {tree: dict(FAKE_ROW), assess: dict(FAKE_ROW),
                       loop_block: dict(FAKE_ROW), appeal_seq: dict(FAKE_ROW, deficit=7.0)}
        tmp_dir = tempfile.mkdtemp()
        tmp_out = Path(tmp_dir) / 'out.csv'

        with patch('lab.exp_claims_degrade.build_claims_tree', return_value=tree), \
             patch('lab.exp_claims_degrade.pm4py.read_xes', return_value=fake_log()), \
             patch('lab.exp_claims_degrade.degrade_target_subprocess', side_effect=degrade_stub), \
             patch('lab.exp_claims_degrade.compute_metrics', return_value=(dict(FAKE_METRICS), FakeDv())), \
             patch('lab.exp_claims_degrade.voidmass_table_pn', return_value=fake_table):
            df = run_claims_degradation(target_names=['appeal_seq'], n_drops=[0],
                                         out_csv=str(tmp_out),
                                         ground_truth_csv=fake_ground_truth_csv(tmp_dir))

        self.assertAlmostEqual(df.iloc[0]['voidmass_deficit'], 7.0)


class ComputeErrorTest(unittest.TestCase):
    '''A compute_metrics/voidmass_table_pn failure at one level should be
    recorded as an error row, not crash the whole sweep.'''

    def test_error_at_one_level_does_not_block_others(self):
        tree, assess, loop_block, appeal_seq = make_tree()
        fake_table = {tree: dict(FAKE_ROW), assess: dict(FAKE_ROW),
                      loop_block: dict(FAKE_ROW), appeal_seq: dict(FAKE_ROW)}
        tmp_dir = tempfile.mkdtemp()
        tmp_out = Path(tmp_dir) / 'out.csv'

        def flaky_compute_metrics(log, tree_arg, slpn_path, return_dv=False):
            if log is FAILING_LOG:
                raise RuntimeError('boom')
            return dict(FAKE_METRICS), FakeDv()

        def degrade(log, target_activities, n_drop_cases, exclude_cases=None):
            return (FAILING_LOG if n_drop_cases == 1 else log), set()

        with patch('lab.exp_claims_degrade.build_claims_tree', return_value=tree), \
             patch('lab.exp_claims_degrade.pm4py.read_xes', return_value=fake_log()), \
             patch('lab.exp_claims_degrade.degrade_target_subprocess', side_effect=degrade), \
             patch('lab.exp_claims_degrade.compute_metrics', side_effect=flaky_compute_metrics), \
             patch('lab.exp_claims_degrade.voidmass_table_pn', return_value=fake_table):
            df = run_claims_degradation(target_names=['assess'], n_drops=[0, 1],
                                         out_csv=str(tmp_out),
                                         ground_truth_csv=fake_ground_truth_csv(tmp_dir))

        self.assertEqual(len(df), 2)
        ok_row = df[df['n_drop_cases'] == 0].iloc[0]
        err_row = df[df['n_drop_cases'] == 1].iloc[0]
        self.assertEqual(ok_row['status'], 'ok')
        self.assertIn('RuntimeError: boom', err_row['status'])


class MergeWriteTargetDimensionTest(unittest.TestCase):
    def test_different_targets_coexist(self):
        tree, assess, loop_block, appeal_seq = make_tree()
        fake_table = {tree: dict(FAKE_ROW), assess: dict(FAKE_ROW),
                      loop_block: dict(FAKE_ROW), appeal_seq: dict(FAKE_ROW)}
        tmp_dir = tempfile.mkdtemp()
        tmp_out = Path(tmp_dir) / 'out.csv'
        ground_truth_csv = fake_ground_truth_csv(tmp_dir)

        with patch('lab.exp_claims_degrade.build_claims_tree', return_value=tree), \
             patch('lab.exp_claims_degrade.pm4py.read_xes', return_value=fake_log()), \
             patch('lab.exp_claims_degrade.degrade_target_subprocess', side_effect=degrade_stub), \
             patch('lab.exp_claims_degrade.compute_metrics', return_value=(dict(FAKE_METRICS), FakeDv())), \
             patch('lab.exp_claims_degrade.voidmass_table_pn', return_value=fake_table):
            run_claims_degradation(target_names=['assess'], n_drops=[0], out_csv=str(tmp_out),
                                    ground_truth_csv=ground_truth_csv)
            df = run_claims_degradation(target_names=['loop_block'], n_drops=[0], out_csv=str(tmp_out),
                                         ground_truth_csv=ground_truth_csv)

        self.assertEqual(set(df['target']), {'assess', 'loop_block'})
        self.assertEqual(len(df), 2)


class ExcludeDeviatedCasesTest(unittest.TestCase):
    '''The deviated-case exclusion must actually reach
    degrade_target_subprocess, not just be computed and dropped.'''

    def test_deviated_case_is_passed_as_exclude_cases(self):
        tree, assess, loop_block, appeal_seq = make_tree()
        fake_table = {tree: dict(FAKE_ROW), assess: dict(FAKE_ROW),
                      loop_block: dict(FAKE_ROW), appeal_seq: dict(FAKE_ROW)}
        tmp_dir = tempfile.mkdtemp()
        tmp_out = Path(tmp_dir) / 'out.csv'
        ground_truth_csv = fake_ground_truth_csv(tmp_dir, deviated_cases=['c2'])

        seen_exclude_cases = []

        def spy_degrade(log, target_activities, n_drop_cases, exclude_cases=None):
            seen_exclude_cases.append(exclude_cases)
            return degrade_stub(log, target_activities, n_drop_cases)

        with patch('lab.exp_claims_degrade.build_claims_tree', return_value=tree), \
             patch('lab.exp_claims_degrade.pm4py.read_xes', return_value=fake_log()), \
             patch('lab.exp_claims_degrade.degrade_target_subprocess', side_effect=spy_degrade), \
             patch('lab.exp_claims_degrade.compute_metrics', return_value=(dict(FAKE_METRICS), FakeDv())), \
             patch('lab.exp_claims_degrade.voidmass_table_pn', return_value=fake_table):
            run_claims_degradation(target_names=['assess'], n_drops=[1], out_csv=str(tmp_out),
                                    ground_truth_csv=ground_truth_csv)

        self.assertTrue(seen_exclude_cases)
        for exclude_cases in seen_exclude_cases:
            self.assertEqual(exclude_cases, {'c2'})


if __name__ == '__main__':
    unittest.main()
