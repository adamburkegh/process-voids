import io
import tempfile
import unittest
from contextlib import redirect_stdout
from datetime import datetime
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from skipalignments import Activity, Tau, Xor

from lab.discovery import DiscoveryCombo, DiscoveryResult
from lab.exp_disco_degrade import (
    run_disco_degrade, main, _node_rows, _classical_metrics, CLASSICAL_METRIC_KEYS,
    PER_NODE_METRIC_KEYS, ALIGNED_DURATION_METRIC_KEYS, NODE_ROW_COLUMNS,
)
from process_voids.coveragemass import TREE_METRIC_KEYS, log_to_traces
from process_voids.voidmass_pn import VoidmassPnResult

FAKE_METRICS = {'weight_coverage': 0.5, 'weight_voidage': 0.5, 'skipprob': 0.1,
                 'salign_coverage': 0.7}
FAKE_CLASSICAL_METRICS = {k: 0.42 for k in CLASSICAL_METRIC_KEYS}
FAKE_ALIGNED_DURATION_METRICS = {k: 0.55 for k in ALIGNED_DURATION_METRIC_KEYS}
FAKE_DV = object()  # compute_metrics(..., return_dv=True)'s second value - opaque here,
                    # since _classical_metrics/_aligned_duration_metrics are both
                    # mocked in every test below


def degrade_stub(log, level):
    if level == 0.0:
        return log, set()
    return f'{log}_degraded', {f'dropped_at_{level}'}


# Every test below also mocks build_id_net (a cheap conversion in
# reality, but 'FAKE_TREE' isn't a real ProcessTree so it can't run for
# real here) and _classical_metrics itself (bypasses voidmass_table_pn/
# coverage_by_alignment_pn - the expensive classical-alignment
# computation, out of scope for these orchestration tests - see
# tests/process_voids/test_voidmass_pn_prototype.py for that).


class ZeroLevelDedupTest(unittest.TestCase):
    # Level 0.0 drops nothing regardless of dimension, so it's the same
    # (log, tree) computation under every dim - run_disco_degrade should
    # call compute_metrics for it once, not once per dimension.

    def setUp(self):
        self.tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        self.combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult('FAKE_TREE'))}
        self.degradations = {'activity': degrade_stub, 'trace': degrade_stub}

    def test_compute_metrics_called_once_for_level_zero(self):
        with patch('lab.exp_disco_degrade.compute_metrics',
                   return_value=(dict(FAKE_METRICS), FAKE_DV)) as mock_cm, \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {}, {}, {})), \
             patch('lab.exp_disco_degrade._aligned_duration_metrics',
                   return_value=dict(FAKE_ALIGNED_DURATION_METRICS)), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            run_disco_degrade(['fake_log.xes'], combos=self.combos,
                               degradations=self.degradations, levels=[0.0, 0.5],
                               out_csv=str(self.tmp_out))
            # 1 call for the shared level-0.0 point + 2 dims x 1 non-zero level
            self.assertEqual(mock_cm.call_count, 3)

    def test_both_dims_get_a_row_at_level_zero_with_matching_metrics(self):
        with patch('lab.exp_disco_degrade.compute_metrics',
                   return_value=(dict(FAKE_METRICS), FAKE_DV)), \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {}, {}, {})), \
             patch('lab.exp_disco_degrade._aligned_duration_metrics',
                   return_value=dict(FAKE_ALIGNED_DURATION_METRICS)), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=self.combos,
                                    degradations=self.degradations, levels=[0.0],
                                    out_csv=str(self.tmp_out))
            self.assertEqual(len(df), 2)
            self.assertEqual(set(df['degradation_dim']), {'activity', 'trace'})
            for _, row in df.iterrows():
                self.assertEqual(row['status'], 'ok')
                self.assertEqual(row['weight_coverage'], FAKE_METRICS['weight_coverage'])
                self.assertEqual(row['salign_coverage'], FAKE_METRICS['salign_coverage'])
                self.assertEqual(row['voidmass_deficit_lower'],
                                  FAKE_CLASSICAL_METRICS['voidmass_deficit_lower'])

    def test_nonzero_levels_still_computed_per_dim_per_level(self):
        with patch('lab.exp_disco_degrade.compute_metrics',
                   return_value=(dict(FAKE_METRICS), FAKE_DV)) as mock_cm, \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {}, {}, {})), \
             patch('lab.exp_disco_degrade._aligned_duration_metrics',
                   return_value=dict(FAKE_ALIGNED_DURATION_METRICS)), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=self.combos,
                                    degradations=self.degradations, levels=[0.5, 1.0],
                                    out_csv=str(self.tmp_out))
            # no level-0.0 shared call here, so it's a plain 2 dims x 2 levels
            self.assertEqual(mock_cm.call_count, 4)
            self.assertEqual(len(df), 4)


class NotImplementedComboTest(unittest.TestCase):
    # A combo whose discover() raises NotImplementedError should be
    # skipped entirely - no compute_metrics call at all, including at
    # level 0.0, and every row reports status='not_implemented'.

    def test_skipped_without_computing_anything(self):
        def not_implemented(log):
            raise NotImplementedError('no discovery yet')

        combos = {'stub': DiscoveryCombo('stub', not_implemented)}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        with patch('lab.exp_disco_degrade.compute_metrics') as mock_cm, \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {}, {}, {})), \
             patch('lab.exp_disco_degrade._aligned_duration_metrics',
                   return_value=dict(FAKE_ALIGNED_DURATION_METRICS)), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                    degradations=degradations, levels=[0.0, 0.5],
                                    out_csv=str(tmp_out))
            mock_cm.assert_not_called()
            self.assertTrue((df['status'] == 'not_implemented').all())
            self.assertTrue(df['weight_coverage'].isna().all())
            self.assertTrue(df['voidmass_deficit_lower'].isna().all())
            self.assertTrue(df['timed_out_count'].isna().all())
            self.assertTrue(df['timed_out_weight'].isna().all())


class ComputeMetricsErrorTest(unittest.TestCase):
    # A compute_metrics failure (for the shared level-0.0 call or a
    # regular per-level call) should be recorded as an error row rather
    # than crashing the sweep, and should not be mistaken for 'ok'.

    def test_level_zero_error_is_recorded_and_not_crashing(self):
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult('FAKE_TREE'))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        with patch('lab.exp_disco_degrade.compute_metrics',
                   side_effect=RuntimeError('boom')), \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {}, {}, {})), \
             patch('lab.exp_disco_degrade._aligned_duration_metrics',
                   return_value=dict(FAKE_ALIGNED_DURATION_METRICS)), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                    degradations=degradations, levels=[0.0],
                                    out_csv=str(tmp_out))
            self.assertEqual(len(df), 1)
            self.assertIn('RuntimeError: boom', df.iloc[0]['status'])
            self.assertIsNone(df.iloc[0]['weight_coverage'])
            self.assertIsNone(df.iloc[0]['voidmass_deficit_lower'])

            # Every cell errored, so node_rows never got populated - the
            # written _nodes CSV must still carry a real header (see
            # EmptyNodeCsvHasHeaderTest for the full round-trip check;
            # this just confirms THIS scenario is the one that triggers
            # it) rather than a bare, columnless empty file.
            self.assertTrue(node_df.empty)
            self.assertGreater(len(node_df.columns), 0)

    def test_nonzero_level_error_does_not_block_other_levels(self):
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult('FAKE_TREE'))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        def flaky(log, tree, slpn_path, ppt_weights=None, return_dv=False):
            if log == 'FAKE_LOG':
                return dict(FAKE_METRICS), FAKE_DV
            raise RuntimeError('boom')

        with patch('lab.exp_disco_degrade.compute_metrics', side_effect=flaky), \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {}, {}, {})), \
             patch('lab.exp_disco_degrade._aligned_duration_metrics',
                   return_value=dict(FAKE_ALIGNED_DURATION_METRICS)), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                    degradations=degradations, levels=[0.0, 0.5],
                                    out_csv=str(tmp_out))
            self.assertEqual(len(df), 2)
            zero_row = df[df['degradation_level'] == 0.0].iloc[0]
            nonzero_row = df[df['degradation_level'] == 0.5].iloc[0]
            self.assertEqual(zero_row['status'], 'ok')
            self.assertIn('RuntimeError: boom', nonzero_row['status'])


class SharedZeroLevelNodeRowsWeightStabilityTest(unittest.TestCase):
    """
    Regression test for a bug found in the 2026-09-10 runs: mass_by_weight/voidage_
    by_weight read tree.weight/child.weight directly off the shared,
    mutable ProcessTree object - transfer_pt_weights (inside
    compute_metrics) overwrites those attributes on EVERY cell's call,
    not just the shared level-0.0 one. A second dim reusing the level-
    0.0 result used to call _node_rows again, reading whatever weight
    state an intervening, unrelated nonzero-level cell of the FIRST dim
    had already left on the tree - silently corrupting weight_coverage/
    weight_voidage for every dim after the first. Fixed by computing
    _node_rows for level 0.0 exactly once, right when the weight state
    is fresh, and reusing that same snapshot (dim stamped in after the
    fact) for every dim instead of recomputing.
    """

    def setUp(self):
        self.a = Activity(None, 'a', 100000)
        self.a.id = '1'
        self.b = Activity(None, 'b', 100000)
        self.b.id = '2'
        self.tree = Xor(None, [self.a, self.b])
        self.tree.id = '3'
        self.a.set_parent(self.tree)
        self.b.set_parent(self.tree)

        class FakeDv:
            def __init__(self, skip_probs):
                self.skip_probs = skip_probs

        # a fully covered, b fully void - the Xor's own weight_coverage
        # is then exactly a.weight / (a.weight + b.weight), so a weight
        # swap between calls is directly visible in the number.
        self.dv = FakeDv({self.tree: 0.0, self.a: 0.0, self.b: 1.0})

        fake_vm_row = {'deficit_lower': 0.0, 'deficit_upper': 0.0, 'movecount': 0.0,
                       'movecount_bound': 0.0,
                       'voidmass_subprocess_lower': 0.0, 'voidmass_subprocess_upper': 0.0,
                       'voidmass_process_lower': 0.0, 'voidmass_process_upper': 0.0}
        self.vm_table = {self.tree: dict(fake_vm_row), self.a: dict(fake_vm_row),
                          self.b: dict(fake_vm_row)}

        self.call_count = 0

        def fake_compute_metrics(log, tree, slpn_path, ppt_weights=None, return_dv=False):
            self.call_count += 1
            if self.call_count == 1:
                # the shared level-0.0 call, against the undegraded log
                self.a.weight, self.b.weight = 3, 1
            else:
                # a later, unrelated nonzero-level cell (first dim) -
                # simulates transfer_pt_weights re-estimating weights
                # from a DEGRADED log, mutating the SAME shared tree
                self.a.weight, self.b.weight = 1, 3
            return dict(FAKE_METRICS), self.dv

        self.fake_compute_metrics = fake_compute_metrics

    def test_both_dims_level_zero_node_rows_agree(self):
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(self.tree))}
        degradations = {'activity': degrade_stub, 'trace': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        with patch('lab.exp_disco_degrade.compute_metrics',
                   side_effect=self.fake_compute_metrics), \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), self.vm_table, {}, {})), \
             patch('lab.exp_disco_degrade._aligned_duration_metrics',
                   return_value=dict(FAKE_ALIGNED_DURATION_METRICS)), \
             patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=0.0), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.0), \
             patch('lab.exp_disco_degrade.voidsat', return_value=0.0), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            df, node_df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                    degradations=degradations, levels=[0.0, 0.5],
                                    out_csv=str(tmp_out))

        # sanity: the mutation actually happened (zero-level call plus
        # at least one nonzero-level call) - otherwise this test would
        # pass vacuously without exercising the bug at all.
        self.assertGreaterEqual(self.call_count, 2)

        root_rows = node_df[(node_df['node_id'] == '3') & (node_df['degradation_level'] == 0.0)]
        self.assertEqual(len(root_rows), 2)  # one per dim
        weight_coverages = set(root_rows['weight_coverage'])
        self.assertEqual(
            len(weight_coverages), 1,
            f"weight_coverage diverged across dims at level 0.0: "
            f"{root_rows[['degradation_dim', 'weight_coverage']].to_dict('records')}")
        # and it must be the snapshot from the FIRST (zero-level) call
        # (a=3,b=1 -> 0.75), not whatever a later cell's mutation left
        # behind (a=1,b=3 -> 0.25).
        self.assertAlmostEqual(weight_coverages.pop(), 0.75, places=6)


class ClassicalMetricsTimeoutDiagnosticsTest(unittest.TestCase):
    """
    _classical_metrics carries voidmass_table_pn's per-cell timeout
    diagnostics into the root row - the count AND the summed probability
    weight, since 0.03% of a log timing out is fine and 20% is not, and
    the count alone can't tell those apart - plus voidmass_movecount_
    bound alongside the observed voidmass_movecount, rather than
    overloading one column whose meaning would depend on whether a
    timeout happened.
    """

    def test_diagnostics_and_bound_denominator_reach_the_metrics(self):
        tree = Activity(None, 'a', 100000)
        tree.id = 'a'
        root_row = {'deficit_lower': 0.0, 'deficit_upper': 2.0,
                    'movecount': 1.0, 'movecount_bound': 3.0,
                    'voidmass_subprocess_lower': 0.0, 'voidmass_subprocess_upper': 2 / 3,
                    'voidmass_process_lower': 0.0, 'voidmass_process_upper': 2 / 3}
        result = VoidmassPnResult(table={tree: root_row}, skip_dict={},
                                  timed_out_count=2, timed_out_weight=0.25)
        log = pd.DataFrame({'case:concept:name': ['c1'], 'concept:name': ['a'],
                            'time:timestamp': [pd.Timestamp('2026-01-01')]})

        class FakeDv:
            skip_probs = {tree: 0.0}

        with patch('lab.exp_disco_degrade.voidmass_table_pn', return_value=result), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.5):
            metrics, *_rest = _classical_metrics(tree, log, 'NET', 'IM', 'FM', {}, set(), [],
                                                 FakeDv())

        self.assertEqual(metrics['timed_out_count'], 2)
        self.assertEqual(metrics['timed_out_weight'], 0.25)
        self.assertEqual(metrics['voidmass_movecount'], 1.0)
        self.assertEqual(metrics['voidmass_movecount_bound'], 3.0)
        self.assertEqual(metrics['voidmass_deficit_upper'], 2.0)


class EmptyNodeCsvHasHeaderTest(unittest.TestCase):
    """
    When every cell in a run errors, node_rows never gets populated -
    pandas' default pd.DataFrame([]) has NO columns at all, and writing
    that produces a bare, headerless file that raises pandas.errors.
    EmptyDataError in any downstream reader expecting an empty-but-
    columned frame instead of a crash. NODE_ROW_COLUMNS fixes this by
    giving the DataFrame its real columns even with zero rows.
    """

    def test_written_csv_round_trips_as_an_empty_but_columned_frame(self):
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult('FAKE_TREE'))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        with patch('lab.exp_disco_degrade.compute_metrics',
                   side_effect=RuntimeError('boom')), \
             patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('lab.exp_disco_degrade.build_id_net', return_value=('NET', 'IM', 'FM', {}, set(), [])), \
             patch('lab.exp_disco_degrade._classical_metrics',
                   return_value=(dict(FAKE_CLASSICAL_METRICS), {}, {}, {})), \
             patch('lab.exp_disco_degrade._aligned_duration_metrics',
                   return_value=dict(FAKE_ALIGNED_DURATION_METRICS)), \
             patch('lab.exp_disco_degrade.mandatory_node_count', return_value=1), \
             patch('lab.exp_disco_degrade.total_node_count', return_value=1):
            _df, node_df = run_disco_degrade(['fake_log.xes'], combos=combos,
                                    degradations=degradations, levels=[0.0],
                                    out_csv=str(tmp_out))

        self.assertTrue(node_df.empty)
        self.assertEqual(list(node_df.columns), NODE_ROW_COLUMNS)

        node_out_csv = tmp_out.with_name('out_nodes.csv')
        self.assertTrue(node_out_csv.exists())
        # The real regression: pd.read_csv on a genuinely headerless
        # empty file raises EmptyDataError - this must not.
        reread = pd.read_csv(node_out_csv)
        self.assertTrue(reread.empty)
        self.assertEqual(list(reread.columns), NODE_ROW_COLUMNS)


class NodeRowsTest(unittest.TestCase):
    """
    _node_rows against a small REAL tree (Xor(Tau, a)) - mass_by_weight/
    voidage_by_weight/mandatory_node_count/total_node_count run for
    real; coverage_by_alignment/coverage_by_alignment_pn are mocked
    since their own correctness is covered elsewhere (test_coveragemass,
    test_voidmass_pn_prototype) - this test is about the row-assembly
    schema (one row per vm_table node, right ids, right columns), not
    re-deriving the metric math.
    """

    def setUp(self):
        self.a = Activity(None, 'a', 100000)
        self.a.id = '1'
        self.tau = Tau(None, 'tau', 0)
        self.tau.id = '2'
        self.choice = Xor(None, [self.tau, self.a])
        self.choice.id = '3'
        self.a.set_parent(self.choice)
        self.tau.set_parent(self.choice)
        # mass_by_weight/infer_operator_weights need .weight already set
        # on every node (including Tau) - normally done by
        # transfer_pt_weights inside compute_metrics before _node_rows
        # ever sees the tree; set directly here since this test bypasses
        # that pipeline.
        self.a.weight = 1
        self.tau.weight = 0
        self.choice.weight = 1

        class FakeDv:
            def __init__(self, skip_probs):
                self.skip_probs = skip_probs

        self.dv = FakeDv({self.choice: 0.1, self.a: 0.2, self.tau: 0.0})
        self.vm_table = {
            self.choice: {'deficit_lower': 1.0, 'deficit_upper': 1.0, 'movecount': 2.0,
                           'movecount_bound': 2.0,
                           'voidmass_subprocess_lower': 0.5, 'voidmass_subprocess_upper': 0.5,
                           'voidmass_process_lower': 0.5, 'voidmass_process_upper': 0.5},
            self.a: {'deficit_lower': 0.5, 'deficit_upper': 0.5, 'movecount': 1.0,
                     'movecount_bound': 1.5,
                     'voidmass_subprocess_lower': 0.5, 'voidmass_subprocess_upper': 0.5,
                     'voidmass_process_lower': 0.5, 'voidmass_process_upper': 0.5},
            self.tau: {'deficit_lower': 0.0, 'deficit_upper': 0.0, 'movecount': 0.0,
                       'movecount_bound': 0.0,
                       'voidmass_subprocess_lower': 0.0, 'voidmass_subprocess_upper': 0.0,
                       'voidmass_process_lower': 1.0, 'voidmass_process_upper': 1.0},
        }
        # coverage_by_alignment_pn/voidsat are mocked in every test
        # below, so these just need to exist to satisfy _node_rows'
        # signature - their actual contents are never read.
        self.variant_probs = {}
        self.skip_dict = {}
        self.log = 'FAKE_LOG'

    def test_one_row_per_vm_table_node(self):
        with patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=0.77), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.88), \
             patch('lab.exp_disco_degrade.voidsat', return_value=0.33):
            rows = _node_rows('mylog', 'mycombo', 'trace', 0.5, self.dv, self.vm_table,
                               self.variant_probs, self.skip_dict, self.log)
        self.assertEqual({r['node_id'] for r in rows}, {'1', '2', '3'})

    def test_row_carries_cell_identity_and_node_identity(self):
        with patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=0.77), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.88), \
             patch('lab.exp_disco_degrade.voidsat', return_value=0.33):
            rows = _node_rows('mylog', 'mycombo', 'trace', 0.5, self.dv, self.vm_table,
                               self.variant_probs, self.skip_dict, self.log)
        row = next(r for r in rows if r['node_id'] == '1')
        self.assertEqual(row['log'], 'mylog')
        self.assertEqual(row['combo'], 'mycombo')
        self.assertEqual(row['degradation_dim'], 'trace')
        self.assertEqual(row['degradation_level'], 0.5)
        self.assertEqual(row['node_type'], 'Activity')
        self.assertEqual(row['alphabet'], 'a')

    def test_row_has_every_metric_column_with_correct_values(self):
        with patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=0.77), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.88), \
             patch('lab.exp_disco_degrade.voidsat', return_value=0.33):
            rows = _node_rows('mylog', 'mycombo', 'trace', 0.5, self.dv, self.vm_table,
                               self.variant_probs, self.skip_dict, self.log)
        row = next(r for r in rows if r['node_id'] == '1')

        for key in (PER_NODE_METRIC_KEYS + CLASSICAL_METRIC_KEYS
                    + ALIGNED_DURATION_METRIC_KEYS + TREE_METRIC_KEYS):
            self.assertIn(key, row)

        self.assertEqual(row['skipprob'], 0.2)
        self.assertEqual(row['voidmass_deficit_lower'], 0.5)
        self.assertEqual(row['voidmass_deficit_upper'], 0.5)
        self.assertEqual(row['voidmass_movecount'], 1.0)
        self.assertEqual(row['voidmass_movecount_bound'], 1.5)
        self.assertEqual(row['salign_coverage'], 0.77)
        self.assertEqual(row['alignment_coverage_pn_lower'], 0.88)
        self.assertEqual(row['alignment_coverage_pn_upper'], 0.88)
        self.assertEqual(row['voidsat'], 0.33)
        # a's own subtree is just itself, no silent alternative from its
        # own perspective (mandatory_node_count/total_node_count are
        # evaluated AT that node, not from an outside ancestor's view).
        self.assertEqual(row['total_node_count'], 1)

    def test_tau_node_gets_a_row_too(self):
        with patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=0.0), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.0), \
             patch('lab.exp_disco_degrade.voidsat', return_value=0.0):
            rows = _node_rows('mylog', 'mycombo', 'trace', 0.5, self.dv, self.vm_table,
                               self.variant_probs, self.skip_dict, self.log)
        tau_row = next(r for r in rows if r['node_id'] == '2')
        self.assertEqual(tau_row['node_type'], 'Tau')
        # Tau is excluded from mandatory/total counts by definition.
        self.assertEqual(tau_row['mandatory_node_count'], 0)
        self.assertEqual(tau_row['total_node_count'], 0)


class NodeRowsBuildsTracesOnceTest(unittest.TestCase):
    """
    voidsat's trace list (log_to_traces: a group-by, sort and conversion
    over the whole log) depends only on the log, never on the node being
    scored, so one _node_rows call must build it exactly once - not once
    per node. voidsat runs for real here; only the alignment-based
    metrics are mocked.
    """

    def setUp(self):
        self.a = Activity(None, 'a', 100000)
        self.a.id = '1'
        self.tau = Tau(None, 'tau', 0)
        self.tau.id = '2'
        self.choice = Xor(None, [self.tau, self.a])
        self.choice.id = '3'
        self.a.set_parent(self.choice)
        self.tau.set_parent(self.choice)
        self.a.weight = 1
        self.tau.weight = 0
        self.choice.weight = 1

        class FakeDv:
            def __init__(self, skip_probs):
                self.skip_probs = skip_probs
                self.skip_dict_backup = {}

        self.dv = FakeDv({self.choice: 0.1, self.a: 0.2, self.tau: 0.0})
        row = {'deficit_lower': 0.0, 'deficit_upper': 0.0,
               'movecount': 0.0, 'movecount_bound': 0.0,
               'voidmass_subprocess_lower': 0.0, 'voidmass_subprocess_upper': 0.0,
               'voidmass_process_lower': 0.0, 'voidmass_process_upper': 0.0}
        self.vm_table = {node: dict(row) for node in (self.choice, self.a, self.tau)}
        self.log = [[{'concept:name': 'a', 'time:timestamp': datetime(2026, 1, 1)}]]

    def test_log_to_traces_runs_once_per_node_rows_call(self):
        with patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=0.0), \
             patch('lab.exp_disco_degrade.coverage_by_alignment_pn', return_value=0.0), \
             patch('process_voids.coveragemass.log_to_traces', wraps=log_to_traces) as counted:
            _node_rows('mylog', 'mycombo', 'trace', 0.5, self.dv, self.vm_table, {}, {},
                       self.log)
        self.assertEqual(counted.call_count, 1)


class DryRunTest(unittest.TestCase):
    """--dry-run prints the resolved Experiment and exits without
    calling compute_metrics or run_disco_degrade at all - main()'s
    describe-then-return path must never fall through into computing
    anything."""

    def test_dry_run_with_run_name_prints_and_computes_nothing(self):
        with patch('sys.argv', ['exp_disco_degrade', '--run', 'smoke', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('lab.exp_disco_degrade.compute_metrics') as mock_cm, \
             patch('lab.exp_disco_degrade.build_id_net') as mock_build_id_net:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_cm.assert_not_called()
        mock_build_id_net.assert_not_called()
        self.assertIn('Experiment: smoke', output)
        self.assertIn('rtfm', output)
        self.assertIn('inductive_noise20', output)
        self.assertIn('cells:', output)

    def test_dry_run_with_ad_hoc_logs_prints_and_computes_nothing(self):
        with patch('sys.argv', ['exp_disco_degrade', 'fake_log.xes',
                                 '--combos', 'inductive_noise20', '--levels', '0.0', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('lab.exp_disco_degrade.compute_metrics') as mock_cm, \
             patch('lab.exp_disco_degrade.build_id_net') as mock_build_id_net:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_cm.assert_not_called()
        mock_build_id_net.assert_not_called()
        self.assertIn('Experiment: ad hoc', output)
        self.assertIn('fake_log', output)

    def test_dry_run_accepts_claims_fixture_combo_and_degradation_names(self):
        """--combos/--degradations must also resolve claims_known/
        appeal_seq etc - these aren't in ALL_COMBOS/ALL_DEGRADATIONS
        (real-discovery only), they come from lab.claims_fixture's
        registrations, merged in only for this ad hoc lookup."""
        with patch('sys.argv', ['exp_disco_degrade', 'data/claims.xes',
                                 '--combos', 'claims_known',
                                 '--degradations', 'appeal_seq', 'loop_block',
                                 '--levels', '0.0', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('lab.exp_disco_degrade.compute_metrics') as mock_cm, \
             patch('lab.exp_disco_degrade.build_id_net') as mock_build_id_net:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_cm.assert_not_called()
        mock_build_id_net.assert_not_called()
        self.assertIn('claims_known', output)
        self.assertIn('appeal_seq', output)
        self.assertIn('loop_block', output)

    def test_unknown_combo_name_errors_out(self):
        with patch('sys.argv', ['exp_disco_degrade', 'fake_log.xes',
                                 '--combos', 'not_a_real_combo', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'):
            with self.assertRaises(SystemExit):
                main()


if __name__ == '__main__':
    unittest.main()
