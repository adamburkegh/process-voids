import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from skipalignments import Activity, Tau, Xor

from lab.discovery import DiscoveryCombo, DiscoveryResult
from lab.exp_disco_degrade import (
    run_disco_degrade, main, ALL_METRICS, CLASSICAL_METRIC_KEYS,
    PER_NODE_METRIC_KEYS, ALIGNED_DURATION_METRIC_KEYS, NODE_ROW_COLUMNS,
)
from process_voids.coveragemass import TREE_METRIC_KEYS
from process_voids.voidmass_pn import VoidmassPnResult


class FakeDv:
    def __init__(self, skip_probs):
        self.skip_probs = skip_probs


def _fake_log():
    """
    A minimal real DataFrame, not a bare string - the classical stage's
    own _variant_probs helper (unmocked - only voidmass_table_pn itself
    is patched) does a real groupby/sort over whatever log it's given.
    """
    return pd.DataFrame([
        {'case:concept:name': 'c1', 'concept:name': 'a', 'time:timestamp': pd.Timestamp('2026-01-01')},
    ])


def degrade_stub(log, level):
    if level == 0.0:
        return log, set()
    return log, {f'dropped_at_{level}'}


def _fake_classical_row(**overrides):
    row = {'deficit_lower': 0.0, 'deficit_upper': 0.0, 'movecount': 0.0, 'movecount_bound': 0.0,
           'voidmass_subprocess_lower': 0.0, 'voidmass_subprocess_upper': 0.0,
           'voidmass_process_lower': 0.0, 'voidmass_process_upper': 0.0}
    row.update(overrides)
    return row


class FakePipelineMixin:
    """
    Patches the two expensive pipeline stages CellContext's per-cell
    'dv' and 'classical' stages call (pvoid.skipprob - the ebi
    subprocess; voidmass_table_pn - the classical-alignment search),
    plus build_id_net, so run_disco_degrade runs against a real (but
    tiny) tree/log without any real ebi or alignment-search computation.
    mass_by_weight/voidage_by_weight/mandatory_node_count/total_node_count
    run for real (cheap, pure tree functions - need real ProcessTree
    nodes with .weight already set, transfer_pt_weights being mocked
    out); coverage_by_alignment/coverage_by_alignment_pn/voidsat are
    patched too, since their real computation needs a live Aligner.
    """
    def _patch(self, target, **kwargs):
        patcher = patch(target, **kwargs)
        self.addCleanup(patcher.stop)
        return patcher.start()

    def patch_pipeline(self, skip_probs, vm_table, timed_out_count=0, timed_out_weight=0.0,
                       skip_dict=None, read_xes_return=None,
                       salign_coverage=0.0, alignment_coverage_pn=0.0, voidsat_value=0.0):
        self._patch('lab.exp_disco_degrade.pm4py.read_xes',
                    return_value=read_xes_return if read_xes_return is not None else _fake_log())
        self._patch('lab.exp_disco_degrade.build_id_net',
                    return_value=('NET', 'IM', 'FM', {}, set(), []))
        self._patch('process_voids.metric_context.pvoid.skipprob',
                    return_value=FakeDv(skip_probs))
        self._patch('process_voids.metric_context.slpn_importer.read_slpn', return_value='FAKE_SLPN')
        self._patch('process_voids.metric_context.transfer_pt_weights')
        self.voidmass_table_pn_mock = self._patch(
            'process_voids.metric_context.voidmass_table_pn',
            return_value=VoidmassPnResult(table=vm_table, skip_dict=skip_dict or {},
                                          timed_out_count=timed_out_count,
                                          timed_out_weight=timed_out_weight))
        self._patch('lab.exp_disco_degrade.coverage_by_alignment', return_value=salign_coverage)
        self._patch('lab.exp_disco_degrade.coverage_by_alignment_pn',
                    return_value=alignment_coverage_pn)
        self._patch('lab.exp_disco_degrade.voidsat', return_value=voidsat_value)


def _single_activity_tree():
    a = Activity(None, 'a', 100000)
    a.id = '1'
    a.weight = 1
    return a


def _xor_tau_activity_tree():
    a = Activity(None, 'a', 100000)
    a.id = '1'
    a.weight = 1
    tau = Tau(None, 'tau', 0)
    tau.id = '2'
    tau.weight = 0
    choice = Xor(None, [tau, a])
    choice.id = '3'
    choice.weight = 1
    a.set_parent(choice)
    tau.set_parent(choice)
    return choice, tau, a


class ZeroLevelDedupTest(FakePipelineMixin, unittest.TestCase):
    # Level 0.0 drops nothing regardless of dimension, so it's the same
    # (log, tree) computation under every dim - run_disco_degrade should
    # compute it once, not once per dimension.

    def setUp(self):
        self.tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        self.tree = _single_activity_tree()
        self.combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(self.tree))}
        self.degradations = {'activity': degrade_stub, 'trace': degrade_stub}
        self.vm_table = {self.tree: _fake_classical_row()}

    def test_dv_stage_called_once_for_level_zero(self):
        self.patch_pipeline({self.tree: 0.1}, self.vm_table)
        with patch('process_voids.metric_context.pvoid.skipprob',
                   return_value=FakeDv({self.tree: 0.1})) as mock_skipprob:
            run_disco_degrade(['fake_log.xes'], combos=self.combos,
                              degradations=self.degradations, levels=[0.0, 0.5],
                              out_csv=str(self.tmp_out))
            # 1 call for the shared level-0.0 point + 2 dims x 1 non-zero level
            self.assertEqual(mock_skipprob.call_count, 3)

    def test_both_dims_get_a_row_at_level_zero_with_matching_metrics(self):
        self.patch_pipeline({self.tree: 0.1}, self.vm_table, salign_coverage=0.7)
        df, node_df, timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=self.combos, degradations=self.degradations, levels=[0.0],
            out_csv=str(self.tmp_out))
        self.assertEqual(len(df), 2)
        self.assertEqual(set(df['degradation_dim']), {'activity', 'trace'})
        for _, row in df.iterrows():
            self.assertEqual(row['status'], 'ok')
            self.assertEqual(row['skipprob'], 0.1)
            self.assertEqual(row['salign_coverage'], 0.7)

    def test_nonzero_levels_still_computed_per_dim_per_level(self):
        self.patch_pipeline({self.tree: 0.1}, self.vm_table)
        with patch('process_voids.metric_context.pvoid.skipprob',
                   return_value=FakeDv({self.tree: 0.1})) as mock_skipprob:
            df, node_df, timings_df = run_disco_degrade(
                ['fake_log.xes'], combos=self.combos, degradations=self.degradations,
                levels=[0.5, 1.0], out_csv=str(self.tmp_out))
            # no level-0.0 shared call here, so it's a plain 2 dims x 2 levels
            self.assertEqual(mock_skipprob.call_count, 4)
            self.assertEqual(len(df), 4)


class NotImplementedComboTest(unittest.TestCase):
    # A combo whose discover() raises NotImplementedError should be
    # skipped entirely - no metric computation at all, including at
    # level 0.0, and every row reports status='not_implemented'.

    def test_skipped_without_computing_anything(self):
        def not_implemented(log):
            raise NotImplementedError('no discovery yet')

        combos = {'stub': DiscoveryCombo('stub', not_implemented)}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        with patch('lab.exp_disco_degrade.pm4py.read_xes', return_value='FAKE_LOG'), \
             patch('process_voids.metric_context.pvoid.skipprob') as mock_skipprob:
            df, node_df, timings_df = run_disco_degrade(
                ['fake_log.xes'], combos=combos, degradations=degradations,
                levels=[0.0, 0.5], out_csv=str(tmp_out))
            mock_skipprob.assert_not_called()
            self.assertTrue((df['status'] == 'not_implemented').all())
            self.assertTrue(df['weight_coverage'].isna().all())
            self.assertTrue(df['voidmass_deficit_lower'].isna().all())
            self.assertTrue(df['timed_out_count'].isna().all())
            self.assertTrue(df['timed_out_weight'].isna().all())


class ComputeErrorTest(FakePipelineMixin, unittest.TestCase):
    # A dv-stage failure (for the shared level-0.0 call or a regular
    # per-level call) should be recorded as an error row rather than
    # crashing the sweep, and should not be mistaken for 'ok'.

    def setUp(self):
        self.tree = _single_activity_tree()
        self.combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(self.tree))}
        self.vm_table = {self.tree: _fake_classical_row()}

    def test_level_zero_error_is_recorded_and_not_crashing(self):
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        self.patch_pipeline({self.tree: 0.1}, self.vm_table)
        self._patch('process_voids.metric_context.pvoid.skipprob', side_effect=RuntimeError('boom'))

        df, node_df, timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=self.combos, degradations=degradations, levels=[0.0],
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
        base_log = _fake_log()

        def degrade(log, level):
            if level == 0.0:
                return log, set()
            return _fake_log(), {f'dropped_at_{level}'}  # a DIFFERENT log object

        degradations = {'activity': degrade}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        self.patch_pipeline({self.tree: 0.1}, self.vm_table, read_xes_return=base_log)

        def flaky(log, tree, slpn_path, ppt_weights=None):
            if log is base_log:
                return FakeDv({self.tree: 0.1})
            raise RuntimeError('boom')

        self._patch('process_voids.metric_context.pvoid.skipprob', side_effect=flaky)

        df, node_df, timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=self.combos, degradations=degradations, levels=[0.0, 0.5],
            out_csv=str(tmp_out))
        self.assertEqual(len(df), 2)
        zero_row = df[df['degradation_level'] == 0.0].iloc[0]
        nonzero_row = df[df['degradation_level'] == 0.5].iloc[0]
        self.assertEqual(zero_row['status'], 'ok')
        self.assertIn('RuntimeError: boom', nonzero_row['status'])


class SharedZeroLevelNodeRowsWeightStabilityTest(FakePipelineMixin, unittest.TestCase):
    """
    Level-0.0 node rows must agree across degradation dims.
    mass_by_weight/voidage_by_weight read tree.weight/child.weight
    directly off the shared, mutable ProcessTree object -
    transfer_pt_weights (inside the 'dv' stage) overwrites those
    attributes on EVERY cell's call, not just the shared level-0.0 one.
    So the per-node rows for level 0.0 are computed exactly once, right
    when the weight state is fresh, and that snapshot is reused (dim
    stamped in after the fact) for every dim. Recomputing them for a
    second dim would read whatever weight state an intervening,
    unrelated nonzero-level cell of the FIRST dim had left on the tree,
    silently corrupting weight_coverage/weight_voidage for every dim
    after the first.
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

        # a fully covered, b fully void - the Xor's own weight_coverage
        # is then exactly a.weight / (a.weight + b.weight), so a weight
        # swap between calls is directly visible in the number.
        self.skip_probs = {self.tree: 0.0, self.a: 0.0, self.b: 1.0}
        self.vm_table = {node: _fake_classical_row() for node in (self.tree, self.a, self.b)}
        self.call_count = 0

        def fake_skipprob(log, tree, slpn_path, ppt_weights=None):
            self.call_count += 1
            if self.call_count == 1:
                # the shared level-0.0 call, against the undegraded log
                self.a.weight, self.b.weight = 3, 1
            else:
                # a later, unrelated nonzero-level cell (first dim) -
                # simulates transfer_pt_weights re-estimating weights
                # from a DEGRADED log, mutating the SAME shared tree
                self.a.weight, self.b.weight = 1, 3
            return FakeDv(self.skip_probs)

        self.fake_skipprob = fake_skipprob

    def test_both_dims_level_zero_node_rows_agree(self):
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(self.tree))}
        degradations = {'activity': degrade_stub, 'trace': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        self.patch_pipeline(self.skip_probs, self.vm_table)
        self._patch('process_voids.metric_context.pvoid.skipprob',
                   side_effect=self.fake_skipprob)

        df, node_df, timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=combos, degradations=degradations, levels=[0.0, 0.5],
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


class ClassicalMetricsTimeoutDiagnosticsTest(FakePipelineMixin, unittest.TestCase):
    """
    run_disco_degrade carries voidmass_table_pn's per-cell timeout
    diagnostics into the root row - the count AND the summed probability
    weight, since 0.03% of a log timing out is fine and 20% is not, and
    the count alone can't tell those apart - plus voidmass_movecount_
    bound alongside the observed voidmass_movecount, rather than
    overloading one column whose meaning would depend on whether a
    timeout happened.
    """

    def test_diagnostics_and_bound_denominator_reach_the_root_row(self):
        tree = _single_activity_tree()
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(tree))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        vm_table = {tree: _fake_classical_row(deficit_lower=0.0, deficit_upper=2.0,
                                              movecount=1.0, movecount_bound=3.0,
                                              voidmass_subprocess_upper=2 / 3,
                                              voidmass_process_upper=2 / 3)}

        self.patch_pipeline({tree: 0.0}, vm_table, timed_out_count=2, timed_out_weight=0.25,
                           alignment_coverage_pn=0.5)

        df, node_df, timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=combos, degradations=degradations, levels=[0.0],
            out_csv=str(tmp_out))

        row = df.iloc[0]
        self.assertEqual(row['timed_out_count'], 2)
        self.assertEqual(row['timed_out_weight'], 0.25)
        self.assertEqual(row['voidmass_movecount'], 1.0)
        self.assertEqual(row['voidmass_movecount_bound'], 3.0)
        self.assertEqual(row['voidmass_deficit_upper'], 2.0)


class EmptyNodeCsvHasHeaderTest(FakePipelineMixin, unittest.TestCase):
    """
    When every cell in a run errors, node_rows never gets populated -
    pandas' default pd.DataFrame([]) has NO columns at all, and writing
    that produces a bare, headerless file that raises pandas.errors.
    EmptyDataError in any downstream reader expecting an empty-but-
    columned frame instead of a crash. NODE_ROW_COLUMNS fixes this by
    giving the DataFrame its real columns even with zero rows.
    """

    def test_written_csv_round_trips_as_an_empty_but_columned_frame(self):
        tree = _single_activity_tree()
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(tree))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

        self.patch_pipeline({tree: 0.0}, {tree: _fake_classical_row()})
        self._patch('process_voids.metric_context.pvoid.skipprob', side_effect=RuntimeError('boom'))

        _df, node_df, timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=combos, degradations=degradations, levels=[0.0],
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


class NodeRowsTest(FakePipelineMixin, unittest.TestCase):
    """
    Full run_disco_degrade against a small REAL tree (Xor(Tau, a)) -
    mass_by_weight/voidage_by_weight/mandatory_node_count/
    total_node_count run for real; coverage_by_alignment/
    coverage_by_alignment_pn/voidsat are mocked since their own
    correctness is covered elsewhere (test_coveragemass,
    test_voidmass_pn_prototype) - this is about the row-assembly schema
    (one row per tree node, right ids, right columns), not re-deriving
    the metric math.
    """

    def setUp(self):
        self.choice, self.tau, self.a = _xor_tau_activity_tree()
        self.skip_probs = {self.choice: 0.1, self.a: 0.2, self.tau: 0.0}
        self.vm_table = {
            self.choice: _fake_classical_row(
                deficit_lower=1.0, deficit_upper=1.0, movecount=2.0, movecount_bound=2.0,
                voidmass_subprocess_lower=0.5, voidmass_subprocess_upper=0.5,
                voidmass_process_lower=0.5, voidmass_process_upper=0.5),
            self.a: _fake_classical_row(
                deficit_lower=0.5, deficit_upper=0.5, movecount=1.0, movecount_bound=1.5,
                voidmass_subprocess_lower=0.5, voidmass_subprocess_upper=0.5,
                voidmass_process_lower=0.5, voidmass_process_upper=0.5),
            self.tau: _fake_classical_row(voidmass_process_lower=1.0, voidmass_process_upper=1.0),
        }
        self.combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(self.choice))}
        self.degradations = {'trace': degrade_stub}
        self.tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'

    def _run(self, **pipeline_kwargs):
        self.patch_pipeline(self.skip_probs, self.vm_table, **pipeline_kwargs)
        self._patch('process_voids.metric_context.pvoid.skipprob',
                   return_value=FakeDv(self.skip_probs))
        return run_disco_degrade(['fake_log.xes'], combos=self.combos,
                                 degradations=self.degradations, levels=[0.5],
                                 out_csv=str(self.tmp_out))

    def test_one_row_per_vm_table_node(self):
        _df, node_df, _timings_df = self._run()
        self.assertEqual(set(node_df['node_id']), {'1', '2', '3'})

    def test_row_carries_cell_identity_and_node_identity(self):
        _df, node_df, _timings_df = self._run()
        row = node_df[node_df['node_id'] == '1'].iloc[0]
        self.assertEqual(row['log'], 'fake_log')
        self.assertEqual(row['combo'], 'fake')
        self.assertEqual(row['degradation_dim'], 'trace')
        self.assertEqual(row['degradation_level'], 0.5)
        self.assertEqual(row['node_type'], 'Activity')
        self.assertEqual(row['alphabet'], 'a')

    def test_row_has_every_metric_column_with_correct_values(self):
        _df, node_df, _timings_df = self._run(
            salign_coverage=0.77, alignment_coverage_pn=0.88, voidsat_value=0.33)
        row = node_df[node_df['node_id'] == '1'].iloc[0]

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
        _df, node_df, _timings_df = self._run()
        tau_row = node_df[node_df['node_id'] == '2'].iloc[0]
        self.assertEqual(tau_row['node_type'], 'Tau')
        # Tau is excluded from mandatory/total counts by definition.
        self.assertEqual(tau_row['mandatory_node_count'], 0)
        self.assertEqual(tau_row['total_node_count'], 0)


class MetricsExclusionTest(FakePipelineMixin, unittest.TestCase):
    """
    run_disco_degrade(..., metrics=...) lets a caller score a SUBSET of
    ALL_METRICS - for a metric that's known-expensive (voidsat on a
    high-case-count log) or known-wrong (pending an upstream fix) and
    not worth paying for on a given run, re-collecting it later being
    cheap. An excluded metric contributes no timing row and its CSV
    column is entirely absent from the scored values (NaN once written
    to a fixed-schema CSV), never a KeyError or an inconsistent column
    set between rows - see the null_metric_values construction in
    run_disco_degrade, which is derived from the SAME metrics argument
    actually used to score, not the frozen ALL_METRICS default.
    """

    def setUp(self):
        self.tree = _single_activity_tree()
        self.combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(self.tree))}
        self.degradations = {'activity': degrade_stub}
        self.tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        self.metrics_without_voidsat = [m for m in ALL_METRICS if m.id != 'voidsat']

    def test_excluded_metric_gets_no_timing_row(self):
        self.patch_pipeline({self.tree: 0.1}, {self.tree: _fake_classical_row()},
                            voidsat_value=0.99)
        _df, _node_df, timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=self.combos, degradations=self.degradations,
            levels=[0.5], metrics=self.metrics_without_voidsat, out_csv=str(self.tmp_out))
        self.assertNotIn('voidsat', set(timings_df['metric_or_stage']))
        # a metric that WAS scored is still there, so this isn't just an
        # empty/broken timings frame
        self.assertIn('skipprob', set(timings_df['metric_or_stage']))

    def test_excluded_metric_column_is_absent_from_scored_values_not_crashed_around(self):
        self.patch_pipeline({self.tree: 0.1}, {self.tree: _fake_classical_row()},
                            voidsat_value=0.99)
        df, node_df, _timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=self.combos, degradations=self.degradations,
            levels=[0.5], metrics=self.metrics_without_voidsat, out_csv=str(self.tmp_out))
        # root df has no fixed schema (plain pd.DataFrame(rows)) - a key
        # absent from every row means the column doesn't exist at all,
        # not a NaN-filled one.
        self.assertNotIn('voidsat', df.columns)
        # node_df DOES have a fixed schema (NODE_ROW_COLUMNS), so the
        # excluded metric still gets its column - just entirely empty.
        self.assertIn('voidsat', node_df.columns)
        self.assertTrue(node_df['voidsat'].isna().all())
        # a metric that WAS scored has its real value, not also nulled
        self.assertEqual(df.iloc[0]['skipprob'], 0.1)

    def test_error_row_null_fallback_matches_the_excluded_metric_set(self):
        # a cell-wide failure (dv/classical stage) must fall back to
        # None for every SCORED metric only - if the null fallback still
        # referenced the frozen ALL_METRICS default, an error row would
        # carry a 'voidsat' key a successful row in the same run does
        # not, which is exactly the column-inconsistency this feature
        # exists to avoid.
        self.patch_pipeline({self.tree: 0.1}, {self.tree: _fake_classical_row()})
        self._patch('process_voids.metric_context.pvoid.skipprob',
                    side_effect=RuntimeError('boom'))
        df, _node_df, _timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=self.combos, degradations=self.degradations,
            levels=[0.5], metrics=self.metrics_without_voidsat, out_csv=str(self.tmp_out))
        row = df.iloc[0]
        self.assertTrue(str(row['status']).startswith('error'))
        # root df has no fixed schema - a metric excluded from BOTH the
        # success and error null-fallback paths is absent from every
        # row, so the column never exists at all. The real point of this
        # test is what did NOT happen: no KeyError from a null fallback
        # that still expected a 'voidsat' key this run never scores.
        self.assertNotIn('voidsat', df.columns)

    def test_default_metrics_argument_still_scores_everything(self):
        # no metrics= passed - existing callers (and the CLI with no
        # --exclude-metric) are unaffected.
        self.patch_pipeline({self.tree: 0.1}, {self.tree: _fake_classical_row()},
                            voidsat_value=0.42)
        df, _node_df, timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=self.combos, degradations=self.degradations,
            levels=[0.5], out_csv=str(self.tmp_out))
        self.assertEqual(df.iloc[0]['voidsat'], 0.42)
        self.assertIn('voidsat', set(timings_df['metric_or_stage']))


class TimingRowsTest(FakePipelineMixin, unittest.TestCase):
    """
    run_disco_degrade's third return value: one long-form timing row
    per stage/metric actually computed, stamped with the cell's own
    identity - see lab.timing.TimingListener and
    process_voids.metric_context.CellContext.
    """

    def test_timings_carry_cell_identity_and_cover_the_dv_and_classical_stages(self):
        tree = _single_activity_tree()
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(tree))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        self.patch_pipeline({tree: 0.1}, {tree: _fake_classical_row()})

        _df, _node_df, timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=combos, degradations=degradations, levels=[0.5],
            out_csv=str(tmp_out))

        self.assertFalse(timings_df.empty)
        row = timings_df.iloc[0]
        for col in ('log', 'combo', 'degradation_dim', 'degradation_level',
                    'metric_or_stage', 'seconds', 'status'):
            self.assertIn(col, timings_df.columns)
        self.assertEqual(set(timings_df['status']), {'ok'})
        stages_seen = set(timings_df['metric_or_stage'])
        self.assertIn('dv', stages_seen)
        self.assertIn('classical', stages_seen)
        self.assertIn('skipprob', stages_seen)

    def test_a_failed_dv_stage_still_produces_one_timing_row_marked_error(self):
        tree = _single_activity_tree()
        combos = {'fake': DiscoveryCombo('fake', lambda log: DiscoveryResult(tree))}
        degradations = {'activity': degrade_stub}
        tmp_out = Path(tempfile.mkdtemp()) / 'out.csv'
        self.patch_pipeline({tree: 0.1}, {tree: _fake_classical_row()})
        self._patch('process_voids.metric_context.pvoid.skipprob', side_effect=RuntimeError('boom'))

        _df, _node_df, timings_df = run_disco_degrade(
            ['fake_log.xes'], combos=combos, degradations=degradations, levels=[0.5],
            out_csv=str(tmp_out))

        dv_rows = timings_df[timings_df['metric_or_stage'] == 'dv']
        self.assertEqual(len(dv_rows), 1)
        self.assertEqual(dv_rows.iloc[0]['status'], 'error')


class DryRunTest(unittest.TestCase):
    """--dry-run prints the resolved Experiment and exits without
    running any real computation at all - main()'s describe-then-return
    path must never fall through into computing anything."""

    def test_dry_run_with_run_name_prints_and_computes_nothing(self):
        with patch('sys.argv', ['exp_disco_degrade', '--run', 'smoke', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('process_voids.metric_context.pvoid.skipprob') as mock_skipprob, \
             patch('lab.exp_disco_degrade.build_id_net') as mock_build_id_net:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_skipprob.assert_not_called()
        mock_build_id_net.assert_not_called()
        self.assertIn('Experiment: smoke', output)
        self.assertIn('rtfm', output)
        self.assertIn('inductive_noise20', output)
        self.assertIn('cells:', output)

    def test_dry_run_with_ad_hoc_logs_prints_and_computes_nothing(self):
        with patch('sys.argv', ['exp_disco_degrade', 'fake_log.xes',
                                 '--combos', 'inductive_noise20', '--levels', '0.0', '--dry-run']), \
             patch('lab.exp_disco_degrade.configure'), \
             patch('process_voids.metric_context.pvoid.skipprob') as mock_skipprob, \
             patch('lab.exp_disco_degrade.build_id_net') as mock_build_id_net:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_skipprob.assert_not_called()
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
             patch('process_voids.metric_context.pvoid.skipprob') as mock_skipprob, \
             patch('lab.exp_disco_degrade.build_id_net') as mock_build_id_net:
            buf = io.StringIO()
            with redirect_stdout(buf):
                main()
            output = buf.getvalue()

        mock_skipprob.assert_not_called()
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
