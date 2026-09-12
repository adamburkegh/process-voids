import tempfile
import unittest
from pathlib import Path

from process_voids.metric_context import ProcessMetric, CellContext, METRIC_ERROR, score_all


class RecordingListener:
    def __init__(self):
        self.events = []

    def __call__(self, event, ctx, id_, node, **extra):
        self.events.append((event, id_, node, extra))


def _make_ctx(stages=None, listeners=()):
    ctx = CellContext(log='fake_log', tree='fake_tree', listeners=listeners)
    if stages:
        ctx.STAGES = dict(ctx.STAGES, **stages)
    return ctx


class StageMemoisationTest(unittest.TestCase):
    def test_stage_computed_once_and_reused(self):
        calls = []

        def factory(ctx):
            calls.append(1)
            return object()

        ctx = _make_ctx(stages={'thing': factory})
        first = ctx.stage('thing')
        second = ctx.stage('thing')
        self.assertIs(first, second)
        self.assertEqual(len(calls), 1)

    def test_different_contexts_do_not_share_stage_cache(self):
        def factory(ctx):
            return object()

        ctx1 = _make_ctx(stages={'thing': factory})
        ctx2 = _make_ctx(stages={'thing': factory})
        self.assertIsNot(ctx1.stage('thing'), ctx2.stage('thing'))

    def test_stage_factory_receives_the_context(self):
        seen = []
        ctx = _make_ctx(stages={'thing': lambda c: seen.append(c)})
        ctx.stage('thing')
        self.assertEqual(seen, [ctx])

    def test_a_failed_stage_is_memoised_as_an_error_not_recomputed(self):
        calls = []

        def flaky(ctx):
            calls.append(1)
            raise RuntimeError('stage boom')

        ctx = _make_ctx(stages={'flaky': flaky})
        with self.assertRaises(RuntimeError):
            ctx.stage('flaky')
        # second access re-raises the SAME cached failure, not a fresh computation -
        # an expensive stage (eg the ebi-backed 'dv' stage) must run at most once
        # per cell even if several metrics need it and it fails.
        with self.assertRaises(RuntimeError) as second:
            ctx.stage('flaky')
        self.assertEqual(len(calls), 1)
        self.assertEqual(str(second.exception), 'stage boom')


class StageLifecycleEventsTest(unittest.TestCase):
    def test_events_fire_once_around_the_computing_call_only(self):
        listener = RecordingListener()
        ctx = _make_ctx(stages={'thing': lambda c: 'value'}, listeners=[listener])
        ctx.stage('thing')
        ctx.stage('thing')  # memoised - no new events
        events = [(e, id_) for e, id_, node, extra in listener.events]
        self.assertEqual(events, [('stage_started', 'thing'), ('stage_finished', 'thing')])

    def test_finished_event_carries_elapsed_seconds(self):
        listener = RecordingListener()
        ctx = _make_ctx(stages={'thing': lambda c: 'value'}, listeners=[listener])
        ctx.stage('thing')
        finished = [extra for e, id_, node, extra in listener.events if e == 'stage_finished']
        self.assertEqual(len(finished), 1)
        self.assertIn('elapsed_s', finished[0])
        self.assertGreaterEqual(finished[0]['elapsed_s'], 0.0)

    def test_a_failing_stage_emits_stage_failed_not_stage_finished(self):
        listener = RecordingListener()
        boom = RuntimeError('stage boom')

        def flaky(ctx):
            raise boom

        ctx = _make_ctx(stages={'flaky': flaky}, listeners=[listener])
        with self.assertRaises(RuntimeError):
            ctx.stage('flaky')
        event_names = [e for e, id_, node, extra in listener.events]
        self.assertEqual(event_names, ['stage_started', 'stage_failed'])
        failed_extra = listener.events[1][3]
        self.assertIs(failed_extra['exception'], boom)
        self.assertIn('elapsed_s', failed_extra)

    def test_a_memoised_failed_stage_emits_no_further_events_on_later_access(self):
        listener = RecordingListener()
        ctx = _make_ctx(stages={'flaky': lambda c: 1 / 0}, listeners=[listener])
        with self.assertRaises(ZeroDivisionError):
            ctx.stage('flaky')
        with self.assertRaises(ZeroDivisionError):
            ctx.stage('flaky')
        self.assertEqual(len(listener.events), 2)  # just the one started/failed pair


class ProcessMetricScoringTest(unittest.TestCase):
    def test_score_calls_compute_with_context_and_node(self):
        seen = []
        metric = ProcessMetric(id='m', scope='node', needs=(),
                         compute=lambda ctx, node: seen.append((ctx, node)) or 'ok')
        ctx = _make_ctx()
        result = ctx.score(metric, node='node1')
        self.assertEqual(result, 'ok')
        self.assertEqual(seen, [(ctx, 'node1')])

    def test_root_scope_defaults_node_to_the_tree(self):
        seen = []
        metric = ProcessMetric(id='m', scope='root', needs=(),
                         compute=lambda ctx, node: seen.append(node))
        ctx = _make_ctx()
        ctx.score(metric)
        self.assertEqual(seen, [ctx.tree])

    def test_compute_can_read_a_needed_stage_via_ctx_stage(self):
        metric = ProcessMetric(id='m', scope='node', needs=('thing',),
                         compute=lambda ctx, node: ctx.stage('thing'))
        ctx = _make_ctx(stages={'thing': lambda c: 42})
        self.assertEqual(ctx.score(metric, node='n'), 42)

    def test_metric_events_fire_in_order(self):
        listener = RecordingListener()
        metric = ProcessMetric(id='m', scope='node', needs=(), compute=lambda ctx, node: 'ok')
        ctx = _make_ctx(listeners=[listener])
        ctx.score(metric, node='n')
        events = [(e, id_, node) for e, id_, node, extra in listener.events]
        self.assertEqual(events, [('metric_started', 'm', 'n'), ('metric_finished', 'm', 'n')])

    def test_metric_finished_event_carries_elapsed_seconds(self):
        listener = RecordingListener()
        metric = ProcessMetric(id='m', scope='node', needs=(), compute=lambda ctx, node: 'ok')
        ctx = _make_ctx(listeners=[listener])
        ctx.score(metric, node='n')
        finished = [extra for e, id_, node, extra in listener.events if e == 'metric_finished']
        self.assertEqual(len(finished), 1)
        self.assertIn('elapsed_s', finished[0])


class ProcessMetricErrorIsolationTest(unittest.TestCase):
    def test_a_raising_metric_returns_the_error_sentinel_not_raise(self):
        metric = ProcessMetric(id='m', scope='node', needs=(),
                         compute=lambda ctx, node: 1 / 0)
        ctx = _make_ctx()
        self.assertIs(ctx.score(metric, node='n'), METRIC_ERROR)

    def test_a_raising_metric_emits_metric_failed_with_the_exception(self):
        listener = RecordingListener()
        boom = ValueError('boom')

        def compute(ctx, node):
            raise boom

        metric = ProcessMetric(id='m', scope='node', needs=(), compute=compute)
        ctx = _make_ctx(listeners=[listener])
        ctx.score(metric, node='n')
        failed = [(id_, node, extra) for e, id_, node, extra in listener.events
                  if e == 'metric_failed']
        self.assertEqual(len(failed), 1)
        id_, node, extra = failed[0]
        self.assertEqual((id_, node), ('m', 'n'))
        self.assertIs(extra['exception'], boom)
        self.assertIn('elapsed_s', extra)

    def test_no_metric_started_finished_pair_on_failure(self):
        listener = RecordingListener()
        metric = ProcessMetric(id='m', scope='node', needs=(),
                         compute=lambda ctx, node: 1 / 0)
        ctx = _make_ctx(listeners=[listener])
        ctx.score(metric, node='n')
        event_names = [e for e, id_, node, extra in listener.events]
        self.assertEqual(event_names, ['metric_started', 'metric_failed'])

    def test_error_sentinel_is_distinct_from_a_genuine_none_result(self):
        none_metric = ProcessMetric(id='m', scope='node', needs=(), compute=lambda ctx, node: None)
        ctx = _make_ctx()
        result = ctx.score(none_metric, node='n')
        self.assertIsNone(result)
        self.assertIsNot(result, METRIC_ERROR)

    def test_a_failing_stage_reached_via_compute_runs_once_even_for_two_metrics(self):
        """
        The realistic failure mode: a metric's compute() calls ctx.stage(),
        not stage() called directly. An expensive, failing stage (eg the
        ebi-backed 'dv' stage) must run at most once per cell even though
        every metric needing it fails independently - not be retried once
        per metric.
        """
        calls = []

        def flaky_dv(ctx):
            calls.append(1)
            raise RuntimeError('ebi boom')

        metric_a = ProcessMetric(id='a', scope='node', needs=('dv',),
                          compute=lambda ctx, node: ctx.stage('dv'))
        metric_b = ProcessMetric(id='b', scope='node', needs=('dv',),
                          compute=lambda ctx, node: ctx.stage('dv'))
        listener = RecordingListener()
        ctx = _make_ctx(stages={'dv': flaky_dv}, listeners=[listener])

        result_a = ctx.score(metric_a, node='n')
        result_b = ctx.score(metric_b, node='n')

        self.assertIs(result_a, METRIC_ERROR)
        self.assertIs(result_b, METRIC_ERROR)
        self.assertEqual(len(calls), 1)
        metric_failed_ids = [id_ for e, id_, node, extra in listener.events if e == 'metric_failed']
        self.assertEqual(metric_failed_ids, ['a', 'b'])
        stage_events = [e for e, id_, node, extra in listener.events if id_ == 'dv']
        self.assertEqual(stage_events, ['stage_started', 'stage_failed'])


class WeightMutationHazardTest(unittest.TestCase):
    """
    Mirrors process_voids.coveragemass's tree.weight mutation hazard
    (transfer_pt_weights overwrites .weight in place on the one tree
    object every cell reuses) with a synthetic shared mutable object,
    proving the contract that protects it: one CellContext per cell,
    fully scored before the next cell's context is built.
    """

    def test_cell_n_weight_dependent_value_survives_building_cell_n_plus_1(self):
        shared = {'weight': None}

        def mutate_stage(ctx):
            shared['weight'] = ctx.log  # stands in for transfer_pt_weights, keyed by this cell's log
            return shared

        weight_metric = ProcessMetric(id='w', scope='root', needs=('mutate',),
                                compute=lambda c, node: c.stage('mutate')['weight'])

        ctx1 = _make_ctx(stages={'mutate': mutate_stage})
        ctx1.log = 'log1'
        value1 = ctx1.score(weight_metric)  # scored to completion before ctx2 exists

        ctx2 = _make_ctx(stages={'mutate': mutate_stage})
        ctx2.log = 'log2'
        ctx2.score(weight_metric)  # mutates the SAME shared object

        self.assertEqual(value1, 'log1')
        self.assertEqual(shared['weight'], 'log2')


class ScoreAllTest(unittest.TestCase):
    """
    score_all(ctx, metrics, node) - the row-assembly counterpart to
    CellContext.score: scores every declared ProcessMetric against one
    node and returns {id: value}, keeping ProcessMetric itself a single-
    quantity declaration rather than a bundle.
    """

    def test_returns_one_value_per_metric_keyed_by_id(self):
        metrics = [
            ProcessMetric(id='a', scope='node', needs=(), compute=lambda ctx, node: node + 1),
            ProcessMetric(id='b', scope='node', needs=(), compute=lambda ctx, node: node * 2),
        ]
        ctx = _make_ctx()
        self.assertEqual(score_all(ctx, metrics, node=5), {'a': 6, 'b': 10})

    def test_a_failing_metric_contributes_its_error_sentinel_not_raise(self):
        metrics = [
            ProcessMetric(id='ok', scope='node', needs=(), compute=lambda ctx, node: 1),
            ProcessMetric(id='bad', scope='node', needs=(), compute=lambda ctx, node: 1 / 0),
        ]
        ctx = _make_ctx()
        result = score_all(ctx, metrics, node='n')
        self.assertEqual(result['ok'], 1)
        self.assertIs(result['bad'], METRIC_ERROR)

    def test_each_metric_fires_its_own_started_finished_events(self):
        listener = RecordingListener()
        metrics = [
            ProcessMetric(id='a', scope='node', needs=(), compute=lambda ctx, node: 1),
            ProcessMetric(id='b', scope='node', needs=(), compute=lambda ctx, node: 2),
        ]
        ctx = _make_ctx(listeners=[listener])
        score_all(ctx, metrics, node='n')
        ids = [id_ for e, id_, node, extra in listener.events]
        self.assertEqual(ids, ['a', 'a', 'b', 'b'])

    def test_root_level_use_passes_the_tree_as_node(self):
        seen = []
        metrics = [ProcessMetric(id='a', scope='root', needs=(),
                                 compute=lambda ctx, node: seen.append(node))]
        ctx = _make_ctx()
        score_all(ctx, metrics, node=ctx.tree)
        self.assertEqual(seen, [ctx.tree])


class ClassicalStageIntegrationTest(unittest.TestCase):
    """
    The 'classical' stage reproduces process_voids.voidmass_pn.
    voidmass_table_pn's own values for the same (tree, log) - proves it
    runs the SAME pipeline exp_disco_degrade._classical_metrics does
    today, not a reimplementation.
    """

    def test_classical_stage_matches_a_direct_voidmass_table_pn_call(self):
        from lab.fixtures import build_running_example_log, build_running_example_tree
        from process_voids.voidmass_pn import build_id_net, voidmass_table_pn

        log = build_running_example_log()
        tree = build_running_example_tree()
        net, im, fm, activity_to_id, tau_ids, id_loop_list = build_id_net(tree)

        def variant_probs(log):
            n_cases = log['case:concept:name'].nunique()
            variants = {}
            for _case, group in log.groupby('case:concept:name', sort=False):
                variant = tuple(group.sort_values('time:timestamp')['concept:name'])
                variants[variant] = variants.get(variant, 0) + 1
            return {v: c / n_cases for v, c in variants.items()}

        expected = voidmass_table_pn(tree, variant_probs(log), net, im, fm,
                                     activity_to_id, tau_ids, id_loop_list=id_loop_list)

        ctx = CellContext(log=log, tree=tree,
                          classical_net=(net, im, fm, activity_to_id, tau_ids, id_loop_list))
        result, actual_variant_probs = ctx.stage('classical')

        self.assertEqual(result.table[tree]['deficit_lower'], expected.table[tree]['deficit_lower'])
        self.assertEqual(result.table[tree]['movecount'], expected.table[tree]['movecount'])
        self.assertEqual(actual_variant_probs, variant_probs(log))


class RealStageIntegrationTest(unittest.TestCase):
    """
    Builds a real CellContext against the running-example fixture and
    checks it reproduces lab.metrics.compute_metrics' values for the same
    (log, tree) - proves the machinery runs the SAME pipeline, not just
    that it runs something.
    """

    def test_skipprob_and_weight_coverage_match_compute_metrics(self):
        from lab.fixtures import build_running_example_log, build_running_example_tree
        from lab.metrics import compute_metrics
        from process_voids.coveragemass import mass_by_weight

        log = build_running_example_log()
        # Temp dir for the intermediate SLPNs: nothing here creates var/lab/.
        with tempfile.TemporaryDirectory() as tmp:
            expected_metrics = compute_metrics(
                log, build_running_example_tree(),
                str(Path(tmp) / 'test_metric_context_expected.slpn'))

            ctx = CellContext(log=log, tree=build_running_example_tree(),
                              slpn_path=str(Path(tmp) / 'test_metric_context_ctx.slpn'))
            skipprob_metric = ProcessMetric(id='skipprob', scope='root', needs=('dv',),
                                     compute=lambda c, node: c.stage('dv').skip_probs[node])
            weight_coverage_metric = ProcessMetric(
                id='weight_coverage', scope='root', needs=('dv',),
                compute=lambda c, node: mass_by_weight(node, c.stage('dv').skip_probs))

            self.assertAlmostEqual(ctx.score(skipprob_metric), expected_metrics['skipprob'])
            self.assertAlmostEqual(ctx.score(weight_coverage_metric),
                                   expected_metrics['weight_coverage'])


if __name__ == '__main__':
    unittest.main()
