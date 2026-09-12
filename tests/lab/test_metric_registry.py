import ast
import inspect
import unittest

import lab.metric_registry as metric_registry_module
from lab.exp_disco_degrade import (
    CLASSICAL_METRIC_KEYS, PER_NODE_METRIC_KEYS, ALIGNED_DURATION_METRIC_KEYS,
)
from lab.exp_surprise import (
    NODE_METRIC_KEYS as SURPRISE_NODE_METRIC_KEYS,
    NODE_METRIC_BASELINE_KEYS, SUMMARY_METRIC_KEYS, SUMMARY_METRIC_BASELINE_KEYS,
)
from lab.exp_voidmass import (
    NODE_METRIC_KEYS as VOIDMASS_NODE_METRIC_KEYS, SUMMARY_METRIC_KEYS as VOIDMASS_SUMMARY_METRIC_KEYS,
)
from lab.metric_registry import METRICS, STATUSES, format_registry
from lab.metrics import METRIC_KEYS
from process_voids.coveragemass import TREE_METRIC_KEYS


def _all_emitted_ids():
    return (set(METRIC_KEYS) | set(CLASSICAL_METRIC_KEYS)
            | set(SURPRISE_NODE_METRIC_KEYS) | set(NODE_METRIC_BASELINE_KEYS)
            | set(SUMMARY_METRIC_KEYS) | set(SUMMARY_METRIC_BASELINE_KEYS)
            | set(TREE_METRIC_KEYS) | set(PER_NODE_METRIC_KEYS)
            | set(ALIGNED_DURATION_METRIC_KEYS) | set(VOIDMASS_NODE_METRIC_KEYS)
            | set(VOIDMASS_SUMMARY_METRIC_KEYS))


class DriftTest(unittest.TestCase):
    """
    Every metric id a script actually builds its output from (the
    *_KEYS constants those scripts construct their result dicts via
    zip() from - see lab.metric_registry's module docstring) must have
    exactly one registry entry, tagged with that script - and the
    registry must contain nothing else. Catches a metric renamed in its
    producing script without an accompanying registry update.
    """

    def test_compute_metrics_keys_match_registry(self):
        skip_alignment_sources = {
            'process_voids.coveragemass.mass_by_weight',
            'process_voids.coveragemass.voidage_by_weight',
            'process_voids.coveragemass.coverage_by_alignment',
            'dv.skip_probs (direct lookup, no computation of its own)',
            'process_voids.voidsalign.voidsalign',
            'lab.metrics.mean_leaf_skipprob',
        }
        skip_alignment_ids = {mid for mid, m in METRICS.items()
                               if m.source in skip_alignment_sources and m.status == 'live'
                               and 'exp_disco_degrade' in m.scripts}
        self.assertEqual(set(METRIC_KEYS), skip_alignment_ids)

    def test_per_node_metric_keys_match_registry(self):
        # weight_coverage/weight_voidage/skipprob/salign_coverage/
        # voidsalign are the same ids/sources as METRIC_KEYS above (same
        # functions/lookups, evaluated at an arbitrary node instead of
        # only the root) - mean_leaf_skipprob is the one METRIC_KEYS id
        # that's root-only, not emitted per-node (always the same
        # whole-tree average regardless of node, so a per-node column
        # would be meaningless).
        per_node_sources = {
            'process_voids.coveragemass.mass_by_weight',
            'process_voids.coveragemass.voidage_by_weight',
            'process_voids.coveragemass.coverage_by_alignment',
            'dv.skip_probs (direct lookup, no computation of its own)',
            'process_voids.voidsalign.voidsalign',
        }
        per_node_ids = {mid for mid, m in METRICS.items() if m.source in per_node_sources
                        and m.status == 'live' and 'exp_disco_degrade' in m.scripts}
        self.assertEqual(set(PER_NODE_METRIC_KEYS), per_node_ids)

    def test_classical_metric_keys_match_registry(self):
        classical_ids = {mid for mid, m in METRICS.items()
                          if m.source.startswith('process_voids.voidmass_pn') and m.status == 'live'
                          and 'exp_disco_degrade' in m.scripts}
        self.assertEqual(set(CLASSICAL_METRIC_KEYS), classical_ids)

    def test_aligned_duration_metric_keys_match_registry(self):
        aligned_duration_ids = {mid for mid, m in METRICS.items()
                                 if m.source == 'process_voids.coveragemass.voidsat'
                                 and m.status == 'live'}
        self.assertEqual(set(ALIGNED_DURATION_METRIC_KEYS), aligned_duration_ids)

    def test_tree_metric_keys_match_registry(self):
        tree_ids = {mid for mid, m in METRICS.items()
                    if m.source in ('process_voids.coveragemass.mandatory_node_count',
                                     'process_voids.coveragemass.total_node_count')
                    and m.status == 'live'}
        self.assertEqual(set(TREE_METRIC_KEYS), tree_ids)

    def test_surprise_node_keys_match_registry(self):
        node_ids = {mid for mid, m in METRICS.items()
                    if m.source in ('process_voids.surprise.surprise_totals',
                                     'process_voids.surprise.predecessor_totals')
                    and m.status == 'live'}
        self.assertEqual(set(SURPRISE_NODE_METRIC_KEYS) | set(NODE_METRIC_BASELINE_KEYS), node_ids)

    def test_surprise_summary_keys_match_registry(self):
        summary_ids = {mid for mid, m in METRICS.items()
                       if m.source == 'lab.exp_surprise._compute_variant' and m.status == 'live'}
        self.assertEqual(set(SUMMARY_METRIC_KEYS) | set(SUMMARY_METRIC_BASELINE_KEYS), summary_ids)

    def test_exp_voidmass_keys_match_registry(self):
        voidmass_ids = {mid for mid, m in METRICS.items()
                        if m.status == 'live' and m.scripts == ('exp_voidmass',)}
        self.assertEqual(set(VOIDMASS_NODE_METRIC_KEYS) | set(VOIDMASS_SUMMARY_METRIC_KEYS),
                         voidmass_ids)

    def test_live_and_evaluation_ids_are_all_emitted(self):
        live_or_evaluation = {mid for mid, m in METRICS.items()
                              if m.status in ('live', 'evaluation')}
        self.assertEqual(live_or_evaluation, _all_emitted_ids())

    def test_retired_ids_are_never_emitted(self):
        retired = {mid for mid, m in METRICS.items() if m.status == 'retired'}
        self.assertEqual(retired & _all_emitted_ids(), set())
        self.assertTrue(retired, 'expected at least one retired id to exist')

    def test_product_only_ids_are_never_emitted_by_lab(self):
        product_only = {mid for mid, m in METRICS.items() if m.status == 'product-only'}
        self.assertEqual(product_only & _all_emitted_ids(), set())
        self.assertTrue(product_only, 'expected at least one product-only id to exist')

    def test_every_status_is_a_known_status(self):
        for metric_id, metric in METRICS.items():
            self.assertIn(metric.status, STATUSES, metric_id)

    def test_superseded_by_targets_exist(self):
        for metric_id, metric in METRICS.items():
            if metric.superseded_by is not None:
                self.assertIn(metric.superseded_by, METRICS,
                              f'{metric_id}.superseded_by={metric.superseded_by!r} is not a registered id')

    def test_every_metric_declares_a_nonempty_script_list_unless_product_only(self):
        for metric_id, metric in METRICS.items():
            if metric.status == 'product-only':
                continue
            self.assertTrue(metric.scripts, f'{metric_id} declares no scripts')

    def test_format_registry_mentions_every_id(self):
        text = format_registry()
        for metric_id in METRICS:
            self.assertIn(metric_id, text)


class NoDuplicateKeysTest(unittest.TestCase):
    """
    Parses metric_registry.py's own SOURCE (not the imported METRICS
    dict) for a duplicated key in the METRICS dict literal. A duplicate
    can't be caught by inspecting METRICS itself - Python silently keeps
    the last of a repeated dict key at parse time, and every drift test
    above compares id SETS, which can't see a key that was never
    missing. A squash merge can leave two adjacent, identical entries -
    not a conflict, since both insertions match - so the registry can
    carry a dead entry while every other check passes.
    """

    def test_metrics_dict_literal_has_no_duplicate_keys(self):
        source_path = inspect.getsourcefile(metric_registry_module)
        with open(source_path, encoding='utf-8') as f:
            tree = ast.parse(f.read(), filename=source_path)

        metrics_dict = next(
            (node.value for node in ast.walk(tree)
             if isinstance(node, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id == 'METRICS' for t in node.targets)),
            None)
        self.assertIsNotNone(metrics_dict, "couldn't find a 'METRICS = {...}' assignment")
        self.assertIsInstance(metrics_dict, ast.Dict)

        keys = [k.value for k in metrics_dict.keys]
        seen = set()
        duplicates = {k for k in keys if k in seen or seen.add(k)}
        self.assertEqual(duplicates, set(),
                          f'duplicate key(s) in METRICS: {sorted(duplicates)}')


if __name__ == '__main__':
    unittest.main()
