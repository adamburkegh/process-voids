import unittest

from lab.exp_disco_degrade import CLASSICAL_METRIC_KEYS
from lab.exp_surprise import NODE_METRIC_KEYS, SUMMARY_METRIC_KEYS
from lab.metric_registry import METRICS, format_registry
from lab.metrics import METRIC_KEYS


class DriftTest(unittest.TestCase):
    """
    Every metric id a script actually builds its output from (the
    *_KEYS constants those scripts construct their result dicts via
    zip() from - see lab.metric_registry's module docstring) must have
    exactly one registry entry, tagged with that script - and the
    registry must contain nothing else. Catches a metric renamed in its
    producing script without an accompanying registry update, the same
    silent-drift failure mode that left plot_dose_response's METRICS
    list stale earlier in this project's history.
    """

    def _ids_for(self, script_name):
        return {metric_id for metric_id, metric in METRICS.items()
                if script_name in metric.scripts}

    def test_compute_metrics_keys_match_registry(self):
        registered = self._ids_for('exp_disco_degrade') & self._ids_for('exp_claims_degrade')
        skip_alignment_ids = {mid for mid, m in METRICS.items()
                               if m.source.startswith('process_voids.coveragemass')
                               or m.source == 'lab.metrics.mean_skipprob'}
        self.assertEqual(set(METRIC_KEYS), skip_alignment_ids)

    def test_classical_metric_keys_match_registry(self):
        classical_ids = {mid for mid, m in METRICS.items()
                          if m.source.startswith('process_voids.voidmass_pn')}
        self.assertEqual(set(CLASSICAL_METRIC_KEYS), classical_ids)

    def test_surprise_node_keys_match_registry(self):
        node_ids = {mid for mid, m in METRICS.items()
                    if m.source in ('process_voids.surprise.surprise_totals',
                                     'process_voids.surprise.predecessor_totals')}
        self.assertEqual(set(NODE_METRIC_KEYS), node_ids)

    def test_surprise_summary_keys_match_registry(self):
        summary_ids = {mid for mid, m in METRICS.items()
                       if m.source == 'lab.exp_surprise._compute_variant'}
        self.assertEqual(set(SUMMARY_METRIC_KEYS), summary_ids)

    def test_every_registered_id_is_emitted_by_at_least_one_script(self):
        all_emitted = (set(METRIC_KEYS) | set(CLASSICAL_METRIC_KEYS)
                       | set(NODE_METRIC_KEYS) | set(SUMMARY_METRIC_KEYS))
        self.assertEqual(set(METRICS), all_emitted)

    def test_every_metric_declares_a_nonempty_script_list(self):
        for metric_id, metric in METRICS.items():
            self.assertTrue(metric.scripts, f'{metric_id} declares no scripts')

    def test_format_registry_mentions_every_id(self):
        text = format_registry()
        for metric_id in METRICS:
            self.assertIn(metric_id, text)


if __name__ == '__main__':
    unittest.main()
