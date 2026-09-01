"""
Metric computation for a single (log, model) run, built on top of the
existing skip-probability pipeline in process_voids.pvoid.

Captures the metrics currently available:
  - weight_coverage: coveragemass.mass_by_weight, weighted by tree/leaf
    occurrence weights
  - skipprob: mean skip probability across Activity leaves
  - duration_coverage: coveragemass.coverage_by_duration, estimating
    coverage from the implied duration of activities present in the
    model but missing from the log

More metrics are expected to land here later.
"""

from pathlib import Path

from skipalignments import Activity

from process_voids import pvoid, slpn_importer
from process_voids.coveragemass import mass_by_weight, transfer_pt_weights, \
    infer_operator_weights, coverage_by_duration, log_to_traces, dur


def mean_skipprob(tree, skip_probs):
    values = [prob for node, prob in skip_probs.items()
              if isinstance(node, Activity)]
    return sum(values) / len(values) if values else 0.0


def compute_metrics(log, tree, slpn_path, has_estimator=True):
    """
    Run the skip-alignment pipeline for (log, tree) and return the
    weight-coverage and skipprob summary metrics.

    has_estimator: False when the discovery method already assigned leaf
    weights (eg a stochastic miner); the separate occurrence-based
    weight-fitting pass is then skipped.
    """
    Path(slpn_path).parent.mkdir(parents=True, exist_ok=True)
    dv = pvoid.skipprob(log, tree, slpn_path)
    if has_estimator:
        slpn = slpn_importer.read_slpn(slpn_path)
        transfer_pt_weights(tree, slpn)
    else:
        infer_operator_weights(tree)
    traces = log_to_traces(log)
    total_dur = sum(dur(sigma) for sigma in traces)
    return {
        'weight_coverage': mass_by_weight(tree, dv.skip_probs),
        'skipprob': mean_skipprob(tree, dv.skip_probs),
        'duration_coverage': coverage_by_duration(
            tree, traces, dv.skip_probs, total_dur),
    }
