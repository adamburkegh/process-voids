"""
Registry of discovery algorithms used by the experiment sweeps in this
package - each combo pairs a name with a discovery method
(log -> DiscoveryResult).
"""

from dataclasses import dataclass
from typing import Callable

import pm4py_config as pm4py

from process_voids.tree import from_pm4py


@dataclass
class DiscoveryResult:
    tree: object  # ProcessTree
    ppt_weights: object = None  # (weights, loop_taus) pair, toothpaste only


@dataclass
class DiscoveryCombo:
    name: str
    discover: Callable  # (log, **kwargs) -> DiscoveryResult


def discover_inductive(log, noise_threshold=0.0):
    """
    noise_threshold direction (verified against pm4py's own IMf filter,
    __filter_dfg_noise in algo/discovery/inductive/variants/imf.py - not
    documented in pm4py's own docstring, so pinned here instead of
    re-deriving it from their source every time this comes up): a
    directly-follows edge survives only if its frequency exceeds
    noise_threshold * (the strongest competing edge from the same source
    activity). So HIGHER noise_threshold is STRICTER/MORE AGGRESSIVE
    filtering - fewer edges survive, not more. At 0.8, only edges within
    the top ~20% relative-strength band of their source activity survive;
    everything weaker gets excluded from cut-finding. This is a relative,
    per-edge frequency cutoff, not a fraction of log volume/traces/cases
    kept or discarded - there's no direct "80% of the log" reading of it.
    Whatever gets filtered out of cut-finding this way is NOT excluded
    from the resulting tree's replay, though: Inductive still wraps it in
    permissive/optional structure (Xor(Tau, ...)) rather than genuinely
    dropping it, so guaranteed fitness holds at every noise_threshold.
    """
    pt_pm4py = pm4py.discover_process_tree_inductive(
        log, noise_threshold=noise_threshold)
    return DiscoveryResult(from_pm4py(pt_pm4py))


def discover_inductive_noise80(log):
    """Aggressive filtering (see discover_inductive's docstring for the
    direction) - only near-dominant edges (>80% as strong as the
    strongest competing edge from their source) survive cut-finding."""
    return discover_inductive(log, noise_threshold=0.8)


def discover_inductive_noise20(log):
    """Mild filtering (see discover_inductive's docstring for the
    direction) - edges need only exceed 20% of their source's strongest
    competing edge to survive cut-finding, so most of the DFG stays in
    play; closer to noise_threshold=0.0 (exact fit) than to
    discover_inductive_noise80's aggressive end."""
    return discover_inductive(log, noise_threshold=0.2)


def discover_toothpaste(log):
    from lab.toothpaste_bridge import discover
    return discover(log)


def discover_toothpaste_noise10(log):
    from lab.toothpaste_bridge import discover
    return discover(log, noise=0.1)


def _not_implemented(name):
    def discover(log, **kwargs):
        raise NotImplementedError(
            f"{name} discovery is not wired up yet - invocation TBD")
    return discover


COMBOS = {
    'inductive': DiscoveryCombo('inductive', discover_inductive),
    # inductive_noise80 (aggressive filtering) is off the default roster -
    # discover_inductive_noise80 itself is untouched/still directly
    # importable, same pattern as degradation.py's 'activity' entry.
    'inductive_noise20': DiscoveryCombo('inductive_noise20', discover_inductive_noise20),
    'indulpet': DiscoveryCombo('indulpet', _not_implemented('indulpet')),
    'toothpaste': DiscoveryCombo('toothpaste', discover_toothpaste),
    'toothpaste_noise10': DiscoveryCombo('toothpaste_noise10', discover_toothpaste_noise10),
}
