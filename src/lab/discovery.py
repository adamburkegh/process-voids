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
    pt_pm4py = pm4py.discover_process_tree_inductive(
        log, noise_threshold=noise_threshold)
    return DiscoveryResult(from_pm4py(pt_pm4py))


def discover_toothpaste(log):
    from lab.toothpaste_bridge import discover
    return discover(log)


def _not_implemented(name):
    def discover(log, **kwargs):
        raise NotImplementedError(
            f"{name} discovery is not wired up yet - invocation TBD")
    return discover


COMBOS = {
    'inductive': DiscoveryCombo('inductive', discover_inductive),
    'indulpet': DiscoveryCombo('indulpet', _not_implemented('indulpet')),
    'toothpaste': DiscoveryCombo('toothpaste', discover_toothpaste),
}
