"""
Registry of discovery-algorithm + estimator combinations used by the
experiment sweeps in this package.

Each combo pairs a discovery method (log -> ProcessTree) with whether a
separate occurrence-based weight-estimation pass is needed afterwards.
Miners that are themselves stochastic (e.g. toothpaste) assign leaf
weights directly and skip that pass.
"""

from dataclasses import dataclass
from typing import Callable

import pm4py_config as pm4py

from process_voids.tree import from_pm4py


@dataclass
class DiscoveryCombo:
    name: str
    discover: Callable  # (log, **kwargs) -> ProcessTree
    has_estimator: bool  # False if discovery already assigns leaf weights


def discover_inductive(log, noise_threshold=0.0):
    pt_pm4py = pm4py.discover_process_tree_inductive(
        log, noise_threshold=noise_threshold)
    return from_pm4py(pt_pm4py)


def _not_implemented(name):
    def discover(log, **kwargs):
        raise NotImplementedError(
            f"{name} discovery is not wired up yet - invocation TBD")
    return discover


COMBOS = {
    'inductive': DiscoveryCombo('inductive', discover_inductive, has_estimator=True),
    'indulpet': DiscoveryCombo('indulpet', _not_implemented('indulpet'), has_estimator=True),
    'toothpaste': DiscoveryCombo('toothpaste', _not_implemented('toothpaste'), has_estimator=False),
}
