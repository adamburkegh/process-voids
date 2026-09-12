"""
Registry of discovery algorithms used by the experiment sweeps in this
package - each combo pairs a name with a discovery method
(log -> DiscoveryResult).
"""

import pickle
from dataclasses import dataclass
from pathlib import Path
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


TREE_CACHE_DIR = Path('var/lab/tree_cache')

# Cache files hold a (tree, ppt_weights) pair. The suffix distinguishes
# them from files written when they held a bare tree, which unpickle
# without error but unpack into the wrong shape.
_CACHE_SUFFIX = 'pair'


def discover_cached(log_name, combo_name, combo, base_log):
    """
    (tree, ppt_weights) for a (log, combo) pair, from
    var/lab/tree_cache/ if it has been discovered before, otherwise
    discovered now and cached.

    Shared by every experiment script: discovery depends only on the log
    and the combo, never on degradation dim/level, so it is paid once -
    which matters most for toothpaste's external-subprocess path. It
    also fixes the tree across processes, where pm4py's Inductive cut
    selection is otherwise hash-seed dependent near a tie (rtfm draws a
    12- or 13-node tree from the same log), so two runs of nominally the
    same experiment can score different models. The cached tree is an
    artefact of a run's provenance: the tie-break it froze is arbitrary,
    not authoritative.

    ppt_weights is toothpaste's fixed PPT weights, which only
    exp_disco_degrade consumes (threaded into pvoid.skipprob); callers
    on the occurrence-weighted combos unpack and discard it. Keyed by
    log NAME, not log content - delete the cache file (or the whole
    var/lab/tree_cache/ dir) if the underlying log changes.

    A failing discovery raises and caches nothing, leaving the caller's
    own status handling (exp_disco_degrade's discover_status) intact.
    """
    cache_path = TREE_CACHE_DIR / f'{log_name}__{combo_name}__{_CACHE_SUFFIX}.pkl'
    if cache_path.exists():
        with open(cache_path, 'rb') as f:
            return pickle.load(f)
    result = combo.discover(base_log)
    pair = (result.tree, result.ppt_weights)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, 'wb') as f:
        pickle.dump(pair, f)
    return pair


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
    # Vanilla inductive (noise_threshold=0.0) and inductive_noise80
    # (aggressive filtering) are both off the default roster -
    # discover_inductive/discover_inductive_noise80 stay directly
    # importable, same pattern as degradation.py's step-wise
    # degrade_activity_wise.
    #
    # Vanilla is off because it makes a poor test: discovered at exact
    # fit from the same log it is then checked against, it absorbs
    # dropped activities at no cost, staying flat-zero on the voidmass
    # metrics until degradation level 0.6 on rtfm where toothpaste
    # responds from 0.2.
    # inductive_noise20 gives a more responsive signal for the same
    # per-cell cost.
    'inductive_noise20': DiscoveryCombo('inductive_noise20', discover_inductive_noise20),
    # indulpet is off the roster until its invocation is worked out -
    # the _not_implemented stub factory is kept for whenever that
    # happens. On the roster it only ever contributed 20 rows of
    # status='not_implemented' per run (2 dims x 10 levels) at zero
    # compute, diluting every results CSV for no signal.
    'toothpaste': DiscoveryCombo('toothpaste', discover_toothpaste),
    'toothpaste_noise10': DiscoveryCombo('toothpaste_noise10', discover_toothpaste_noise10),
}
