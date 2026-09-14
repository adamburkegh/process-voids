"""
Metric computation for a single (log, model) run, built on top of the
existing skip-probability pipeline in process_voids.pvoid.

METRIC_KEYS: weight_coverage, weight_voidage, skipprob, salign_coverage,
voidsalign2, mean_leaf_skipprob - see lab.metric_registry for what each id
means and where it's computed. duration_coverage is product-only, not on
this roster.
"""

from pathlib import Path

from skipalignments import Activity

from process_voids import pvoid, slpn_importer
from process_voids.coveragemass import mass_by_weight, transfer_pt_weights, \
    coverage_by_alignment, voidage_by_weight
from process_voids.voidsalign2 import voidsalign2


METRIC_KEYS = ('weight_coverage', 'weight_voidage', 'skipprob', 'salign_coverage',
               'voidsalign2', 'mean_leaf_skipprob')


def mean_leaf_skipprob(tree, skip_probs):
    """
    Mean of skip_probs[leaf] over every Activity leaf in `skip_probs` -
    see lab.metric_registry for how this differs from skipprob. Ignores
    its own `tree` argument for scoping: averages over every Activity
    anywhere in skip_probs regardless of which subtree was passed, so
    calling this on a non-root node returns the same value as the root.
    """
    values = [prob for node, prob in skip_probs.items()
              if isinstance(node, Activity)]
    return sum(values) / len(values) if values else 0.0


def compute_metrics(log, tree, slpn_path, ppt_weights=None, return_dv=False):
    """
    Run the skip-alignment pipeline for (log, tree) and return the
    METRIC_KEYS metrics, scored at the tree root.

    ppt_weights: the (weights, loop_taus) pair from a toothpaste
    discovery - passed through to pvoid.skipprob so DerivationPipeline
    uses DiscoverySource.TOOTHPASTE (exact PPT weights) instead of
    estimating occurrence-based weights from the log. Either way,
    compute() writes a resolved SLPN to slpn_path (estimated for
    OCCURANCE, exactly compiled from the PPT for TOOTHPASTE) with the
    same transitions/label/weight shape, so the same transfer_pt_weights
    read-back works unconditionally for both.

    return_dv=True also returns the computed DerivationPipeline (as
    (metrics, dv)) - for a caller that additionally needs dv.skip_probs
    itself (eg process_voids.voidmass_pn.coverage_by_alignment_pn reuses
    it unchanged rather than recomputing anything skip-alignments
    already gives us for free) without paying for the expensive pipeline
    twice.
    """
    Path(slpn_path).parent.mkdir(parents=True, exist_ok=True)
    dv = pvoid.skipprob(log, tree, slpn_path, ppt_weights=ppt_weights)
    slpn = slpn_importer.read_slpn(slpn_path)
    transfer_pt_weights(tree, slpn)
    values = (
        mass_by_weight(tree, dv.skip_probs),
        voidage_by_weight(tree, dv.skip_probs),
        dv.skip_probs[tree],
        coverage_by_alignment(tree, dv),
        voidsalign2(tree, dv.skip_dict_backup, dv.pl, dv.skip_probs),
        mean_leaf_skipprob(tree, dv.skip_probs),
    )
    metrics = dict(zip(METRIC_KEYS, values))
    return (metrics, dv) if return_dv else metrics
