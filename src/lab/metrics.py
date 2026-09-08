"""
Metric computation for a single (log, model) run, built on top of the
existing skip-probability pipeline in process_voids.pvoid.

Captures the metrics currently available:
  - weight_coverage: coveragemass.mass_by_weight, a structural weighted
    average over the tree driven by .weight (occurrence frequency read
    off an estimated SLPN) and skip_probs (at the leaf only, 1 -
    skip_probs[leaf]). Not built on skip-alignments' alignment
    machinery at all - no executions()/Aligner involved anywhere in it,
    only a dict lookup for skip_probs, the same input skipprob itself
    uses directly. Known incorrect independent of that (see session
    notes) - the bug is in the weight-interpolation logic itself, not
    which alignment engine supplies skip_probs, so this keeps the
    unqualified weight_coverage name; a fix replaces it in place rather
    than sitting alongside a differently-sourced version.
  - weight_voidage: coveragemass.voidage_by_weight - the same tree
    aggregation as weight_coverage, but skip_probs[leaf] directly
    instead of its complement. Exactly 1 - weight_coverage at every
    node (averaging distributes linearly over the complement), so this
    is a convenience read-out, not a separately-derived quantity.
  - skipprob: mean skip probability across Activity leaves - skip-
    alignments' own dv.skip_probs output, unmodified beyond averaging.
  - salign_coverage: coveragemass.coverage_by_alignment (zero
    convention - matches the paper; favours flagging a possible void
    over silently absorbing a subprocess that ran but wasn't recorded).
    Named _salign (skip-alignment): unlike weight_coverage, this one
    genuinely is built on skip-alignments' machinery - alignment_mass ->
    executions() -> Aligner.align2's lumped normal form, the same
    representation voidmass's deficit had to move off of. Kept under
    this name until a classical-alignment replacement lands under the
    unqualified alignment_coverage name.

duration_coverage (coveragemass.coverage_by_duration) has been dropped
from this roster - not computed here for now.

More metrics are expected to land here later.
"""

from pathlib import Path

from skipalignments import Activity

from process_voids import pvoid, slpn_importer
from process_voids.coveragemass import mass_by_weight, transfer_pt_weights, \
    coverage_by_alignment, voidage_by_weight


METRIC_KEYS = ('weight_coverage', 'weight_voidage', 'skipprob', 'salign_coverage')


def mean_skipprob(tree, skip_probs):
    values = [prob for node, prob in skip_probs.items()
              if isinstance(node, Activity)]
    return sum(values) / len(values) if values else 0.0


def compute_metrics(log, tree, slpn_path, ppt_weights=None, return_dv=False):
    """
    Run the skip-alignment pipeline for (log, tree) and return the
    weight-coverage and skipprob summary metrics.

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
    twice. Defaults to False so existing callers are unaffected.
    """
    Path(slpn_path).parent.mkdir(parents=True, exist_ok=True)
    dv = pvoid.skipprob(log, tree, slpn_path, ppt_weights=ppt_weights)
    slpn = slpn_importer.read_slpn(slpn_path)
    transfer_pt_weights(tree, slpn)
    values = (
        mass_by_weight(tree, dv.skip_probs),
        voidage_by_weight(tree, dv.skip_probs),
        mean_skipprob(tree, dv.skip_probs),
        coverage_by_alignment(tree, dv),
    )
    metrics = dict(zip(METRIC_KEYS, values))
    return (metrics, dv) if return_dv else metrics
