"""
Single source of truth for every analytical metric id this package's
experiment scripts (lab.exp_disco_degrade, lab.exp_surprise) can emit
into a results CSV. exp_claims_degrade.py was retired once
lab.claims_fixture's CLAIMS_COMBOS/CLAIMS_DEGRADATIONS let the claims
fixture run through exp_disco_degrade directly - same schema, no
separate bespoke script or registry entries needed for it.

Deliberately scoped to analytical quantities only - not administrative/
bookkeeping columns each script also writes (log, combo, degradation_dim,
degradation_level, status, dropped_count, elapsed_s, node_id, node_type,
alphabet, n_events, and the surprise diagnostics
ambiguous_event_count/unattributable_event_count/out_of_alphabet_event_count,
and exp_disco_degrade's timeout diagnostics timed_out_count/
timed_out_weight).

This is the executable source of truth, not documentation of it: every
id here is one of the *_KEYS constants each producing script actually
builds its output dict from (METRIC_KEYS, CLASSICAL_METRIC_KEYS,
NODE_METRIC_KEYS, SUMMARY_METRIC_KEYS, TREE_METRIC_KEYS) - see tests/lab/
test_metric_registry.py, which imports those constants directly and
asserts this registry's ids match them exactly. A metric renamed in its
producing script without a matching registry update fails that test,
the same silent-drift failure mode that left plot_dose_response's
METRICS list stale earlier in this project's history.

Every exp_surprise metric is registered twice: an unsuffixed id ('self'
- the deployable mode, tail distribution estimated from the log being
scored) and a '_baseline'-suffixed id (scored against the *undegraded*
log's distribution instead - a benchmark-only oracle comparison,
isolating real degradation signal from the self-estimator degrading
along with the log; unavailable under real missing-activity scenarios,
where no undegraded reference exists). These were originally a single
id plus a 'distribution' column, but self and baseline are different
quantities with different deployability - self is a dead end on its
own (all signal vanishes without a baseline to compare against, hence
this split), baseline is a genuine standing benchmark metric - so they
get separate ids like every other metric here, not a shared one plus a
modifier column. See lab.exp_surprise's module docstring and
_compute_cell.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Metric:
    id: str
    description: str
    source: str    # module.function that computes it
    scripts: tuple  # lab.exp_* scripts that emit it


METRICS = {
    'weight_coverage': Metric(
        id='weight_coverage',
        description="Structural weighted average of (1 - skip_prob) over the "
                    "tree, driven by discovered occurrence-frequency weights "
                    "and skip_probs at the leaves. Not built on skip-alignments' "
                    "alignment machinery - only borrows skip_probs as a leaf "
                    "input, the same as skipprob does directly.",
        source='process_voids.coveragemass.mass_by_weight',
        scripts=('exp_disco_degrade',),
    ),
    'weight_voidage': Metric(
        id='weight_voidage',
        description='1 - weight_coverage at every node (averaging distributes '
                    'linearly over the complement - a convenience read-out, '
                    'not a separately-derived quantity).',
        source='process_voids.coveragemass.voidage_by_weight',
        scripts=('exp_disco_degrade',),
    ),
    'skipprob': Metric(
        id='skipprob',
        description="Skip-alignments' own published skip probability - "
                    'dv.skip_probs[node], unmodified, evaluated at whichever node '
                    'is being scored (the root in the root-level CSV, any node in '
                    "the per-node CSV - same id, same computation, both places). "
                    "Until a naming fix, this id was wired to mean_leaf_skipprob "
                    "(lab.metrics) by mistake - a different, unrelated statistic "
                    "that happened to get this term's name. A result CSV written "
                    "before that fix has the wrong quantity under this column.",
        source='dv.skip_probs (direct lookup, no computation of its own)',
        scripts=('exp_disco_degrade',),
    ),
    'mean_leaf_skipprob': Metric(
        id='mean_leaf_skipprob',
        description='Mean of skip_probs[leaf] over EVERY Activity leaf in the '
                    "whole tree, regardless of which node is passed in (see that "
                    "function's own docstring on the tree argument it ignores for "
                    "scoping) - NOT skip-alignments' own \"skipprob\" (see that "
                    'id\'s entry - this used to be wired to the \'skipprob\' name '
                    'by mistake). A separate, home-grown blended statistic, kept '
                    'under an honest name rather than dropped outright since '
                    "it's not yet established whether it's actually informative. "
                    'Root-level only - not meaningful per-node, since it always '
                    'returns the same whole-tree average regardless of the node '
                    'argument.',
        source='lab.metrics.mean_leaf_skipprob',
        scripts=('exp_disco_degrade',),
    ),
    'salign_coverage': Metric(
        id='salign_coverage',
        description="Coverage by alignment correspondence over skip-alignments' "
                    'lumped normal form (zero convention: a unit with no valid '
                    'execution contributes 0). Genuinely alignment-machinery-'
                    'based, unlike weight_coverage.',
        source='process_voids.coveragemass.coverage_by_alignment',
        scripts=('exp_disco_degrade',),
    ),
    'mandatory_node_count': Metric(
        id='mandatory_node_count',
        description='Number of non-Tau nodes in the scored tree with no silent '
                    'alternative (see coveragemass.py\'s Mandatory Node Count '
                    'section) - a diagnostic, not a void/coverage metric. '
                    'Inductive Miner at noise_threshold=0.0 wraps nearly every '
                    'leaf in Xor(Tau, activity), making this near-zero and every '
                    'void metric uninformative by construction on such a tree. '
                    'Scored at the tree root in the root-level CSV, and at every '
                    'node (including claims-fixture ablation targets like '
                    'appeal_seq) in the per-node CSV.',
        source='process_voids.coveragemass.mandatory_node_count',
        scripts=('exp_disco_degrade',),
    ),
    'total_node_count': Metric(
        id='total_node_count',
        description='Number of non-Tau nodes in the scored tree - the '
                    'denominator for reading mandatory_node_count as a fraction '
                    'rather than a bare count only meaningful relative to a '
                    'specific tree\'s size. Same scoping as mandatory_node_count.',
        source='process_voids.coveragemass.total_node_count',
        scripts=('exp_disco_degrade',),
    ),
    'voidmass_deficit_lower': Metric(
        id='voidmass_deficit_lower',
        description='Pooled voidmass deficit at the scored node, from classical '
                    "(non-lumped) Petri-net alignments - the formally-correct "
                    "replacement for skip-alignments' lumped normal form. LOWER "
                    "bound: a variant whose align_variant_all search timed out "
                    "to zero alignments (a real, observed failure mode, not "
                    "hypothetical - see voidmass_table_pn) is credited 0 deficit "
                    "there, as if it fit perfectly. Equal to voidmass_deficit_"
                    "upper whenever no variant times out - the two only diverge "
                    "on a cell that actually hit this case, making the gap "
                    "itself a visible signal rather than a hidden assumption.",
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade',),
    ),
    'voidmass_deficit_upper': Metric(
        id='voidmass_deficit_upper',
        description="Same as voidmass_deficit_lower, but a timed-out variant is "
                    "credited its largest possible deficit instead: w * X_max, "
                    "X_max = voidmass_pn.timed_out_movecount_bound (2|sigma| + "
                    "the model's cheapest complete path length) - every move of "
                    "the longest path an optimal alignment could take counted "
                    "as a model move. Derived from the aligner's cost model; "
                    "see that function for the proof.",
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade',),
    ),
    'voidmass_movecount': Metric(
        id='voidmass_movecount',
        description='Pooled non-silent move count at the scored node, OBSERVED '
                    'from completed variants only - a measurement, not a '
                    'substitution. A timed-out variant contributes nothing '
                    'here; see voidmass_movecount_bound for the denominator the '
                    'lower/upper bounds divide by.',
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade',),
    ),
    'voidmass_movecount_bound': Metric(
        id='voidmass_movecount_bound',
        description="voidmass_movecount plus w * X_max for every timed-out "
                    "variant (voidmass_pn.timed_out_movecount_bound) - the "
                    "shared denominator of voidmass_subprocess_lower/upper, "
                    "reported so the bounds' arithmetic is reproducible. Its "
                    "own column rather than overloading voidmass_movecount, "
                    "whose meaning would otherwise depend on whether a timeout "
                    "happened. Equal to voidmass_movecount when nothing timed "
                    "out.",
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade',),
    ),
    'voidmass_subprocess_lower': Metric(
        id='voidmass_subprocess_lower',
        description='voidmass_deficit_lower / voidmass_movecount_bound at the '
                    "scored node - voidmass as a fraction of that node's own "
                    'moves. A valid lower bound on the no-timeout value.',
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade',),
    ),
    'voidmass_subprocess_upper': Metric(
        id='voidmass_subprocess_upper',
        description='voidmass_deficit_upper / voidmass_movecount_bound at the '
                    "scored node - voidmass as a fraction of that node's own "
                    'moves. A valid upper bound on the no-timeout value.',
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade',),
    ),
    'voidmass_process_lower': Metric(
        id='voidmass_process_lower',
        description="voidmass_deficit_lower / the ROOT's voidmass_movecount_bound "
                    '- missing moves under the scored node as a fraction of the '
                    "whole model's moves (additive over any cut through the "
                    'tree). A valid lower bound on the no-timeout value.',
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade',),
    ),
    'voidmass_process_upper': Metric(
        id='voidmass_process_upper',
        description="voidmass_deficit_upper / the ROOT's voidmass_movecount_bound "
                    '- see voidmass_process_lower. A valid upper bound on the '
                    'no-timeout value: X_max is the same at every node, so '
                    'numerator and denominator share it.',
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade',),
    ),
    'alignment_coverage_pn_lower': Metric(
        id='alignment_coverage_pn_lower',
        description="\\covermove (defn:move-coverage), the classical-alignment "
                    "analogue of salign_coverage: (1 - skip_prob) * a per-"
                    "execution match/movecount ratio, averaged uniformly across "
                    "executions within an alignment and across an alignment's "
                    "tied alternatives, then across a variant's own weight - "
                    "verified term-by-term against the formal definition, not "
                    "pooled (see voidmass_process/voidmass_subprocess for the "
                    "pooled quantities, which this id used to be computed from "
                    "by mistake - convenient since that table was already "
                    "built, but not what the definition specifies). Reuses "
                    "skip-alignments' own skip_probs unchanged rather than "
                    "deriving a separate estimate. LOWER bound: a variant whose "
                    "alignment search timed out (no alignments at all - "
                    "previously silently excluded from the weighted average "
                    "entirely) is treated as contributing a ratio of 0 (as if "
                    "it matched nothing), the SMALLER of the two coverage "
                    "readings.",
        source='process_voids.voidmass_pn.coverage_by_alignment_pn',
        scripts=('exp_disco_degrade',),
    ),
    'alignment_coverage_pn_upper': Metric(
        id='alignment_coverage_pn_upper',
        description="Same as alignment_coverage_pn_lower, but a timed-out "
                    "variant is treated as contributing a ratio of 1 (as if it "
                    "matched perfectly) instead of being excluded - the LARGER "
                    "of the two coverage readings. Equal to alignment_coverage_"
                    "pn_lower whenever no variant times out.",
        source='process_voids.voidmass_pn.coverage_by_alignment_pn',
        scripts=('exp_disco_degrade',),
    ),
    'voidsat': Metric(
        id='voidsat',
        description="\\voidsat (defn:aligned-duration), voidat computed over "
                    "skip-alignments' own lumped optimal alignments: "
                    "skip_prob * a real-elapsed-time mass estimate (Definition "
                    "[Move Durations]/[Aligned Duration Mass]), averaged over "
                    "every actual trace instance (not deduplicated variants - "
                    "duration is per-instance) and every tied optimal "
                    "alignment. A lumped skip move over an entirely-missing "
                    "subprocess shares its time gap with whatever real event "
                    "immediately follows it, rather than claiming the whole "
                    "gap - see coveragemass.py's own Coverage By Aligned "
                    "Duration section for the worked example this was "
                    "verified against.",
        source='process_voids.coveragemass.voidsat',
        scripts=('exp_disco_degrade',),
    ),
    'containment_bits': Metric(
        id='containment_bits',
        description='Interval-surprise bits charged to a tree node under '
                    'containment attribution: every node whose alphabet '
                    'contains the event activity (leaf and all ancestors). '
                    'Puts the signal on whatever ran late, which for a '
                    "missing subprocess is the wrong node. Self-distribution "
                    "(tail estimated from the log being scored) - see "
                    "containment_bits_baseline for the oracle-comparison "
                    'counterpart.',
        source='process_voids.surprise.surprise_totals',
        scripts=('exp_surprise',),
    ),
    'containment_bits_baseline': Metric(
        id='containment_bits_baseline',
        description='containment_bits scored against the undegraded log\'s '
                    'interval distribution instead of the (possibly degraded) '
                    'scored log\'s own - a benchmark-only oracle comparison, '
                    'unavailable under real missing-activity scenarios.',
        source='process_voids.surprise.surprise_totals',
        scripts=('exp_surprise',),
    ),
    'predecessor_bits': Metric(
        id='predecessor_bits',
        description='Interval-surprise bits charged to a tree node under '
                    'process-tree-aware predecessor attribution: the '
                    'subprocess the model says should run immediately before '
                    'the event, then its ancestors. A structural '
                    'approximation of true (alignment-resolved) predecessor '
                    'attribution - Xor/And exits are genuinely ambiguous and '
                    'get charged to the composite node rather than a branch. '
                    'Self-distribution - see predecessor_bits_baseline for '
                    'the oracle-comparison counterpart.',
        source='process_voids.surprise.predecessor_totals',
        scripts=('exp_surprise',),
    ),
    'predecessor_bits_baseline': Metric(
        id='predecessor_bits_baseline',
        description='predecessor_bits scored against the undegraded log\'s '
                    'interval distribution instead of the scored log\'s own - '
                    'a benchmark-only oracle comparison, unavailable under '
                    'real missing-activity scenarios.',
        source='process_voids.surprise.predecessor_totals',
        scripts=('exp_surprise',),
    ),
    'headline_bits': Metric(
        id='headline_bits',
        description="Root-level total interval-surprise bits (containment's "
                    "total at the tree root) - the whole log's surprise, "
                    'self-distribution.',
        source='lab.exp_surprise._compute_variant',
        scripts=('exp_surprise',),
    ),
    'headline_bits_baseline': Metric(
        id='headline_bits_baseline',
        description='headline_bits scored against the undegraded log\'s '
                    'interval distribution - a benchmark-only oracle '
                    'comparison.',
        source='lab.exp_surprise._compute_variant',
        scripts=('exp_surprise',),
    ),
    'bits_per_event': Metric(
        id='bits_per_event',
        description='headline_bits / n_events - mean surprise per event, '
                    'self-distribution.',
        source='lab.exp_surprise._compute_variant',
        scripts=('exp_surprise',),
    ),
    'bits_per_event_baseline': Metric(
        id='bits_per_event_baseline',
        description='bits_per_event scored against the undegraded log\'s '
                    'interval distribution - a benchmark-only oracle '
                    'comparison.',
        source='lab.exp_surprise._compute_variant',
        scripts=('exp_surprise',),
    ),
}


def format_registry():
    """
    Plain-text cross-reference table (id, source function, producing
    scripts, description) - not a substitute for a proper table with
    paper symbols/definition refs, just a quick sanity-check dump of
    what this registry currently holds and where each entry comes from.
    """
    lines = []
    for metric_id in sorted(METRICS):
        metric = METRICS[metric_id]
        lines.append(metric_id)
        lines.append(f'    source:      {metric.source}')
        lines.append(f'    scripts:     {", ".join(metric.scripts)}')
        lines.append(f'    description: {metric.description}')
    return '\n'.join(lines)


if __name__ == '__main__':
    print(format_registry())
