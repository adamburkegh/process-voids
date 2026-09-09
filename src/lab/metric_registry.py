"""
Single source of truth for every analytical metric id this package's
experiment scripts (lab.exp_disco_degrade, lab.exp_claims_degrade,
lab.exp_surprise) can emit into a results CSV.

Deliberately scoped to analytical quantities only - not administrative/
bookkeeping columns each script also writes (log, combo, degradation_dim,
degradation_level, status, dropped_count, elapsed_s, node_id, node_type,
alphabet, target, n_drop_cases, n_events, and the surprise diagnostics
ambiguous_event_count/unattributable_event_count/out_of_alphabet_event_count).

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
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
    ),
    'weight_voidage': Metric(
        id='weight_voidage',
        description='1 - weight_coverage at every node (averaging distributes '
                    'linearly over the complement - a convenience read-out, '
                    'not a separately-derived quantity).',
        source='process_voids.coveragemass.voidage_by_weight',
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
    ),
    'skipprob': Metric(
        id='skipprob',
        description='Mean skip probability across Activity leaves - '
                    "skip-alignments' own dv.skip_probs output, unmodified "
                    'beyond averaging.',
        source='lab.metrics.mean_skipprob',
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
    ),
    'salign_coverage': Metric(
        id='salign_coverage',
        description="Coverage by alignment correspondence over skip-alignments' "
                    'lumped normal form (zero convention: a unit with no valid '
                    'execution contributes 0). Genuinely alignment-machinery-'
                    'based, unlike weight_coverage.',
        source='process_voids.coveragemass.coverage_by_alignment',
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
    ),
    'mandatory_node_count': Metric(
        id='mandatory_node_count',
        description='Number of non-Tau nodes in the scored tree with no silent '
                    'alternative (see coveragemass.py\'s Mandatory Node Count '
                    'section) - a diagnostic, not a void/coverage metric. '
                    'Inductive Miner at noise_threshold=0.0 wraps nearly every '
                    'leaf in Xor(Tau, activity), making this near-zero and every '
                    'void metric uninformative by construction on such a tree. '
                    'Scored at the tree root in exp_disco_degrade, at the target '
                    'node in exp_claims_degrade - same scoping as every other '
                    'per-node metric in each script.',
        source='process_voids.coveragemass.mandatory_node_count',
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
    ),
    'total_node_count': Metric(
        id='total_node_count',
        description='Number of non-Tau nodes in the scored tree - the '
                    'denominator for reading mandatory_node_count as a fraction '
                    'rather than a bare count only meaningful relative to a '
                    'specific tree\'s size. Same scoping as mandatory_node_count.',
        source='process_voids.coveragemass.total_node_count',
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
    ),
    'voidmass_deficit': Metric(
        id='voidmass_deficit',
        description='Pooled voidmass deficit at the scored node, from '
                    'classical (non-lumped) Petri-net alignments - the '
                    "formally-correct replacement for skip-alignments' lumped "
                    'normal form.',
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
    ),
    'voidmass_movecount': Metric(
        id='voidmass_movecount',
        description='Pooled non-silent move count at the scored node - the '
                    'denominator of voidmass_subprocess.',
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
    ),
    'voidmass_subprocess': Metric(
        id='voidmass_subprocess',
        description='voidmass_deficit / voidmass_movecount at the scored node - '
                    "voidmass as a fraction of that node's own moves.",
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
    ),
    'voidmass_process': Metric(
        id='voidmass_process',
        description='1 - voidmass_subprocess (alignment_mass_pooled) at the '
                    'scored node.',
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
    ),
    'alignment_coverage_pn': Metric(
        id='alignment_coverage_pn',
        description='(1 - skip_prob) * voidmass_process at the scored node - '
                    'the classical-alignment analogue of salign_coverage '
                    "(paper's \\covermove), reusing skip-alignments' own "
                    'skip_probs unchanged rather than deriving a separate '
                    'estimate.',
        source='process_voids.voidmass_pn.coverage_by_alignment_pn',
        scripts=('exp_disco_degrade', 'exp_claims_degrade'),
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
