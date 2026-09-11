"""
Append-only historical dictionary of every CSV column id this package's
experiment scripts (lab.exp_disco_degrade, lab.exp_surprise,
lab.exp_voidmass) have ever emitted, plus product-only ids computed by
process_voids for pvoid's own output. Entries are never deleted, only
appended or given status='retired' - a result CSV from any point in
this project's history should still have every one of its columns
findable here. The claims fixture runs through exp_disco_degrade
(lab.claims_fixture's CLAIMS_COMBOS/CLAIMS_DEGRADATIONS), so its
columns are exp_disco_degrade's, with no entries of their own.

Each Metric carries a status ('live': currently emitted; 'evaluation':
currently emitted, but only to evaluate other metrics - oracle
baselines, diagnostics - not a metric in its own right; 'product-only':
computed by process_voids, never emitted by a lab script; 'retired': no
current script emits it), an optional superseded_by (the id that
replaced it, for a clean rename/merge), and a history dict (commit or
version -> prior meaning) for an id whose CURRENT meaning changed
in-place without a rename - as opposed to superseded_by, which is for
an old id abandoned in favour of a new one.

Deliberately scoped to analytical quantities only - not administrative/
bookkeeping columns each script also writes (log, combo, degradation_dim,
degradation_level, status, dropped_count, elapsed_s, node_id, node_type,
alphabet, n_events, and the surprise diagnostics
ambiguous_event_count/unattributable_event_count/out_of_alphabet_event_count,
and exp_disco_degrade's timeout diagnostics timed_out_count/
timed_out_weight).

This is the executable source of truth, not documentation of it: every
live and evaluation id here is in one of the *_KEYS constants its
producing script actually builds its output dict from (in lab.metrics,
lab.exp_disco_degrade, lab.exp_surprise, lab.exp_voidmass and
process_voids.coveragemass) - see tests/lab/test_metric_registry.py,
which imports those constants directly and asserts this registry's ids
match them exactly. A metric renamed in its producing script without a
matching registry update fails that test.

Every exp_surprise metric is registered twice: an unsuffixed id ('self'
- the deployable mode, tail distribution estimated from the log being
scored) and a '_baseline'-suffixed id (scored against the *undegraded*
log's distribution instead - a benchmark-only oracle comparison,
isolating real degradation signal from the self-estimator degrading
along with the log; unavailable under real missing-activity scenarios,
where no undegraded reference exists). Self and baseline are different
quantities with different deployability - self is a dead end on its
own (all signal vanishes without a baseline to compare against),
baseline is a genuine standing benchmark metric - so they get separate
ids like every other metric here, not a shared one plus a modifier
column. See lab.exp_surprise's module docstring and _compute_cell.
"""

from dataclasses import dataclass, field


STATUSES = {'live', 'evaluation', 'product-only', 'retired'}


@dataclass(frozen=True)
class Metric:
    id: str
    description: str
    source: str    # module.function that computes it
    scripts: tuple  # lab.exp_* scripts that emit it - () for product-only ids
    status: str = 'live'          # one of STATUSES
    superseded_by: str = None     # id of the metric that replaces this one, if retired
    history: dict = field(default_factory=dict)  # commit/version -> prior meaning


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
                    "Not mean_leaf_skipprob (lab.metrics), a different, "
                    "unrelated statistic - see this id's history for result "
                    "CSVs where this column held that one instead.",
        source='dv.skip_probs (direct lookup, no computation of its own)',
        scripts=('exp_disco_degrade',),
        history={'f2aa44c': 'Blended mean of skip_probs[leaf] over every Activity '
                            'leaf in the tree (lab.metrics.mean_skipprob), not '
                            "dv.skip_probs[node] itself - the id was wired to the "
                            "wrong quantity. Renamed to mean_leaf_skipprob and "
                            "skipprob repointed to the correct lookup."},
    ),
    'mean_leaf_skipprob': Metric(
        id='mean_leaf_skipprob',
        description='Mean of skip_probs[leaf] over EVERY Activity leaf in the '
                    "whole tree, regardless of which node is passed in (see that "
                    "function's own docstring on the tree argument it ignores for "
                    "scoping) - NOT skip-alignments' own \"skipprob\" (see that "
                    "id's entry and history). A separate, home-grown blended "
                    'statistic, kept under an honest name rather than dropped '
                    'outright since '
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
                    "pooled (see voidmass_subprocess_lower/_upper for the "
                    "pooled quantity, and this id's history). Reuses "
                    "skip-alignments' own skip_probs unchanged rather than "
                    "deriving a separate estimate. LOWER bound: a variant whose "
                    "alignment search timed out (no alignments at all) is "
                    "treated as contributing a ratio of 0 (as if "
                    "it matched nothing), the SMALLER of the two coverage "
                    "readings.",
        source='process_voids.voidmass_pn.coverage_by_alignment_pn',
        scripts=('exp_disco_degrade',),
        history={'ab77b03': 'Pooled deficit/movecount ratio (1 - skip_prob) * '
                            'voidmass_process at the scored node, not the '
                            'per-execution average the definition specifies - '
                            'computed from voidmass_process by mistake, since '
                            'that table was already built.'},
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
        history={'e255c45': "Shared a single id plus a 'distribution' column "
                            "with the baseline comparison, rather than a "
                            "separate id per distribution."},
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
        history={'e255c45': "Shared a single id plus a 'distribution' column "
                            "with the baseline comparison, rather than a "
                            "separate id per distribution."},
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
        history={'e255c45': "Shared a single id plus a 'distribution' column "
                            "with the baseline comparison, rather than a "
                            "separate id per distribution."},
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
        history={'e255c45': "Shared a single id plus a 'distribution' column "
                            "with the baseline comparison, rather than a "
                            "separate id per distribution."},
    ),
    'bits_per_event_baseline': Metric(
        id='bits_per_event_baseline',
        description='bits_per_event scored against the undegraded log\'s '
                    'interval distribution - a benchmark-only oracle '
                    'comparison.',
        source='lab.exp_surprise._compute_variant',
        scripts=('exp_surprise',),
    ),

    # Retired ids - no current script emits these. Kept so a reader of an
    # old result CSV, or a git-blame trail through this file, can still
    # find out what the column meant.
    'voidmass_deficit': Metric(
        id='voidmass_deficit',
        description='Pooled voidmass deficit at the scored node, from classical '
                    "(non-lumped) Petri-net alignments - exp_disco_degrade's "
                    "predecessor to voidmass_deficit_lower/voidmass_deficit_upper, "
                    "before a timed-out variant's bound was split into a "
                    "lower/upper pair.",
        source='process_voids.voidmass_pn.voidmass_table_pn',
        scripts=('exp_disco_degrade',),
        status='retired',
        history={'4116fa3': 'Retired and split into voidmass_deficit_lower/'
                            'voidmass_deficit_upper, to give a timed-out '
                            "variant's deficit an honest bracket instead of "
                            'the ZeroDivisionError this unsplit id hit.'},
    ),
    'node_skip_prob': Metric(
        id='node_skip_prob',
        description='Per-node dv.skip_probs[node] in the per-node CSV, kept '
                    "under a distinct id from the root CSV's skipprob "
                    "(itself wired to the wrong quantity at the time) so the "
                    "per-node column wouldn't silently inherit that mistake.",
        source='dv.skip_probs (direct lookup, no computation of its own)',
        scripts=('exp_disco_degrade',),
        status='retired',
        superseded_by='skipprob',
        history={'f2aa44c': "Retired once skipprob itself was fixed to mean "
                            "dv.skip_probs[node] directly - node_skip_prob "
                            "and skipprob became the same id/computation at "
                            "the root and at every node alike, so the split "
                            "was no longer needed."},
    ),
    'alignment_coverage': Metric(
        id='alignment_coverage',
        description='Coverage by alignment correspondence over skip-alignments\' '
                    "lumped normal form - salign_coverage's name before the "
                    "classical-alignment replacement (alignment_coverage_pn) "
                    "claimed the unqualified name.",
        source='process_voids.coveragemass.coverage_by_alignment',
        scripts=('exp_disco_degrade',),
        status='retired',
        superseded_by='salign_coverage',
        history={'fa85701': 'Renamed to salign_coverage (skip-alignment) once '
                            'alignment_coverage_pn (classical-alignment) was '
                            'added, freeing the unqualified name for the '
                            'eventual classical-alignment replacement.'},
    ),

    # Product-only ids - computed by process_voids for pvoid's own output,
    # never emitted into a lab results CSV.
    'duration_coverage': Metric(
        id='duration_coverage',
        description='Coverage by aligned duration: (1 - skip_prob) weighted by '
                    "the scored subprocess's share of total real elapsed time, "
                    "rather than by move count - pvoid's own demo output "
                    '(show_tree_coverage_by_duration), not a lab metric.',
        source='process_voids.coveragemass.coverage_by_duration',
        scripts=(),
        status='product-only',
    ),

    # exp_voidmass.py's own ids - a target-subprocess dose-response sweep
    # using the SAME lumped skip-alignment machinery as salign_coverage/
    # voidsat (coveragemass.voidmass_table), not the classical Petri-net
    # alignments exp_disco_degrade uses. voidmass_subprocess/voidmass_process
    # here are lumped, unsuffixed, and unrelated to exp_disco_degrade's own
    # (_lower/_upper-suffixed) classical ids of the same short name -
    # see their history entries below.
    'skip_prob': Metric(
        id='skip_prob',
        description="dv.skip_probs[node] at the scored node - the same "
                    "underlying lookup as skipprob, under exp_voidmass.py's "
                    "own column name (with an underscore, unlike skipprob).",
        source='dv.skip_probs (direct lookup, no computation of its own)',
        scripts=('exp_voidmass',),
    ),
    'deficit': Metric(
        id='deficit',
        description="Pooled voidmass deficit at the scored node, from "
                    "skip-alignments' lumped normal form - the numerator "
                    "voidmass_subprocess/voidmass_process divide by their "
                    "respective denominators.",
        source='process_voids.coveragemass.voidmass_table',
        scripts=('exp_voidmass',),
    ),
    'movecount': Metric(
        id='movecount',
        description="Pooled non-silent move count at the scored node, from "
                    "skip-alignments' lumped normal form - "
                    "voidmass_subprocess's own denominator.",
        source='process_voids.coveragemass.voidmass_table',
        scripts=('exp_voidmass',),
    ),
    'voidmass_subprocess': Metric(
        id='voidmass_subprocess',
        description="Pooled voidmass deficit as a fraction of the scored "
                    "node's own move count, from skip-alignments' lumped "
                    'normal form (coveragemass.voidmass_table) - the target '
                    "subprocess's own share of missing mass under an "
                    'ablation sweep.',
        source='process_voids.coveragemass.voidmass_table',
        scripts=('exp_voidmass',),
        history={'v0.4.1': "exp_disco_degrade.py used this same unsuffixed id "
                           "for a DIFFERENT (classical Petri-net alignment) "
                           "quantity through v0.4.1 - the two scripts' CSVs "
                           "disagreed on what the column meant. Resolved in "
                           "v0.4.2 when exp_disco_degrade's own use split "
                           "into voidmass_subprocess_lower/_upper; "
                           "exp_voidmass.py's lumped use of the plain id is "
                           "unaffected and continues unchanged."},
    ),
    'voidmass_process': Metric(
        id='voidmass_process',
        description="Pooled voidmass deficit as a fraction of the ROOT's "
                    'move count, from skip-alignments\' lumped normal form '
                    '(coveragemass.voidmass_table) - see voidmass_subprocess.',
        source='process_voids.coveragemass.voidmass_table',
        scripts=('exp_voidmass',),
        history={'v0.4.1': "exp_disco_degrade.py used this same unsuffixed id "
                           "for a DIFFERENT (classical Petri-net alignment) "
                           "quantity through v0.4.1 - see voidmass_subprocess's "
                           "history entry."},
    ),
    'voidage_subprocess': Metric(
        id='voidage_subprocess',
        description='voidmass_subprocess * skip_prob(node) - see '
                    'coveragemass.voidmass_table.',
        source='process_voids.coveragemass.voidmass_table',
        scripts=('exp_voidmass',),
    ),
    'voidage_process': Metric(
        id='voidage_process',
        description='voidmass_process * skip_prob(node) - see '
                    'coveragemass.voidmass_table.',
        source='process_voids.coveragemass.voidmass_table',
        scripts=('exp_voidmass',),
    ),
    'target_voidmass_subprocess': Metric(
        id='target_voidmass_subprocess',
        description="Summary-row copy of the ablation target node's own "
                    'voidmass_subprocess, for the dose-response curve.',
        source='process_voids.coveragemass.voidmass_table',
        scripts=('exp_voidmass',),
    ),
    'target_voidmass_process': Metric(
        id='target_voidmass_process',
        description="Summary-row copy of the ablation target node's own "
                    'voidmass_process, for the dose-response curve.',
        source='process_voids.coveragemass.voidmass_table',
        scripts=('exp_voidmass',),
    ),
    'target_voidage_subprocess': Metric(
        id='target_voidage_subprocess',
        description="Summary-row copy of the ablation target node's own "
                    'voidage_subprocess, for the dose-response curve.',
        source='process_voids.coveragemass.voidmass_table',
        scripts=('exp_voidmass',),
    ),
    'target_voidage_process': Metric(
        id='target_voidage_process',
        description="Summary-row copy of the ablation target node's own "
                    'voidage_process, for the dose-response curve.',
        source='process_voids.coveragemass.voidmass_table',
        scripts=('exp_voidmass',),
    ),
    'target_rank_voidmass_process': Metric(
        id='target_rank_voidmass_process',
        description='1-indexed rank of the target node among all nodes in '
                    'the tree, sorted by voidmass_process descending - is '
                    'the ablated subprocess actually the biggest void?',
        source='lab.exp_voidmass._rank_descending',
        scripts=('exp_voidmass',),
    ),
    'target_rank_voidage_process': Metric(
        id='target_rank_voidage_process',
        description='Same as target_rank_voidmass_process, ranked by '
                    'voidage_process instead.',
        source='lab.exp_voidmass._rank_descending',
        scripts=('exp_voidmass',),
    ),
    'n_optimal_alignments': Metric(
        id='n_optimal_alignments',
        description='Total optimal alignments found across all trace '
                    'variants (the sum of |Gamma_sigma|) - a '
                    'diagnostic on the alignment search itself, not a '
                    'void/coverage metric.',
        source='lab.exp_voidmass._n_optimal_alignments',
        scripts=('exp_voidmass',),
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
