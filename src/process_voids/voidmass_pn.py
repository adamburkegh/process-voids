'''
Prototype: voidmass deficit computed on CLASSICAL (standard) Petri-net
alignments, via skipalignments.alignall.align_pn_all, rather than
skip-alignments' Aligner.align2/executions().

Why this exists: coveragemass.executions() is built on skip-alignments'
normal form, which lumps an entirely-unwitnessed subtree into a single
Skip/TauPath move on its coarsest ancestor. The voidmass brief's
deficit is defined to count every constituent activity individually -
skipprob is the one allowed to lift to the highest block, not deficit -
so building deficit on the lumped representation is wrong. Classical
alignments have no Skip(subtree) construct at all: every missing leaf
is necessarily its own model-move, so no lumping can occur.

This module is a prototype (see tests/process_voids/test_voidmass_pn_prototype.py).
voidmass_deficit/subprocess/process (terms_by_node, below) don't need
the Definition [Executions] grouping (loop iterations collapsing into
one traversal) - Bär et al.'s conditions (2)/(3) exist to keep one
traversal's moves separate from another's, which only matters when you
AVERAGE per traversal, and summing over the whole alignment gives the
same total regardless of how moves are partitioned into traversals,
since it's the same set of moves either way. coverage_by_alignment_pn
(\\covermove, defn:move-coverage) DOES average per traversal, so it
DOES need that grouping - see _to_alignment_mass_path, which translates
a classical alignment into coveragemass.alignment_mass's own path
shape and reuses that machinery (executions()/matchcount/movecount/the
averaging structure itself) entirely unmodified, verified term-by-term
against the formal definition. An earlier version of this function
reused voidmass_table_pn's pooled sums instead of building this
translation - convenient, but not what the definition specifies; see
that function's own docstring for where the pooled quantity still
lives, honestly labelled.

If this approach is adopted, this needs folding into coveragemass.py;
if not, delete both files.

Requires skip-alignments' tau_ids fix AND the id_loop_list/cycle-guard
fix (build_petri_net's 6-tuple return - net, im, fm, activity_to_id,
tau_ids, id_loop_list; align_pn_all's tau_ids parameter) -
process-voids is pinned to a local editable install of the sibling
skip-alignments repo for this (see pyproject.toml) until these fixes
are released. The cycle-guard fix matters beyond correctness: id_loop_list
was always [] before (a genuine gap, not just an unused default) - the
A* search couldn't cycle-detect on any tau-skippable loop, which is the
suspected cause of the near-timeout clustering found on rtfm (42/48
variants at ~100-104s regardless of complexity - see session notes).
Passing the real id_loop_list from build_id_net is expected to fix that,
not just theoretically close the gap.
'''

import logging
import time
from types import SimpleNamespace

from skipalignments import Skip, TauPath
from skipalignments.probabilities import EbiOccurance
from skipalignments.alignall import align_pn_all

from process_voids.coveragemass import alignment_mass, _variant_key, min_activity_count_by_node

logger = logging.getLogger(__name__)


def build_id_net(tree):
    '''
    (net, im, fm, activity_to_id, tau_ids, id_loop_list) for tree, where
    net's transition labels are internal ids, not activity names -
    align_variant() handles the renaming this requires on the trace
    side. tau_ids is the set of ids that are genuine Tau leaves (e.g.
    the tau branch of an Xor(activity, Tau)) - align_variant() passes it
    straight to align_pn_all so those branches are correctly zero-cost
    rather than priced as a deviation. id_loop_list is the set of loop
    ids the A* search's cycle guard needs to actually engage on a tau-
    skippable loop - MUST be threaded through to align_pn_all/
    align_variant_all/voidmass_table_pn, not left as [] (their default),
    or the guard silently does nothing (see module docstring).
    '''
    return EbiOccurance().build_petri_net(tree)


def align_variant_all(activities, net, im, fm, activity_to_id, tau_id_set, id_loop_list=None, timeout=100):
    '''
    Every tied optimal classical alignment for the trace `activities`
    against net/im/fm. Returns a list of move lists (one per tied
    optimal alignment), each a list of pm4py PetriNet.Transition with
    .label a (trace_side, model_side) pair using '>>' for the skipped side.
    '''
    if id_loop_list is None:
        id_loop_list = []
    renamed = [activity_to_id.get(a, a) for a in activities]
    result = align_pn_all(renamed, net, im, fm, id_loop_list, timeout=timeout, tau_ids=tau_id_set)
    _time_end, (opt_agns, _code, _first_time) = result[0]
    return [agn['alignment'] for agn in opt_agns]


def align_variant(activities, net, im, fm, activity_to_id, tau_id_set, id_loop_list=None, timeout=100):
    '''The first tied optimal alignment only - see align_variant_all.'''
    return align_variant_all(activities, net, im, fm, activity_to_id, tau_id_set,
                              id_loop_list=id_loop_list, timeout=timeout)[0]


def classify_move(t, id_to_activity, tau_id_set):
    '''
    (kind, activity) for one alignment move, kind in 'sync'/'log'/'model'/'tau'.
    activity is None for 'log' (nothing on the model side) and for 'tau'
    (a silent/helper transition, or a genuine Tau leaf in tau_id_set -
    neither is a deviation, so neither counts towards deficit).
    '''
    trace_side, model_side = t.label
    if model_side == '>>':
        return 'log', None
    if model_side is None or model_side in tau_id_set:
        return 'tau', None
    activity = id_to_activity.get(model_side)
    if activity is None:
        return 'tau', None
    if trace_side == '>>':
        return 'model', activity
    return 'sync', activity


def deficit_by_node(alignment, tree, activity_to_id, tau_id_set):
    '''
    {node: deficit} for every node in tree, where deficit(node) is the
    count of 'model' moves in `alignment` whose activity is one of
    node's leaf labels. deficit(execution) = movecount - matchcount
    simplifies to exactly this count, since movecount excludes tau
    moves and matchcount is the sync-move count: (sync + model) - sync.
    '''
    return {node: d for node, (d, _m) in
            terms_by_node(alignment, tree, activity_to_id, tau_id_set).items()}


def terms_by_node(alignment, tree, activity_to_id, tau_id_set):
    '''
    {node: (deficit, movecount)} for every node in tree, over one
    alignment - deficit = count of 'model' moves under node's leaves,
    movecount = count of 'sync'+'model' moves under node's leaves (tau
    moves excluded from both, matching coveragemass.movecount/matchcount).
    '''
    id_to_activity = {v: k for k, v in activity_to_id.items()}
    model_activities = []
    sync_activities = []
    for t in alignment:
        kind, activity = classify_move(t, id_to_activity, tau_id_set)
        if kind == 'model':
            model_activities.append(activity)
        elif kind == 'sync':
            sync_activities.append(activity)

    table = {}

    def _walk(node):
        leaves = set(node.get_leaf_labels())
        d = sum(1 for a in model_activities if a in leaves)
        m = d + sum(1 for a in sync_activities if a in leaves)
        table[node] = (d, m)
        for child in node.children:
            _walk(child)

    _walk(tree)
    return table


def _leaf_nodes_by_name(tree):
    '''{leaf.name: leaf} for every Activity/Tau leaf in tree - the
    name->node resolution _to_alignment_mass_path needs to wrap a
    classical alignment move's activity name back into the real tree
    node coveragemass.py's Skip/TauPath wrappers expect. Same
    duplicate-label caveat as terms_by_node's leaf-label-set membership
    check elsewhere in this module: not solved, just not any worse here
    than it already is throughout this codebase.'''
    nodes = {}

    def _walk(node):
        if not node.children:
            nodes[node.name] = node
        for child in node.children:
            _walk(child)

    _walk(tree)
    return nodes


def _to_alignment_mass_path(alignment, id_to_activity, nodes_by_name, tau_id_set):
    '''
    Translates one classical alignment (list of pm4py Transitions,
    .label = (trace_side, model_side)) into the (log_elem, model_elem)
    path shape coveragemass.executions()/alignment_mass already expect -
    reusing that machinery entirely unchanged for classical alignments,
    rather than reimplementing execution-grouping/averaging for a
    second move representation. See defn:move-coverage - this module's
    docstring on why classical (non-lumped) alignments matter for the
    formal definition, and coveragemass.alignment_mass's own docstring
    for the wrapper convention being reproduced here:
      - a bare leaf node: a synchronous move (log event matched a real
        leaf), same as skip-alignments' own bare-leaf convention
      - Skip(leaf, leaf.skip_cost): a required activity present in the
        model but missing from the log (movecount-counted, non-silent)
      - TauPath(leaf): a silent move through a genuine Tau leaf
        (movecount-excluded) - resolved via tau_id_set + activity_to_id,
        which already maps genuine Tau leaves by name (confirmed against
        build_id_net's real output, not just its docstring)
      - '>>' : a pure log move, OR a structural/helper silent transition
        with no tree correspondence at all (build_petri_net inserts
        these for net routing - label=None, never in tau_id_set). Both
        cases are unattributable to any subprocess, so both become a gap
        in executions()' consecutive-run grouping - the same treatment
        an unrelated move already gets there, not a special case.

    id_to_activity/nodes_by_name are precomputed ONCE by the caller
    (voidmass_table_pn) and passed in, not rebuilt per call - this is
    called once per deduped alignment, and rebuilding a full tree walk
    (_leaf_nodes_by_name) every time measurably slowed a real sweep
    (rtfm: voidmass_table_pn 105s -> 307s) before this was hoisted out.
    '''
    path = []
    for t in alignment:
        trace_side, model_side = t.label
        if model_side == '>>':
            path.append((trace_side, '>>'))
            continue
        name = id_to_activity.get(model_side)
        node = nodes_by_name.get(name) if name is not None else None
        if node is None:
            path.append((trace_side, '>>'))
            continue
        if model_side in tau_id_set:
            path.append((trace_side, TauPath(node)))
        elif trace_side == '>>':
            path.append((trace_side, Skip(node, node.skip_cost)))
        else:
            path.append((trace_side, node))
    return path


def sum_safe_signature(alignment, id_to_activity, tau_id_set):
    '''
    Canonical key for "these two tied optimal alignments tell the same
    causal story": the ORDERED list of (kind, activity) for 'sync'/
    'model' moves only - 'log' and 'tau' moves are dropped entirely.

    Why this is safe for voidmass specifically, and NOT safe to reuse
    unmodified for anything else: align_pn_all's all-optimal search
    returns every total ordering of a move sequence, including orderings
    that differ ONLY in where two commuting no-ops fall relative to each
    other (a silent/tau move and an unrelated log move don't constrain
    one another, so the search enumerates both orderings as if they were
    distinct alignments - see session notes, claim_25 in the claims
    fixture: 9 raw alignments, 3 real causal stories, a 6:1:2 raw split
    instead of the correct 1:1:1). Dropping log/tau collapses exactly
    those spurious duplicates, because dropping them is exactly what
    voidmass needs anyway - deficit/movecount only ever look at sync and
    model moves (see terms_by_node), and their SUM is invariant to how
    those two move-kinds are interleaved with log/tau. A duration-
    weighted metric (coverage_by_alignment-style) would NOT be safe with
    this signature: a move's duration depends on which block/traversal
    it falls in, so its position relative to a log move is real
    information there, not a commuting no-op - that metric would need a
    finer signature that keeps log moves (though still not tau moves,
    which never carry duration).

    This is a cheaper, per-metric version of what Bär et al.'s normal
    form does at the alignment level: canonicalise away accidental
    enumeration order before counting/weighting, rather than let raw
    search output privilege whichever story has more free slots to
    permute.
    '''
    return tuple(
        (kind, activity) for t in alignment
        for kind, activity in [classify_move(t, id_to_activity, tau_id_set)]
        if kind in ('sync', 'model')
    )


# TODO(\covat): a duration-weighted metric on this classical-alignment path
# needs its OWN signature function, not this one. It must keep 'log' moves
# in the key (their position relative to sync/model moves determines which
# block/traversal a duration gets attributed to - see the docstring above)
# but can still drop 'tau' moves, since silent transitions never carry a
# duration. Get this wrong and durations will get attributed to the wrong
# traversal exactly the way deficit counts got attributed to the wrong
# activity here before the dedup fix - same failure mode, worse to debug
# because it'd show up as a plausible-looking wrong number, not a
# visible asymmetry between two leaves.


def dedupe_alignments(alignments, id_to_activity, tau_id_set, signature_fn=sum_safe_signature):
    '''
    One representative alignment per distinct signature_fn value, in
    first-seen order - collapses combinatorial reorderings of commuting
    moves that align_pn_all's all-optimal search returns as separate
    "optimal alignments" but which represent the same causal story.
    '''
    seen = {}
    for alignment in alignments:
        sig = signature_fn(alignment, id_to_activity, tau_id_set)
        seen.setdefault(sig, alignment)
    return list(seen.values())


def voidmass_table_pn(tree, variant_probs, net, im, fm, activity_to_id, tau_id_set,
                       id_loop_list=None, timeout=100):
    '''
    (table, skip_dict) - table is the classical-alignment analogue of
    coveragemass.voidmass_table: one pass over the tree, {node: {
    'deficit_lower', 'deficit_upper', 'movecount', 'voidmass_subprocess_
    lower', 'voidmass_subprocess_upper', 'voidmass_process_lower',
    'voidmass_process_upper'}}, aggregated (SUMMED, not averaged - see
    voidmass-brief.md) over every variant and every distinct causal
    story among that variant's tied optimal alignments (see dedupe_
    alignments/sum_safe_signature - raw alignment count is NOT used,
    since align_pn_all's all-optimal search can return the same causal
    story multiple times under different commuting-move orders).

    lower/upper: align_variant_all can legitimately return ZERO
    alignments for a variant (a per-variant timeout, not a hang - see
    the near_timeout logging below) - a real, observed failure mode
    (labnotes.md finding C), not a hypothetical one. There is no
    principled single point estimate for such a variant's deficit, so
    two conservative bounds are reported instead of guessing: LOWER
    treats it as a perfect fit (deficit 0), UPPER treats it as the
    worst possible fit (deficit = movecount, every expected move a
    model move - the metric can only be overstated, never understated).
    Both bounds use the SAME movecount contribution either way - the
    model's own minimum executable length for that node (coveragemass.
    min_activity_count), computable without any alignment succeeding -
    so only deficit, not movecount, needs two accumulators; the
    denominator's own meaning doesn't change between bounds, only how
    void the timed-out variant is assumed to be. A variant that never
    times out contributes identically to both bounds, so this is the
    same cost as a single accumulator except for the (hopefully rare)
    timed-out variants themselves, which need no alignment search at
    all for their synthetic contribution.

    skip_dict is the SAME deduped alignments, translated (via
    _to_alignment_mass_path) into coveragemass.alignment_mass's expected
    input shape - built in the same pass so the expensive
    align_variant_all call is only ever made once per variant, not once
    for the pooled sums here and again for coverage_by_alignment_pn's
    (correctly non-pooled - see defn:move-coverage) mass term. A timed-
    out variant still gets an entry (an empty list, not simply absent)
    so alignment_mass's own timed_out_ratio handling can recognise it.

    Weighting choice, named explicitly rather than left implicit:
    UNIFORM OVER DISTINCT CAUSAL SIGNATURES - a variant's probability is
    split equally across its deduped alignments, treating every
    surviving causal story as equally likely. This is a cheap
    approximation, not the ICPM-consistent answer: the "correct"
    stochastic weighting would split by each story's own path
    probability under the model's stochastic language (e.g. in the
    claims fixture, the tau branch alone carries prior weight 0.75, so
    a stochastic weighting would favour the "skip entirely" story far
    more than uniform-over-signatures does). Revisit if that gap turns
    out to matter empirically.
    '''
    deficit_sum_lower = {}
    deficit_sum_upper = {}
    movecount_sum = {}
    skip_dict = {}

    def _add(node, d_lower, d_upper, m, share):
        deficit_sum_lower[node] = deficit_sum_lower.get(node, 0.0) + share * d_lower
        deficit_sum_upper[node] = deficit_sum_upper.get(node, 0.0) + share * d_upper
        movecount_sum[node] = movecount_sum.get(node, 0.0) + share * m

    id_to_activity = {v: k for k, v in activity_to_id.items()}
    nodes_by_name = _leaf_nodes_by_name(tree)
    activity_counts = min_activity_count_by_node(tree)
    n_variants = len(variant_probs)
    total_started = time.monotonic()
    n_near_timeout = 0
    n_timed_out = 0
    for i, (variant, weight) in enumerate(variant_probs.items(), start=1):
        started = time.monotonic()
        raw_alignments = align_variant_all(list(variant), net, im, fm, activity_to_id,
                                            tau_id_set, id_loop_list=id_loop_list, timeout=timeout)
        elapsed = time.monotonic() - started
        alignments = dedupe_alignments(raw_alignments, id_to_activity, tau_id_set)
        # Flag anything running close to the per-variant timeout - the
        # likely signature of the id_loop_list gap (a tau-skippable loop
        # the A* search can't cycle-detect, so it burns the full budget
        # rather than resolving quickly or hanging outright) rather than
        # genuine alignment cost. See session notes, the rtfm baseline
        # run's ~75 minutes unaccounted for outside skip-alignments.
        near_timeout = elapsed > 0.9 * timeout
        n_near_timeout += near_timeout
        log = logger.warning if near_timeout else logger.debug
        log('variant %d/%d (%d activities, weight=%.4g): %.2fs, %d raw -> %d deduped '
            'alignments%s', i, n_variants, len(variant), weight, elapsed,
            len(raw_alignments), len(alignments),
            ' [near timeout - possible id_loop_list gap]' if near_timeout else '')
        if not alignments:
            n_timed_out += 1
            logger.warning('variant %d/%d: 0 alignments (timed out) - contributing the '
                            'conservative lower/upper deficit bound instead of crashing',
                            i, n_variants)
            for node, count in activity_counts.items():
                _add(node, 0.0, count, count, weight)
        else:
            share = weight / len(alignments)
            for alignment in alignments:
                for node, (d, m) in terms_by_node(alignment, tree, activity_to_id, tau_id_set).items():
                    _add(node, d, d, m, share)
        skip_dict[_variant_key(variant)] = [
            SimpleNamespace(path=_to_alignment_mass_path(alignment, id_to_activity,
                                                           nodes_by_name, tau_id_set))
            for alignment in alignments
        ]

    logger.info('voidmass_table_pn: %d variants in %.1fs (%d near-timeout, %d timed out)',
                n_variants, time.monotonic() - total_started, n_near_timeout, n_timed_out)

    root_movecount = movecount_sum.get(tree, 0.0)
    table = {}

    def _walk(node):
        d_lower = deficit_sum_lower.get(node, 0.0)
        d_upper = deficit_sum_upper.get(node, 0.0)
        m = movecount_sum.get(node, 0.0)
        v1_lower = d_lower / m if m else 0.0
        v1_upper = d_upper / m if m else 0.0
        table[node] = {
            'deficit_lower': d_lower,
            'deficit_upper': d_upper,
            'movecount': m,
            'voidmass_subprocess_lower': v1_lower,
            'voidmass_subprocess_upper': v1_upper,
            'voidmass_process_lower': d_lower / root_movecount if root_movecount else 0.0,
            'voidmass_process_upper': d_upper / root_movecount if root_movecount else 0.0,
            # pooled match ratio (matchcount/movecount summed across every
            # variant/execution, not averaged per-execution the way
            # coveragemass.alignment_mass does) - matchcount = movecount -
            # deficit, so this is exactly 1 - voidmass_subprocess_*.
            # Attribution of moves to this node is still via terms_by_
            # node's leaf-label membership (classical alignments have no
            # lumping to disambiguate, unlike skip-alignments'
            # executions(), so plain leaf-set membership is sufficient -
            # no _mandatorily_implies-style inference needed here).
            'alignment_mass_pooled_lower': 1 - v1_lower,
            'alignment_mass_pooled_upper': 1 - v1_upper,
        }
        for child in node.children:
            _walk(child)

    _walk(tree)
    return table, skip_dict


def coverage_by_alignment_pn(node, skip_prob, skip_dict, variant_probs, convention='zero',
                              executions_cache=None, timed_out_ratio=None):
    '''
    \\covermove (defn:move-coverage), computed on classical (non-lumped)
    alignments: (1 - skip_prob) * alignment_mass(node, skip_dict,
    variant_probs, convention). skip_dict is voidmass_table_pn's second
    return value (already translated into coveragemass.alignment_mass's
    expected path shape - see _to_alignment_mass_path) - reuses
    coveragemass.alignment_mass entirely unmodified, verified term-by-
    term against the formal definition (see session notes): no pooling,
    per-execution match/movecount ratios averaged uniformly within an
    alignment, alignments averaged uniformly within a variant, both
    zero-denominator conventions matching the definition's own "treated
    as zero" clause exactly.

    Deliberately reuses skip-alignments' own skip_prob unchanged
    (dv.skip_probs[node]) rather than inventing a classical-alignment
    replacement for it - skip_prob answers a per-variant "was this node
    skipped at all" question via node_reached, not a move-count, so it
    was never subject to the lumping bug voidmass's deficit (or this
    metric's own prior pooled implementation) had to move off of.

    Prior implementations of this id used voidmass_table_pn's POOLED
    alignment_mass_pooled (matchcount/movecount summed across every
    execution before dividing once) - convenient since that table was
    already being computed, but not what \\covermove's own definition
    specifies. That pooled quantity is still available, honestly
    labelled, as 1 - voidmass_process (equivalently
    table[node]['alignment_mass_pooled']) - nothing was lost, this id
    just no longer claims to be that.

    executions_cache: optional, see coveragemass.alignment_mass /
    coveragemass.make_executions_cache - pass one shared cache across
    every node's call in a per-node report row (see
    lab.exp_disco_degrade._node_rows) to avoid re-walking each
    alignment's path once per node.

    timed_out_ratio: passed straight through to alignment_mass - a
    variant whose align_variant_all search timed out to zero
    alignments (present in skip_dict as an empty list, not absent -
    see voidmass_table_pn) is excluded entirely by default (None), or
    treated as contributing this synthetic per-variant ratio instead.
    See lab.exp_disco_degrade's _lower/_upper wiring for why this
    exists: align_variant_all returning nothing is a real, latent
    failure mode on any log with one slow-enough variant, previously an
    unhandled ZeroDivisionError in voidmass_table_pn's own pooled sums
    (labnotes.md finding C) rather than a defined result here.
    '''
    return (1 - skip_prob) * alignment_mass(node, skip_dict, variant_probs, convention,
                                             executions_cache, timed_out_ratio)
