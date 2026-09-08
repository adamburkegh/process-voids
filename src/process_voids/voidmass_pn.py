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
It has NOT been through the same design scrutiny as coveragemass.py -
in particular it does not implement the Definition [Executions] grouping
(loop iterations collapsing into one traversal). That omission is narrower
than it looks: Bär et al.'s conditions (2)/(3) exist to keep one
traversal's moves separate from another's, which only matters when you
AVERAGE per traversal - which alignment_mass does and voidmass doesn't.
Summing over the whole alignment gives the same total regardless of how
moves are partitioned into traversals, since it's the same set of moves
either way - so traversal-separation is genuinely unnecessary for a pure
sum. It is NOT a licence to drop move-to-subprocess attribution, though:
voidmass_process's own divisor (root movecount) and voidmass_subprocess's
divisor (the node's own movecount) both still need to know which moves
belong to which subprocess - terms_by_node below provides that via leaf-
label membership (a purely structural check, no traversal concept
needed), which is sufficient for that narrower job. A future coverage_by_
alignment-style, duration-weighted metric on this classical-alignment
path WOULD need the full executions() machinery back, since \\covat
weights by duration per traversal, not just a flat sum.

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

from skipalignments.probabilities import EbiOccurance
from skipalignments.alignall import align_pn_all

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
    Classical-alignment analogue of coveragemass.voidmass_table: one
    pass over the tree, {node: {'deficit', 'movecount', 'voidmass_subprocess',
    'voidmass_process'}}, aggregated (SUMMED, not averaged - see
    voidmass-brief.md) over every variant and every distinct causal
    story among that variant's tied optimal alignments (see
    dedupe_alignments/sum_safe_signature - raw alignment count is NOT
    used, since align_pn_all's all-optimal search can return the same
    causal story multiple times under different commuting-move orders).

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
    deficit_sum = {}
    movecount_sum = {}

    def _add(node, d, m, share):
        deficit_sum[node] = deficit_sum.get(node, 0.0) + share * d
        movecount_sum[node] = movecount_sum.get(node, 0.0) + share * m

    id_to_activity = {v: k for k, v in activity_to_id.items()}
    n_variants = len(variant_probs)
    total_started = time.monotonic()
    n_near_timeout = 0
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
        share = weight / len(alignments)
        for alignment in alignments:
            for node, (d, m) in terms_by_node(alignment, tree, activity_to_id, tau_id_set).items():
                _add(node, d, m, share)

    logger.info('voidmass_table_pn: %d variants in %.1fs (%d near-timeout)',
                n_variants, time.monotonic() - total_started, n_near_timeout)

    root_movecount = movecount_sum.get(tree, 0.0)
    table = {}

    def _walk(node):
        d = deficit_sum.get(node, 0.0)
        m = movecount_sum.get(node, 0.0)
        v1 = d / m if m else 0.0
        table[node] = {
            'deficit': d,
            'movecount': m,
            'voidmass_subprocess': v1,
            'voidmass_process': d / root_movecount if root_movecount else 0.0,
            # pooled match ratio (matchcount/movecount summed across every
            # variant/execution, not averaged per-execution the way
            # coveragemass.alignment_mass does) - matchcount = movecount -
            # deficit, so this is exactly 1 - voidmass_subprocess. Attribution
            # of moves to this node is still via terms_by_node's leaf-label
            # membership (classical alignments have no lumping to
            # disambiguate, unlike skip-alignments' executions(), so plain
            # leaf-set membership is sufficient - no _mandatorily_implies-
            # style inference needed here).
            'alignment_mass_pooled': 1 - v1,
        }
        for child in node.children:
            _walk(child)

    _walk(tree)
    return table


def coverage_by_alignment_pn(node, skip_prob, table):
    '''
    Classical-alignment analogue of coveragemass.coverage_by_alignment:
    (1 - skip_prob) * alignment_mass_pooled. Deliberately reuses skip-
    alignments' own skip_prob unchanged (dv.skip_probs[node]) rather
    than inventing a classical-alignment replacement for it - skip_prob
    answers a per-variant "was this node skipped at all" question via
    node_reached, not a move-count, so it was never subject to the
    lumping bug voidmass's deficit had to move off of. Only the
    alignment_mass side (a move-count ratio) needed replacing - see
    session notes.
    '''
    return (1 - skip_prob) * table[node]['alignment_mass_pooled']
