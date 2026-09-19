'''
Prototype: voidmass deficit computed on CLASSICAL (standard) Petri-net
alignments, via skipalignments.alignall.align_pn_all, rather than
skip-alignments' Aligner.align_normal_form/executions().

Why this exists: coveragemass.executions() is built on skip-alignments'
normal form, which lumps an entirely-unwitnessed subtree into a single
Skip/TauPath move on its coarsest ancestor. Voidmass's deficit is
defined to count every constituent activity individually - skipprob is
the one allowed to lift to the highest block, not deficit - so building
deficit on the lumped representation is wrong. Classical alignments
have no Skip(subtree) construct at all: every missing leaf is
necessarily its own model-move, so no lumping can occur.

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
a classical alignment into coveragemass's own path shape and reuses
that machinery (executions()/matchcount/movecount and the averaging
structure) rather than reimplementing it for a second move
representation. The pooled ratio voidmass_table_pn reports
(alignment_mass_pooled_lower/_upper) is a different quantity - see
coverage_by_alignment_pn.

If this approach is adopted, this needs folding into coveragemass.py;
if not, delete both files.

Relies on skip-alignments' build_petri_net returning a 6-tuple (net,
im, fm, activity_to_id, tau_ids, id_loop_list) and on align_pn_all's
tau_ids parameter - see pyproject.toml for the pinned skip-alignments
version. id_loop_list must reach align_pn_all (see build_id_net): left
as [], the A* search can't cycle-detect on a tau-skippable loop.
'''

import copy
import logging
import time
from dataclasses import dataclass
from types import SimpleNamespace

from skipalignments import Skip, TauPath
from skipalignments.probabilities import EbiOccurance
from skipalignments.alignall import align_pn_all

from process_voids.coveragemass import observed_alignment_mass, _variant_key, min_activity_count

logger = logging.getLogger(__name__)

# Seconds a single variant's optimal-alignment search may run before
# align_pn_all abandons it and the caller falls back to bounds. The
# default this module falls back to when a caller states no budget of
# its own, so that importing it standalone works - not a claim about
# what any particular experiment should allow, which is the caller's
# policy to set (lab.params.CLASSICAL_ALIGNMENT_TIMEOUT is the lab's).
DEFAULT_ALIGNMENT_TIMEOUT = 100


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

    Every leaf's transition also records which leaf it is, in
    transition.properties[LEAF_ID], so a move can be credited to the leaf
    that fired - see leaf_by_transition. The labels alone cannot say:
    leaves sharing an activity label are given the same transition
    label, which is what lets a trace event synchronise with any of
    them. So the net is built from a copy of the tree whose leaves are
    each labelled with their own node id, which makes every transition
    label name one leaf; that id is recorded, and the label then reset to
    the one the tree itself produces. The result is the net a plain build
    gives, transition for transition, plus the record.
    '''
    by_id = {}
    for leaf in _leaves(tree):
        if leaf.id in by_id:
            raise ValueError(f'leaf id {leaf.id!r} is not unique in the tree; moves '
                             'cannot be attributed to leaves')
        by_id[leaf.id] = leaf

    net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)
    relabelled = copy.deepcopy(tree)
    for leaf in _leaves(relabelled):
        leaf.name = leaf.id
    id_net, id_im, id_fm, id_to_own_id, _id_taus, id_loops = \
        EbiOccurance().build_petri_net(relabelled)
    if sorted(map(str, id_loops)) != sorted(map(str, id_loop_list)):
        raise AssertionError('relabelling leaves changed id_loop_list; transition '
                             'identity cannot be recorded this way')

    leaf_id_of_label = {label: own for own, label in id_to_own_id.items()}
    for transition in id_net.transitions:
        leaf_id = leaf_id_of_label.get(transition.label)
        if leaf_id is None:
            continue
        leaf = by_id[leaf_id]
        transition.properties[LEAF_ID] = leaf_id
        transition.label = activity_to_id[leaf.name]
    return id_net, id_im, id_fm, activity_to_id, tau_ids, id_loop_list


LEAF_ID = 'process_voids_leaf_id'


def _leaves(tree):
    if not tree.children:
        return [tree]
    return [leaf for child in tree.children for leaf in _leaves(child)]


def leaf_by_transition(net, tree):
    '''
    {net transition name: tree leaf} for every leaf's transition in a net
    from build_id_net. A classical alignment move names the net
    transition that fired, so this is how it is credited to one leaf
    rather than to every leaf sharing its label.
    '''
    by_id = {leaf.id: leaf for leaf in _leaves(tree)}
    return {transition.name: by_id[transition.properties[LEAF_ID]]
            for transition in net.transitions if LEAF_ID in transition.properties}


def align_variant_all(activities, net, im, fm, activity_to_id, tau_id_set,
                       id_loop_list=None, timeout=DEFAULT_ALIGNMENT_TIMEOUT):
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


def align_variant(activities, net, im, fm, activity_to_id, tau_id_set,
                   id_loop_list=None, timeout=DEFAULT_ALIGNMENT_TIMEOUT):
    '''The first tied optimal alignment only - see align_variant_all.'''
    return align_variant_all(activities, net, im, fm, activity_to_id, tau_id_set,
                              id_loop_list=id_loop_list, timeout=timeout)[0]


def classify_move(t, leaves, tau_id_set):
    '''
    (kind, leaf) for one alignment move, kind in 'sync'/'log'/'model'/'tau'.
    leaf is the tree leaf whose transition fired, found through `leaves`
    (leaf_by_transition's map) by the move's model-side transition name -
    the one leaf, even where several share its label. None for 'log'
    (nothing on the model side) and for 'tau' (a silent/helper transition,
    or a genuine Tau leaf in tau_id_set - neither is a deviation, so
    neither counts towards deficit).
    '''
    trace_side, model_side = t.label
    if model_side == '>>':
        return 'log', None
    if model_side is None or model_side in tau_id_set:
        return 'tau', None
    leaf = leaves.get(t.name[1])
    if leaf is None:
        return 'tau', None
    if trace_side == '>>':
        return 'model', leaf
    return 'sync', leaf


def deficit_by_node(alignment, tree, leaves, tau_id_set):
    '''
    {node: deficit} for every node in tree, where deficit(node) is the
    count of 'model' moves in `alignment` on a leaf of node's subtree.
    deficit(execution) = movecount - matchcount simplifies to exactly
    this count, since movecount excludes tau moves and matchcount is the
    sync-move count: (sync + model) - sync.
    '''
    return {node: d for node, (d, _m) in
            terms_by_node(alignment, tree, leaves, tau_id_set).items()}


def terms_by_node(alignment, tree, leaves, tau_id_set):
    '''
    {node: (deficit, movecount)} for every node in tree, over one
    alignment - deficit = count of 'model' moves on a leaf of node's
    subtree, movecount = count of 'sync'+'model' moves there (tau moves
    excluded from both, matching coveragemass.movecount/matchcount).

    Each move counts at the one leaf whose transition fired and that
    leaf's ancestors - Definition [Void by Process-Relative Alignment
    Moves] attributes a move to exactly one execution. Classical
    alignments put every move on a leaf, so the values over any
    antichain covering the labelled leaves sum to the root's.
    '''
    per_leaf = {}
    for t in alignment:
        kind, leaf = classify_move(t, leaves, tau_id_set)
        if kind in ('model', 'sync'):
            d, m = per_leaf.get(id(leaf), (0, 0))
            per_leaf[id(leaf)] = (d + (kind == 'model'), m + 1)

    table = {}

    def _walk(node):
        if not node.children:
            table[node] = per_leaf.get(id(node), (0, 0))
            return table[node]
        d = m = 0
        for child in node.children:
            child_d, child_m = _walk(child)
            d += child_d
            m += child_m
        table[node] = (d, m)
        return table[node]

    _walk(tree)
    return table


def _to_alignment_mass_path(alignment, leaves, tau_id_set):
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
        (movecount-excluded), recognised by its label being in
        tau_id_set
      - '>>' : a pure log move, OR a structural/helper silent transition
        with no tree correspondence at all (build_petri_net inserts
        these for net routing - label=None, never in tau_id_set). Both
        cases are unattributable to any subprocess, so both become a gap
        in executions()' consecutive-run grouping - the same treatment
        an unrelated move already gets there, not a special case.

    Each move is resolved to the one leaf whose transition fired, through
    `leaves` (leaf_by_transition's map, built once by the caller), so
    leaves sharing a label are told apart.
    '''
    path = []
    for t in alignment:
        trace_side, model_side = t.label
        if model_side == '>>':
            path.append((trace_side, '>>'))
            continue
        node = leaves.get(t.name[1])
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


def sum_safe_signature(alignment, leaves, tau_id_set):
    '''
    Canonical key for "these two tied optimal alignments tell the same
    causal story": the ORDERED list of (kind, leaf) for 'sync'/'model'
    moves only - 'log' and 'tau' moves are dropped entirely.

    Keyed by leaf, not by activity label. Two alignments differing only in
    which of several same-labelled leaves fired are different alignments
    in the definition's Gamma_sigma and credit different leaves, so they
    are different stories; merging them would give the variant one story
    where it has two, and mis-weight every node including the root.

    Why this is safe for voidmass specifically, and NOT safe to reuse
    unmodified for anything else: align_pn_all's all-optimal search
    returns every total ordering of a move sequence, including orderings
    that differ ONLY in where two commuting no-ops fall relative to each
    other (a silent/tau move and an unrelated log move don't constrain
    one another, so the search enumerates both orderings as if they were
    distinct alignments - e.g. claim_25 in the claims fixture: 9 raw
    alignments, 3 real causal stories, a 6:1:2 raw split instead of the
    correct 1:1:1). Dropping log/tau collapses exactly
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
        (kind, leaf.id) for t in alignment
        for kind, leaf in [classify_move(t, leaves, tau_id_set)]
        if kind in ('sync', 'model')
    )


# TODO(\covat): a duration-weighted metric on this classical-alignment path
# needs its OWN signature function, not this one. It must keep 'log' moves
# in the key (their position relative to sync/model moves determines which
# block/traversal a duration gets attributed to - see the docstring above)
# but can still drop 'tau' moves, since silent transitions never carry a
# duration. Get this wrong and durations will get attributed to the wrong
# traversal - the same failure mode dedupe_alignments guards against for
# deficit counts, but worse to debug because it'd show up as a
# plausible-looking wrong number, not a visible asymmetry between two
# leaves.


def dedupe_alignments(alignments, leaves, tau_id_set, signature_fn=sum_safe_signature):
    '''
    One representative alignment per distinct signature_fn value, in
    first-seen order - collapses combinatorial reorderings of commuting
    moves that align_pn_all's all-optimal search returns as separate
    "optimal alignments" but which represent the same causal story.
    '''
    seen = {}
    for alignment in alignments:
        sig = signature_fn(alignment, leaves, tau_id_set)
        seen.setdefault(sig, alignment)
    return list(seen.values())


@dataclass
class VoidmassPnResult:
    '''
    voidmass_table_pn's result.

    table: {node: row} - see voidmass_table_pn for the row keys.
    skip_dict: the same deduped alignments, translated for
        coveragemass.alignment_mass. A timed-out variant maps to an
        empty list, not absent, so alignment_mass's timed_out_ratio can
        recognise it.
    timed_out_count / timed_out_weight: how many variants' alignment
        search timed out to zero alignments, and their summed
        probability. The weight is what matters - 0.03% of a log timing
        out is immaterial, 20% is not - which the count alone can't show.
    '''
    table: dict
    skip_dict: dict
    timed_out_count: int
    timed_out_weight: float


def timed_out_movecount_bound(trace_length, tree):
    '''
    X_max = 2|sigma| + min_activity_count(tree): an upper bound on the
    movecount of ANY optimal alignment of a length-|sigma| trace against
    tree, at every node - what a timed-out variant is bounded by.

    Rests on align_pn_all's cost model: log move = labelled model move =
    100000, sync = 0, every silent/tau transition = 0 (pinned
    behaviourally by CostModelPinTest - if the costs change, the bound
    must be re-derived). An alignment pairs sigma with one model path;
    every event is consumed exactly once, so sync + log = |sigma|. The
    alignment "every event a log move, plus the cheapest complete model
    path" always exists - min_activity_count(tree) is exactly the non-
    silent length of a cheapest traversal to the final state - and costs
    100000 * (|sigma| + C_root). No optimal alignment costs more, so
    log + model <= |sigma| + C_root, and

        movecount = sync + model = (|sigma| - log) + model
                  <= 2|sigma| + C_root - 2*log  <=  2|sigma| + C_root.

    A node's movecount counts a subset of the whole alignment's moves
    (terms_by_node credits each move to one leaf and its ancestors, and
    a node's leaves are a subset of the root's), so the same X bounds
    every node.

    No bound from the tree's shape alone can do this job: a loop lets
    model moves outnumber any leaf count (loop(seq(a,b), tau) on <a,a,a>
    has an optimal alignment with movecount 6 > |sigma| + 2). The bound
    is deliberately node-independent - see voidmass_table_pn on why
    voidmass_process needs one shared X.
    '''
    return 2 * trace_length + min_activity_count(tree)


def voidmass_table_pn(tree, variant_probs, net, im, fm, activity_to_id, tau_id_set,
                       id_loop_list=None, timeout=DEFAULT_ALIGNMENT_TIMEOUT):
    '''
    VoidmassPnResult whose table is the classical-alignment analogue of
    coveragemass.voidmass_table: one pass over the tree, aggregated
    (SUMMED, not averaged - see coveragemass's Voidmass / Voidage section)
    over every variant and every distinct causal story among that
    variant's tied optimal alignments (see dedupe_alignments/
    sum_safe_signature - raw alignment count is NOT used, since
    align_pn_all's all-optimal search can return the same causal story
    multiple times under different commuting-move orders).

    Timed-out variants: align_variant_all can return ZERO alignments for
    a variant (a per-variant timeout, seen in real runs). Its real
    contribution (d, m) is unknown, but 0 <= d <= m <= w * X_max, with
    X_max = timed_out_movecount_bound(|sigma|, tree), the same X at
    every node. The pooled ratio (D0 + d) / (M0 + m) increases in d and,
    at d = m, in m (since D0 <= M0), so

        upper: (d, m) = (w * X_max, w * X_max)
        lower: (d, m) = (0,         w * X_max)

    bracket the value the cell would have had without the timeout. That
    holds for voidmass_process too - its numerator is the node's deficit
    and its denominator the ROOT's movecount - but only because X is
    the same for every node: a larger X underneath than on top is not a
    bound. Several timed-out variants: the same argument, summed. Bounds
    within reporting precision of each other mean the timeouts were
    immaterial for the cell; a wide interval means little can be said.
    A timed-out variant never raises, whatever policy sits on top.

    Row keys, with W = sum of w * X_max over timed-out variants:
      deficit_lower / deficit_upper      D0 / D0 + W
      movecount                          M0 - OBSERVED from completed
                                         variants only, a measurement
                                         rather than a substitution
      movecount_bound                    M0 + W - what both bounds divide by
      voidmass_subprocess_lower/_upper   deficit_* / movecount_bound
      voidmass_process_lower/_upper      deficit_* / root movecount_bound
      alignment_mass_pooled_lower/_upper 1 - voidmass_subprocess_upper/
                                         _lower (so _lower is the
                                         smaller number, as elsewhere)

    skip_dict is the SAME deduped alignments, translated (via
    _to_alignment_mass_path) into coveragemass.alignment_mass's expected
    input shape - built in the same pass so the expensive
    align_variant_all call is only ever made once per variant, not once
    for the pooled sums here and again for coverage_by_alignment_pn's
    (correctly non-pooled - see defn:move-coverage) mass term.

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
    skip_dict = {}
    timed_out_count = 0
    timed_out_weight = 0.0
    timed_out_bound = 0.0

    def _add(node, d, m, share):
        deficit_sum[node] = deficit_sum.get(node, 0.0) + share * d
        movecount_sum[node] = movecount_sum.get(node, 0.0) + share * m

    leaves = leaf_by_transition(net, tree)
    n_variants = len(variant_probs)
    total_started = time.monotonic()
    n_near_timeout = 0
    for i, (variant, weight) in enumerate(variant_probs.items(), start=1):
        started = time.monotonic()
        raw_alignments = align_variant_all(list(variant), net, im, fm, activity_to_id,
                                            tau_id_set, id_loop_list=id_loop_list, timeout=timeout)
        elapsed = time.monotonic() - started
        alignments = dedupe_alignments(raw_alignments, leaves, tau_id_set)
        # Flag anything running close to the per-variant timeout - the
        # likely signature of the id_loop_list gap (a tau-skippable loop
        # the A* search can't cycle-detect, so it burns the full budget
        # rather than resolving quickly or hanging outright) rather than
        # genuine alignment cost.
        near_timeout = elapsed > 0.9 * timeout
        n_near_timeout += near_timeout
        log = logger.warning if near_timeout else logger.debug
        log('variant %d/%d (%d activities, weight=%.4g): %.2fs, %d raw -> %d deduped '
            'alignments%s', i, n_variants, len(variant), weight, elapsed,
            len(raw_alignments), len(alignments),
            ' [near timeout - possible id_loop_list gap]' if near_timeout else '')
        if not alignments:
            timed_out_count += 1
            timed_out_weight += weight
            timed_out_bound += weight * timed_out_movecount_bound(len(variant), tree)
            logger.warning('variant %d/%d (weight=%.4g): 0 alignments (timed out) - '
                            'bounded, not estimated', i, n_variants, weight)
        else:
            share = weight / len(alignments)
            for alignment in alignments:
                for node, (d, m) in terms_by_node(alignment, tree, leaves, tau_id_set).items():
                    _add(node, d, m, share)
        skip_dict[_variant_key(variant)] = [
            SimpleNamespace(path=_to_alignment_mass_path(alignment, leaves, tau_id_set))
            for alignment in alignments
        ]

    logger.info('voidmass_table_pn: %d variants in %.1fs (%d near-timeout, %d timed out, '
                'timed-out weight %.4g)', n_variants, time.monotonic() - total_started,
                n_near_timeout, timed_out_count, timed_out_weight)

    root_bound = movecount_sum.get(tree, 0.0) + timed_out_bound
    table = {}

    def _walk(node):
        d = deficit_sum.get(node, 0.0)
        m = movecount_sum.get(node, 0.0)
        d_upper = d + timed_out_bound
        m_bound = m + timed_out_bound
        subprocess_lower = d / m_bound if m_bound else 0.0
        subprocess_upper = d_upper / m_bound if m_bound else 0.0
        table[node] = {
            'deficit_lower': d,
            'deficit_upper': d_upper,
            'movecount': m,
            'movecount_bound': m_bound,
            'voidmass_subprocess_lower': subprocess_lower,
            'voidmass_subprocess_upper': subprocess_upper,
            'voidmass_process_lower': d / root_bound if root_bound else 0.0,
            'voidmass_process_upper': d_upper / root_bound if root_bound else 0.0,
            # pooled match ratio (matchcount/movecount summed across every
            # variant/execution, not averaged per-execution the way
            # coveragemass.alignment_mass does) - matchcount = movecount -
            # deficit, so this is 1 - voidmass_subprocess. A move is
            # attributed to this node through the one leaf whose
            # transition fired (terms_by_node).
            'alignment_mass_pooled_lower': 1 - subprocess_upper,
            'alignment_mass_pooled_upper': 1 - subprocess_lower,
        }
        for child in node.children:
            _walk(child)

    _walk(tree)
    return VoidmassPnResult(table, skip_dict, timed_out_count, timed_out_weight)


def coverage_by_alignment_pn(node, skip_prob, skip_dict, variant_probs,
                              executions_cache=None, timed_out_ratio=None):
    '''
    \\covermove (defn:move-coverage), computed on classical (non-lumped)
    alignments: (1 - skip_prob) * coveragemass.observed_alignment_mass(
    node, skip_dict, variant_probs). skip_dict is voidmass_table_pn's
    result.skip_dict, already translated into that function's expected
    path shape - see _to_alignment_mass_path.

    The mass is conditioned on observation: an execution with no
    synchronous move is excluded, and so is a trace whose alignments
    hold no observed execution of node. That keeps the two factors
    independent - the skip probability answers whether node was recorded
    at all, the mass how completely it was recorded where it was - so
    coverage falls linearly in a submodel's absence rather than
    quadratically. See observed_alignment_mass for the definition term
    by term.

    Deliberately reuses skip-alignments' own skip_prob unchanged
    (dv.skip_probs[node]) rather than inventing a classical-alignment
    replacement for it - skip_prob answers a per-variant "was this node
    skipped at all" question via node_reached, not a move-count, so the
    normal form's lumping that rules skip-alignments out for voidmass's
    deficit doesn't affect it.

    The pooled match ratio - matchcount/movecount summed across every
    execution before dividing once - is a different quantity, reported by
    voidmass_table_pn as table[node]['alignment_mass_pooled_lower'/
    'alignment_mass_pooled_upper'] (equivalently 1 - voidmass_subprocess).

    executions_cache: optional, see coveragemass.observed_alignment_mass /
    coveragemass.make_executions_cache - pass one shared cache across
    every node's call in a per-node report row (see
    lab.exp_disco_degrade._node_rows) to avoid re-walking each
    alignment's path once per node.

    timed_out_ratio: passed straight through - a variant whose
    align_variant_all search timed out to zero alignments (present in
    skip_dict as an empty list, not absent - see voidmass_table_pn)
    counts as one observed unit at this synthetic ratio, so 0.0 and 1.0
    bracket the value the cell would have had; the default None drops it
    and renormalises over the variants that are left. See
    lab.exp_disco_degrade's _lower/_upper wiring for why this exists:
    align_variant_all returning nothing is a real, latent failure mode
    on any log with one slow-enough variant.
    '''
    return (1 - skip_prob) * observed_alignment_mass(node, skip_dict, variant_probs,
                                                      executions_cache, timed_out_ratio)
