
import bisect

from skipalignments import *


'''
Currently ignores silents. Assumes labels are activity ids
'''
def update_activity_weights(pt:ProcessTree,slpn):
    leaves = pt.get_leafs()
    for leaf in leaves:
        for tran in slpn.transitions:
            # print(f'uaw( {tran} ... {leaf.id} )')
            # inefficient due to list instead of dict
            if tran['label'] == leaf.id:
                leaf.weight = tran['weight']

def infer_operator_weights(pt:ProcessTree):
    if isinstance(pt,Activity) or isinstance(pt,Tau):
        return
    for child in pt.children:
        infer_operator_weights(child)
    if isinstance(pt,Xor):
        pt.weight = sum([ child.weight for child in pt.children  ])
    if isinstance(pt,And) or isinstance(pt,Sequence) or isinstance(pt,Loop):
        pt.weight = sum([ child.weight for child in pt.children  ]) \
                    / len(pt.children)


def transfer_pt_weights(pt:ProcessTree, slpn):
    update_activity_weights(pt,slpn)
    infer_operator_weights(pt)


'''
=====================================================================================
Mandatory Node Count

Diagnostic, not a void/coverage metric - no log, alignment, or skip_probs
needed, purely a property of the discovered tree's structure. Counts how
much of a tree the void/coverage metrics can actually speak about.

A node has no silent alternative - is "mandatory" - if nowhere on its
path to the root can it be silently skipped: no ancestor Xor has a Tau
sibling along that path (the model could choose to do nothing instead
of this node's subtree), and no ancestor Loop position on that path is
the redo-child (a Loop can always iterate zero extra times, regardless
of whether Tau appears explicitly there).

Inductive Miner at noise_threshold=0.0 wraps nearly every leaf in
Xor(Tau, activity) to guarantee it perfectly replays its own source log
(see lab.discovery's discover_inductive docstring) - on such a tree
almost nothing is mandatory, skip_prob(root) is trivially 0, and every
void metric built on skip_prob is uninformative by construction, not
because there is genuinely nothing missing. mandatory_node_count /
total_node_count is a cheap sanity check for that failure mode: a low
ratio is a warning sign that a void metric's near-zero reading reflects
the model's own permissiveness, not real coverage.

Both counts exclude Tau nodes themselves - a Tau leaf represents "do
nothing", not a thing a void/coverage metric could ever meaningfully be
asked about.
'''

TREE_METRIC_KEYS = ('mandatory_node_count', 'total_node_count')


def has_silent_alternative(pt: ProcessTree):
    '''
    True if `pt` can be silently skipped: somewhere on its path to the
    root, either an ancestor Xor has a Tau sibling (the model can choose
    to do nothing instead of pt's subtree), or an ancestor Loop position
    on that path is the redo-child (structurally optional regardless of
    Tau - a Loop can always iterate zero extra times).

    Walks the FULL ancestor chain, not just pt's immediate parent - a
    silent alternative anywhere above pt makes pt itself skippable even
    when pt's own immediate parent is a plain Sequence/And.
    '''
    current = pt
    while current.parent is not None:
        parent = current.parent
        if isinstance(parent, Xor) and any(isinstance(sibling, Tau)
                                            for sibling in parent.children):
            return True
        if isinstance(parent, Loop) and parent.children[1] is current:
            return True
        current = parent
    return False


def mandatory_node_count(pt: ProcessTree):
    '''Number of non-Tau nodes (leaf and internal) in `pt` with no
    silent alternative - see this section's module docstring.'''
    count = 0

    def _walk(node):
        nonlocal count
        if not isinstance(node, Tau) and not has_silent_alternative(node):
            count += 1
        for child in node.children:
            _walk(child)

    _walk(pt)
    return count


def total_node_count(pt: ProcessTree):
    '''Number of non-Tau nodes (leaf and internal) in `pt` - the
    denominator for reading mandatory_node_count as a fraction rather
    than a bare count only meaningful relative to a specific tree's
    size.'''
    count = 0

    def _walk(node):
        nonlocal count
        if not isinstance(node, Tau):
            count += 1
        for child in node.children:
            _walk(child)

    _walk(pt)
    return count


'''
Metric which indicates how much of the process is tracked by the data.

Pre: Tree has weights
'''
def mass_by_weight(pt:ProcessTree, skip_probs:dict):
    if isinstance(pt,Activity) or isinstance(pt,Tau):
        return 1 - skip_probs[pt]
    child_coverage = []
    total_weight = sum( [child.weight for child in pt.children] )
    if isinstance(pt,Xor):
        return sum([ mass_by_weight(child,skip_probs)* child.weight \
                        / total_weight \
                        for child in pt.children  ])
    if isinstance(pt,And) or isinstance(pt,Sequence) or isinstance(pt,Loop):
        return sum([ mass_by_weight(child,skip_probs) \
                          for child in pt.children  ]) \
                          / len(pt.children)
    raise ValueError('Unrecognised process tree node')


'''
Weight-averaged void mass: same tree aggregation as mass_by_weight
(Xor weighted by relative child weight, And/Sequence/Loop averaged
uniformly across children) but skip_probs[leaf] directly rather than
its complement, ie skipprob * mass instead of (1-skipprob) * mass.

Both aggregations are (weighted) averages, and averaging distributes
linearly over the complement (avg(1-x_i) == 1-avg(x_i)) - so by
induction voidage_by_weight(pt) == 1 - mass_by_weight(pt, skip_probs)
for every node, leaf or not. No separate tree-walk needed.
'''
def voidage_by_weight(pt:ProcessTree, skip_probs:dict):
    return 1 - mass_by_weight(pt, skip_probs)


'''
=====================================================================================
Coverage by Duration

Uses the implied duration of the activities present in the model but missing from the
log to estimate coverage.

    dur(sigma)           : duration of a trace sigma from its timestamps
    mi(m, sigma)          : trace indexes (excl. index 1) whose activity label is part
                            of the alphabet act(m) of a (sub)model m
    sdur(m, sigma)        : summed / averaged duration attributable to m
    cov_dt(m, m_sub, L)   : coverage by duration of submodel m_sub w.r.t. model m,
                            computed over log L

A trace sigma is expected to be an ordered sequence of events, each of which provides
    - event['concept:name']    the activity label, and
    - event['time:timestamp']  a datetime timestamp.
This matches the pm4py / XES event representation already used elsewhere in this
project (see e.g. probabilities.py, disco.py, pvoid.py), so a pm4py event log (or
any list of such traces) can be passed directly as `log`.

skip_probs is the same dict already used by mass_by_weight: ProcessTree node ->
P_skip(node). Here it supplies P_skip(m, m_sub, L) for the submodel pt.

Ported from Josh Gong's coverage-by-duration work.
'''

def dur(sigma):
    '''
    dur(sigma) = pi_time( sigma[|sigma|] ) - pi_time( sigma[1] )
    dur(<>)    = 0
    '''
    if len(sigma) == 0:
        return 0
    return (sigma[-1]['time:timestamp'] - sigma[0]['time:timestamp']).total_seconds()


def mi(activities, sigma):
    '''
    mi(m,sigma) = { i | 2 <= i <= |sigma|  and  pi_act(sigma[i]) in act(m) }

    0-indexed here as positions 1 .. len(sigma)-1 (i.e. excluding the first event,
    which has no predecessor to measure an inter-event duration against).
    '''
    return [ i for i in range(1, len(sigma))
                if sigma[i]['concept:name'] in activities ]


def sdur(activities, sigma, has_submodels):
    '''
    sdur(m,sigma) = sum_{i in mi(m,sigma)} ( pi_time(sigma[i]) - pi_time(sigma[i-1]) )   , sub(m) != {}
    sdur(m,sigma) = ( ... same sum ... ) / |mi(m,sigma)|                                 , sub(m) == {} and |mi(m,sigma)| > 0
    sdur(m,sigma) = 0                                                                    , otherwise

    `has_submodels` corresponds to sub(m) != {}, i.e. whether m is a composite
    (operator) node with children, as opposed to a leaf/atomic activity.
    '''
    indices = mi(activities, sigma)
    if len(indices) == 0:
        return 0
    total = sum([ (sigma[i]['time:timestamp'] - sigma[i-1]['time:timestamp']).total_seconds()
                    for i in indices ])
    if has_submodels:
        return total
    else:
        return total / len(indices)


def log_to_traces(log):
    '''
    Normalises a pm4py event log into a list of traces, each trace an ordered
    list of events (dict-like, exposing ['concept:name'] / ['time:timestamp']),
    as expected by coverage_by_duration.

    Works whether pm4py.read_xes returned a pandas DataFrame (modern pm4py
    default) or a classic pm4py EventLog / list of Trace objects (each Trace
    is already dict-like per event and iterable in order).
    '''
    if hasattr(log, 'groupby'):
        # pandas DataFrame representation
        traces = []
        for _, group in log.groupby('case:concept:name'):
            group = group.sort_values('time:timestamp')
            traces.append(group.to_dict('records'))
        return traces
    else:
        # classic pm4py EventLog: already a list of ordered Trace objects
        return list(log)


def coverage_by_duration(pt:ProcessTree, log, skip_probs:dict, total_dur=None):
    '''
    cov_dt(m, m_sub, L) = (1 - P_skip(m, m_sub, L)) *  sum_{sigma in L} sdur(m_sub,sigma)/sum_{sigma in L} dur(sigma)

    pt is treated as the submodel m_sub (pt in gsub(m)) whose coverage is computed.
    P_skip(m, m_sub, L) is taken from skip_probs[pt], mirroring how mass_by_weight
    looked up skip probabilities per node.
    log is an iterable of traces sigma (each a sequence of events as described above).
    Requires log to contain at least one trace of non-zero duration.

    total_dur (optional): the denominator sum_{sigma in L} dur(sigma), precomputed.
    Pass this in when calling coverage_by_duration for many nodes over the same
    log (e.g. once per tree node) to avoid recomputing it every time; if omitted
    it is computed from `log` as usual.
    '''
    activities = set(pt.get_leaf_labels())
    has_submodels = not isinstance(pt, LeafNode)

    if total_dur is None:
        total_dur = sum([ dur(sigma) for sigma in log ])
    if total_dur == 0:
        raise ValueError('Log L must contain at least one trace of non-zero duration')

    total_sdur = sum([ sdur(activities, sigma, has_submodels) for sigma in log ])

    p_skip = skip_probs[pt]

    return (1 - p_skip) * total_sdur / total_dur


'''
=====================================================================================
Coverage by Alignment Correspondence

For each trace variant, examines its skip-alignment(s) against the full
model, extracts the executions of the submodel pt within each alignment,
and averages how much of each execution's non-silent moves were actually
synchronous (matched a real log event) rather than model-only.

An execution is a maximal run of a skip-alignment's moves (by original
position) whose model element falls inside pt's subtree - except when pt
itself is a Loop, whose entire set of moves (all iterations) counts as
one execution (Definition [Executions]: a loop execution captures all of
its iterations). Log moves never appear in an execution.

matchcount(execution) counts synchronous moves (log event matched to a
real, unwrapped leaf). movecount(execution) counts every non-silent move
(synchronous, plus non-silent Skip - ie a required activity inserted with
real cost); TauPath (silent, cost 0) moves are excluded from movecount.

Two conventions for combining ratios across executions/alignments/trace
variants when a unit has no valid (movecount>0) execution at all:
  'zero'         - such a unit contributes 0 (the literal definition:
                   1/|Gamma| and 1/|P| are treated as zero when the
                   denominator would be zero)
  'renormalised' - such units are excluded and the remaining weight is
                   renormalised, so a submodel that is rarely exercised
                   but well-corroborated when it is is not conflated with
                   one that is poorly recorded. Defaults to mass=1 if no
                   unit anywhere ever contributes (never exercised at all).

Ported from a design worked out against a reference implementation; see
tests/test_coverage_by_alignment*.py for the worked examples this was
checked against.
'''

def _classify_move(model_elem):
    '''
    Returns (node, kind) for a skip-alignment path's model-side element,
    or (None, None) for a pure log move ('>>').
    kind is one of 'sync' (real, matched execution), 'skip' (non-silent
    model move), 'tau' (silent model move).
    '''
    if model_elem == '>>':
        return None, None
    if isinstance(model_elem, TauPath):
        return model_elem.node, 'tau'
    if isinstance(model_elem, Skip):
        return model_elem.node, 'skip'
    return model_elem, 'sync'


def _mandatorily_implies(ancestor_node:ProcessTree, pt:ProcessTree):
    '''
    True if executing/skipping ancestor_node (a strict ancestor of pt)
    unambiguously implies pt was executed/skipped too - ie every step
    between ancestor_node and pt is a mandatory position (a Sequence/And
    child, or a Loop's do-child), never an Xor branch or a Loop's
    redo-child, both of which are genuinely ambiguous/optional.
    '''
    current = pt
    while current is not ancestor_node:
        parent = current.parent
        if parent is None:
            return False
        if isinstance(parent, Xor):
            return False
        if isinstance(parent, Loop) and parent.children[1] is current:
            return False
        current = parent
    return True


def _root_of(pt:ProcessTree):
    node = pt
    while node.parent is not None:
        node = node.parent
    return node


def _ancestor_chains(tree:ProcessTree):
    '''
    {node: [ancestor-or-self, ...]} for every node in tree - the set a
    move classified at `node` is relevant to under executions()'s first
    (pt.contains_tree(node)) branch, for every pt at once. Built
    top-down in one pass (each node's chain = its parent's chain plus
    itself), so this is O(n) total rather than one upward walk per
    node.
    '''
    chains = {}

    def _walk(node, parent_chain):
        chain = parent_chain + [node]
        chains[node] = chain
        for child in node.children:
            _walk(child, chain)

    _walk(tree, [])
    return chains


def _implied_descendant_sets(tree:ProcessTree):
    '''
    {node: {node itself, plus every descendant it mandatorily implies}}
    for every node in tree - the forward/downward direction of
    _mandatorily_implies, precomputed once per tree rather than
    re-walked (upward, from pt to a candidate ancestor) per move per
    queried node. Purely structural - does not depend on any alignment
    - so this only needs recomputing when the tree itself changes.
    '''
    implied = {}

    def _walk(node):
        result = {node}
        if isinstance(node, Xor):
            for child in node.children:
                _walk(child)
        else:
            for idx, child in enumerate(node.children):
                if isinstance(node, Loop) and idx == 1:
                    _walk(child)
                    continue
                result |= _walk(child)
        implied[node] = result
        return result

    _walk(tree)
    return implied


def _group_executions(relevant, pt:ProcessTree):
    '''
    Groups an already-filtered, path-order list of (i, log_elem,
    model_elem) triples into pt's executions, per Definition
    [Executions]: one execution for the whole list if pt is a Loop
    (all iterations count as one execution), otherwise one execution
    per maximal run of consecutive original path positions.
    '''
    if not relevant:
        return []
    if isinstance(pt, Loop):
        return [[(log_elem, model_elem) for _, log_elem, model_elem in relevant]]
    groups = []
    current = [relevant[0]]
    for prev, curr in zip(relevant, relevant[1:]):
        if curr[0] == prev[0] + 1:
            current.append(curr)
        else:
            groups.append(current)
            current = [curr]
    groups.append(current)
    return [[(log_elem, model_elem) for _, log_elem, model_elem in group]
            for group in groups]


def executions_by_node(path, tree:ProcessTree, ancestor_chains=None, implied_sets=None):
    '''
    {node: executions(path, node)} for EVERY node in tree, from one
    pass over `path` - batches what calling executions(path, node) once
    per node would otherwise redo (re-walking the same path from
    scratch each time). Pass precomputed ancestor_chains/implied_sets
    (this function's own _ancestor_chains/_implied_descendant_sets
    outputs) when calling repeatedly against the same tree, to also
    skip re-deriving those structural, alignment-independent maps - see
    make_executions_cache/alignment_mass.

    executions() (below) is a thin facade over this for the single-node
    case the existing test suite and simpler callers use; this is the
    real implementation both share, not a second maintained copy of the
    partitioning logic.
    '''
    if ancestor_chains is None:
        ancestor_chains = _ancestor_chains(tree)
    if implied_sets is None:
        implied_sets = _implied_descendant_sets(tree)

    relevant_by_node = {}
    for i, (log_elem, model_elem) in enumerate(path):
        node, kind = _classify_move(model_elem)
        if node is None:
            continue
        for pt in ancestor_chains[node]:
            relevant_by_node.setdefault(pt, []).append((i, log_elem, model_elem))
        if kind in ('skip', 'tau'):
            for pt in implied_sets[node]:
                if pt is node:
                    continue
                relevant_by_node.setdefault(pt, []).append((i, log_elem, model_elem))

    result = {pt: _group_executions(relevant, pt)
              for pt, relevant in relevant_by_node.items()}
    for pt in ancestor_chains:
        result.setdefault(pt, [])
    return result


def executions(path, pt:ProcessTree):
    '''
    Partitions a skip-alignment path (a list of (log_elem, model_elem)
    pairs) into the executions of submodel pt, per Definition
    [Executions]. Each execution is itself a list of (log_elem,
    model_elem) pairs, restricted to pt's subtree.

    A skip-alignment normal form lumps an entirely-unwitnessed subtree
    into one Skip/TauPath on its coarsest ancestor, rather than naming
    every descendant leaf. Where pt sits only on mandatory positions
    beneath such a lump (see _mandatorily_implies), that lump is also
    pt's own execution - pt inherits the ancestor's fate rather than
    being reported as vacuous/never-reached.

    A facade over executions_by_node - computes every node's executions
    from one pass over `path` and returns just pt's, rather than
    walking `path` again for a single node. Fine for this function's
    own callers (tests, one-off lookups); a caller that needs many
    nodes' executions over the same path (alignment_mass, once per tree
    node in a report row) should call executions_by_node directly with
    a shared cache instead - see make_executions_cache.
    '''
    tree = _root_of(pt)
    return executions_by_node(path, tree).get(pt, [])


def matchcount(execution):
    return sum(1 for log_elem, model_elem in execution
               if _classify_move(model_elem)[1] == 'sync' and log_elem != '>>')


def movecount(execution):
    return sum(1 for _, model_elem in execution
               if _classify_move(model_elem)[1] in ('sync', 'skip'))


def _variant_key(variant):
    return ', '.join(variant)


def make_executions_cache(tree:ProcessTree):
    '''
    Structural maps for `tree` (ancestor_chains, implied_sets) plus a
    per-path memo, shared across many alignment_mass calls against
    different nodes of the SAME tree/skip_dict (see
    lab.exp_disco_degrade._node_rows, which calls alignment_mass once
    per node via coverage_by_alignment/coverage_by_alignment_pn) - pass
    the SAME cache object to every one of those calls so both the
    structural maps and each alignment's batched executions_by_node
    result are computed once per report row, not once per node. Safe to
    share between coverage_by_alignment's skip-alignment paths and
    coverage_by_alignment_pn's classical-alignment paths at once (both
    score the same tree) - the per-path memo is keyed by path object
    identity, so entries from either source just coexist.
    '''
    return {
        'tree': tree,
        'ancestor_chains': _ancestor_chains(tree),
        'implied_sets': _implied_descendant_sets(tree),
        'by_path': {},
    }


def _executions_for(path, pt, cache):
    by_node = cache['by_path'].get(id(path))
    if by_node is None:
        by_node = executions_by_node(path, cache['tree'], cache['ancestor_chains'],
                                      cache['implied_sets'])
        cache['by_path'][id(path)] = by_node
    return by_node.get(pt, [])


def alignment_mass(pt:ProcessTree, skip_dict:dict, variant_probs:dict,
                    convention='zero', executions_cache=None):
    '''
    The mass term of Coverage by Alignment Correspondence: the (1-P_skip)
    factor is not applied here (see coverage_by_alignment), so this can
    be computed and tested without a real skip-probability estimation
    (no ebi dependency) - skip_dict is a variant-key -> list of skip-
    alignment State (as produced by Aligner.align2/DerivationPipeline.
    compute_skip_alignments), variant_probs is a variant tuple -> weight.

    convention: 'zero' or 'renormalised', see module docstring above.

    executions_cache: optional, from make_executions_cache(tree) where
    tree is pt's tree - when a caller is going to call alignment_mass
    for many nodes of the same tree over the same skip_dict (the usual
    per-node report-row case), passing a shared cache avoids re-walking
    each alignment's path once per node. Omit for a one-off single-node
    call (falls back to executions(), itself a facade doing the
    equivalent single-tree-walk work with no cache to share).
    '''
    if convention not in ('zero', 'renormalised'):
        raise ValueError("convention must be 'zero' or 'renormalised'")

    variant_terms = []
    for variant, weight in variant_probs.items():
        states = skip_dict.get(_variant_key(variant), [])
        if not states:
            continue
        alignment_values = []
        for state in states:
            execs_all = (_executions_for(state.path, pt, executions_cache)
                         if executions_cache is not None
                         else executions(state.path, pt))
            execs = [e for e in execs_all if movecount(e) > 0]
            if not execs:
                alignment_values.append(None)
                continue
            ratios = [matchcount(e) / movecount(e) for e in execs]
            alignment_values.append(sum(ratios) / len(ratios))
        variant_terms.append((weight, alignment_values))

    if convention == 'zero':
        total = 0.0
        for weight, alignment_values in variant_terms:
            if not alignment_values:
                continue
            per_variant = sum(v if v is not None else 0.0
                               for v in alignment_values) / len(alignment_values)
            total += weight * per_variant
        return total

    weighted_sum = 0.0
    weight_total = 0.0
    for weight, alignment_values in variant_terms:
        valid = [v for v in alignment_values if v is not None]
        if not valid:
            continue
        weighted_sum += weight * (sum(valid) / len(valid))
        weight_total += weight
    if weight_total == 0:
        return 1.0
    return weighted_sum / weight_total


def coverage_by_alignment(pt:ProcessTree, dv, convention='zero', executions_cache=None):
    '''
    Coverage by Alignment Correspondence: (1 - P_skip(pt)) * alignment_mass(...).
    dv is a computed DerivationPipeline (dv.pl for variant weights,
    dv.skip_dict_backup for the per-variant skip-alignments, dv.skip_probs
    for P_skip) - the full, ebi-backed pipeline. See alignment_mass for
    the ebi-free mass computation this wraps, and for executions_cache.
    '''
    mass = alignment_mass(pt, dv.skip_dict_backup, dv.pl, convention, executions_cache)
    return (1 - dv.skip_probs[pt]) * mass


'''
=====================================================================================
Voidmass / Voidage

Per voidmass-brief.md: where alignment_mass averages a match/movecount
ratio per execution and then averages those ratios up (smoothing out
severity), voidmass sums deficit (= movecount - matchcount, the model
moves) and movecount separately across every execution of every optimal
alignment of every trace variant, and divides once at the end. This
keeps severity: a subprocess entirely missing across many traces reads
as more void than the same subprocess missing in only one.

Two divisor choices give two variants, sharing everything except the
divisor:
  voidmass_subprocess (variant 1) - sum(deficit)/sum(movecount) over
      pt's own executions. Scale-free: a 1-activity and a 20-activity
      subprocess both entirely missing both score 1.0.
  voidmass_process (variant 2) - sum(deficit) over pt / sum(movecount)
      over the WHOLE model. Size-preserving and additive: subprocess
      values sum over any antichain through the decomposition, so the
      root equals the sum of everything below it - a decomposable root
      headline, which neither alignment_mass nor skipprob provide.

A further factor of skip_prob(pt) (from the same DerivationPipeline
already used by coverage_by_alignment) gives voidage, variants 3 and 4
- voidmass_table's skip_probs parameter. These lose the additivity
variant 2 has, since skip_prob varies per node and a product of
per-node quantities does not sum over a cut the way voidmass_process
does - confirmed numerically rather than assumed (per the brief), see
tests/test_voidmass.py.

Reference: voidmass-brief.md's worked example on the payment running
example; see tests/test_voidmass.py, checked against the real aligner's
output on that same fixture (not the brief's own hand-picked
alignments, which may not match what the real aligner actually finds -
see test_coverage_by_alignment.py's own TiesAcrossOptimalAlignmentsTest
for a precedent of that divergence).
'''

def deficit(execution):
    '''Model moves in this execution: movecount - matchcount.'''
    return movecount(execution) - matchcount(execution)


def voidmass_terms(pt:ProcessTree, skip_dict:dict, variant_probs:dict):
    '''
    (deficit_sum, movecount_sum) for pt: summed - NOT averaged - over
    every execution of every optimal alignment of every trace variant,
    weighted by variant probability. Multiple tied optimal alignments
    within one variant share that variant's weight equally (same
    convention as alignment_mass's per-state averaging).
    '''
    deficit_sum = 0.0
    movecount_sum = 0.0
    for variant, weight in variant_probs.items():
        states = skip_dict.get(_variant_key(variant), [])
        if not states:
            continue
        share = weight / len(states)
        for state in states:
            for e in executions(state.path, pt):
                deficit_sum += share * deficit(e)
                movecount_sum += share * movecount(e)
    return deficit_sum, movecount_sum


def voidmass_table(tree:ProcessTree, skip_dict:dict, variant_probs:dict, skip_probs:dict=None):
    '''
    node -> {deficit, movecount, voidmass_subprocess, voidmass_process}
    (plus voidage_subprocess, voidage_process when skip_probs is given -
    see below) for every node in tree, computed in one pass (well, one
    pass per node - see voidmass_terms; performance not yet a concern
    here, no real log has exercised this path at scale).

    Pass skip_probs (a computed DerivationPipeline's dv.skip_probs, the
    same P_skip already used by coverage_by_alignment) to also get
    voidage_subprocess and voidage_process (variants 3 and 4 -
    voidmass * skip_prob(node)). Omit for just variants 1/2, which need
    no DerivationPipeline/skip_prob at all.

    Variants 3/4 lose the additivity variant 2 has - but this is
    conditional, not universal: it only breaks when two or more siblings
    each have nonzero deficit AND different skip_probs (a product of
    per-node quantities doesn't sum over a cut in general). It can
    coincidentally hold, e.g. when only one child ever contributes
    deficit. See tests/test_voidmass.py's VoidageTest (a degenerate case
    where it happens to hold, documented as such) and
    VoidageAdditivityLossTest (a constructed non-degenerate case
    confirming the brief's claim genuinely does hold in general) - per
    the brief's own instruction not to assume it either way.
    '''
    _, root_movecount = voidmass_terms(tree, skip_dict, variant_probs)
    table = {}

    def _walk(node):
        d, m = voidmass_terms(node, skip_dict, variant_probs)
        v1 = d / m if m else 0.0
        v2 = d / root_movecount if root_movecount else 0.0
        row = {
            'deficit': d,
            'movecount': m,
            'voidmass_subprocess': v1,
            'voidmass_process': v2,
        }
        if skip_probs is not None:
            sp = skip_probs[node]
            row['voidage_subprocess'] = sp * v1
            row['voidage_process'] = sp * v2
        table[node] = row
        for child in node.children:
            _walk(child)

    _walk(tree)
    return table


'''
=====================================================================================
Coverage By Aligned Duration

Uses the real elapsed time an alignment's moves are responsible for (per
Definition [Move Durations]) to estimate a submodel's mass, rather than a
move count or structural weight. \\voidsat is \\voidat computed over
skip-alignments' own lumped optimal alignments - see voidsat below; a
classical-alignment \\covat/\\voidat (not yet requested) would instead need
paths translated the way voidmass_pn._to_alignment_mass_path does for
\\covermove.

consumes/nxt/block/mdur implement Definition [Move Durations] exactly,
including its two corrections found this session against the original
draft: the move-weights reference was removed (a lumped skip move over an
entirely-absent subprocess is simply the only non-silent move in its own
block, so it already takes the whole gap - no separate weighting needed),
and mdur's zero-guard is on nxt(agn,j) being the FIRST element of
consumes(agn) (no preceding consumed event to measure a gap from), not on
the literal alignment position j=1 - several leading model-only moves
before the first-ever consumed log event all resolve to that same first
element and must all be zeroed, not just position 1.

adur/admass implement Definition [Aligned Duration Mass]. admass iterates
every REAL TRACE in the log, not deduplicated variants like
alignment_mass/voidmass_terms elsewhere in this module - duration is a
per-instance quantity: two traces sharing one activity-sequence variant
(and hence the same alignment(s), reused here exactly like everywhere
else) can still have entirely different real elapsed time. covat/voidat
implement Definition [Coverage and Voidage by Aligned Duration].
'''


def consumes(path):
    '''
    Positions (0-indexed) in `path` that consume a log event: log or
    synchronous moves, ie log_elem is a real event rather than the
    nomove symbol '>>'. Definition [Move Durations]' consumes(agn) - the
    k-th smallest position there (1-indexed) corresponds to sigma[k];
    here, the k-th element of this list (0-indexed) corresponds to
    trace[k] where trace is 0-indexed too.
    '''
    return [j for j, (log_elem, _model_elem) in enumerate(path) if log_elem != '>>']


def nxt(path, j, consumed=None):
    '''
    Smallest position >= j in `path` that consumes a log event, or None
    if there is none - Definition [Move Durations]' nxt(agn,j). Pass
    consumed (consumes(path), precomputed) when calling this repeatedly
    over the same path, as block/mdur both do.
    '''
    if consumed is None:
        consumed = consumes(path)
    idx = bisect.bisect_left(consumed, j)
    return consumed[idx] if idx < len(consumed) else None


def block(path, j, consumed=None):
    '''
    Positions sharing the same nxt(path,j) as j, excluding pure silent
    (tau) moves - Definition [Move Durations]' block(agn,j). Computed as
    the contiguous run (prev_consumed, target] rather than scanning
    nxt() for every position in `path`: every position in that range
    shares the same next-consumed position (target) by construction,
    since none of them is itself a consuming position before target.
    '''
    if consumed is None:
        consumed = consumes(path)
    target = nxt(path, j, consumed)
    if target is None:
        return []
    idx = bisect.bisect_left(consumed, target)
    prev = consumed[idx - 1] if idx > 0 else -1
    return [p for p in range(prev + 1, target + 1)
            if _classify_move(path[p][1])[1] != 'tau']


def mdur(path, trace, j, consumed=None):
    '''
    Definition [Move Durations]' mdur(agn,sigma,j) - `trace` is a real
    event sequence (each event exposing ['time:timestamp'], the same
    shape coverage_by_duration's log/trace already use), 0-indexed.
    Zero where path[j] is a silent (tau) move, where nxt(path,j) is
    undefined (nothing left to consume after j), or where nxt(path,j)
    is the FIRST element of consumes(path) (no preceding consumed event
    to measure a gap from - not merely j==0, see this section's module
    docstring on this specific correction).
    '''
    if consumed is None:
        consumed = consumes(path)
    _log_elem, model_elem = path[j]
    if _classify_move(model_elem)[1] == 'tau':
        return 0
    target = nxt(path, j, consumed)
    if target is None:
        return 0
    k = bisect.bisect_left(consumed, target)
    if k == 0:
        return 0
    gap = (trace[k]['time:timestamp'] - trace[k - 1]['time:timestamp']).total_seconds()
    blk = block(path, j, consumed)
    return gap / len(blk) if blk else 0


def _relevant_positions_by_node(path, tree, ancestor_chains=None, implied_sets=None):
    '''
    {node: [positions]} for every node in tree - the union of
    exec(m,node,path)'s members, for every node at once, WITHOUT
    grouping into discrete executions (unlike executions_by_node).
    adur has no per-execution averaging (a flat sum over every relevant
    position, see this section's module docstring), so the grouped
    structure executions_by_node builds is unneeded work here; this
    shares only the structural ancestor_chains/implied_sets maps with
    it (see make_executions_cache), not its grouping step.
    '''
    if ancestor_chains is None:
        ancestor_chains = _ancestor_chains(tree)
    if implied_sets is None:
        implied_sets = _implied_descendant_sets(tree)

    positions_by_node = {}
    for i, (_log_elem, model_elem) in enumerate(path):
        node, kind = _classify_move(model_elem)
        if node is None:
            continue
        for pt in ancestor_chains[node]:
            positions_by_node.setdefault(pt, []).append(i)
        if kind in ('skip', 'tau'):
            for pt in implied_sets[node]:
                if pt is node:
                    continue
                positions_by_node.setdefault(pt, []).append(i)
    return positions_by_node


def make_aligned_duration_cache(tree):
    '''
    Structural maps for `tree` (ancestor_chains, implied_sets) plus a
    per-path memo of _relevant_positions_by_node's own output - shared
    across many adur/admass/covat/voidat/voidsat calls against
    different nodes of the SAME tree/log (see
    lab.exp_disco_degrade._node_rows, which calls voidsat once per
    node), same purpose and shape as coveragemass.make_executions_cache
    for alignment_mass. Pass the SAME cache object to every one of
    those calls so the structural maps and each alignment's relevant-
    positions computation (_relevant_positions_by_node already computes
    every node's positions from one pass over a path) are done once per
    report row, not once per node.
    '''
    return {
        'tree': tree,
        'ancestor_chains': _ancestor_chains(tree),
        'implied_sets': _implied_descendant_sets(tree),
        'by_path': {},
    }


def _relevant_positions_for(path, cache):
    by_node = cache['by_path'].get(id(path))
    if by_node is None:
        by_node = _relevant_positions_by_node(path, cache['tree'], cache['ancestor_chains'],
                                                cache['implied_sets'])
        cache['by_path'][id(path)] = by_node
    return by_node


def adur(pt, tree, alignments, trace, cache=None):
    '''
    Definition [Aligned Duration Mass]'s adur(m,msub,sigma) for one
    trace paired with its variant's Gamma_sigma (`alignments` - a list
    of (log_elem, model_elem) paths, eg skip-alignments' own
    State.path objects for voidsat). (1/|Gamma_sigma|) times the sum,
    over every alignment and every one of pt's relevant positions in
    it, of mdur - equivalently, grouping into exec(m,pt,path)'s
    discrete executions first and summing mdur over every position in
    every execution, since adur has no per-execution weighting to make
    that grouping matter (see _relevant_positions_by_node).

    cache: optional, from make_aligned_duration_cache(tree) - pass a
    shared cache when calling this for many nodes of the same tree over
    the same alignments (the usual per-node report-row case) to avoid
    re-deriving each alignment's relevant positions once per node.
    '''
    if not alignments:
        return 0.0
    total = 0.0
    for path in alignments:
        consumed = consumes(path)
        if cache is not None:
            positions_by_node = _relevant_positions_for(path, cache)
        else:
            positions_by_node = _relevant_positions_by_node(path, tree)
        for j in positions_by_node.get(pt, []):
            total += mdur(path, trace, j, consumed)
    return total / len(alignments)


def admass(pt, tree, log, alignments_by_variant, cache=None):
    '''
    Definition [Aligned Duration Mass]'s admass(m,msub,L). `log` is a
    real event log/list of traces (log_to_traces(log)'s shape, same as
    coverage_by_duration's own `log` parameter - each event exposing
    ['concept:name']/['time:timestamp']); `alignments_by_variant` maps
    a variant key (_variant_key of that trace's activity sequence) to
    its list of optimal-alignment paths, eg
    {k: [s.path for s in v] for k, v in dv.skip_dict_backup.items()}
    for voidsat (skip-alignments' own State.path objects already use
    this module's (log_elem, model_elem) wrapper shape directly - no
    translation needed, unlike voidmass_pn.py's classical-alignment
    path).

    Iterates every trace in `log`, not deduplicated variants like
    alignment_mass/voidmass_terms elsewhere in this module - duration
    is a per-instance quantity, so two traces sharing one variant (and
    hence the same alignment(s)) can still have different real elapsed
    time. A trace whose adur is 0 is excluded from both the sum and the
    count (Definition [Aligned Duration Mass]'s L' - "traces in which
    msub has observable duration"), not just given a zero numerator -
    0.0 overall if no trace has any observable duration for pt.

    cache: optional, see adur - shared across the traces here too, so
    the SAME alignment reused by several traces of one variant only has
    its relevant positions computed once.
    '''
    total = 0.0
    n = 0
    for trace in log_to_traces(log):
        activities = tuple(e['concept:name'] for e in trace)
        alignments = alignments_by_variant.get(_variant_key(activities), [])
        d = adur(pt, tree, alignments, trace, cache)
        if d > 0:
            # d > 0 implies duration(trace) > 0: d is a sum of mdur
            # terms, each a non-negative fraction of some internal gap
            # of trace, and those gaps telescope to at most trace's own
            # total duration (see coverage_by_duration's Bounded Mass
            # lemma for the same argument) - so a positive d can never
            # be paired with a zero denominator here.
            total += d / dur(trace)
            n += 1
    return total / n if n else 0.0


def covat(pt, tree, skip_prob, log, alignments_by_variant, cache=None):
    '''Definition [Coverage and Voidage by Aligned Duration]'s covat(m,msub,L).'''
    return (1 - skip_prob) * admass(pt, tree, log, alignments_by_variant, cache)


def voidat(pt, tree, skip_prob, log, alignments_by_variant, cache=None):
    '''Definition [Coverage and Voidage by Aligned Duration]'s voidat(m,msub,L).'''
    return skip_prob * admass(pt, tree, log, alignments_by_variant, cache)


def voidsat(pt, tree, dv, log, cache=None):
    '''
    \\voidsat - voidat computed over skip-alignments' own lumped optimal
    alignments (dv.skip_dict_backup), reusing their State.path objects
    directly (already in this module's wrapper shape - no translation
    needed, unlike voidmass_pn.py's classical-alignment path). dv is a
    computed DerivationPipeline, same as coverage_by_alignment/
    coverage_by_duration's own dv parameter.
    '''
    alignments_by_variant = {k: [s.path for s in v] for k, v in dv.skip_dict_backup.items()}
    return voidat(pt, tree, dv.skip_probs[pt], log, alignments_by_variant, cache)

