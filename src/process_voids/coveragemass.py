
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
    '''
    relevant = []
    for i, (log_elem, model_elem) in enumerate(path):
        node, kind = _classify_move(model_elem)
        if node is None:
            continue
        if pt.contains_tree(node):
            relevant.append((i, log_elem, model_elem))
        elif (kind in ('skip', 'tau') and node.contains_tree(pt)
                and _mandatorily_implies(node, pt)):
            relevant.append((i, log_elem, model_elem))
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


def matchcount(execution):
    return sum(1 for log_elem, model_elem in execution
               if _classify_move(model_elem)[1] == 'sync' and log_elem != '>>')


def movecount(execution):
    return sum(1 for _, model_elem in execution
               if _classify_move(model_elem)[1] in ('sync', 'skip'))


def _variant_key(variant):
    return ', '.join(variant)


def alignment_mass(pt:ProcessTree, skip_dict:dict, variant_probs:dict,
                    convention='zero'):
    '''
    The mass term of Coverage by Alignment Correspondence: the (1-P_skip)
    factor is not applied here (see coverage_by_alignment), so this can
    be computed and tested without a real skip-probability estimation
    (no ebi dependency) - skip_dict is a variant-key -> list of skip-
    alignment State (as produced by Aligner.align2/DerivationPipeline.
    compute_skip_alignments), variant_probs is a variant tuple -> weight.

    convention: 'zero' or 'renormalised', see module docstring above.
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
            execs = [e for e in executions(state.path, pt) if movecount(e) > 0]
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


def coverage_by_alignment(pt:ProcessTree, dv, convention='zero'):
    '''
    Coverage by Alignment Correspondence: (1 - P_skip(pt)) * alignment_mass(...).
    dv is a computed DerivationPipeline (dv.pl for variant weights,
    dv.skip_dict_backup for the per-variant skip-alignments, dv.skip_probs
    for P_skip) - the full, ebi-backed pipeline. See alignment_mass for
    the ebi-free mass computation this wraps.
    '''
    mass = alignment_mass(pt, dv.skip_dict_backup, dv.pl, convention)
    return (1 - dv.skip_probs[pt]) * mass

