
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


