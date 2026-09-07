'''
Interval surprise: a coverage-void metric needing no alignments and no
tunable parameters, based on the empirical tail probability of each
event's inter-event interval.

    P(D_a >= d) = |{d' in D_a : d' >= d}| / (n + 1)      (plus-one empirical p-value)
    surprise    = -log2 P(D_a >= d)                       (bits)

Two attribution schemes for charging an event's surprise to tree nodes,
both additive (a composite node's total is the sum of its children's -
see surprise_totals/predecessor_totals below):

  containment - charge to every node whose alphabet contains the event's
                activity (leaf and all ancestors). Puts the signal on the
                node that ran late, which for a *missing* subprocess is
                the wrong node - the gap inflates the *following*
                activity's interval, so the signal appears downstream.
  predecessor - charge to the subprocess the model says should have run
                immediately before the event, then its ancestors. Moves
                the signal onto the actually-missing subprocess.

See docs/void-entropy-brief-v2.md style writeup (not checked in here) for
the full design rationale and open questions.
'''

import math
from bisect import bisect_left
from collections import defaultdict
from functools import reduce

from skipalignments import Activity, Tau, Sequence, Xor, And, Loop, ProcessTree


def event_intervals(trace):
    '''
    [(i, activity, delta_seconds), ...] for i >= 1 - the interval since
    the previous event in the trace. The first event (i == 0) owns no
    interval and is excluded, matching the paper's convention that
    waiting time is charged to the following activity.
    '''
    return [(i, trace[i]['concept:name'],
              (trace[i]['time:timestamp'] - trace[i - 1]['time:timestamp']).total_seconds())
            for i in range(1, len(trace))]


def observed_intervals(traces):
    '''
    activity -> [delta, ...] over the whole log (traces, not variants),
    sorted ascending - tail_probability requires this, so it can locate
    the tail by bisection (O(log n)) rather than a linear scan.
    '''
    obs = defaultdict(list)
    for trace in traces:
        for _, activity, delta in event_intervals(trace):
            obs[activity].append(delta)
    return {activity: sorted(deltas) for activity, deltas in obs.items()}


def tail_probability(vals, d):
    '''
    P(D >= d), plus-one denominator so the maximum observation is finite.
    vals must be sorted ascending (observed_intervals returns its lists
    pre-sorted for exactly this reason) - this is called once per event,
    and "cheap, no alignments" is the metric's whole argument, so it
    locates the tail by bisection rather than an O(n) linear scan.

    Floored at 1/(n+1): the same probability the in-sample maximum gets,
    rather than letting it hit exactly 0 (and -log2(0) blow up) for a
    value more extreme than anything vals ever observed. Unreachable
    when vals is an activity's own self-estimated distribution (its own
    value is always already in the sample), but real when scoring
    against a different reference distribution - e.g. event_surprise's
    baseline-obs path, where a degraded log's event can legitimately
    exceed the undegraded log's historical range for that activity.
    '''
    n = len(vals)
    count = n - bisect_left(vals, d)
    return max(count, 1) / (n + 1)


def event_surprise(traces, obs=None):
    '''
    [(activity, delta, p, bits), ...] for every event with a predecessor
    in its trace. obs (activity -> observed deltas) is computed from
    `traces` if not supplied - pass it in explicitly to reuse the same
    distribution across many calls (e.g. once per ablation level).
    '''
    if obs is None:
        obs = observed_intervals(traces)
    rows = []
    for trace in traces:
        for _, activity, delta in event_intervals(trace):
            p = tail_probability(obs[activity], delta)
            rows.append((activity, delta, p, -math.log2(p)))
    return rows


def _by_activity(rows):
    totals = defaultdict(float)
    for activity, _, _, bits in rows:
        totals[activity] += bits
    return totals


def _leaves_by_name(tree):
    by_name = defaultdict(list)
    for leaf in tree.get_leafs():
        if isinstance(leaf, Activity):
            by_name[leaf.name].append(leaf)
    return by_name


def _lca_of(nodes):
    return reduce(ProcessTree.get_lca, nodes)


def surprise_totals(tree, rows):
    '''
    node -> total bits under containment attribution (see module
    docstring). Additive bottom-up sum over the tree, keyed by node
    identity (matching the skip_probs dict convention elsewhere in this
    project).

    A duplicated activity label (the same name on more than one leaf)
    can't be split among its leaves without knowing which tree position
    produced which event - its full total is instead charged to the
    leaves' lowest common ancestor (and propagated up from there); the
    leaves themselves, and any node strictly between the LCA and a
    specific leaf, get 0 from it.

    Returns (totals, out_of_alphabet_event_count). An activity with NO
    leaf anywhere in the tree (eg noise-pruned out of a discovered
    model) would otherwise have its bits silently excluded from every
    node's total, including the root - breaking the invariant that the
    root always holds the log's full surprise total regardless of the
    model. Those bits are charged directly to the root instead, and the
    affected event count is returned as a diagnostic (analogous to
    predecessor_totals' unattributable_event_count, but for "missing
    from the model" rather than "no predecessor").
    '''
    by_activity = _by_activity(rows)
    leaves_by_name = _leaves_by_name(tree)
    duplicated_names = {name for name, leaves in leaves_by_name.items() if len(leaves) > 1}
    out_of_alphabet_names = set(by_activity) - set(leaves_by_name)
    totals = {}

    def _walk(node):
        if isinstance(node, Activity):
            val = 0.0 if node.name in duplicated_names else by_activity.get(node.name, 0.0)
        elif isinstance(node, Tau):
            val = 0.0
        else:
            val = sum(_walk(child) for child in node.children)
        totals[node] = val
        return val

    _walk(tree)

    for name in duplicated_names:
        lca = _lca_of(leaves_by_name[name])
        bits = by_activity.get(name, 0.0)
        node = lca
        while node is not None:
            totals[node] = totals.get(node, 0.0) + bits
            node = node.parent

    for name in out_of_alphabet_names:
        totals[tree] = totals.get(tree, 0.0) + by_activity[name]
    out_of_alphabet_event_count = sum(1 for activity, _, _, _ in rows
                                      if activity in out_of_alphabet_names)

    return totals, out_of_alphabet_event_count


def _exit_owner(node):
    '''
    (owner_node, ambiguous) - the node that structurally executes last
    within `node`'s subtree, i.e. what a following sibling's predecessor
    should be.

    Sequence: its last child, recursively. Loop: the do-child
    (children[0]) always executes last (do (redo do)*). Xor/And:
    genuinely ambiguous which branch's exit is real - charged to the
    composite node itself rather than guessing a branch. NOTE: this is
    known to be too coarse when a branch is a silent Tau (the true
    predecessor can reach back past the Xor/And entirely) - flagged for
    review, not resolved here.
    '''
    if isinstance(node, (Activity, Tau)):
        return node, False
    if isinstance(node, Sequence):
        return _exit_owner(node.children[-1])
    if isinstance(node, Loop):
        return _exit_owner(node.children[0])
    if isinstance(node, (Xor, And)):
        return node, True
    raise ValueError(f'Unrecognised process tree node: {node}')


def compute_predecessors(tree):
    '''
    leaf -> (predecessor_node_or_None, ambiguous) for every Activity/Tau
    leaf in the tree, derived structurally:

      Sequence child: the preceding sibling's exit owner (or, for the
      first child, the Sequence's own predecessor).
      Xor/And child: every child's predecessor is the composite node's
      own predecessor (entering a choice/parallel is never ambiguous).
      Loop do-child and redo-child: both allocate to the body (do-child)
      node unconditionally - the body precedes both the redo and every
      subsequent iteration, so there is no first-iteration special case
      to track (a purely structural rule, no per-trace state needed).

    predecessor_node may be an internal node (e.g. an Xor), not
    necessarily a leaf, when the exit was ambiguous.
    '''
    predecessors = {}

    def _assign(node, pred, ambiguous):
        if isinstance(node, (Activity, Tau)):
            predecessors[node] = (pred, ambiguous)
            return
        if isinstance(node, Sequence):
            cur_pred, cur_amb = pred, ambiguous
            for child in node.children:
                _assign(child, cur_pred, cur_amb)
                cur_pred, cur_amb = _exit_owner(child)
            return
        if isinstance(node, (Xor, And)):
            for child in node.children:
                _assign(child, pred, ambiguous)
            return
        if isinstance(node, Loop):
            do = node.children[0]
            redo = node.children[1]
            _assign(do, do, False)
            _assign(redo, do, False)
            return
        raise ValueError(f'Unrecognised process tree node: {node}')

    _assign(tree, None, False)
    return predecessors


def predecessor_totals(tree, rows, predecessors):
    '''
    (node -> total bits, ambiguous_event_count, unattributable_event_count)
    under predecessor attribution: each event's surprise is charged to
    the subprocess the model says should have run immediately before it,
    then propagated to every ancestor.

    ambiguous_event_count is how many events were charged via an Xor/And
    exit-ambiguity (ie the predecessor could not be pinned to one
    branch) - a diagnostic for how often that arises, not a correctness
    signal by itself.

    unattributable_event_count is how many events had no predecessor at
    all (nothing structurally precedes them - only possible for a leaf
    reachable from the root with no preceding sibling anywhere on its
    path). Their bits are charged directly to the root rather than
    dropped, so the root total always agrees with surprise_totals'
    (containment always attributes root the full total).

    A duplicated activity label (the same name on more than one leaf)
    generally has a DIFFERENT predecessor per leaf - unlike containment,
    there's no single converging ancestor to punt to - so its surprise
    is shared out evenly: bits/k charged through each of the k leaves'
    own predecessor chains. ambiguous_event_count and
    unattributable_event_count become fractional in this case (e.g. one
    of two duplicate leaves resolving via an Xor-exit ambiguity counts
    as half an ambiguous event), which is the honest reading of "how
    many events, in expectation, hit this" once an event's true leaf is
    unknown and shared uniformly across the candidates.

    Performance: the per-leaf walk to the root is the same for every
    event of a given activity, regardless of that event's own bits - so
    it's computed once per distinct activity (a per-unit-bit "pattern"
    of node -> share) and reused for every row, rather than repeated per
    event. Matters a lot when a small number of activities account for
    many events (e.g. most of a log's activities dropped by
    degradation) combined with a tree that gives that activity many
    duplicate leaves (e.g. a stochastically mined tree) - the naive
    per-event walk is O(events x duplicates x tree depth), this is
    O(activities x duplicates x tree depth) + O(events x pattern size).
    '''
    leaves_by_name = _leaves_by_name(tree)
    totals = defaultdict(float)
    ambiguous_event_count = 0.0
    unattributable_event_count = 0.0
    patterns = {}

    def _pattern_for(activity):
        cached = patterns.get(activity)
        if cached is not None:
            return cached
        leaves = leaves_by_name.get(activity, [])
        node_shares = defaultdict(float)
        ambiguous_share = 0.0
        unattributable_share = 0.0
        if leaves:
            share = 1.0 / len(leaves)
            for leaf in leaves:
                pred, ambiguous = predecessors.get(leaf, (None, False))
                if pred is None:
                    unattributable_share += share
                    node_shares[tree] += share
                    continue
                if ambiguous:
                    ambiguous_share += share
                node = pred
                while node is not None:
                    node_shares[node] += share
                    node = node.parent
        pattern = (dict(node_shares), ambiguous_share, unattributable_share)
        patterns[activity] = pattern
        return pattern

    for activity, _, _, bits in rows:
        node_shares, ambiguous_share, unattributable_share = _pattern_for(activity)
        for node, share in node_shares.items():
            totals[node] += share * bits
        ambiguous_event_count += ambiguous_share
        unattributable_event_count += unattributable_share
    return dict(totals), ambiguous_event_count, unattributable_event_count
