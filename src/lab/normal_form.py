'''
=====================================================================================
Normal-form move precedence

A skip alignment is in NORMAL FORM when moves that could be transposed
follow a fixed precedence - log, then synchronous, then skip (Definition
of Skip Alignment, Bar et al. BPM 2025). skip-alignments returns that one
representative of each class of transposable alignments, so the
precedence, rather than the evidence, decides which member a metric
reads.

Two adjacent moves are INDEPENDENT when transposing them yields another
valid skip alignment of the same trace against the same model:

    log / skip                  always - a log move belongs to no
                                execution, and a skip consumes no event
    sync / skip, skip / skip    when both orderings of the model path are
                                in the skip language, which for a process
                                tree means the two nodes sit in different
                                branches of a common And
    log / sync, log / log       when the two carry the same activity

Two moves that both consume a log event and carry different activities
are never independent: the trace pins their order.

This module reorders an alignment the aligner has already returned:
nothing here reaches into skip-alignments' own search or cost function.

Invariance under the precedence is a property of skip-alignments, whose
rule it is, so that library is this code's strategic home. It lives here
because process-voids cannot modify a sibling repo, which left
post-processing the returned alignment as the only route available - a
boundary, not a judgement that the harness is where it belongs.

The concurrency test is structural: a transposition of two model-side
moves is admitted only where the nodes' lowest common ancestor is an And
and they descend through different children of it. Sequence and Loop fix
the order of their children, and two moves in different branches of an
Xor cannot both appear in one alignment, so And is the only operator
that makes both orderings valid.
'''

from skipalignments import And, Skip, TauPath

LOG_FIRST = ('log', 'sync', 'skip')
SKIP_FIRST = ('skip', 'sync', 'log')

PRECEDENCES = {'log_first': LOG_FIRST, 'skip_first': SKIP_FIRST}


def move_kind(move):
    '''
    'log', 'sync' or 'skip' for one (log_elem, model_elem) pair. A
    silent TauPath counts as a skip: skip alignments have no silent move
    type, so it is a deviation like any other (see coveragemass.block).
    '''
    log_elem, model_elem = move
    if model_elem == '>>':
        return 'log'
    if isinstance(model_elem, (Skip, TauPath)):
        return 'skip'
    return 'sync'


def move_node(move):
    '''The tree node a move is classified at, or None for a log move.'''
    _log_elem, model_elem = move
    if model_elem == '>>':
        return None
    if isinstance(model_elem, (Skip, TauPath)):
        return model_elem.node
    return model_elem


def _activity(move):
    '''The activity label a move consumes, or None where it consumes none.'''
    log_elem, _model_elem = move
    if log_elem == '>>':
        return None
    return log_elem.rstrip('0123456789')


def _ancestors(node):
    chain = []
    while node is not None:
        chain.append(node)
        node = node.parent
    return chain


def _concurrent(first, second):
    '''
    Whether two nodes sit in different branches of a common And, so both
    orderings of their moves are in the skip language.
    '''
    if first is None or second is None or first is second:
        return False
    first_chain = _ancestors(first)
    second_chain = _ancestors(second)
    second_index = {id(node): i for i, node in enumerate(second_chain)}
    for i, node in enumerate(first_chain):
        j = second_index.get(id(node))
        if j is None:
            continue
        # node is the lowest common ancestor; the two moves reach it
        # through first_chain[i-1] and second_chain[j-1].
        return isinstance(node, And) and i > 0 and j > 0 and first_chain[i - 1] is not second_chain[j - 1]
    return False


def independent(first, second, tree=None):
    '''
    Whether the adjacent moves `first`, `second` may be transposed - the
    independence table in this module's docstring. `tree` is accepted for
    call-site clarity and future use; independence is decided from the
    moves' own nodes, which carry their parents.
    '''
    first_kind, second_kind = move_kind(first), move_kind(second)
    first_activity, second_activity = _activity(first), _activity(second)

    if first_activity is not None and second_activity is not None:
        # both consume an event: only the same label commutes
        return first_activity == second_activity
    if first_kind == 'log' or second_kind == 'log':
        # a log move against a skip - the remaining case, always free
        return True
    return _concurrent(move_node(first), move_node(second))


def reorder(path, tree, precedence=LOG_FIRST):
    '''
    `path` rewritten so that independent adjacent moves follow
    `precedence` (LOG_FIRST reproduces skip-alignments' own normal form;
    SKIP_FIRST is its reverse).

    A bubble pass over independent adjacent pairs only, so the result is
    reachable from `path` by valid transpositions and is still an optimal
    alignment of the same trace. Repeated until no admitted swap is left,
    which terminates because every swap strictly reduces the number of
    out-of-precedence independent pairs.
    '''
    rank = {kind: i for i, kind in enumerate(precedence)}
    moves = list(path)
    swapped = True
    while swapped:
        swapped = False
        for i in range(len(moves) - 1):
            first, second = moves[i], moves[i + 1]
            if rank[move_kind(first)] <= rank[move_kind(second)]:
                continue
            if not independent(first, second):
                continue
            moves[i], moves[i + 1] = second, first
            swapped = True
    return moves
