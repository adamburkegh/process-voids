'''
=====================================================================================
Coverage and Void by Skip Alignment Correspondence

Definitions [Coverage by Skip Alignment Correspondence] and [Void by Skip
Alignment Correspondence].

Skip alignments lift a deviation to the highest applicable block, so a
subprocess absent from a trace is explained by a single skip move
however many activities it contains. Counting moves uniformly over skip
alignments (coveragemass.matchcount/movecount) therefore underweights
exactly the case of interest. smovecount counts each skip move as the
subprocess it stands for: |leaves(msub') \\ {silent}|, the number of
labelled activities in the subprocess skipped.

The weight follows the model's own articulation of a subprocess: a
region described as three activities carries three times the weight of
one described as a single activity. aligncost(<>, msub') - the cost of a
cost-minimal execution against the empty trace - is not used, because a
subprocess with a traversal performing no labelled activity (any
optional block, or a loop whose redo-child need not run) would weigh
less than its activities, down to nothing, and its absence would be
understated or invisible in the mass.

Everything else is Coverage by Alignment Correspondence's arithmetic,
conditioned on observation at both levels: an execution with no
synchronous move is excluded, a trace whose alignments hold no such
execution leaves both the average and the normaliser obscount, and each
observing alignment carries 1/|Upsilon_sigma| of its trace's weight.
That is coveragemass.observed_alignment_mass, called with the
skip-weighted ratio. Coverage is 0 where obscount is 0, so a subprocess
observed nowhere reads void 1.

skip_dict is a variant-key -> list of skip-alignment State (the shape
coveragemass.alignment_mass takes, eg dv.skip_dict_backup);
variant_probs is a variant tuple -> weight (eg dv.pl); skip_probs is a
ProcessTree node -> P_skip (eg dv.skip_probs).
'''

from skipalignments import Skip, Tau

from process_voids.coveragemass import matchcount, observed_alignment_mass


def smatchcount(execution):
    '''Synchronous moves in `execution` - coveragemass.matchcount, named
    per the definition's skip-weighted move counts.'''
    return matchcount(execution)


def labelled_leaf_count(node):
    '''|leaves(node) \\ {silent}|: the labelled activities in node's
    subtree, counting node itself when it is a leaf.'''
    return sum(1 for leaf in node.get_leafs() if not isinstance(leaf, Tau))


def smovecount(execution):
    '''
    smatchcount(execution) plus, for every skip move, the number of
    labelled activities in the subprocess it skips. A lumped skip over a
    large subtree counts for every activity in it; a skip over a wholly
    silent subprocess counts 0.
    '''
    total = smatchcount(execution)
    for _, model_elem in execution:
        if isinstance(model_elem, Skip):
            total += labelled_leaf_count(model_elem.node)
    return total


def _skip_weighted_ratio(execution):
    return smatchcount(execution) / smovecount(execution)


def observed_skip_weighted_mass(pt, skip_dict, variant_probs,
                                 executions_cache=None, timed_out_ratio=None):
    '''
    The mass term of Definition [Coverage by Skip Alignment
    Correspondence]: observed_alignment_mass's conditioning and
    weighting, over smatchcount/smovecount. 0 where obscount is 0.

    executions_cache: optional, from coveragemass.make_executions_cache
    - pass one shared cache across every node scored in a report row.
    timed_out_ratio: see coveragemass.alignment_mass.
    '''
    return observed_alignment_mass(pt, skip_dict, variant_probs,
                                   executions_cache=executions_cache,
                                   timed_out_ratio=timed_out_ratio,
                                   ratio=_skip_weighted_ratio)


def coversalign3(pt, skip_dict, variant_probs, skip_probs,
                  executions_cache=None, timed_out_ratio=None):
    '''\\coversalign (defn:salign-coverage): (1 - skipprob(pt)) times the
    observed skip-weighted mass.'''
    mass = observed_skip_weighted_mass(pt, skip_dict, variant_probs,
                                       executions_cache, timed_out_ratio)
    return (1 - skip_probs[pt]) * mass


def voidsalign3(pt, skip_dict, variant_probs, skip_probs,
                 executions_cache=None, timed_out_ratio=None):
    '''\\voidsalign (defn:salign-void): 1 - coversalign3.'''
    return 1 - coversalign3(pt, skip_dict, variant_probs, skip_probs,
                            executions_cache, timed_out_ratio)
