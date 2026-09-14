'''
=====================================================================================
Coverage and Void by Skip-Weighted Alignment Correspondence

Definitions [Coverage by Skip-Weighted Alignment Correspondence] and
[Void by Skip-Weighted Alignment Correspondence].

Skip alignments lift a deviation to the highest applicable block, so a
subprocess absent from a trace is explained by a single skip move
however many activities it contains. Counting moves uniformly over skip
alignments (coveragemass.matchcount/movecount) therefore underweights
exactly the case of interest. smovecount restores the missing magnitude
by counting each skip move as the execution it stands for:
aligncost(<>, msub'), the least number of labelled activities any
traversal of msub' performs (coveragemass.min_activity_count). A skip
over a wholly silent subprocess contributes nothing, so no separate
exclusion of silent moves is needed.

Everything else is Coverage by Alignment Correspondence's arithmetic
unchanged, conditioned on observation at both levels: an execution with
no synchronous move is excluded, a trace whose alignments hold no such
execution leaves both the average and the W that normalises it, and each
observing alignment carries 1/|Upsilon_sigma| of its trace's weight
rather than 1/|O_sigma|. That is coveragemass.observed_alignment_mass,
which this calls with the skip-weighted ratio in place of
matchcount/movecount - the only term that differs. Coverage is 0 where
W is 0, so a subprocess observed nowhere reads void 1 rather than
collapsing to 0 the way an unconditioned ratio does.

skip_dict is a variant-key -> list of skip-alignment State (the shape
coveragemass.alignment_mass takes, eg dv.skip_dict_backup);
variant_probs is a variant tuple -> weight (eg dv.pl); skip_probs is a
ProcessTree node -> P_skip (eg dv.skip_probs).
'''

from skipalignments import Skip

from process_voids.coveragemass import matchcount, min_activity_count, observed_alignment_mass


def smatchcount(execution):
    '''
    Synchronous moves in `execution` - identical to
    coveragemass.matchcount, named per the definition's own
    skip-weighted move counts.
    '''
    return matchcount(execution)


def smovecount(execution):
    '''
    smatchcount(execution) plus aligncost(<>, msub') for every skip move
    over a subprocess msub' - the minimum number of labelled activities
    any traversal of that subprocess performs. Unlike
    coveragemass.movecount, a lumped skip move over a large subtree
    counts for its whole minimum size rather than as a single move; a
    skip over a wholly silent subprocess counts 0.
    '''
    total = smatchcount(execution)
    for _, model_elem in execution:
        if isinstance(model_elem, Skip):
            total += min_activity_count(model_elem.node)
    return total


def _skip_weighted_ratio(execution):
    return smatchcount(execution) / smovecount(execution)


def observed_skip_weighted_mass(pt, skip_dict, variant_probs,
                                 executions_cache=None, timed_out_ratio=None):
    '''
    The mass term of Definition [Coverage by Skip-Weighted Alignment
    Correspondence] - observed_alignment_mass's conditioning and
    weighting, over smatchcount/smovecount. 0 where W is 0.

    executions_cache: optional, from
    coveragemass.make_executions_cache(tree) - pass one shared cache
    across every node scored in a report row. timed_out_ratio: see
    coveragemass.alignment_mass.
    '''
    return observed_alignment_mass(pt, skip_dict, variant_probs,
                                   executions_cache=executions_cache,
                                   timed_out_ratio=timed_out_ratio,
                                   ratio=_skip_weighted_ratio)


def coversalign2(pt, skip_dict, variant_probs, skip_probs,
                  executions_cache=None, timed_out_ratio=None):
    '''
    \\coversalign (defn:salign-coverage): (1 - skipprob(pt)) times the
    observed skip-weighted mass. The two factors answer separate
    questions - whether pt was recorded at all, and how completely where
    it was - so coverage falls linearly in a subprocess's absence rather
    than quadratically.
    '''
    mass = observed_skip_weighted_mass(pt, skip_dict, variant_probs,
                                       executions_cache, timed_out_ratio)
    return (1 - skip_probs[pt]) * mass


def voidsalign2(pt, skip_dict, variant_probs, skip_probs,
                 executions_cache=None, timed_out_ratio=None):
    '''\\voidsalign (defn:salign-void): 1 - coversalign2.'''
    return 1 - coversalign2(pt, skip_dict, variant_probs, skip_probs,
                            executions_cache, timed_out_ratio)
