'''
=====================================================================================
Voidage by Skip-Weighted Alignment Moves

Skip alignments lift a deviation to the highest applicable block, so a
subprocess absent from a trace is explained by a single skip move
regardless of how many activities it contains. Counting moves uniformly
over skip alignments (as coveragemass.matchcount/movecount do)
underweights exactly that case. smovecount restores the missing
magnitude by weighting each skip move by the minimum executable length
of the subprocess it stands for (coveragemass.min_activity_count),
rather than expanding skip alignments to full alignments.

voidsalign is skipprob(msub) times a SIZE share - smovetotal(msub) over
smovetotal(the whole model) - not a match/movecount completeness ratio.
A completeness ratio is zero on a wholly-skipped execution regardless of
its skip weight (0/w = 0 for any w), so skipprob times that ratio is an
inverted U: it peaks at partial recording and returns to zero exactly
where a subprocess is entirely missing, the opposite of what a voidage
metric for missing subprocesses needs. smovetotal is a size, not a
completeness rate, so it does not collapse to zero under skipping the
way a ratio does - a skipped execution still contributes its
smovecount, so the two factors (how often skipped, how large) stay
independent evidence instead of the second cancelling under the first.

skip_dict is a variant-key -> list of skip-alignment State (the same
shape coveragemass.alignment_mass takes, eg dv.skip_dict_backup);
variant_probs is a variant tuple -> weight (eg dv.pl); skip_probs is a
ProcessTree node -> P_skip (eg dv.skip_probs).
'''

from skipalignments import Skip

from process_voids.coveragemass import _executions_for, executions, matchcount, min_activity_count


def _variant_key(variant):
    return ', '.join(variant)


def smatchcount(execution):
    '''Number of synchronous moves in `execution` - identical to
    coveragemass.matchcount, named per the paper's skip-weighted move
    counts.'''
    return matchcount(execution)


def smovecount(execution):
    '''
    smatchcount(execution) plus, for every skip move over a subprocess
    node, min_activity_count(node) - the minimum number of labelled
    activities any traversal of that subprocess performs. Unlike
    coveragemass.movecount, a lumped skip move over a large subtree
    counts for its whole minimum size rather than as a single move. A
    skip over a wholly silent subprocess contributes 0 (its minimum
    traversal performs no labelled activity), so no separate silent-
    exclusion filter is needed here, unlike movecount.
    '''
    total = smatchcount(execution)
    for _, model_elem in execution:
        if isinstance(model_elem, Skip):
            total += min_activity_count(model_elem.node)
    return total


def smovetotal(pt, skip_dict, variant_probs, executions_cache=None):
    '''
    Summed (not averaged) smovecount over every execution of `pt`, across
    every optimal skip alignment of every trace variant, weighted by
    variant probability - multiple tied optimal alignments within one
    variant share that variant's weight equally, same convention as
    coveragemass.voidmass_terms's own movecount_sum.

    executions_cache: optional, from coveragemass.make_executions_cache(tree)
    where tree is pt's tree - pass a shared cache when scoring many nodes
    of the same tree over the same skip_dict (the usual per-node report-
    row case) to avoid re-walking each alignment's path once per node;
    the SAME cache passed to alignment_mass for the same tree/skip_dict
    is shared automatically, since both key their per-path memo by path
    object identity. Omit for a one-off single-node call.
    '''
    total = 0.0
    for variant, weight in variant_probs.items():
        states = skip_dict.get(_variant_key(variant), [])
        if not states:
            continue
        share = weight / len(states)
        for state in states:
            execs = (_executions_for(state.path, pt, executions_cache)
                     if executions_cache is not None
                     else executions(state.path, pt))
            for e in execs:
                total += share * smovecount(e)
    return total


def voidsalign(pt, tree, skip_dict, variant_probs, skip_probs, executions_cache=None):
    '''
    Voidage by Skip-Weighted Alignment Moves: skipprob(pt) times pt's
    share of the whole model's skip-weighted move volume -
    smovetotal(pt) / smovetotal(tree). `tree` is the root of the same
    process tree pt belongs to - the denominator is the WHOLE model's
    expected moves, not pt's own, so this is a share of the process
    rather than a per-subprocess rate (and so does not decompose
    additively over a cut the way voidmass_process does, since the
    product with skipprob does not distribute over a sum the way the
    size term alone would).

    voidsalign = 0 where smovetotal(tree) = 0 (the whole model has no
    skip-weighted moves at all to share).
    '''
    root_total = smovetotal(tree, skip_dict, variant_probs, executions_cache)
    if root_total == 0:
        return 0.0
    return skip_probs[pt] * smovetotal(pt, skip_dict, variant_probs, executions_cache) / root_total
