'''
=====================================================================================
Coverage and Void by Aligned Duration

Definition [Coverage and Void by Aligned Duration] (defn:aligned-duration).

How much of a subprocess's elapsed time was actually recorded. Each
interval between consecutive consumed events is divided equally among
the moves preceding it (coveragemass' consumes/nxt/block/mdur), and a
subprocess's share splits by whether the move that earned it was
synchronous: obsdur for the moves that match a real event, misdur for
the moves that do not.

adratio is that split read as a rate - observed over observed plus
missing - so it answers "how much of this subprocess's time is
accounted for", not "how large a share of the process it occupies".
That is the difference from the retired voidsat, whose mass divided by
the WHOLE TRACE's duration and so reported a share. A share has no
reading at which a subprocess is wholly void: voidat multiplied its mass
by skip_prob, so an entirely missing subprocess, whose mass is 0, read
void 0 - the reverse of the truth. Here void is 1 - covsat, so it reads
1.

Conditioned on observation, exactly as coveragemass.
observed_alignment_mass is: an alignment counts only where the
subprocess has a synchronous move in it (O_sigma), a trace whose
alignments never observe the subprocess leaves both the sum and the
obscount that normalises it, and each observing alignment carries
1/|Upsilon_sigma| of its trace rather than 1/|O_sigma|, so a variant
whose tied alignments disagree about whether the subprocess was observed
contributes only the share its observing alignments carry. Mass is 0
where obscount is 0.

It does NOT call observed_alignment_mass, despite that identical
conditioning, because two things differ and neither is a parameter of
it. It iterates REAL TRACES rather than deduplicated weighted variants,
since duration is per-instance - two traces sharing one activity
sequence share an alignment but not their timings. And adratio is a
ratio of sums over all of a subprocess's positions, where
observed_alignment_mass averages a per-execution ratio.

An observed subprocess can still have no measurable duration: one
recorded at the very start of a trace has no interval bounding it, so
every mdur is 0 and obsdur + misdur is 0. adratio is 1 there, not 0 - it
was observed, and the log simply cannot time it.

alignments_by_variant maps a variant key (coveragemass._variant_key of a
trace's activity sequence) to that variant's optimal skip-alignment
paths, eg {k: [s.path for s in v] for k, v in dv.skip_dict_backup.items()};
log is a real event log in coveragemass.log_to_traces' shape; skip_probs
is a ProcessTree node -> P_skip (eg dv.skip_probs).
'''

from process_voids.coveragemass import (
    _split_duration, _variant_key, has_synchronous_move, log_to_traces)


def adratio(pt, tree, path, trace, cache=None):
    '''
    Definition [Coverage and Void by Aligned Duration]'s
    adratio(m,msub,sigma,agn) - obsdur / (obsdur + misdur).

    1 where that denominator is 0. Callers reach this only for an
    alignment in O_sigma, where the subprocess has a synchronous move,
    so a zero denominator means its moves were unmeasurable rather than
    unobserved - see the module docstring.
    '''
    observed, missing = _split_duration(pt, tree, path, trace, cache)
    total = observed + missing
    return 1.0 if total == 0 else observed / total


def admass(pt, tree, log, alignments_by_variant, cache=None):
    '''
    Definition [Coverage and Void by Aligned Duration]'s admass(m,msub,L):

        (1/obscount) * sum_sigma sum_{agn in O_sigma}
            (1/|Upsilon_sigma|) * adratio(m,msub,sigma,agn)

    where obscount = sum_sigma |O_sigma|/|Upsilon_sigma|. Zero where
    obscount is zero - a subprocess observed nowhere.

    Iterates every trace in `log`, not deduplicated variants - see the
    module docstring. A cache built with this same log
    (coveragemass.make_aligned_duration_cache(tree, log)) supplies the
    trace list; otherwise it is built here.
    '''
    if cache is not None and cache.get('traces_log') is log and cache.get('traces') is not None:
        traces = cache['traces']
    else:
        traces = log_to_traces(log)

    weighted_sum = 0.0
    obscount = 0.0
    for trace in traces:
        activities = tuple(e['concept:name'] for e in trace)
        paths = alignments_by_variant.get(_variant_key(activities), [])
        if not paths:
            continue
        share = 1.0 / len(paths)
        for path in paths:
            if not has_synchronous_move(pt, tree, path, cache):
                continue
            weighted_sum += share * adratio(pt, tree, path, trace, cache)
            obscount += share
    return weighted_sum / obscount if obscount else 0.0


def covsat2(pt, tree, log, alignments_by_variant, skip_probs, cache=None):
    '''Definition [Coverage and Void by Aligned Duration]'s covsat(m,msub,L).'''
    return (1 - skip_probs[pt]) * admass(pt, tree, log, alignments_by_variant, cache)


def voidsat2(pt, tree, log, alignments_by_variant, skip_probs, cache=None):
    '''Definition [Coverage and Void by Aligned Duration]'s voidsat(m,msub,L).'''
    return 1 - covsat2(pt, tree, log, alignments_by_variant, skip_probs, cache)
