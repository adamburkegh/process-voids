'''
Synthetic claims-handling fixture (see claims-fixture-spec.md, provided
by the user): a mid-sized (~30 trace) log/tree pair generated from a
stated weighted process tree with per-activity duration distributions,
so ground truth is known on both control flow and timing.

Model: seq(register, and(assess, loop(request_docs, receive_docs)),
           xor(seq(lodge_appeal, decide_appeal), tau), close)

Bigger than the six-trace payment running example (lab.fixtures), small
enough that the alignment pipeline runs in seconds. Unlike the payment
example, this log has a non-zero baseline deficit by construction (a
few injected deviating traces - see claims-fixture-spec.md's "Baseline
deficit" section) and one deliberately optional subprocess (the appeal)
as a specificity negative control.

Run this module directly to (re)write the checked-in fixtures at
data/claims.xes and data/claims_ground_truth.csv.
'''

import math
import random
from pathlib import Path

import pandas as pd
from skipalignments import Activity, And, Loop, Sequence, Tau, Xor

from process_voids import dtlog

ACTIVITY_COST = 100000
SEED = 42

CLAIMS_XES = 'data/claims.xes'
CLAIMS_GROUND_TRUTH_CSV = 'data/claims_ground_truth.csv'

ACTIVITY_ALPHABET = {
    'register', 'assess', 'request_docs', 'receive_docs',
    'lodge_appeal', 'decide_appeal', 'close',
}

LOOP_CONTINUE_PROB = 0.4
APPEAL_PROB = 0.25
EXTRA_ACTIVITY_LABEL = 'call_customer'

# activity -> (median_hours, tight_sigma, lognormal_sigma). "constant"
# duration_mode uses tight_sigma for everything, so ablation math stays
# exactly reasonable ("a missing loop block costs exactly N hours");
# "lognormal" mode uses the wider sigma, notably receive_docs's heavy
# tail - the "legitimately slow occurrence" false-positive test.
_DURATIONS = {
    'register': (0.5, 0.15, 0.15),
    'assess': (2.0, 0.15, 0.15),
    'request_docs': (1.0, 0.15, 0.15),
    'receive_docs': (24.0, 0.15, 0.7),
    'lodge_appeal': (1.0, 0.15, 0.15),
    'decide_appeal': (4.0, 0.15, 0.15),
    'close': (1.0, 0.15, 0.15),
}


def build_claims_tree():
    register = Activity(None, 'register', ACTIVITY_COST)
    register.id = '1'
    assess = Activity(None, 'assess', ACTIVITY_COST)
    assess.id = '2'
    request_docs = Activity(None, 'request_docs', ACTIVITY_COST)
    request_docs.id = '3'
    receive_docs = Activity(None, 'receive_docs', ACTIVITY_COST)
    receive_docs.id = '4'
    lodge_appeal = Activity(None, 'lodge_appeal', ACTIVITY_COST)
    lodge_appeal.id = '5'
    decide_appeal = Activity(None, 'decide_appeal', ACTIVITY_COST)
    decide_appeal.id = '6'
    tau = Tau(None, 'skip-appeal', 0)
    tau.id = '7'
    close = Activity(None, 'close', ACTIVITY_COST)
    close.id = '8'

    doc_loop = Loop(None, [request_docs, receive_docs])
    doc_loop.id = '9'
    request_docs.set_parent(doc_loop)
    receive_docs.set_parent(doc_loop)

    parallel = And(None, [assess, doc_loop])
    parallel.id = '10'
    assess.set_parent(parallel)
    doc_loop.set_parent(parallel)

    appeal_seq = Sequence(None, [lodge_appeal, decide_appeal])
    appeal_seq.id = '11'
    lodge_appeal.set_parent(appeal_seq)
    decide_appeal.set_parent(appeal_seq)

    appeal_choice = Xor(None, [appeal_seq, tau])
    appeal_choice.id = '12'
    appeal_seq.set_parent(appeal_choice)
    tau.set_parent(appeal_choice)

    tree = Sequence(None, [register, parallel, appeal_choice, close])
    tree.id = '13'
    register.set_parent(tree)
    parallel.set_parent(tree)
    appeal_choice.set_parent(tree)
    close.set_parent(tree)
    return tree


def _sample_duration(activity, rng, duration_mode):
    median, tight_sigma, lognormal_sigma = _DURATIONS[activity]
    sigma = tight_sigma if duration_mode == 'constant' else lognormal_sigma
    return rng.lognormvariate(math.log(median), sigma)


def _sample_loop_iterations(rng, p_continue=LOOP_CONTINUE_PROB):
    '''Number of times request_docs (the do-child) fires - do always
    fires once, then continues with p_continue after each iteration.'''
    n = 1
    while rng.random() < p_continue:
        n += 1
    return n


def _doc_loop_events(n_iterations):
    '''request_docs, receive_docs, request_docs, ..., request_docs -
    n_iterations request_docs, n_iterations - 1 receive_docs.'''
    events = ['request_docs']
    for _ in range(n_iterations - 1):
        events.append('receive_docs')
        events.append('request_docs')
    return events

def _interleave(assess_activity, doc_events, rng):
    '''Uniformly random insertion position for the single-event assess
    branch among doc_events' fixed internal order - the only degree of
    freedom in a valid interleaving here, since assess is a singleton.'''
    position = rng.randint(0, len(doc_events))
    return doc_events[:position] + [assess_activity] + doc_events[position:]


def _sample_trace(rng, duration_mode):
    n_iterations = _sample_loop_iterations(rng)
    doc_events = _doc_loop_events(n_iterations)
    parallel_events = _interleave('assess', doc_events, rng)
    appeal = rng.random() < APPEAL_PROB

    activities = ['register'] + parallel_events
    if appeal:
        activities += ['lodge_appeal', 'decide_appeal']
    activities += ['close']

    durations = [_sample_duration(a, rng, duration_mode) for a in activities]
    return activities, durations, n_iterations, appeal


def _apply_deviation(activities, durations, kind, rng):
    activities = list(activities)
    durations = list(durations)
    if kind == 'reordered':
        # Deliberately swap lodge_appeal/decide_appeal, not a random
        # adjacent pair: appeal_seq is the one strict Sequence in the
        # model, so this is the only place a reordering is guaranteed
        # to be a genuine deviation rather than a free reinterleaving
        # (a random position could land in the order-tolerant
        # parallel/loop region and cost nothing - see session notes).
        # Callers must only pass a trace that already has an appeal.
        i = activities.index('lodge_appeal')
        j = activities.index('decide_appeal')
        activities[i], activities[j] = activities[j], activities[i]
        durations[i], durations[j] = durations[j], durations[i]
        label = 'reordered:lodge_appeal<->decide_appeal'
    elif kind == 'extra_activity':
        i = rng.randint(1, len(activities) - 1)
        activities.insert(i, EXTRA_ACTIVITY_LABEL)
        durations.insert(i, rng.lognormvariate(math.log(1.0), 0.15))
        label = f'extra_activity:{EXTRA_ACTIVITY_LABEL}'
    else:
        raise ValueError(f'unknown deviation kind {kind!r}')
    return activities, durations, label


def generate_claims_log(n_traces=30, duration_mode='constant', seed=SEED,
                         n_deviated=2):
    '''
    (log, ground_truth) - log is a pm4py-format event log DataFrame
    (see process_voids.dtlog), ground_truth a DataFrame with one row per
    case: which path was taken (loop_iterations, appeal), and whether
    (and how) the trace was deviated from the model.
    '''
    if n_deviated not in (0, 1, 2, 3):
        raise ValueError('n_deviated should be 0-3 per the fixture spec')

    rng = random.Random(seed)
    case_names = [f'claim_{i + 1}' for i in range(n_traces)]

    samples = [list(_sample_trace(rng, duration_mode)) for _ in range(n_traces)]
    deviation_labels = [''] * n_traces

    # Deliberate deviation targets, not blind random indices (see
    # _apply_deviation's docstring): 'reordered' must land on a trace
    # that already has an appeal, since it swaps the appeal pair -
    # picking a case at random and hoping it qualifies would make the
    # fixture's baseline deficit depend on incidental rng-stream luck.
    if n_deviated >= 1:
        appeal_indices = [i for i, s in enumerate(samples) if s[3]]
        if not appeal_indices:
            raise RuntimeError(
                'no traces with an appeal to deviate - increase n_traces or APPEAL_PROB')
        i = rng.choice(appeal_indices)
        activities, durations, deviation_labels[i] = _apply_deviation(
            samples[i][0], samples[i][1], 'reordered', rng)
        samples[i][0], samples[i][1] = activities, durations

    remaining_kinds = (['extra_activity'] * (n_deviated - 1))
    for kind in remaining_kinds:
        candidates = [i for i in range(n_traces) if deviation_labels[i] == '']
        i = rng.choice(candidates)
        activities, durations, deviation_labels[i] = _apply_deviation(
            samples[i][0], samples[i][1], kind, rng)
        samples[i][0], samples[i][1] = activities, durations

    trace_strings = []
    ground_truth_rows = []
    for i, case in enumerate(case_names):
        activities, durations, n_iterations, appeal = samples[i]
        deviation_label = deviation_labels[i]

        offsets = []
        t = 0.0
        for d in durations:
            t += d
            offsets.append(t)
        trace_strings.append(' '.join(f'{a}:{o}' for a, o in zip(activities, offsets)))

        ground_truth_rows.append({
            'case': case,
            'loop_iterations': n_iterations,
            'appeal': appeal,
            'deviation': deviation_label,
        })

    log = dtlog.convert_timed(*trace_strings, names=case_names, time_unit='hours')
    ground_truth = pd.DataFrame(ground_truth_rows)
    return log, ground_truth


def summarize_ground_truth(ground_truth):
    '''
    Human-readable summary of a generate_claims_log ground-truth table:
    loop-iteration distribution, appeal count, and every deviated row -
    the "did the fixture come out the way I expect" check, without
    needing a one-off script every time the fixture is regenerated.
    '''
    lines = [
        f'n_traces: {len(ground_truth)}',
        f'appeal: {int(ground_truth["appeal"].sum())} of {len(ground_truth)}',
        'loop_iterations distribution:',
        ground_truth['loop_iterations'].value_counts().sort_index().to_string(),
        'deviated traces:',
        ground_truth.loc[ground_truth['deviation'] != '',
                          ['case', 'loop_iterations', 'appeal', 'deviation']].to_string(index=False),
    ]
    return '\n'.join(lines)


if __name__ == '__main__':
    # 60, not the function's own default of 30: gives the appeal_seq
    # subprocess a bigger eligible-case pool for degradation sweeps
    # (exp_claims_degrade.py), so a handful of ablated cases is a
    # smaller fraction of it - exp_claims_degrade also deliberately
    # excludes the deviated case from ablation eligibility either way
    # (see degrade_target_subprocess's exclude_cases), so this isn't
    # load-bearing for correctness, just headroom.
    log, ground_truth = generate_claims_log(n_traces=60)
    dtlog.write_xes(log, CLAIMS_XES)
    Path(CLAIMS_GROUND_TRUTH_CSV).parent.mkdir(parents=True, exist_ok=True)
    ground_truth.to_csv(CLAIMS_GROUND_TRUTH_CSV, index=False)
    print(f'Wrote {CLAIMS_XES} and {CLAIMS_GROUND_TRUTH_CSV}')
    print()
    print(summarize_ground_truth(ground_truth))
