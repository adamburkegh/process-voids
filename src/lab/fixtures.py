"""
The invoice-payment running example: an order is placed (o), then
approved (a), possibly after escalation to a manager (e); payment is
then scheduled (s) - unless the payment is urgent, in which case
scheduling is skipped - and executed (p).

Model:  N = seq( o, loop(a, e), xor(s, tau), p )

Approval, with optional escalation, executes in the order-inventory
system; urgent payments skip scheduling.

Log L_ex: six cases, of which sigma3 and sigma5 escalated in reality
but the order-inventory events (e) were not captured in the extract.
sigma6 is an urgent payment that skips scheduling entirely (the only
trace omitting s).

This is the paper's worked example. It doubles as a quick smoke test for
the harness (see exp_premodel.py) and a fixture for the coverage metric
unit tests.

NOTE: the log below fits a loop body with a silent skip (sigma3 and
sigma4 have no `a` at all), not this model's mandatory `a` - don't rely
on this pair for anything beyond a smoke test.

payment_partial is a second, larger fixture over the same story, in the
shape inductive_noise20 actually discovers from it:

    N = seq( o, xor( tau, loop(a, e) ), s, p )

Approval is optional under a silent branch, and the loop is reached only
when that branch is not taken. It keeps all six traces above and adds
two where a subprocess was traversed but recorded INCOMPLETELY, which
none of the six do: every absence there is a whole traversal, so a
subprocess either runs and is fully recorded or does not run at all.

sigma7 is the discriminating case. A loop's language is a (e a)*, so the
fragment <a, e> is not a complete traversal - it needs a closing a, and
the alignment pairs that closing a with no log event. The loop therefore
has one traversal holding two synchronous moves and one model move: it
ran, so its skip probability is 0, yet a third of its moves went
unobserved. sigma8 is the same variant with different timings, so the
duration-based metrics have something to separate that the move-based
ones cannot see - the same role sigma3/sigma4 play above.

partial_sequence is a third fixture, and the only one here that models
nothing - it is named for what it tests. seq(o, seq(x, y, z), p) over
three traces recording the inner sequence completely, without its middle
step, and without two of its three. Each trace has exactly one optimal
alignment, so a partial term reads its full magnitude, which
payment_partial cannot offer: there, completing a loop traversal and
discarding an escalation cost the same, so the partial explanation is
always one of two ties.

Between them: payment_approval for the paper's worked example,
payment_partial for what the metrics do on a realistic discovered model,
and partial_sequence for the definitional property on its own.

Run this module directly to (re)write all three checked-in XES copies
from the trace strings below.
"""

import pandas as pd
from skipalignments import Activity, Tau, Sequence, Xor, Loop

from process_voids import dtlog

ACTIVITY_COST = 100000

RUNNING_EXAMPLE_XES = 'data/payment_approval.xes'

# sigma1 .. sigma6, each "activity:hour-offset" per the paper's a:t notation
_TRACES = [
    "o:0 a:1 s:2 p:3",
    "o:0 a:1 e:5 a:9 s:10 p:11",
    "o:0 s:10 p:11",
    "o:0 s:1 p:2",
    "o:0 a:1 s:10 p:11",
    "o:0 a:1 p:2",
]
_CASE_NAMES = [f'sigma{i + 1}' for i in range(len(_TRACES))]

PAYMENT_PARTIAL_XES = 'data/payment_partial.xes'
PARTIAL_SEQUENCE_XES = 'data/partial_sequence.xes'

# The six above, plus two traversals of the loop that stop after `e`:
# its closing `a` is unrecorded, so the traversal is partial rather than
# absent. sigma7 and sigma8 are the same variant with different timings.
_PARTIAL_TRACES = _TRACES + [
    "o:0 a:1 e:5 s:10 p:11",
    "o:0 a:1 e:5 s:6 p:7",
]
_PARTIAL_CASE_NAMES = [f'sigma{i + 1}' for i in range(len(_PARTIAL_TRACES))]


def build_running_example_log() -> pd.DataFrame:
    return dtlog.convert_timed(*_TRACES, names=_CASE_NAMES, time_unit='hours')


def build_payment_partial_log() -> pd.DataFrame:
    return dtlog.convert_timed(*_PARTIAL_TRACES, names=_PARTIAL_CASE_NAMES,
                               time_unit='hours')


# One subprocess recorded completely, then missing its middle step, then
# missing two of three. Named for what it tests rather than for a
# domain: it models nothing, unlike the two payment fixtures.
_PARTIAL_SEQUENCE_TRACES = [
    "o:0 x:1 y:2 z:3 p:4",
    "o:0 x:1 z:3 p:4",
    "o:0 x:1 p:4",
]
_PARTIAL_SEQUENCE_CASE_NAMES = ['complete', 'one_missing', 'two_missing']


def build_partial_sequence_log() -> pd.DataFrame:
    return dtlog.convert_timed(*_PARTIAL_SEQUENCE_TRACES,
                               names=_PARTIAL_SEQUENCE_CASE_NAMES, time_unit='hours')


def build_running_example_tree():
    o = Activity(None, 'o', ACTIVITY_COST)
    o.id = '1'
    a = Activity(None, 'a', ACTIVITY_COST)
    a.id = '2'
    e = Activity(None, 'e', ACTIVITY_COST)
    e.id = '3'
    s = Activity(None, 's', ACTIVITY_COST)
    s.id = '4'
    tau = Tau(None, 'skip-schedule', 0)
    tau.id = '5'
    p = Activity(None, 'p', ACTIVITY_COST)
    p.id = '6'

    approval = Loop(None, [a, e])
    approval.id = '7'
    a.set_parent(approval)
    e.set_parent(approval)

    schedule_choice = Xor(None, [s, tau])
    schedule_choice.id = '8'
    s.set_parent(schedule_choice)
    tau.set_parent(schedule_choice)

    tree = Sequence(None, [o, approval, schedule_choice, p])
    tree.id = '9'
    o.set_parent(tree)
    approval.set_parent(tree)
    schedule_choice.set_parent(tree)
    p.set_parent(tree)
    return tree


def build_payment_partial_tree():
    """
    seq( o, xor( tau, loop(a, e) ), s, p ) - the shape
    inductive_noise20 discovers from this log, hand-built so a test
    against it does not depend on discovery's own cross-process
    tie-break. Unlike build_running_example_tree, approval sits under a
    silent branch (so skipping it costs nothing) and s is mandatory.
    """
    o = Activity(None, 'o', ACTIVITY_COST)
    o.id = '1'
    tau = Tau(None, 'skip-approval', 0)
    tau.id = '2'
    a = Activity(None, 'a', ACTIVITY_COST)
    a.id = '3'
    e = Activity(None, 'e', ACTIVITY_COST)
    e.id = '4'
    s = Activity(None, 's', ACTIVITY_COST)
    s.id = '5'
    p = Activity(None, 'p', ACTIVITY_COST)
    p.id = '6'

    approval = Loop(None, [a, e])
    approval.id = '7'
    a.set_parent(approval)
    e.set_parent(approval)

    approval_choice = Xor(None, [tau, approval])
    approval_choice.id = '8'
    tau.set_parent(approval_choice)
    approval.set_parent(approval_choice)

    tree = Sequence(None, [o, approval_choice, s, p])
    tree.id = '9'
    o.set_parent(tree)
    approval_choice.set_parent(tree)
    s.set_parent(tree)
    p.set_parent(tree)
    return tree


def build_partial_sequence_tree():
    """
    seq( o, seq(x, y, z), p ) - five nodes, no choice and no loop.

    The inner sequence is the node under test: it is traversed in every
    trace, so its skip probability is 0 throughout, and the only thing
    separating the three traces is how much of it was recorded. Each
    trace has exactly ONE optimal alignment, so a partial term reads its
    full magnitude here - unlike payment_partial's loop, where
    completing a traversal and discarding an escalation cost the same
    and the partial explanation is one of two ties.
    """
    o = Activity(None, 'o', ACTIVITY_COST)
    o.id = '1'
    x = Activity(None, 'x', ACTIVITY_COST)
    x.id = '2'
    y = Activity(None, 'y', ACTIVITY_COST)
    y.id = '3'
    z = Activity(None, 'z', ACTIVITY_COST)
    z.id = '4'
    p = Activity(None, 'p', ACTIVITY_COST)
    p.id = '5'

    inner = Sequence(None, [x, y, z])
    inner.id = '6'
    for node in (x, y, z):
        node.set_parent(inner)

    tree = Sequence(None, [o, inner, p])
    tree.id = '7'
    for node in (o, inner, p):
        node.set_parent(tree)
    return tree


if __name__ == '__main__':
    dtlog.write_xes(build_running_example_log(), RUNNING_EXAMPLE_XES)
    print(f'Wrote {RUNNING_EXAMPLE_XES}')
    dtlog.write_xes(build_payment_partial_log(), PAYMENT_PARTIAL_XES)
    print(f'Wrote {PAYMENT_PARTIAL_XES}')
    dtlog.write_xes(build_partial_sequence_log(), PARTIAL_SEQUENCE_XES)
    print(f'Wrote {PARTIAL_SEQUENCE_XES}')
