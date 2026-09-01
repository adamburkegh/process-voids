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
the harness (see exp_premodel.py) and is bookmarked for the coverage
metric unit tests still to be written.

NOTE: the log below predates this model revision (loop body is now the
mandatory `a`, not a choice with a silent skip) - see the bookmarked
"running example not explanatory" item before relying on this pair for
anything beyond a smoke test.

Run this module directly to (re)write the checked-in XES copy at
data/payment_approval.xes from the trace strings below.
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


def build_running_example_log() -> pd.DataFrame:
    return dtlog.convert_timed(*_TRACES, names=_CASE_NAMES, time_unit='hours')


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


if __name__ == '__main__':
    dtlog.write_xes(build_running_example_log(), RUNNING_EXAMPLE_XES)
    print(f'Wrote {RUNNING_EXAMPLE_XES}')
