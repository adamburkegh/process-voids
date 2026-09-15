"""
Full catalog of parameters available to the experiment harness.

A named run in lab/runs.py selects a subset of these (or supplies its
own values) to run as a Cartesian product via run_disco_degrade(). This
module is the single place new logs/levels get registered.
"""

from lab.degradation import DEGRADATIONS
from lab.discovery import COMBOS

ALL_LOGS = {
    'payment_approval': 'data/payment_approval.xes',
    'payment_partial': 'data/payment_partial.xes',
    'partial_sequence': 'data/partial_sequence.xes',
    'rtfm': 'C:/working/data/rtfm.xes',
    'sepsis': 'C:/working/data/sepsis.xes',
    'bpi2013_incidents': 'C:/working/data/BPI_Challenge_2013_incidents.xes',
    'bpic2020_rfp': 'C:/working/data/BPIC2020_rfp.xes',
}

ALL_COMBOS = COMBOS

ALL_DEGRADATIONS = DEGRADATIONS

# 1.0 (full degradation) deliberately excluded - the degraded log is
# empty there, so every metric falls through its own (disagreeing) "no
# data" default rather than measuring anything - see plots.py's
# _exclude_degenerate, which also filters level=1.0 out defensively
# for any CSV that has it anyway (eg from an explicit --levels
# override).
ALL_LEVELS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


# Seconds this lab allows one variant's alignment search before it is
# abandoned and the cell falls back to bounds - a budget, chosen for the
# logs and hardware here, not a property of any search. Applied to the
# classical Petri-net search (process_voids.voidmass_pn, whose own
# DEFAULT_ALIGNMENT_TIMEOUT is only the fallback for a caller that
# states nothing) and to skip-alignments' align_sk_all alike, since the
# question it answers - how long are we willing to wait - is the same
# either way.
CLASSICAL_ALIGNMENT_TIMEOUT = 100
