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
    'rtfm': 'C:/working/data/rtfm.xes',
    'sepsis': 'C:/working/data/sepsis.xes',
    'bpi2013_incidents': 'C:/working/data/BPI_Challenge_2013_incidents.xes',
}

ALL_COMBOS = COMBOS

ALL_DEGRADATIONS = DEGRADATIONS

# 1.0 (full degradation) deliberately excluded - the degraded log is
# empty there, so every metric falls through its own (disagreeing) "no
# data" default rather than measuring anything - see plots.py's
# _exclude_degenerate, which still filters level=1.0 out defensively
# for any CSV that has it anyway (old results, or an explicit --levels
# override) even though the default roster no longer produces it.
ALL_LEVELS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
