"""
Full catalog of parameters available to the experiment harness.

A named run in lab/runs.py selects a subset of these (or supplies its
own values) to run as a Cartesian product via run_experiment1(). This
module is the single place new logs/levels get registered.
"""

from lab.degradation import DEGRADATIONS
from lab.discovery import COMBOS

ALL_LOGS = {
    'rtfm_fine_appeal': 'logs/rtfm_fine_appeal.xes.gz',
    'payment_approval': 'data/payment_approval.xes',
}

ALL_COMBOS = COMBOS

ALL_DEGRADATIONS = DEGRADATIONS

ALL_LEVELS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
