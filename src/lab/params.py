"""
Full catalog of parameters available to the experiment harness.

A named run in lab/runs.py selects a subset of these (or supplies its
own values) to run as a Cartesian product via run_disco_degrade(). This
module is the single place new logs/levels get registered.
"""

import os
from pathlib import Path

from lab.degradation import DEGRADATIONS
from lab.discovery import COMBOS
from process_voids.config import value


class ExternalLog(os.PathLike):
    """
    A log kept outside the repository, registered by filename. Its
    directory is [paths] data_dir in pvoid.toml, read only when the path
    is actually used - os.fspath, or Path() around it - so building
    ALL_LOGS, and lab.runs' named runs from it, needs no configuration.

    Resolves with forward slashes. pm4py's reader is typed for a str, so
    a caller opening one passes os.fspath(log_path) rather than the
    object.
    """

    def __init__(self, filename):
        self.filename = filename

    def __fspath__(self):
        return (Path(value('paths', 'data_dir')) / self.filename).as_posix()

    def __repr__(self):
        return f'ExternalLog({self.filename!r})'


ALL_LOGS = {
    'payment_approval': 'data/payment_approval.xes',
    'payment_partial': 'data/payment_partial.xes',
    'partial_sequence': 'data/partial_sequence.xes',
    'rtfm': ExternalLog('rtfm.xes'),
    'sepsis': ExternalLog('sepsis.xes'),
    'bpic2020_rfp': ExternalLog('BPIC2020_rfp.xes'),
    'bpi2013_closed_problems': ExternalLog('BPI_Challenge_2013_closed_problems.xes'),
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
