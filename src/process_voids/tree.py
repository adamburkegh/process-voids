"""
Conversion from a pm4py process tree to this project's skip-alignments
ProcessTree, using standard cost constants.
"""

from skipalignments import ProcessTree, update_pair_taus

MM_COST = 100000
TAU_COST = 0
SYNCH_COST = 0


def from_pm4py(pt_pm4py) -> ProcessTree:
    """
    Convert a pm4py process tree using this project's standard cost
    constants, with paired taus resolved.
    """
    pt = ProcessTree.from_pm4py(pt_pm4py, MM_COST, TAU_COST, SYNCH_COST)
    update_pair_taus(pt)
    return pt
