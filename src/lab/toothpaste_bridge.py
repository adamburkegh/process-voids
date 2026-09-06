"""
Bridge to the external toothpaste miner tool (not part of this repo -
installed separately at TOOTHPASTE_DIR). Mines a log with toothpaste's
own dcdt pipeline, then translates its .ptree output into skip-alignments'
ProcessTree representation.

Mirrors the pattern already used for ebi in process_voids/pvoid.py:
rather than requiring the executable on PATH, we import toothpaste's own
Python wrapper module directly from its install location and point its
TBIN constant at the exe - see README for the equivalent ebi setup.
"""

import sys
import tempfile
from pathlib import Path

from skipalignments.ppt import parse_ptree, translate_ppt

from lab.discovery import DiscoveryResult
from process_voids.dtlog import write_xes

TOOTHPASTE_DIR = r'C:\working\tools\toothpaste'
TOOTHPASTE_EXECUTABLE = str(Path(TOOTHPASTE_DIR) / 'toothpaste.exe')

if TOOTHPASTE_DIR not in sys.path:
    sys.path.insert(0, TOOTHPASTE_DIR)

import toothpaste as _toothpaste_tool  # the external tool's own wrapper module

_toothpaste_tool.TBIN = TOOTHPASTE_EXECUTABLE


def discover(log):
    """
    Mine `log` (a pm4py-format DataFrame) with toothpaste, returning a
    DiscoveryResult whose ppt_weights is the (weights, loop_taus) pair
    translate_ppt() returns, needed by DerivationPipeline's
    DiscoverySource.TOOTHPASTE path (toothpaste's weights are exact from
    the PPT, not estimated separately from the log).

    toothpaste.mine() only takes an XES file path (not an in-memory
    log), so the log is written to a temp file first; it also only
    returns the mined PNML/PetriNet, not the .ptree path, so that's
    recovered from the deterministic sibling filename it writes.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        xes_path = str(Path(tmpdir) / 'log.xes')
        write_xes(log, xes_path)
        _, pnfile = _toothpaste_tool.mine(xes_path)
        ptree_path = str(Path(pnfile).with_suffix('.ptree'))
        ppt = parse_ptree(Path(ptree_path).read_text())

    tree, weights, loop_taus = translate_ppt(ppt)
    return DiscoveryResult(tree, ppt_weights=(weights, loop_taus))
