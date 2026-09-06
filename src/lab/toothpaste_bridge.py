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

import os.path
import subprocess
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


def _mine(logfile, noise=0.0):
    """
    Reimplements toothpaste.mine() (see TOOTHPASTE_DIR/toothpaste.py) to
    add the --noise flag, which that wrapper doesn't expose - see
    `toothpaste --help`: -n/--noise=NUM, prune subtrees below this
    threshold, range [0,1], default 0. Returns the .pnml path (the
    .ptree path is its deterministic sibling); unlike mine(), doesn't
    bother reading the PNML back in, since callers here only need the
    .ptree.
    """
    outdir = os.path.dirname(logfile)
    pref = os.path.basename(logfile).split('.')[0]
    dfile = os.path.join(outdir, pref + '.dcdt')
    pnfile = os.path.join(outdir, pref + '.pnml')
    ptfile = os.path.join(outdir, pref + '.ptree')
    _toothpaste_tool.xestodcdt(logfile, dfile)
    subprocess.run([_toothpaste_tool.TBIN,
                    '--logformat=dcdt',
                    '--eventlog', dfile,
                    '--pnetfile', pnfile,
                    '--ptreefile', ptfile,
                    f'--noise={noise}'])
    return pnfile


def discover(log, noise=0.0):
    """
    Mine `log` (a pm4py-format DataFrame) with toothpaste, returning a
    DiscoveryResult whose ppt_weights is the (weights, loop_taus) pair
    translate_ppt() returns, needed by DerivationPipeline's
    DiscoverySource.TOOTHPASTE path (toothpaste's weights are exact from
    the PPT, not estimated separately from the log).

    noise: toothpaste's own pruning threshold (see _mine), default 0 -
    unrelated to the pm4py inductive miner's noise_threshold, just the
    same idea for this different miner.

    toothpaste only takes an XES file path (not an in-memory log), so
    the log is written to a temp file first.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        xes_path = str(Path(tmpdir) / 'log.xes')
        write_xes(log, xes_path)
        pnfile = _mine(xes_path, noise=noise)
        ptree_path = str(Path(pnfile).with_suffix('.ptree'))
        ppt = parse_ptree(Path(ptree_path).read_text())

    tree, weights, loop_taus = translate_ppt(ppt)
    return DiscoveryResult(tree, ppt_weights=(weights, loop_taus))
