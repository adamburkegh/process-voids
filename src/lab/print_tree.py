'''
Prints a discovered process tree as text - a quick look at what a combo
actually discovered on a log, without recomputing skip-probabilities
just to see the tree's shape. The tree itself already has a readable
str() in skip-alignments' own process-tree notation (sequence/xor/
and/loop symbols); this is just the wiring to get one via the shared
discovery cache (lab.discovery.discover_cached), same as every
experiment script uses.
'''

import argparse
import sys
from pathlib import Path

import pm4py_config as pm4py

from lab.discovery import discover_cached
from lab.params import ALL_COMBOS, ALL_LOGS


def discovered_tree_text(log, combo_name):
    '''
    log: a registered name in lab.params.ALL_LOGS, or an XES path
    directly. combo_name: a key in lab.params.ALL_COMBOS.

    Goes through discover_cached, so repeated calls for the same
    (log, combo) are free after the first, and the tree matches
    whatever an experiment run against the same pair would score -
    not a separate, possibly different discovery.
    '''
    log_path = ALL_LOGS.get(log, log)
    log_name = log if log in ALL_LOGS else Path(log_path).stem
    base_log = pm4py.read_xes(log_path)
    tree, _ppt_weights = discover_cached(log_name, combo_name, ALL_COMBOS[combo_name], base_log)
    return str(tree)


def main():
    parser = argparse.ArgumentParser(
        description='Print a discovered process tree as text.')
    parser.add_argument('log', help=f'Registered log name ({list(ALL_LOGS)}) or an XES path')
    parser.add_argument('combo', choices=list(ALL_COMBOS), help='Discovery combo name')
    args = parser.parse_args()

    # skip-alignments' tree notation uses non-ASCII operator symbols
    # (->, x, ^, loop) - the default Windows console codepage can't
    # encode them, so force UTF-8 for this process's stdout rather than
    # let a perfectly fine tree crash on print.
    sys.stdout.reconfigure(encoding='utf-8')
    print(discovered_tree_text(args.log, args.combo))


if __name__ == '__main__':
    main()
