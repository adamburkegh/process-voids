'''
Process voids in an event log against a process tree model.

Three void metrics, each returning its value at every node of the tree -
0 where the node is fully backed by event data, 1 where it is never
observed:

    voidsalign(log, tree)        void by skip alignment correspondence
    voidsat(log, tree)           void by aligned duration
    voidmass_process(log, tree)  void by process-relative alignment moves

`log` is an event log DataFrame (as pm4py reads one), `tree` a
skipalignments process tree (process_voids.tree.from_pm4py converts a
pm4py one). Each call does its own alignment work.

From the command line, one metric at a time, as a tree:

    python -m process_voids.pvoid <log.xes> <model.ptml> [--metric voidsalign|voidsat|voidmass_process]
'''

import argparse
import datetime
import sys
import warnings

sys.stdout.reconfigure(encoding='utf-8')
DEBUG = False

if DEBUG:
    print( f'Importing pm4py {datetime.datetime.now()}')
import pm4py_config as pm4py


from process_voids.coveragemass import *
from skipalignments import (
    DerivationPipeline, DiscoverySource, LeafNode, Activity, Tau,
    Sequence, Xor, And, Loop, probabilities,
)
from process_voids import slpn_importer
from process_voids.config import ebi_executable
from process_voids.tree import from_pm4py
from process_voids.voidmass_pn import build_id_net
from process_voids.voidsalign3 import voidsalign3
from process_voids.voidsat2 import voidsat2

probabilities.EBI_EXECUTABLE = ebi_executable()   # bare 'ebi' on PATH unless pvoid.toml says otherwise

SLPN_PATH = 'var/spmodel.slpn'


def show_skip_outcome(dv):
    if not DEBUG:
        return
    print(dv.print_blinded())
    print('=====')
    print(dv.stats())
    print('=====')
    for d in dv.skip_dict_backup:
        print(f'Trace:     {d}')
        sk = dv.skip_dict_backup[d]
        for state in sk:
            pstr  = ', '.join([name for (name,obj) in state.path])
            pstr2 = ', '.join([str(obj)  for (name,obj) in state.path])
            print(f'    Path:  {pstr}')
            print(f'    Path:  {pstr2}')
            print(f'    Path:  {state.path}')
            # print(f'    State: {state.state}')
            print(f'    Trace: {state.trace}')
            print(f'    Costs: {state.acc_costs}')


def skipprob(log, pt, slpn_path, ppt_weights=None):
    if ppt_weights is not None:
        # Toothpaste's weights are exact from the PPT - no pn_log/
        # estimation pass needed.
        dv = DerivationPipeline(pt, log,
                                pn_method=DiscoverySource.TOOTHPASTE,
                                pn_ppt_weights=ppt_weights,
                                sagn_timeout=600)
    else:
        dv = DerivationPipeline(pt, log, pn_log=log,
                                pn_method=DiscoverySource.OCCURANCE,
                                sagn_timeout=600)
    dv.compute(path='var', slpn_path=slpn_path )
    return dv


def _nodes(tree):
    yield tree
    for child in getattr(tree, 'children', []) or []:
        yield from _nodes(child)


def _context(log, tree):
    # Imported here, not at the top: metric_context imports this module
    # for skipprob.
    from process_voids.metric_context import CellContext
    return CellContext(log=log, tree=tree, slpn_path=SLPN_PATH,
                       classical_net=build_id_net(tree))


def _voidsalign(ctx):
    dv = ctx.stage('dv')
    cache = ctx.stage('executions_cache')
    return {node: voidsalign3(node, dv.skip_dict_backup, dv.pl, dv.skip_probs,
                              executions_cache=cache)
            for node in _nodes(ctx.tree)}


def _voidsat(ctx):
    dv = ctx.stage('dv')
    paths = {k: [s.path for s in v] for k, v in dv.skip_dict_backup.items()}
    cache = ctx.stage('aligned_duration_cache')
    return {node: voidsat2(node, ctx.tree, ctx.log, paths, dv.skip_probs, cache=cache)
            for node in _nodes(ctx.tree)}


def _voidmass_process(ctx):
    result, _variant_probs = ctx.stage('classical')
    if result.timed_out_count:
        warnings.warn(f'{result.timed_out_count} trace variant(s) timed out in classical '
                      f'alignment (weight {result.timed_out_weight:.4g}); voidmass_process '
                      f'reports its lower bound')
    return {node: result.table[node]['voidmass_process_lower'] for node in _nodes(ctx.tree)}


def voidsalign(log, tree):
    '''Void by skip alignment correspondence, {node: value} for every node.'''
    return _voidsalign(_context(log, tree))


def voidsat(log, tree):
    '''Void by aligned duration, {node: value} for every node.'''
    return _voidsat(_context(log, tree))


def voidmass_process(log, tree):
    '''
    Void by process-relative alignment moves, {node: value} for every
    node. Where a classical alignment times out, the value is a lower
    bound, and a warning says so.
    '''
    return _voidmass_process(_context(log, tree))


METRICS = {'voidsalign': voidsalign, 'voidsat': voidsat, 'voidmass_process': voidmass_process}

_BY_CONTEXT = {'voidsalign': _voidsalign, 'voidsat': _voidsat,
               'voidmass_process': _voidmass_process}


_OPERATORS = ((Sequence, '→'), (Xor, '×'), (And, '∧'), (Loop, '↺'))


def show_tree(tree, skip_probs, values):
    '''
    One line per node, indented by depth: skip probability ([ ] where the
    node can be traversed silently), then the metric's value.
    '''
    silent = tree.get_cheapest_execution(0)[1]
    skip = f'[ {skip_probs[tree]} ]' if silent else str(skip_probs[tree])
    fields = f'{skip}, {values[tree]}'
    if isinstance(tree, LeafNode):
        return f'{tree} : {fields}'
    operator = next((op for cls, op in _OPERATORS if isinstance(tree, cls)), 'UNKNOWN')
    line = ' ' * tree.get_distance_to_root() * 2 + f'{operator} : {fields}'
    return '\n'.join([line] + [show_tree(c, skip_probs, values) for c in tree.children])


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='python -m process_voids.pvoid',
        description='Report a void metric at every node of a process tree.')
    parser.add_argument('log', help='XES event log')
    parser.add_argument('model', help='PTML process tree')
    parser.add_argument('--metric', choices=tuple(METRICS), default='voidsalign',
                        help='void metric to report (default: voidsalign)')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    print( f'Started at {datetime.datetime.now()}')
    log = pm4py.read_xes(args.log)
    tree = from_pm4py(pm4py.read_ptml(args.model))
    ctx = _context(log, tree)
    values = _BY_CONTEXT[args.metric](ctx)
    dv = ctx.stage('dv')
    show_skip_outcome(dv)
    print( f'Calculated at {datetime.datetime.now()}')
    print(f'node : skip probability, {args.metric}')
    print(show_tree(tree, dv.skip_probs, values))

if __name__ == '__main__':
    main()
