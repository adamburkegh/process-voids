
import datetime
import sys

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
from process_voids.tree import from_pm4py

probabilities.EBI_EXECUTABLE="ebi.exe"   # resolved via PATH; see README



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


def show_tree_weights(tree,dv):
    if isinstance(tree, LeafNode):
        return tree.__str__() + " : " + str(tree.weight) + \
                ", " + ("[ " if tree.get_cheapest_execution(0)[1] else "") + \
                str(dv.skip_probs[tree]) + \
                (" ]" if tree.get_cheapest_execution(0)[1] else "")
    else:
        if isinstance(tree, Sequence):
            operator = "→"
        elif isinstance(tree, Xor):
            operator = "×"
        elif isinstance(tree, And):
            operator = "∧"
        elif isinstance(tree, Loop):
            operator = "↺"
        else:
            operator = "UNKNOWN"
        operator += " : " + str(tree.weight)
        child_string = "\n".join([show_tree_weights(c,dv) \
                                    for c in tree.children])
        return (" " * tree.get_distance_to_root()*2) + operator + ", " + ("[ " if tree.get_cheapest_execution(0)[1] else "") + str(dv.skip_probs[tree]) + (" ]" if tree.get_cheapest_execution(0)[1] else "") + "\n" + child_string


def show_tree_coverage_by_duration(tree, dv, traces, total_dur=None):
    if total_dur is None:
        total_dur = sum([ dur(sigma) for sigma in traces ])
    cov = coverage_by_duration(tree, traces, dv.skip_probs, total_dur)
    if isinstance(tree, LeafNode):
        return tree.__str__() + " : " + str(cov) + \
                ", " + ("[ " if tree.get_cheapest_execution(0)[1] else "") + \
                str(dv.skip_probs[tree]) + \
                (" ]" if tree.get_cheapest_execution(0)[1] else "")
    else:
        if isinstance(tree, Sequence):
            operator = "→"
        elif isinstance(tree, Xor):
            operator = "×"
        elif isinstance(tree, And):
            operator = "∧"
        elif isinstance(tree, Loop):
            operator = "↺"
        else:
            operator = "UNKNOWN"
        operator += " : " + str(cov)
        child_string = "\n".join([show_tree_coverage_by_duration(c,dv,traces,total_dur) \
                                    for c in tree.children])
        return (" " * tree.get_distance_to_root()*2) + operator + ", " + ("[ " if tree.get_cheapest_execution(0)[1] else "") + str(dv.skip_probs[tree]) + (" ]" if tree.get_cheapest_execution(0)[1] else "") + "\n" + child_string


def skipprob(log, pt, slpn_path ):
    dv = DerivationPipeline(pt, log, pn_log=log, 
                            pn_method=DiscoverySource.OCCURANCE,
                            sagn_timeout=600)
    dv.compute(path='var', slpn_path=slpn_path ) 
    return dv


def main():
    print( f'Started at {datetime.datetime.now()}')
    logx   = pm4py.read_xes( sys.argv[1] )
    modelt = pm4py.read_ptml( sys.argv[2] )
    pt = from_pm4py(modelt)
    slpn_path = 'var/spmodel.slpn'
    dv = skipprob(logx, pt, slpn_path)
    show_skip_outcome(dv)
    print( f'Skip probabilities calculated at {datetime.datetime.now()}')
    slpn = slpn_importer.read_slpn(slpn_path)
    transfer_pt_weights(pt,slpn)
    print(show_tree_weights(pt,dv))
    print( '==========' )
    print( f'Coverage: {mass_by_weight(pt, dv.skip_probs)}' )
    print( '==========' )
    traces = log_to_traces(logx)
    total_dur = sum([ dur(sigma) for sigma in traces ])
    print(show_tree_coverage_by_duration(pt, dv, traces, total_dur))
    print( '==========' )
    print( f'Coverage by Duration: {coverage_by_duration(pt, traces, dv.skip_probs, total_dur)}' )
    print( '==========' )
    print( f'Finished at {datetime.datetime.now()}')

if __name__ == '__main__':
    main()


