
import datetime
import random
import sys

sys.stdout.reconfigure(encoding='utf-8')


print( f'Importing pm4py {datetime.datetime.now()}')
import pm4py


from coveragemass import *
from derivation import DerivationPipeline, EbiWeights
import probabilities
from processtree import *
import slpn_importer

probabilities.EBI_EXECUTABLE="./ebi"   # Path to ebi link

MM_COST = 100000
TAU_COST = 0
SYNCH_COST = 0



def show_skip_outcome(dv):
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


'''
This ... what one can only call filthy hack ... taken from im_models.ipynb
'''
def update_pair_taus(tree:ProcessTree):
    if isinstance(tree, Tau):
        if tree.parent is not None and len(tree.parent.children) == 2:
            other = tree.parent.children[0]
            if other == tree:
                other = tree.parent.children[1]
            if isinstance(other, Activity):
                # set tau
                tree.name = "TAU_" + other.name
            else:
                tree.name = "TAU_" + other.id
        else:
            tree.name = "TAU_" + str(tree.get_distance_to_root()) + str(random.random())
        return
    elif not isinstance(tree, Activity):
        for c in tree.children:
            update_pair_taus(c)
        return
    else:
        return


def skipprob(log, pt, slpn_path ):
    dv = DerivationPipeline(pt, log, pn_log=log, 
                            pn_method=EbiWeights.OCCURANCE,
                            sagn_timeout=600)
    dv.compute(path='var', slpn_path=slpn_path ) 
    return dv


def main():
    print( f'Started at {datetime.datetime.now()}')
    logx   = pm4py.read_xes( sys.argv[1] )
    modelt = pm4py.read_ptml( sys.argv[2] )
    pt = ProcessTree.from_pm4py( modelt, MM_COST, TAU_COST, 
                                             SYNCH_COST )
    update_pair_taus(pt)
    slpn_path = 'var/spmodel.slpn'
    dv = skipprob(logx, pt, slpn_path)
    show_skip_outcome(dv)
    print( f'Skip probabilities calculated at {datetime.datetime.now()}')
    slpn = slpn_importer.read_slpn(slpn_path)
    transfer_pt_weights(pt,slpn)
    print(show_tree_weights(pt,dv))
    print( '==========' )
    print( f'Coverage: {coverage_mass(pt, dv.skip_probs)}' )
    print( '==========' )
    print( f'Finished at {datetime.datetime.now()}')

if __name__ == '__main__':
    main()


