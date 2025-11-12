
from processtree import *


'''
Currently ignores silents. Assumes labels are activity ids
'''
def update_activity_weights(pt:ProcessTree,slpn):
    leaves = pt.get_leafs()
    for leaf in leaves:
        for tran in slpn.transitions:
            # print(f'uaw( {tran} ... {leaf.id} )')
            # inefficient due to list instead of dict
            if tran['label'] == leaf.id:
                leaf.weight = tran['weight']

def infer_operator_weights(pt:ProcessTree):
    if isinstance(pt,Activity) or isinstance(pt,Tau):
        return
    for child in pt.children:
        infer_operator_weights(child)
    if isinstance(pt,Xor):
        pt.weight = sum([ child.weight for child in pt.children  ])
    if isinstance(pt,And) or isinstance(pt,Sequence) or isinstance(pt,Loop):
        pt.weight = sum([ child.weight for child in pt.children  ]) \
                    / len(pt.children)


def transfer_pt_weights(pt:ProcessTree, slpn):
    update_activity_weights(pt,slpn)
    infer_operator_weights(pt)


'''
Metric which indicates how much of the process is tracked by the data.

Pre: Tree has weights
'''
def coverage_mass(pt:ProcessTree, skip_probs:dict):
    if isinstance(pt,Activity) or isinstance(pt,Tau):
        return 1 - skip_probs[pt]
    child_coverage = []
    total_weight = sum( [child.weight for child in pt.children] )
    if isinstance(pt,Xor):
        return sum([ coverage_mass(child,skip_probs)* child.weight \
                        / total_weight \
                        for child in pt.children  ])
    if isinstance(pt,And) or isinstance(pt,Sequence) or isinstance(pt,Loop):
        return sum([ coverage_mass(child,skip_probs) \
                          for child in pt.children  ]) \
                          / len(pt.children) 
    raise ValueError('Unrecognised process tree node')


