
import sys
sys.stdout.reconfigure(encoding='utf-8')


import pm4py

from processtree import *


def printTree(tree):
    print()
    print(tree)
    print()

aj = Activity(None, 'Appeal to Judge', 100000)
aj.id = "1"
sf = Activity(None, 'Send Fine', 100000)
sf.id = "2"
ifn = Activity(None, 'Insert Fine Notification', 100000)
ifn.id = "3"
ap = Activity(None, 'Add penalty', 100000)
ap.id = "4"
rr = Activity(None, 'Receive Result Appeal from Prefecture', 100000)
rr.id = "5"
ida = Activity(None, 'Insert Date Appeal to Prefecture', 100000)
ida.id = "6"
nra = Activity(None, 'Notify Result Appeal to Offender', 100000)
nra.id = "7"



seqfine = Sequence(None, [sf,ifn,ap])
sf.set_parent(seqfine)
ifn.set_parent(seqfine)
ap.set_parent(seqfine)
seqfine.id = "10"

choicefine = Xor(None, [aj,seqfine])
seqfine.set_parent(choicefine)
aj.set_parent(choicefine)
choicefine.id = "11"

tauappeal = Tau(None,'appeal',0)
tauappeal.id = "12"

choiceappeal = Xor(None, [rr,tauappeal] )
rr.set_parent(choiceappeal)
tauappeal.set_parent(choiceappeal)
choiceappeal.id = "13"

seqbureau = Sequence(None, [ida,nra] )
ida.set_parent(seqbureau)
nra.set_parent(seqbureau)
seqbureau.id = "14"

concappeal = And(None, [choiceappeal,seqbureau] )
choiceappeal.set_parent(concappeal)
seqbureau.set_parent(concappeal)
concappeal.id = "15"

tree = Sequence(None, [choicefine,concappeal])
choicefine.set_parent(tree)
concappeal.set_parent(tree)
tree.id = "20"

print(tree)

pm4py.write_ptml( tree.to_pm4py(), 'var/rtfma.ptml' )
print('Wrote paper model')


cj = Activity(None, 'Certify Judgement', 100000)
cj.id = "8"

tree = Sequence(None, [choicefine,cj,concappeal])
choicefine.set_parent(tree)
cj.set_parent(tree)
concappeal.set_parent(tree)
tree.id = "20"

print(tree)


pm4py.write_ptml( tree.to_pm4py(), 'var/rtfm_extra.ptml' )
print('Wrote model with extra Certify activity')

seqfine.set_parent(None)
printTree(seqfine)
pm4py.write_ptml( seqfine.to_pm4py(), 'var/rtfm_fine.ptml' )
print('Wrote fine sequence model')




choicefine = Xor(None, [aj,seqfine])
seqfine.set_parent(choicefine)
aj.set_parent(choicefine)
choicefine.id = "11"

tree = choicefine

printTree(tree)
pm4py.write_ptml( tree.to_pm4py(), 'var/rtfm_appeal.ptml' )
print('Wrote appeal choice model')




seqfine = Sequence(None, [sf,ifn,ap])
sf.set_parent(seqfine)
ifn.set_parent(seqfine)
ap.set_parent(seqfine)
seqfine.id = "10"

choicefine = Xor(None, [aj,seqfine])
seqfine.set_parent(choicefine)
aj.set_parent(choicefine)
choicefine.id = "11"

seqbureau = Sequence(None, [ida,nra] )
ida.set_parent(seqbureau)
nra.set_parent(seqbureau)
seqbureau.id = "14"

tree = Sequence(None, [choicefine,seqbureau])
choicefine.set_parent(tree)
seqbureau.set_parent(tree)
tree.id = "20"

printTree(tree)
pm4py.write_ptml( tree.to_pm4py(), 'var/rtfm_fa2.ptml' )
print('Wrote fa2 model')


seqfine = Sequence(None, [sf,ifn,ap,cj])
sf.set_parent(seqfine)
cj.set_parent(seqfine)
ifn.set_parent(seqfine)
ap.set_parent(seqfine)
seqfine.id = "10"

choicefine = Xor(None, [aj,seqfine])
seqfine.set_parent(choicefine)
aj.set_parent(choicefine)
choicefine.id = "11"

tree = choicefine

printTree(tree)
pm4py.write_ptml( tree.to_pm4py(), 'var/rtfm_acextra.ptml' )
print('Wrote appeal choice + extra model')




seqfine = Sequence(None, [sf,ifn,ap])
sf.set_parent(seqfine)
ifn.set_parent(seqfine)
ap.set_parent(seqfine)
seqfine.id = "10"

choicefine = Xor(None, [aj,seqfine])
seqfine.set_parent(choicefine)
aj.set_parent(choicefine)
choicefine.id = "11"

seqtop = Sequence( None, [choicefine,cj] )
seqtop.id = "30"
choicefine.set_parent(seqtop)
cj.set_parent(seqtop)
tree = seqtop

printTree(tree)
pm4py.write_ptml( tree.to_pm4py(), 'var/rtfm_acextra2.ptml' )
print('Wrote appeal choice + extra model 2')





