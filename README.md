# process-voids
Characterising unobserved activities in process event data.

Background and motivation can be found in [this blog post](https://adamburkeware.net/2025/11/11/pvoid-adsn.html). This was presented at ADSN 2025.

_Burke, A., Wynn, M.T. (2025). Process Voids: Data Science Without Data. Talk. Australian Data Science Network 2025._

The following related paper is under review:

_Burke, A., Wynn, M.T. (2026). Measuring Process Voids Using Skip Alignments. In review._


## Building and Running

In process-voids, perhaps in a venv:

```
pip install -e .
```

Install [ebi](https://bpm.rwth-aachen.de/ebi/).

Configure paths to ebi, data files, and so on in `pvoid.toml`.

To calculate the voids in a PTML process tree model against a XES event log:

```
python -m process_voids.pvoid <log> <model> [--metric voidsalign|voidsat|voidmass_process]
```

`--metric` chooses the void metric (default `voidsalign`). The same three are available from Python as `process_voids.pvoid.voidsalign`, `voidsat` and `voidmass_process`, each taking a log and a process tree and returning a value for every node.

## Sample Output

This uses a filtered version of the [Road Traffic Fines dataset](https://data.4tu.nl/articles/_/12683249/1). Firstly a number of activities are filtered out to make a clearer example. Secondly a process model is used,  based on a discovered inductive miner model, but which introduces the fictional activity _Certify Judgement_ in the middle of the main process sequence. This makes it a compulsory step which is never observed in the log, ie, a process void.

```
$ python -u -m process_voids.pvoid logs/rtfm_fine_appeal.xes.gz models/rtfm_extra.ptml

...

Calculated at 2026-09-20 00:25:35.858709
node : skip probability, voidsalign
→ : 0.0012237967107096171, 0.5984318351636353
  × : 7.419133891686282e-05, 0.09363011182292846
    Act( Appeal to Judge ) : 0.0, 0.0
    → : 0.0, 0.09358220403988748
      Act( Send Fine ) : 0.0, 0.0
      Act( Insert Fine Notification ) : 0.23186426331685353, 0.2318642633168535
      Act( Add penalty ) : 0.24921031127139176, 0.2492103112713917
  Act( Certify Judgement ) : 1.0, 1.0
  ∧ : 0.9680111510731652, 0.9730005057488195
    × : [ 0.7613227893601723 ], 0.7613227893601723
      Tau( TAU_Receive Result Appeal from Prefecture ) : [ 0.0 ], 1.0
      Act( Receive Result Appeal from Prefecture ) : 0.0, 0.0
    → : 0.01998646296292657, 0.17995351910221014
      Act( Insert Date Appeal to Prefecture ) : 0.15238992211297941, 0.15238992211297941
      Act( Notify Result Appeal to Offender ) : 0.7851834092543752, 0.7851834092543752
```

Each node of the process tree shows its skip probability (in brackets where the node can be traversed silently) and its void, from 0 where the node is fully backed by event log data to 1 where it is never observed. The root's void summarises the whole process. _Certify Judgement_ reads 1.
