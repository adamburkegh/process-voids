# process-voids
Characterising unobserved activities in process event data.

Background and motivation can be found in [this blog post](https://adamburkeware.net/2025/11/11/pvoid-adsn.html). This was presented at ADSN 2025; a preprint is forthcoming.

Burke, A., Wynn, M.T. (2025). Process Voids: Data Science Without Data. Talk. Australian Data Science Network 2025.

This project is based on a hard fork of the [skip-probabilities](https://git.rwth-aachen.de/philippbaer/skip-probabilities) repository created by Philipp Bär.

## Building and Running

Install dependencies (perhaps in a dedicated environment)

```
pip install -r requirements.txt
```

Install [ebi](https://bpm.rwth-aachen.de/ebi/) and create a link in the local directory, or set the `probabilities.EBI_EXECUTABLE` constant to the path of the executable.

Calculate skip probabilities and coverage on a XES event log and a PTML process tree model.

```
python pvoid.py <log> <model>
```

## Sample Output

This uses a filtered version of the [Road Traffic Fines dataset](https://data.4tu.nl/articles/_/12683249/1). Firstly a number of activities are filtered out to make a clearer example. Secondaly an process model is used, that is based on a discovered inductive miner model, but introduces the fictional activity Certify Judgement in the middle of the main process sequence. This makes it a compulsory step which is never observed in the log, ie, a process void.

```
$ python -u pvoid.py logs/rtfm_fine_appeal.xes.gz models/rtfm_extra.ptml

...

Skip probabilities calculated at 2025-11-10 16:44:03.078990
→ : 0.4562279269951479, 0.0012237967107096175
  × : 0.8506576142530637, 7.419133891686282e-05
    Act( Appeal to Judge ) : 0.005337205612239992, 0.0
    → : 0.8453204086408238, 0.0
      Act( Send Fine ) : 1.0, 0.0
      Act( Insert Fine Notification ) : 0.7679806129612355, 0.23186426331685353
      Act( Add penalty ) : 0.7679806129612355, 0.24921031127139176
  Act( Certify Judgement ) : 0.001, 1.0
  ∧ : 0.51702616673238, 0.9680111510731652
    × : 1.009606970102032, [ 0.7613227893601723 ]
      Tau( TAU_Receive Result Appeal from Prefecture ) : 1.0, [ 0.0 ]
      Act( Receive Result Appeal from Prefecture ) : 0.009606970102031985, 0.0
    → : 0.024445363362728033, 0.01998646296292657
      Act( Insert Date Appeal to Prefecture ) : 0.04027426505236231, 0.15238992211297941
      Act( Notify Result Appeal to Offender ) : 0.008616461673093751, 0.7851834092543752
==========
Coverage: 0.5354181989351403
```

Output describes:
 *  Weights and skip probabilities of the process tree
 *  A coverage metric summarising how much of the process is backed by event log data. 



