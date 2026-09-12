# Design of process-voids

*Maintained by the house critic. This document distils the ideas the code and
experiments rest on, and gives guidance to anyone, human or bot, working here.
Where this document and executable code disagree, the code is authoritative,
and one of them needs fixing.*

Some things have an executable source of truth and are deliberately not
repeated here:

| Question | Authoritative source |
|---|---|
| Which metrics exist, what every result column has ever meant, status | `src/lab/metric_registry.py`, checked by `tests/lab/test_metric_registry.py` |
| Which experiments can run | `src/lab/params.py`, `src/lab/runs.py` |
| How to work here (environment, tests, source control) | `AGENTS.md` |
| What changed, and why | `CHANGELOG.md`, commit messages |

## 1. The problem

Event logs record what systems wrote down, not what happened. Whole
subprocesses go unrecorded: they run in another system or on paper, or they
leave events that can't be linked to a case. We call these gaps *process
voids*. Classical process mining can't see them. Discovery can't discover what
was never logged, and conformance checking reads a missing activity as a
behavioural deviation.

process-voids measures voids by comparing two witnesses:

- a **reference model** that says what should happen (ideally human-authored,
  or a hybrid of human-authored and discovered structure), and
- an **event log** that says what was recorded.

Where the model expects a subprocess and the log has no trace of it, we have a
candidate void. The metrics quantify how often that happens and how much of
the process is at stake.

## 2. Principles

These come from the project's own results, failures included. Most have
cost at least one candidate metric its life.

1. **No expectation, no void.** Nothing in a log points at what the log
   doesn't contain. Missing data can be detected only against an expectation
   that comes from outside the log: a model or a baseline.

2. **The record can't be its own witness.** Anything estimated from the log
   has already been shaped by the void. Exact-fit discovery absorbs voids as
   optional structure. Occurrence-based weights turn voids into rarity. An
   interval distribution estimated from the log is never surprised by a void
   that has always been there. So the expectation has to come from somewhere
   the void couldn't reach: a domain expert, another system, an earlier
   period, or a deliberately regularised model.

3. **Only obligations can be broken.** Optional structure (a choice with a
   silent branch, a loop's redo child) explains any absence at no cost, so
   voids can be measured only where the model is mandatory.
   `mandatory_node_count` marks the metrics' jurisdiction: on a tree where
   nothing is mandatory, every void metric reads zero by construction.

4. **Time is conserved; only its attribution changes.** Removing an interior
   event merges its interval into the next one, and trace duration is
   unchanged. Once a model says where a void sits, elapsed time can *size*
   it. Time on its own can never *find* one.

5. **"Whether" and "how much" are different questions.** Skip alignments lift
   a deviation to the highest block that explains it. That's the right
   granularity for skip probability ("was it skipped?"), but it throws away
   magnitude: a missing subprocess of eight activities becomes a single move,
   the same as one missing leaf. Every mass term in the project answers one
   question: how do we put the size back? The answers so far are per-activity
   model moves (classical alignments), skip moves weighted by minimum size,
   and elapsed time.

6. **A disagreement isn't a verdict.** A model move admits three readings:
   the work wasn't done (conformance), the model overstates the process
   (model repair), or the work was done but not recorded (void). The
   alignment is identical in all three cases. A void metric measures the
   disagreement, and attributing it to the record is an assumption the
   analyst has to own. Void metrics are therefore screening instruments. They
   should favour sensitivity: the zero convention prefers a false alarm to
   silently absorbing a subprocess that ran unrecorded. They should present
   evidence (where, how often, how much), not conclusions.

## 3. The metric design space

Every void or coverage metric pairs a **source of model information** (none,
structure only, classical alignments, skip alignments) with a **source of
magnitude** (activity alphabet, estimated stochastic weights, uniform move
counts, skip-weighted move counts, interval distributions, attributed
duration). The paper in preparation arranges the candidates on this grid, and
it's the right map. A new metric proposal should say which cell it occupies
and why the neighbouring cells failed.

Most metrics take a *product form*: `coverage = (1 − P(skip)) × mass` or
`void = P(skip) × mass`. The product means something only when the two factors
are independent evidence: how often the subprocess was skipped, and how big it
is. See open question 2 for masses built from alignment match ratios, which
already contain the skip evidence.

The roster itself, with each metric's status, is in the registry. It isn't
copied here.

## 4. What a void metric must do

These are the criteria the experiments and tests apply. A candidate that fails
one should say so in its registry description; the extremes test enforces
that rule for its own criterion.

| Criterion | Meaning | Where it's exercised |
|---|---|---|
| Responsiveness | Rises monotonically as a mandatory subprocess is degraded | Dose-response runs of `exp_disco_degrade` |
| Specificity | No response to ablating an optional subprocess, or to thinning the log (fewer cases, same proportions) | Claims fixture `appeal_seq`; trace-wise degradation |
| Size sensitivity | A missing subprocess of eight activities outweighs one of two | `SizeSensitivityTest` in `test_voidmass.py` |
| Extremes | Zero when nothing is missing, maximal when the subprocess is always missing; each metric declares its promise as a registry `scale` | `tests/lab/test_metric_extremes.py`, with today's exceptions pinned in `KNOWN_FAILURES` |
| Decomposability (desirable) | The root value can be explained from its parts: a sum over any antichain cut, or a stated weighted average | `test_variant2_root_equals_sum_of_children` |
| Honesty under failure | Timeouts produce provable bounds, not guesses or crashes | `_lower`/`_upper` columns, `timed_out_weight` |
| Interpretable units | A reader can say what 0.25 means (a share of expected moves, of elapsed time, ...) | Registry descriptions |
| Cost | Runs on logs the size of Road Traffic Fines inside a sweep | Per-cell timing in run logs |

## 5. Architecture

- **`process_voids`** is the product library. **`process_voids.pvoid`** is
  the product entry point, and the only path that computes a metric and
  writes it somewhere real. It lags the research by design and will be
  revised once the metrics settle. Shared metric machinery belongs in
  `process_voids`, so that pvoid and the lab use one definition of each
  metric.
- **`lab`** is the research harness: discovery combos, degradations,
  fixtures, the metric registry, experiment runners and result CSVs.
  Lab-only concerns (experiment cells, oracle baselines, CSV writing) stay
  here.
- **skip-alignments** is an external engine, pinned in `pyproject.toml`. It
  provides alignments, skip alignments, executions, and skip probabilities
  via ebi. Bugs in it are reported upstream, with regression tests written
  there; process-voids neither patches nor duplicates it.
- **Two alignment representations coexist.** Skip alignments (the lumped
  normal form) feed skip probability and the skip-alignment metrics. A lumped
  skip is an execution of the lumped node and its ancestors only, per
  Definition [Executions]: a descendant beneath it has no execution there, so
  both factors of a skip probability times mass metric count the same traces.
  Classical Petri-net alignments feed voidmass and `alignment_coverage_pn2`,
  because they name every missing activity. Classical alignments are
  translated into the skip-alignment path shape, so a single execution and
  averaging implementation serves both.
- **Experiment design.** For each (log, discovery combo), a reference model
  is discovered once from the undegraded log and then held fixed. Each
  degradation dimension and level degrades the log and scores it against that
  model, which gives a dose-response curve. Fixtures with known ground truth
  (the payment running example and the synthetic claims log) anchor the
  metrics to numbers that can be checked by hand. Results go to `var/`. Runs
  are named in `lab/runs.py`, so what was run is committed alongside the code
  that ran it.
- **One runner.** Experiments are being consolidated into a single
  registry-driven runner. Harness features are designed for every experiment
  at once, and no new `exp_*` entry points are added.

## 6. Guidance

- **A column name is a promise.** The worst failure mode here isn't a crash.
  It's the wrong number under the right name: see the `skipprob` and
  `alignment_coverage_pn` fixes in 0.4.1, and the stale weights in 0.4.2.
  Name each quantity for what it computes (`mean_leaf_skipprob`), cite the
  formal definition it implements, and pin a hand-computed worked example in
  a test.
- **Registry ids are append-only.** Retire ids; never delete or reuse them. A
  changed meaning gets a new id or a `history` entry.
- **Schema as code.** Anything schema-shaped (metric ids, column sets, run
  configurations) is an executable constant with a drift test. Prose,
  including this document, is downstream of it.
- **Reference models are read-only during a sweep.** Per-cell quantities
  (weights, skip probabilities) belong in maps keyed by node, not in
  attributes written onto a shared tree.
- **Pin the reference model across runs.** Discovery isn't reproducible
  between processes by default: a hash-seeded tie-break in pm4py's cut
  selection can change the discovered tree. Any comparison spanning runs
  needs the same tree, either cached as `exp_surprise` and `exp_voidmass` do
  under `var/lab/tree_cache/`, or written to ptml and reloaded.
- **Bound what you can't compute.** When a search times out, report provable
  bounds, and say how much of the log is affected (`timed_out_weight`).
- **Cache per report row, not per node.** Structural maps and per-path
  results are computed once per (tree, log) and shared across nodes. Every
  performance fix so far has come down to this.
- **Docstrings describe the present.** History and rationale belong in
  CHANGELOG or commit messages. Never cite documents that aren't in source
  control: this repository is public.
- **When an implementation narrows or contradicts the agreed design, stop and
  raise it.** A test that pins surprising behaviour is no substitute for that
  conversation.

## 7. Open questions

These are unresolved. They're listed so that nobody mistakes the current code
for a decision.

1. **Which metric is the headline?** The candidates are `voidmass_process`
   (decomposable, interpretable units, expensive to compute), `voidsat`
   (time-based, unaffected by lumping), and a skip-weighted count. This waits
   on the experiments.
2. **Where the product-form double count remains.** Coverage by Alignment
   Correspondence now conditions its mass on observation
   (`alignment_coverage_pn2`), so absence is counted once, by the skip
   probability. `salign_coverage` and `exp_voidmass`'s `voidage_*` family
   still count it twice, and are pinned in the extremes test's
   `KNOWN_FAILURES`. Condition them the same way, or retire them.
3. **How should a skip move share a time gap?** The definition and the tests
   split the gap with the following event, but the Coverage By Aligned
   Duration docstring in `coveragemass` says the skip takes the whole gap. In
   concurrent regions, the normal form's ordering convention, not evidence,
   decides which gap a skip lands in. A move also gets no duration where
   nothing follows it to consume, so an always-missing submodel at the end of
   a process reads 0, exactly as a fully observed one does.
4. **Model provenance.** By principle 2, voids are meaningful only against
   external obligations, yet nothing records which nodes of a hybrid model a
   person asserted and which were discovered.
5. **What does the product show?** pvoid prints weight coverage and coverage
   by duration, and `bpmn_colour` paints skip probabilities onto tasks.
   Neither is yet a subprocess-level void metric.

## Further reading

- Skip alignments and skip probabilities: Bär, Burke, Wynn and Leemans,
  *Skip Probabilities for Subprocesses*. Implementation:
  [skip-alignments](https://github.com/adamburkegh/skip-alignments).
- Probabilistic process trees and the Toothpaste miner (the `toothpaste`
  discovery combo): Burke, Leemans and Wynn, Petri Nets 2021.
- Background: [Process Voids: Data Science Without Data](https://adamburkeware.net/2025/11/11/pvoid-adsn.html)
  (ADSN 2025).
