# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

Now requires Python 3.11.

### Added

* A property-based test of Lemma [Additivity] for Definition [Void by
  Process-Relative Alignment Moves] - `voidmass_process` from classical
  alignments (`voidmass_pn.voidmass_table_pn`): over any antichain of nodes
  covering the labelled leaves, values sum to the root's. `hypothesis`
  generates the tree (sequences, choices, parallel blocks, loops and silent
  leaves), the log (activities dropped, reordered, repeated, and ones the
  model does not have) and the cut (silent leaves may be left out). A
  second property checks the test has teeth: leaving out a labelled leaf
  that carries deficit leaves the sum short, which hypothesis would fail
  as unsatisfiable were the generated logs all perfectly aligned. The
  lemma holds only where leaf labels are distinct: `terms_by_node`
  credits a move to every node whose leaves carry its label, so on
  `seq(a, a)` against `<a>` both leaves claim the one missing `a` and sum
  to 1.0 where the root reads 0.5. That case is pinned at what the code
  reads today. `hypothesis` joins the dependencies for this.

* `lab.runtime_table`: metric x (log, combo) -> total wall-clock seconds,
  read from each log's own `*_timings.csv`. A metric's own row is
  usually near-zero - the real cost sits in shared stages (`classical`,
  `dv`, ...) several metrics depend on (`lab.run`'s `ProcessMetric.
  needs`) - so this sums a metric's own row plus its declared stage
  dependencies, reading as "what this metric alone would have cost", not
  additive across metrics without double-counting a shared stage two of
  them both need. Outputs both Markdown and LaTeX. Same voidsalign2/
  voidsalign3 fallback as `lab.cross_log_plots`, since different logs'
  sweeps can be pinned to different registry versions.
* `lab.current_coverage`: "are we current" for a log - a table of today's
  live `exp_disco_degrade` metrics against each log's own MOST RECENT
  result CSV (by file mtime), one row per metric, one column per log.
  Deliberately narrower than `lab.collection_report` (which answers "has
  this ever been collected, across every result CSV ever written",
  including retired ids and long-stale runs by design): a log can show
  Done there from a sweep that predates a metric's current form, which
  is the wrong view for deciding what to rerun. A metric absent from a
  log's newest file reads as not done here even if an older file for
  the same log once had it.
* `lab.cross_log_plots`: dose-response plots comparing DIFFERENT logs'
  own result CSVs side by side - `lab.plots` compares combos or
  dimensions within one log's own CSV, not across logs. Three axis
  arrangements via `--panel-by`: `metric` (default, panel per metric,
  line per log), `log` (panel per log, line per metric), `log_metric`
  (panel per log/metric pair, line per combo - figures then split by
  degradation dimension only, since combo is a line here rather than a
  figure axis). `--ylim MIN MAX` fixes every panel's y-axis instead of
  autoscaling per panel. A log's CSV missing a metric's column (a sweep
  predating that metric) is skipped for that line, not an error, since
  different logs' sweeps can be pinned to different metric-registry
  versions.
* `lab.check_run`: a quick sanity summary of a result CSV - row count,
  status breakdown, and which logs/combos/degradation dims it covers.
  The "did this run go the way I expected" question after a sweep,
  answered once instead of as a fresh `python -c` snippet each time. A
  column a given CSV doesn't have (eg a `*_timings.csv` has no
  `combo`/`degradation_dim`) is skipped, not an error.

* `pvoid.toml`, a gitignored settings file for machine-specific paths,
  read by `process_voids.config`, with a tracked `pvoid.example.toml`
  documenting every key. It replaces constants in code: `lab.params`'
  large-log paths, `lab.toothpaste_bridge.TOOTHPASTE_DIR`, and
  `process_voids.pvoid`'s `EBI_EXECUTABLE`. Keys are declared in
  `config.SCHEMA`, and a test holds the example file to it; an unknown
  section or key is rejected rather than ignored, so a misspelt key cannot
  read as unset. Lookup tries the running checkout's root, then the main
  checkout's, so one file serves every worktree and a worktree may
  override it. Nothing reads a key until it is needed: a machine without
  the file imports everything and runs the suite, and reading an unset key
  with no default raises `ConfigError` naming the key, what it is for, and
  the paths searched. A large log is registered as
  `lab.params.ExternalLog('<filename>')` and resolves under
  `[paths] data_dir` only when opened or named, with forward slashes.
  `[tools] ebi` is optional and defaults to bare `ebi` on PATH - the
  skip-alignments default - so an installed package with no file still
  finds it. An `[agent]` section holds machine facts that only agents and
  people read - the interpreter to build a worktree venv with, and where
  the papers are - and a test holds that no source module reads it.

* `voidsalign3` (`process_voids.voidsalign3`), Definitions [Coverage by
  Skip Alignment Correspondence] and [Void by Skip Alignment
  Correspondence]: `voidsalign2` with one term changed. A skip move now
  counts for `|leaves(msub) \ {silent}|`, the number of labelled
  activities in the subprocess it skips, rather than `aligncost(<>,
  msub)`. The cost-minimal execution against the empty trace falls below
  a subprocess's own activities wherever a traversal performs fewer of
  them - to nothing for an optional block, so a skipped optional block's
  absence was invisible in the mass. The weights now differ for
  `xor(a, tau)` (0 -> 1), `xor(a, b)` (1 -> 2) and `loop(a, e)` (1 -> 2),
  and agree for a single activity or a sequence. On the payment running
  example only the root moves, since `<o, s, p>` aligns with a skip of
  the approval loop: its void reading goes 0.0833 -> 2/15. Registered
  `scale='void'`, so `test_metric_extremes` holds it to 0, 1 and 0.5.
  The dose-response plots picked it up from the registry with no edit.

* `matchprob`, `1 - skipprob` at the scored node, emitted per node by
  `lab.run` and at the root by `lab.metrics`, registered live with
  `scale='coverage'` - so `test_metric_extremes` holds it to 1, 0 and
  0.5 through the real pipeline. It carries nothing `skipprob` does not;
  it exists so a consumer wanting skip probability coverage-way up reads
  a registered column rather than inverting `skipprob` itself. The first
  such consumer is the dose-response plots, which invert `skipprob` by
  name today.
* `process_voids.util.squash_review`: the whole squash review as one
  command, defaulting to `HEAD` against the index like the other review
  tools. Prints the diff stat, `registry_diff`'s and `test_name_diff`'s
  reports (called as functions, not shelled out to), a scan of added lines
  for anything that must not reach this public repository, and a note
  where `src/` changed without `CHANGELOG.md`. The scan's patterns are a
  module constant, each commented with what it catches: references to the
  private lab notebook, the gitignored paper and report directories, the
  private archive repository, per-machine agent configuration, absolute
  Windows paths, and merge conflict markers at the start of a line.
  Only the tool and its tests are allowed any pattern, since they have to
  spell each one out. The CHANGELOG note is worded as a report - a
  docstring-only change needs no entry. Like the others it reports rather
  than gates, exiting non-zero only on an internal error, and it does not
  run the test suite.
* `partial_sequence`, a third fixture (`lab.fixtures.
  build_partial_sequence_log`/`build_partial_sequence_tree`,
  `data/partial_sequence.xes`): `seq(o, seq(x, y, z), p)` over three
  traces recording the inner sequence completely, without its middle
  step, and without two of three. Named for what it tests rather than a
  domain, since it models nothing. Every trace has exactly one optimal
  alignment, so the inner sequence reads partial 0, 1/3 and 2/3
  undiluted - which `payment_partial` cannot show, because completing a
  loop traversal and discarding an escalation cost the same there, so
  the partial explanation is always one of two ties. The subprocess is
  traversed in all three traces, so its skip probability is 0
  throughout and only completeness varies.
* `payment_partial`, a second running-example fixture beside
  `payment_approval` (`lab.fixtures.build_payment_partial_log`/
  `build_payment_partial_tree`, `data/payment_partial.xes`, registered in
  `lab.params.ALL_LOGS`). Same story, in the shape `inductive_noise20`
  discovers from it - `seq(o, xor(tau, loop(a, e)), s, p)`, hand-built so
  a test against it does not depend on discovery's cross-process
  tie-break - and two extra traces where a subprocess is TRAVERSED BUT
  RECORDED INCOMPLETELY. Every absence in `payment_approval` is a whole
  traversal, so skip probability accounts for all of it and a partial
  term has nothing to fire on; in `sigma7` the loop's closing `a` is
  unrecorded, so the loop ran (skip probability 0) while a third of its
  moves went unobserved. `sigma8` is the same variant with different
  timings, for the duration-based metrics. `payment_approval` is
  unchanged, so numbers pinned against it stay valid.
* `voidsalign2` (`process_voids.voidsalign2`), Definitions [Coverage by
  Skip-Weighted Alignment Correspondence] and [Void by Skip-Weighted
  Alignment Correspondence]: `1 - (1 - skip_prob) * mass`, where the mass
  averages `smatchcount/smovecount` over skip alignments, counting each
  skip move as `aligncost(<>, msub)` - the least number of labelled
  activities any traversal of the skipped subprocess performs - so a
  lumped skip over a large subtree is not underweighted, and a skip over
  a wholly silent subprocess counts 0 without a separate silent-move
  exclusion. The mass is conditioned on observation exactly as
  `alignment_coverage_pn2`'s is (executions with no synchronous move
  excluded, traces observing the node nowhere dropped from both the
  average and the normaliser `W`, coverage 0 where `W` is 0), which it
  gets by calling `coveragemass.observed_alignment_mass` with the
  skip-weighted ratio rather than repeating that arithmetic.
  `observed_alignment_mass` and `_alignment_values` gained a `ratio=`
  parameter for this; their default behaviour is unchanged.

* `voidsat2` (`process_voids.voidsat2`), Definition [Coverage and Void by
  Aligned Duration]: `1 - (1 - skip_prob) * mass`, where the mass averages
  `obsdur / (obsdur + misdur)` - the share of a subprocess's own aligned
  elapsed time that a synchronous move accounts for. A rate of the
  subprocess's time, where the retired `voidsat` divided by the whole
  trace's duration and so reported a share. Conditioned on observation as
  `alignment_coverage_pn2`'s mass is: only alignments where the subprocess
  has a synchronous move count, traces observing it nowhere leave both the
  average and the `obscount` normalising it, and the mass is 0 where
  `obscount` is 0. It does not call `observed_alignment_mass` despite that
  identical conditioning, because it iterates real traces rather than
  deduplicated weighted variants - duration is per-instance - and takes a
  ratio of summed totals rather than a mean of per-execution ratios.
  `coveragemass` gained `obsdur`, `misdur` and `has_synchronous_move`
  beside the existing move-duration primitives. A subprocess observed but
  unmeasurable - recorded at the very start of a trace, where no interval
  bounds it - reads ratio 1 rather than 0.

* `run.sh --seed N` sets `PYTHONHASHSEED` before invoking the command, so a
  first discovery (a tree-cache miss) can be made reproducible. It has to
  live in the wrapper script rather than as a Python CLI flag: the
  interpreter reads `PYTHONHASHSEED` at startup, so by the time a `main()`
  could parse a flag its own hash seed is already fixed. Recognised in
  first position only, so it can't swallow an argument meant for the
  command.

* `lab.run_history`: a run-history CSV (`var/lab/results/run_history.csv`
  by default), one row per run, appended and never overwritten. Holds the
  run timestamp, `out_csv` as the join key back to the results, the
  resolved config (logs, combos, degradations, levels), which metrics were
  scored and which were therefore excluded, both packages' versions with
  git commit and dirty state, and `PYTHONHASHSEED` as the interpreter
  actually saw it - recorded as `unset` when it was not set, since a blank
  cell reads back as NaN and can't be told from a column that was never
  written. `lab.run` appends a row after the result CSVs are safely on
  disk; neither building nor writing the row can fail a completed run.
  Deliberately not an upsert like `_merge_write`: a rerun to the same
  `out_csv` is a second run, not a correction of the first.

  Paths in both this file and the results' own `tree_cache_file` column
  are written with forward slashes rather than the producing machine's
  separator, so a result set stays legible - and joinable against a path
  typed by hand - away from Windows.

* `activity_frequency_gradual`, a third degradation dimension, on the
  default roster alongside `activity_gradual` and `trace`: activities
  drop rarest-first by event count, and the level is the fraction of the
  log's total event volume removed rather than a count of distinct
  labels. `activity_gradual`'s drop order is a seeded shuffle, so its
  dose is one confounded sample - which activities were hit does more
  work than how much was removed, and a structural response can't be
  told from an unlucky draw without repeating across seeds. This one is
  deterministic. Both are kept, so earlier result CSVs stay comparable.
* `lab.collection_report`: which `lab.metric_registry` ids have real
  (non-null) values recorded for which logs, across every root-level
  result CSV in `var/lab/results`, as a Markdown table. Answers "has this
  been collected at all", not what its latest value is or whether the id
  is still live - a retired id still reads as collected from before it
  was retired. Rows are attributed by a CSV's own `log` column rather
  than its filename, so ad hoc and probe runs count like a named sweep.
* `bpic2020_rfp` joins the log catalogue and gains a named run.
* Every per-node result row carries `tree_source` and `tree_cache_file`:
  whether that cell's tree came from the shared discovery cache or was
  discovered by this run, and which cache file. `lab.discovery.
  discover_cached` returns a `CachedDiscovery` (tree, ppt_weights,
  source, cache_path) rather than a bare pair. Run-level bookkeeping
  rather than metrics, so no `lab.metric_registry` entries - but the
  cache is what makes tree identity ambiguous, so a result file has to
  answer it on its own.
* `lab.print_tree`: prints what a combo discovered on a log, in
  skip-alignments' own process-tree notation, without recomputing skip
  probabilities. It reads through `lab.discovery.discover_cached`, so the
  tree shown is the one an experiment against the same `(log, combo)`
  would score.

* Three squash-review reports in `process_voids.util`, beside
  `release_check`, each runnable as a single command and each defaulting
  to `HEAD` against the index (the staged change under review):
  `registry_diff` (ids added/removed, per-id changes to status,
  superseded_by, scripts, source or scale, and whether any
  `lab.metric_registry` history entry was extended, altered or dropped -
  only the last two break the append-only rule); `test_name_diff` (which
  `def test_*` names exist at the base ref and not the target, with names
  that merely changed file reported as moved rather than lost - a count
  can rise while behaviour tests disappear); and `branch_survey` (every
  local branch's worktree, whether its content has already landed, and
  whether it would squash cleanly or in which files it conflicts -
  "commits ahead" means nothing in a squash workflow, where a landed
  branch still reads as ahead). They report rather than gate: each exits
  non-zero only on error, never on findings. `release_check` remains the
  gate and is unchanged. `process_voids.util.gitread` holds the shared
  "read this path at that ref" helper, where the index is a ref whose
  prefix is empty.

* `lab.run`: a consolidated experiment entry point (`cells = log x combo x
  degradation dim x level`), with `--metrics` selection (by id or group -
  `skip_alignment`/`classical`/`aligned_duration`) and `--no-degradation`
  (one row per `(log, combo)`, no dim/level loop at all - the "fixed model,
  as-is log" case). Only the stages a selected metric actually needs are
  triggered per cell, so a selection that needs no classical alignment
  genuinely skips that search rather than computing and discarding it.
  `lab.exp_disco_degrade`'s metric roster and cell loop now live here -
  `run_disco_degrade` is a thin wrapper over `lab.run.run`, keeping its own
  stable CLI/signature and default filenames, not a second copy. The
  wrapper's `--exclude-metric` maps onto `metrics=`, which accepts either a
  list of `ProcessMetric`s (that case) or the CLI's inclusive id/group-name
  shape. Discovery goes through
  `lab.discovery.discover_cached`, so every experiment entering through this
  runner gets the same cached, process-stable tree.

### Changed

* Now requires Python 3.11 or later, previously 3.10 (via tomllib). 

* `voidsalign2` is retired, superseded by `voidsalign3`, which `lab.run`
  and `lab.metrics` now emit in its place. Result CSVs written before
  this carry the `aligncost`-weighted reading under the `voidsalign2`
  column; the registry keeps that id with its description.
  `process_voids.voidsalign2` and its tests stay in place, off the
  roster.
* `lab.plots` reads its panels from `lab.metric_registry` instead of
  keeping its own list: one panel per live metric with a scale that
  `exp_disco_degrade` emits and that is plotted, with a
  `<base>_lower`/`<base>_upper` pair drawn as one banded panel.
  `lab.plots.METRICS` is gone; `metric_panels()` replaces it. The
  hand-written list had gone stale four times - each time a metric was
  renamed or retired the panel stayed, and plotting a current result CSV
  raised a `KeyError` - and the suite never caught it, because its
  fixtures hand-typed the same names. They are now built from
  `metric_panels()` too, and a test renames and retires a metric in a
  synthetic registry to show the plots follow with no edit.
* `lab.metric_registry.Metric` gains `plotted` (default `True`), `False`
  only where another live metric already shows the same information.
  `skipprob` is the one such metric: its panel is gone, and `matchprob`
  shows it coverage-way up instead of `lab.plots` inverting `skipprob` by
  name. `process_voids.util.registry_diff` compares the new field, and
  reads a ref from before it existed as `plotted=True`, so a diff across
  its introduction reports nothing spurious.
* `weight_coverage` is retired. `lab.run` and `lab.metrics` no longer
  emit it, so result CSVs written from here carry no `weight_coverage`
  column; `weight_voidage` stays, and is what the dose-response plots now
  show in its place. The two were always exact complements
  (`weight_voidage = 1 - weight_coverage` at every node), so a value from
  an earlier result file can still be read against a later one as
  `1 - weight_voidage`. The registry keeps the id with its description and
  history, as it does every retired id; there is no `superseded_by`,
  since `weight_voidage` is its mirror rather than a replacement.
* `bpi2013_closed_problems` replaces `bpi2013_incidents` in the log
  catalogue and as a named run, so `--run full` now sweeps closed
  problems rather than incidents.
* `sepsis` dropped from the log catalogue: both combos tried against it
  failed structurally (`ShuffleExplosionError` on `inductive_noise20`,
  a plain `MemoryError` after 55 minutes on `toothpaste`), the same
  treatment `bpi2013_incidents` got for the same reason.

* `coveragemass.block` no longer excludes silent (`TauPath`) moves, and
  `coveragemass.mdur` no longer returns 0 for one. Skip alignments have no
  silent move type - Definition 5 of the skip-alignment paper makes every
  non-synchronous model-side move a skip, and every skip a deviation - so
  the exclusion was dropping real deviations from the block that shares a
  gap, and inflating what the surviving moves were each charged. The
  clause it was written against also described a different thing from what
  the code tested: a zero-cost skip, not a `TauPath`. Every metric built on
  move durations (`adur`, `admass`, `covat`, `voidat`, `voidsat`) changes
  value where an alignment contains a silent move.

* `voidsat` is retired, superseded by `voidsat2`, and the runner scores
  `voidsat2` in its place at the root and per node. The retired form was
  `skip_prob` times a share of the whole trace's duration, so an entirely
  missing subprocess - whose mass is 0 - read void 0, the reverse of the
  truth. It was pinned in `test_metric_extremes`' `KNOWN_FAILURES` at 0.0
  for `always_missing` and 0.5 for `half_missing`; `voidsat2` meets both
  promises and needs no pin. Result CSVs written before this carry the
  share under the `voidsat` column; the registry keeps that id with its
  description, as it does every retired id.

* `voidsalign` is retired, superseded by `voidsalign2`, and the runner
  scores `voidsalign2` in its place at the root and per node. The retired
  form was `skip_prob` times a SIZE share of the whole tree's
  skip-weighted move volume, which is not what the definition specifies:
  it reports a share of the process rather than a rate, so it has no
  reading at which a wholly missing subprocess is 1. Result CSVs written
  before this carry the size share under the `voidsalign` column; the
  registry keeps that id with its description, as it does every retired
  id. `process_voids.voidsalign` and its tests stay in place, now unused
  by the roster.
* The discovered-tree cache is one shared `lab.discovery.discover_cached`,
  used by `exp_disco_degrade`, `exp_surprise` and `exp_voidmass`, instead
  of a copy in each of the latter two. `exp_disco_degrade` now caches its
  discovery too, so a tree is discovered once per `(log, combo)` across
  every script and every process. That also fixes which tree a rerun
  scores: pm4py's Inductive cut selection is hash-seed dependent where a
  `noise_threshold` cut sits near a tie, so rtfm's `inductive_noise20`
  draws a 12- or 13-node tree from the same log depending on the process,
  and two runs of nominally the same experiment could score different
  models with no warning. Caching freezes that tie-break rather than
  resolving it: the cached tree records which model a result was scored
  against, and is arbitrary, not authoritative.
* Cache files now hold a `(tree, ppt_weights)` pair, since
  `exp_disco_degrade` threads toothpaste's fixed PPT weights into
  `pvoid.skipprob`, and all three scripts share these files. They are
  written under a new `<log>__<combo>__pair.pkl` name: a file written
  when the cache held a bare tree unpickles without error but unpacks
  into the wrong shape, so the old files are ignored rather than
  misread. Delete `var/lab/tree_cache/*.pkl` without the suffix at
  leisure; nothing reads them any more.

* `lab/metrics.py`'s module and `mean_leaf_skipprob` docstrings trimmed -
  they re-explained what each metric means, duplicating
  `lab.metric_registry`, which now owns that.

* The alignment search timeout is now two named constants rather than one
  `lab.run` constant reached through `lab.exp_disco_degrade` (which would
  have kept that wrapper alive as an import shim) plus a bare `100`
  repeated across three `process_voids.voidmass_pn` signatures. They are
  different decisions: `voidmass_pn.DEFAULT_ALIGNMENT_TIMEOUT` is the
  fallback that module uses when a caller states no budget, so it imports
  standalone; `lab.params.CLASSICAL_ALIGNMENT_TIMEOUT` is the budget this
  lab allows, chosen for these logs and this hardware, and is what
  `lab.run` and `lab.mass_term_probe` pass. The lab's applies to
  skip-alignments' `align_sk_all` as well as the classical search, so it
  was never a property of `voidmass_pn`'s own search.

* `voidmass_subprocess_lower`/`voidmass_subprocess_upper` dropped from
  `lab.run`'s default roster and retired in `lab.metric_registry` - not a
  candidate for inclusion, and not informative in the dose-response plots.
  `lab.plots` loses its `voidmass_subprocess` panel with them.
  `process_voids.voidmass_pn.voidmass_table_pn` still computes both fields,
  and the registry entries stay, so old result CSVs carrying those columns
  remain readable.

### Fixed

* `process_voids.util.squash_review` given two refs compared their trees,
  so reviewing a branch against a trunk that had moved on since the branch
  was cut reported the trunk's newer work as the branch removing it -
  registry ids as removed, tests as dropped. It now compares from their
  merge base, and says so in the report's first line when that differs
  from the base named. HEAD against the index, the default, is unchanged.
  Its allowance for absolute paths in `lab/params.py` is gone: machine
  paths live in `pvoid.toml`, so one reappearing there is now reported.

* `lab.print_tree` raised `TypeError` on every invocation: it unpacked
  `discover_cached`'s return as a `(tree, ppt_weights)` pair after that
  became a `CachedDiscovery`. Its tests mocked the old tuple, which let
  the unpack pass; they now mock a `CachedDiscovery`.

* `lab.plots`' panel list still named `voidsalign` after it was retired
  in favour of `voidsalign2`, so plotting a result CSV written since
  raised a `KeyError`. The panel now reads `voidsalign2`, and `voidsat2`
  gains a panel of its own.

* `lab.plots`' panel list still named the `alignment_coverage_pn_*`
  columns retired in 0.5.0, so plotting any result CSV written since
  raised a `KeyError`; its own test fixtures used the same retired names,
  which is why the suite stayed green. The panels now read
  `alignment_coverage_pn2_*`, and `voidsalign` gains one.

* `lab.run`'s root-level row always included every `ALL_METRICS` id
  regardless of a restricted `metrics=` selection - `_compute_cell` built it
  from `ALL_METRICS` instead of the `metrics` it was actually passed, so an
  excluded metric's column stayed present (though unpopulated) instead of
  being absent, the one thing `--exclude-metric` is for.

* `lab.run`'s `_null_metric_values` derives its null fallback from the
  metrics actually scored, so an error row can't carry a key no successful
  row has. Previously a second copy of this lived in
  `lab.exp_disco_degrade`.

## [0.5.0] - 2026-09-12

### Added

* `lab.metric_registry`'s `Metric` now carries `status` (`live`,
  `evaluation`, `product-only`, `retired`), `superseded_by`, and `history`
  (commit/version -> prior meaning), so the registry is an append-only
  historical record of every CSV column this package has ever written,
  not just the currently-emitted ones. Restored `voidmass_deficit`,
  `node_skip_prob` and `alignment_coverage` as retired ids; registered
  `duration_coverage` (`process_voids.coveragemass.coverage_by_duration`)
  as product-only; and registered `lab.exp_voidmass`'s previously
  unregistered ids (`skip_prob`, `deficit`, `movecount`,
  `voidmass_subprocess`/`voidmass_process`, `voidage_subprocess`/
  `voidage_process`, `target_voidmass_*`/`target_voidage_*`,
  `target_rank_*`, `n_optimal_alignments`). `lab.exp_voidmass` gained
  `NODE_METRIC_KEYS`/`SUMMARY_METRIC_KEYS` constants so its ids are
  covered by the same registry drift test as the other experiment
  scripts.
* `process_voids.metric_context`: `ProcessMetric` (id, scope, needs, compute) and
  `CellContext`, a per-cell object with named, lazily-computed and
  memoised stages (`dv`, `executions_cache`, `traces`,
  `aligned_duration_cache`, `surprise_self`) and `stage_started`/
  `stage_finished`/`stage_failed`/`metric_started`/`metric_finished`/
  `metric_failed` lifecycle events. `CellContext.score` isolates a
  metric's own exception to a sentinel (`METRIC_ERROR`) rather than
  raising; a stage's failure is memoised and re-raised (not recomputed)
  on every later access within the same cell, so an expensive, failing
  stage runs at most once per cell even when several metrics need it.
  `lab.timing` gains `TimingListener`, collecting one long-form timing
  row per stage/metric actually computed in a cell from these events.
  `CellContext` also gained a `classical` stage (`voidmass_pn.
  voidmass_table_pn`) and a `score_all(ctx, metrics, node)` helper -
  scores a list of `ProcessMetric`s against one node, keeping
  `ProcessMetric` itself a single declared quantity rather than a bundle.
* `release_check` fails on leftover merge conflict markers (`<<<<<<<`/
  `>>>>>>>` at the start of a line) in tracked files.
* `lab.metric_registry`'s `Metric` carries a `scale` (`coverage`, `void`
  or `void_share`): what a metric promises to read at the node it's
  scored at when that submodel is never missing, always missing, or
  missing from half the traces. `tests/lab/test_metric_extremes.py`
  holds every live metric with a scale to that promise on `seq(a, b)`,
  computing each through the experiment's own `ProcessMetric`
  declarations and the real skip-probability pipeline (the extremes
  criterion in `docs/DESIGN.md`).
* `voidsalign` (Voidage by Skip-Weighted Alignment Moves,
  `process_voids.voidsalign`): like `salign_coverage`, it works over
  skip alignments, but each skip move is weighted by the minimum
  number of activities its subprocess performs, so a lumped skip over
  a large subtree isn't counted as a single move. `skip_prob` times a
  size SHARE of the whole tree's skip-weighted moves, not a
  match/movecount completeness ratio - a ratio collapses to zero
  exactly where a subprocess is wholly missing, which a voidage metric
  for missing subprocesses can't afford. Emitted by `exp_disco_degrade`
  at the root and per node.
* `lab.exp_disco_degrade` now runs every metric through `CellContext`/
  `ProcessMetric` (`ALL_METRICS`) instead of hand-assembling each row at
  its call site - the root row is the same per-node scoring at
  `node=tree`, not a separate computation, and the three independently
  hand-typed "every metric is None" fallback dicts are one
  `NULL_METRIC_VALUES` constant. `run_disco_degrade` now also returns
  (and writes) a third, long-form `_timings` CSV with one row per
  stage/metric actually computed per cell.
* `run_disco_degrade(..., metrics=...)` and `lab.exp_disco_degrade`'s
  `--exclude-metric` flag skip scoring the given metrics for a run (e.g.
  `voidsat`, whose per-trace cost is prohibitive on high case-count
  logs). An excluded metric gets no timing row and an empty column in
  every row; unknown ids are rejected.

### Changed

* `\covermove`'s mass is conditioned on observation (Definition
  [Coverage by Alignment Correspondence]): executions with no
  synchronous move are excluded, and so are traces whose alignments hold
  no observed execution of the submodel, with the remaining weights
  renormalised. Absence is then counted once, by the skip probability,
  instead of twice, so coverage falls linearly in a submodel's absence
  rather than as `(1 - p)²` — a submodel missing from half the traces
  reads 0.5 where it read 0.25. The ids change with the definition:
  `alignment_coverage_pn2_lower`/`_upper` replace
  `alignment_coverage_pn_lower`/`_upper`, now retired, so results
  written under the two definitions can't be confused.
  `coveragemass.observed_alignment_mass` is the new mass term;
  `alignment_mass` is unchanged and still backs `salign_coverage`.
* skip-alignments is now a published dependency, `skipalignments==0.3.0`,
  in place of the `v0.2.3+p4` git tag. 0.3.0 is the released form of the
  skip-probability fix recorded under Fixed below, plus alignment-search
  performance work; per-node skip probabilities and the metrics built on
  them are unchanged from `v0.2.3+p4` on the payment, claims and rtfm
  fixtures. The version is pinned exactly rather than to a range while
  the two projects move together, so a result CSV's recorded
  skip-alignments version is unambiguous. Direct-reference dependencies
  block `release_check`, so this is also what makes a process-voids
  release possible.
* Tests call skip-alignments' `Aligner.align_normal_form` instead of the
  deprecated `align2`, clearing its deprecation warnings from the suite.
* `lab.metric_registry` records the skip-probability correction below as a
  `v0.5.0` history entry on every id built on `dv.skip_probs`, alongside
  the executions entry the ids computed from `coveragemass.executions`
  already carry, so a pre-0.5.0 result CSV's columns are still readable
  from the registry alone.
* The `smoke` run uses `payment_approval` instead of `rtfm`.
* `coveragemass.make_aligned_duration_cache` takes an optional `traces=`
  parameter, for a caller that already has a log's `log_to_traces` result
  and wants the cache to reuse it instead of recomputing it.

### Fixed

* Upgraded to skip-alignments 0.3.0, which fixes a
  skip-probability bug present since skip-alignments 0.2.0: a node
  nested inside a subtree that an alignment skips as a single lumped
  move was counted as skipped whenever anything else in that alignment
  synchronised, although it has no execution there at all. Skip
  probabilities of such nested nodes were overstated (e.g. the children
  of a lumped `seq(b, c)` read 1/2 instead of 0); the lumped node's own
  value is unchanged. On the payment running example the approval
  loop's children `a` and `e` go from 1/3 and 2/3 to 0, and the
  scheduling choice's `s` and silent branch from 1/6 and 1 to 0, while
  the loop (1/3) and the choice (1/6) keep theirs. Metrics built on skip
  probabilities (`skipprob`/`skip_prob`, `weight_coverage`/
  `weight_voidage`, `mean_leaf_skipprob`, `salign_coverage`,
  `alignment_coverage_pn_lower`/`_upper`, `voidsat`, the voidage columns
  `voidage_subprocess`/`voidage_process`/`target_voidage_*`/
  `target_rank_voidage_process`, and pvoid's `duration_coverage`) change
  at such nodes, and values aggregated over descendants can change too.
  `weight_coverage`/`weight_voidage` read skip probabilities only at
  leaves, so they no longer register a subtree skipped as one lumped
  move at all. Results written before this upgrade carry the overstated
  values.
* `voidsat` was very slow on high case-count logs: `admass` rebuilt the
  log's trace list (a full group-by, sort and conversion over the log)
  once per tree node instead of once per report row - about 4,270s
  extra per cell on rtfm. The trace list is now built once, with the
  per-row cache, and reused for every node.
* `lab.plots` read `alignment_coverage_pn`/`voidmass_subprocess`/
  `voidmass_process` as single columns, which no longer exist since 0.4.2
  split them into `_lower`/`_upper` bounds. Each is now plotted as a
  midpoint line with a shaded band between its bounds, which collapses
  to a plain line when nothing timed out.
* The skip-alignment metrics counted an entirely unwitnessed subtree's
  absence again at every node beneath it. The aligner records such a
  subtree as one lumped skip move on its coarsest node, and
  `coveragemass.executions` (and voidsat's equivalent) treated that move
  as an execution of every descendant on a mandatory position
  (Sequence/And child, Loop do-child) too. Definition [Executions],
  which the metrics cite, gives such a descendant no execution there: it
  was neither traversed nor skipped, and the lumped node alone carries
  the void. With the skip-alignments upgrade above, skip probabilities
  follow the same rule, so both factors of every skip_prob * mass
  metric count the same traces. Affects nodes under a lumped subtree in
  `salign_coverage`, `voidsalign` and `voidsat` (per-node CSV) and in
  every `exp_voidmass` column except `skip_prob` and
  `n_optimal_alignments`. Root values, and the classical-alignment
  metrics' mass terms, are unchanged.

### Known issues

* The extremes test pins the live metrics that fail it.
  `salign_coverage` and `lab.exp_voidmass`'s `voidage_subprocess`/
  `voidage_process`, with their `target_` copies, read a submodel
  missing from half the traces at half the value their scale promises
  (0.25 rather than 0.5 for the coverage): the mass term and the skip
  probability both count the absence, which is what conditioning fixed
  for `\covermove`. `voidsat` reads 0 for an always-missing submodel
  that is the last activity, since a skip move with no following event
  gets no duration.

## [0.4.2] - 2026-09-11

### Added

* `\voidsat` (Coverage by Aligned Duration, computed over skip
  alignments): a real-elapsed-time mass estimate, attributing each
  alignment move its share of the time between observed events and
  averaging over every actual trace instance rather than deduplicated
  variants (duration is a per-instance quantity, unlike every other
  metric here). Corrects two errors found in the paper's own draft
  definition along the way: a dangling reference to an unused move-
  weighting definition, and a zero-guard that keyed off the wrong
  index, under-zeroing a leading run of model-only moves before the
  first-ever observed event.

### Changed

* Default discovery roster: vanilla `inductive` and `indulpet` are off.
  Vanilla inductive is discovered at exact fit from the log it is then
  checked against, so it absorbs dropped activities at no cost, and at
  ~142s per cell it was most of a sweep's cost; `indulpet` only ever
  produced `not_implemented` rows. Both discovery functions stay
  importable, and the smoke run now uses `inductive_noise20`.
* `voidmass_table_pn` returns a `VoidmassPnResult` (`table`,
  `skip_dict`, `timed_out_count`, `timed_out_weight`) instead of a
  `(table, skip_dict)` tuple.

### Fixed

* `weight_coverage`/`weight_voidage` in the per-node CSV could pick up
  stale weight estimates left behind by an unrelated, already-run cell,
  for a second degradation dimension reusing the shared level-0.0
  result - `mass_by_weight`/`voidage_by_weight` read weights directly
  off the shared, mutable discovered tree, which every cell's own
  weight-estimation step overwrites in place.
* `voidmass_table_pn` crashed with `ZeroDivisionError` when a variant's
  alignment search timed out to zero alignments - a real, observed
  failure mode on any log with one slow-enough variant. Such a
  variant's contribution is unknown, so each cell now reports provable
  bounds instead of a guess or a crash: `voidmass_deficit`/
  `voidmass_subprocess`/`voidmass_process` split into `_lower`/`_upper`
  columns, substituting zero or full deficit against a movecount
  bounded by `2|σ|` plus the model's cheapest path length. The bound is
  derived from the aligner's cost model, since loops make any bound
  from the tree's shape alone unsound. `alignment_coverage_pn` gets the
  same `_lower`/`_upper` split. `voidmass_movecount` stays the observed
  total from completed variants; the bounds' shared denominator is the
  new `voidmass_movecount_bound`. New `timed_out_count`/
  `timed_out_weight` columns report how many variants timed out and
  their summed probability. Where nothing times out, `_lower == _upper`
  everywhere.
* An empty `_nodes` CSV (every cell in a run errored) had no column
  header, raising `EmptyDataError` in any downstream reader expecting
  an empty-but-columned frame.

### Known issues

* `voidsat` is very slow on high case-count logs: `admass` rebuilds the
  log's trace list (a full group-by, sort and conversion) once per tree
  node instead of once per cell - about 4,270s extra per cell on rtfm.
  To be fixed in the planned single experiment runner.

## [0.4.1] - 2026-09-10

### Added

* Interval surprise metric (`process_voids.surprise`): tail-probability
  surprise, with two attribution schemes (containment, predecessor) for
  locating missing-subprocess signal in a process tree.
* `lab.exp_surprise` experiment runner, with a baseline distribution
  variant (scored against the undegraded log) alongside the
  self-estimated one.
* `inductive_noise80` and `toothpaste_noise10` discovery combo variants.
* Gradual (piecewise-linear) activity-wise degradation
  (`degrade_activity_wise_gradual`), which ramps the currently-targeted
  activity's events out continuously rather than dropping them in one
  step. This gives a visible dose-response curve on logs whose activity
  alphabet is too small for the original step-wise degradation to show
  any response between levels. Kept alongside the original, not
  replacing it, so results already gathered with the step version stay
  comparable.
* Classical-alignment voidmass metrics (`process_voids.voidmass_pn`):
  `voidmass_subprocess`/`voidmass_process`, computed via non-lumped
  Petri-net alignments rather than skip-alignments' lumped normal form.
* `lab.claims_fixture`: a synthetic claims/insurance log with known
  ground truth, for ablation testing against a known-correct model.
* `lab.metric_registry`: single source of truth tying every experiment
  script's emitted metric ids to a description, with a drift test
  against each script's own output columns.
* `weight_voidage` metric.
* `mandatory_node_count`/`total_node_count` diagnostic - flags when a
  discovered tree makes every void metric uninformative by construction
  (e.g. Inductive Miner at noise_threshold=0.0, which wraps nearly
  every leaf in an optional choice).
* `--dry-run` flag and a shared per-cell timer across experiment
  scripts.
* Node-averaged dose-response plot (`--average-nodes`) alongside the
  existing root-level one.
* Every experiment run log now opens with the process-voids and
  skip-alignments version it executed against.

### Changed

* `exp_surprise`'s CSV schema: self/baseline are now separate metric
  ids/columns rather than a shared id plus a `distribution` column.
* Default discovery/degradation rosters: dropped `inductive_noise80`
  and degradation level `1.0` (a guaranteed-empty log - wasted compute,
  not just unplottable).
* `exp_claims_degrade` retired - the claims fixture now runs through
  `exp_disco_degrade` directly (registered as a combo/degradation set).
  `exp_disco_degrade` also gained a per-node CSV alongside the
  root-level one, with every metric re-evaluated at every tree node
  rather than only at the root.

### Fixed

* Two silent-correctness gaps in interval surprise: bits from an
  activity pruned entirely out of a discovered model's alphabet were
  dropped from every node's total instead of charged to root; a value
  more extreme than anything in a reference distribution produced a
  `-log2(0)` blowup.
* Stale rows could survive a successful rerun of `exp_surprise`.
* `skipprob` was silently computed as the mean of `skip_probs[leaf]`
  over every Activity leaf in the whole tree (and ignoring its own node
  argument for scoping), not skip-alignments' own published skip
  probability. It's now `skip_probs[node]` directly, at the root and at
  every node in the per-node CSV alike. The old blended average is kept
  under an honest name, `mean_leaf_skipprob`, since it's not yet
  established whether it's informative on its own. Any results CSV
  written before this fix has the wrong quantity under the `skipprob`
  column.
* `alignment_coverage_pn` (the paper's `\covermove`) was silently
  computed from a pooled deficit/movecount ratio (summed across every
  execution before dividing once) rather than the formal definition's
  per-execution averaging - a quantity closer in spirit to
  `voidmass_subprocess` than to the metric it was meant to be. It's now
  computed by translating classical alignments into skip-alignments'
  own execution-partitioning machinery, verified term-by-term against
  the formal definition. This also batches that computation across all
  of a tree's nodes from one pass over each alignment, to offset the
  extra per-node cost the fix itself introduced.

## [0.4.0] - 2026-09-06

### Added

* Wired up the `toothpaste` stochastic process-tree miner as a discovery
  combo (`lab.discovery.discover_toothpaste`), via a new bridge module
  (`lab/toothpaste_bridge.py`) that shells out to the external tool and
  translates its `.ptree` output through skip-alignments'
  `parse_ptree`/`translate_ppt`.

### Changed

* Extracted the shared skip-probability engine (process tree, alignment,
  execution, probability derivation) into its own library,
  [skip-alignments](https://github.com/adamburkegh/skip-alignments), which
  process-voids now depends on instead of vendoring copies of the same code.
* Switched to a `pyproject.toml`-based install (`pip install -e .`) as the
  documented way to install process-voids and its dependencies.
* `DiscoveryCombo.discover` now returns a `DiscoveryResult` dataclass
  (`tree`, `ppt_weights`) instead of a bare tuple.
* Removed `has_estimator` from the discovery-combo interface; weight
  transfer from a discovered model now works the same way regardless of
  discovery source.
* Introduced  `pm4py_config` dependency to configure pm4py
* `requirements.txt` is now a `pip freeze` record of the exact environment
  used for a build, not a hand-maintained duplicate of the dependency
  ranges in `pyproject.toml`.

### Fixed

* `write_xes` no longer crashes on logs with timezone-aware timestamps
  (normalises to naive UTC before writing, matching what `pm4py.read_xes`
  produces for tz-aware sources).
* Two `ZeroDivisionError`s in skip-alignments' `compute()`/`coninciding_agns()`
  when a fully-degraded (empty) log leaves no variants to average over -
  fixed upstream in skip-alignments 0.2.1; process-voids now surfaces this
  boundary case as a clear `ValueError` from `coverage_by_duration` instead.

### Removed

* Removed the duplicated engine modules (`alignall.py`, `alignment.py`,
  `derivation.py`, `execution.py`, `probabilities.py`, `processtree.py`,
  `skips.py`) from this repo; the same code now lives only in
  `skip-alignments`.

## [0.3.1] - 2026-07-10

### Added

* Added a visualisation tool for process-voids metric. Colour BPMN activities according to their metric values. Can be rendered using `bpmn.io`. 

