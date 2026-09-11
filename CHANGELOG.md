# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

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

### Changed

* Tests call skip-alignments' `Aligner.align_normal_form` instead of the
  deprecated `align2`, clearing its deprecation warnings from the suite.

### Fixed

* `voidsat` was very slow on high case-count logs: `admass` rebuilt the
  log's trace list (a full group-by, sort and conversion over the log)
  once per tree node instead of once per report row - about 4,270s
  extra per cell on rtfm. The trace list is now built once, with the
  per-row cache, and reused for every node.

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

