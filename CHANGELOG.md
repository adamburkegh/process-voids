# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]


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

### Fixed

* `write_xes` no longer crashes on logs with timezone-aware timestamps
  (normalises to naive UTC before writing, matching what `pm4py.read_xes`
  produces for tz-aware sources).

### Removed

* Removed the duplicated engine modules (`alignall.py`, `alignment.py`,
  `derivation.py`, `execution.py`, `probabilities.py`, `processtree.py`,
  `skips.py`) from this repo; the same code now lives only in
  `skip-alignments`.

## [0.3.1] - 2026-07-10

### Added

* Added a visualisation tool for process-voids metric. Colour BPMN activities according to their metric values. Can be rendered using `bpmn.io`. 

