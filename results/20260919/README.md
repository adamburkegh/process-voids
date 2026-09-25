# Results snapshot: 2026-09-19

Published result data and plots for the paper's headline logs
(`bpi2013_closed_problems`, `bpic2020_rfp`, `rtfm`) and the synthetic
fixtures (`claims`, `partial_sequence`, `payment_approval`,
`payment_partial`), built by [`lab.publish_results`](../../src/lab/publish_results.py)
from full-size sweep output.

## Contents

- `csv/` - one root-level result CSV and one aggregated timings CSV per
  log, plus `runs.csv`: which source file each log came from and the
  versions, seed and configuration `lab.run_history` recorded for it.
  Per-node data is not published.
- `plots/` - dose-response plots from [`lab.cross_log_plots`](../../src/lab/cross_log_plots.py),
  reading directly from `csv/`. Ten PNGs: `inductive_noise20` and
  `toothpaste_noise10`, each under `trace` and `activity_frequency_gradual`
  degradation, in three panel arrangements (panel per metric, panel per
  log, panel per log/metric pair with ref model as the line). Regenerate
  with:

  ```
  python -m lab.cross_log_plots \
      --log bpi2013_closed_problems=csv/bpi2013_closed_problems.csv \
      --log rtfm=csv/rtfm.csv \
      --log bpic2020_rfp=csv/bpic2020_rfp.csv \
      --combos inductive_noise20 toothpaste_noise10 \
      --degradation-dims trace activity_frequency_gradual \
      --ylim 0 1 --panel-by metric --out-dir plots
  ```

  (`--panel-by log` and `--panel-by log_metric` produce the other two
  arrangements; `log_metric` additionally takes plain `toothpaste` in
  `--combos`, since comparing ref models is its point.)

  The PNGs were generated on 2026-09-20, one day after the `csv/` data
  they read - the folder is dated for when the data was published and
  submitted against, not for when its plots happened to be drawn.

## Provenance

The `csv/` files were generated before the classical per-node attribution
fix ("Credit classical alignment moves to the leaf that fired", 0.5.2's
`CHANGELOG.md`), so they carry the fix's retired column names
(`voidmass_process_lower`/`_upper` rather than `voidmass_process2_*`).
`lab.cross_log_plots` and `lab.trace_variability` read either name, so the
plots above are unaffected regardless.

This has been checked, not assumed: on these same three logs, the fix
moves the classical per-node root value by exactly 0 in every degraded
cell of `rtfm` and `bpic2020_rfp` (all combos) and of
`bpi2013_closed_problems` `inductive_noise20`, and by at most 0.0025 on
`bpi2013_closed_problems` `toothpaste_noise10` (at most 0.0002 on plain
`toothpaste`). No `inductive_noise20` tree used in these results repeats
an activity label, which is what the fix changes the behaviour of.

This folder is a snapshot, not a rolling update - a later publish adds a
new dated folder alongside this one rather than overwriting it.
