"""
Named, reproducible experiment-1 invocations.

Each entry selects a subset of the full parameter catalog in
lab/params.py (or overrides it outright) and is run with:

    python -m lab.exp_disco_degrade --run <name>

Add a new run by adding an entry here rather than passing ad hoc flags,
so what was actually run for a given result is committed alongside the
code that produced it.
"""

from lab.params import ALL_COMBOS, ALL_DEGRADATIONS, ALL_LEVELS, ALL_LOGS

RUNS = {
    'smoke': dict(
        log_paths=[ALL_LOGS['rtfm']],
        combos={'inductive': ALL_COMBOS['inductive']},
        degradations=ALL_DEGRADATIONS,
        levels=[0.0, 0.5],
    ),
    'rtfm': dict(
        log_paths=[ALL_LOGS['rtfm']],
        combos=ALL_COMBOS,
        degradations=ALL_DEGRADATIONS,
        levels=ALL_LEVELS,
    ),
    'bpi2013_incidents': dict(
        log_paths=[ALL_LOGS['bpi2013_incidents']],
        combos=ALL_COMBOS,
        degradations=ALL_DEGRADATIONS,
        levels=ALL_LEVELS,
    ),
    'payment_approval': dict(
        log_paths=[ALL_LOGS['payment_approval']],
        combos=ALL_COMBOS,
        degradations=ALL_DEGRADATIONS,
        levels=ALL_LEVELS,
    ),
    'full': dict(
        log_paths=list(ALL_LOGS.values()),
        combos=ALL_COMBOS,
        degradations=ALL_DEGRADATIONS,
        levels=ALL_LEVELS,
    ),
}
