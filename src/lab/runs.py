"""
Named, reproducible experiment-1 invocations.

Each entry is an Experiment - a typed configuration (logs, combos,
degradations, levels) rather than an ad hoc dict of kwargs, so the same
object backs both a real run and a --dry-run print of what that run
would compute, and the two can never drift apart from each other.

Run with:

    python -m lab.exp_disco_degrade --run <name>
    python -m lab.exp_disco_degrade --run <name> --dry-run

Add a new run by adding an entry here rather than passing ad hoc flags,
so what was actually run for a given result is committed alongside the
code that produced it.
"""

from dataclasses import dataclass
from pathlib import Path

from lab.params import ALL_COMBOS, ALL_DEGRADATIONS, ALL_LEVELS, ALL_LOGS


@dataclass(frozen=True)
class Experiment:
    """
    A reproducible exp_disco_degrade configuration. out_csv is left
    unset (None) here - exp_disco_degrade.main() resolves and
    timestamps it, since that's specific to how a run was invoked
    (--run name vs ad hoc logs), not to the experiment's own shape.
    """
    name: str
    log_paths: list
    combos: dict
    degradations: dict
    levels: list
    out_csv: str = None

    def cell_count(self):
        """
        Upper bound on (log, combo, dim, level) cells this run
        computes - an overcount when 0.0 is among the levels, since
        exp_disco_degrade.run_disco_degrade shares that one computation
        across every degradation dimension for a given (log, combo)
        rather than repeating it per dimension.
        """
        return len(self.log_paths) * len(self.combos) * len(self.degradations) * len(self.levels)

    def describe(self):
        overcounts_zero = 0.0 in self.levels and len(self.degradations) > 1
        note = (' (overcounts: level 0.0 is shared across degradation '
                 'dims, not repeated per dim)' if overcounts_zero else '')
        return '\n'.join([
            f'Experiment: {self.name}',
            f'  logs:         {[Path(p).stem for p in self.log_paths]}',
            f'  combos:       {list(self.combos)}',
            f'  degradations: {list(self.degradations)}',
            f'  levels:       {self.levels}',
            f'  cells:        {self.cell_count()}{note}',
        ])


RUNS = {
    'smoke': Experiment(
        name='smoke',
        log_paths=[ALL_LOGS['payment_approval']],
        combos={'inductive_noise20': ALL_COMBOS['inductive_noise20']},
        degradations=ALL_DEGRADATIONS,
        levels=[0.0, 0.5],
    ),
    'rtfm': Experiment(
        name='rtfm',
        log_paths=[ALL_LOGS['rtfm']],
        combos=ALL_COMBOS,
        degradations=ALL_DEGRADATIONS,
        levels=ALL_LEVELS,
    ),
    'bpi2013_incidents': Experiment(
        name='bpi2013_incidents',
        log_paths=[ALL_LOGS['bpi2013_incidents']],
        combos=ALL_COMBOS,
        degradations=ALL_DEGRADATIONS,
        levels=ALL_LEVELS,
    ),
    'payment_approval': Experiment(
        name='payment_approval',
        log_paths=[ALL_LOGS['payment_approval']],
        combos=ALL_COMBOS,
        degradations=ALL_DEGRADATIONS,
        levels=ALL_LEVELS,
    ),
    'full': Experiment(
        name='full',
        log_paths=list(ALL_LOGS.values()),
        combos=ALL_COMBOS,
        degradations=ALL_DEGRADATIONS,
        levels=ALL_LEVELS,
    ),
}
