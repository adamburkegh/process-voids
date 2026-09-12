"""
One-off probe, not a permanent metric: is skip-alignments' lumped mass-
term computation (align_sk_all) actually cheaper than the classical
(non-lumped) one (align_variant_all, via voidmass_table_pn's own
'classical' stage), as the lumping theory predicts?

Earlier measurements only supported a DEDUCTION, not a direct comparison:
the classical stage's cost was measured repeatedly and cleanly (it's
already an isolated call), but the skip-alignment mass term was only ever
isolated once, via parsing DerivationPipeline.compute()'s own stdout
progress markers - a single uncontrolled data point, at a different code
version than whatever it was being compared against. This probe removes
the stdout-parsing entirely: align_sk_all is directly importable, same
tier of public API as align_pn_all (which voidmass_pn.py already calls
directly for the classical side) - see its own docstring, importable
from skipalignments.alignall.

Design: reuses the real harness's own CellContext/ProcessMetric/
TimingListener (process_voids.metric_context, lab.timing) rather than a
parallel, ad hoc timing mechanism - the mass-term stage this adds
(_dv_mass_term_stage) is exactly as legitimate a stage as the production
'classical' one, just not registered in process_voids.metric_context's
own STAGES/ALL_METRICS, since this is a one-off investigative question,
not a metric any real run needs - registry hygiene, not hack-quarantine
(see DiagnosticCellContext).

A fresh CellContext per repetition, not one reused across reps - .stage()/
.score() memoise per-context, so a second call on the SAME context would
just return the first measurement, not a new one. Safe to reuse the same
(log, tree, classical_net) objects across repetitions: neither
align_sk_all nor align_variant_all writes back into the tree (unlike the
production 'dv' stage's transfer_pt_weights, which is why THAT one needs
a fresh context per real cell) - see process_voids.metric_context's own
module docstring for that distinction.

Usage:
    python -m lab.mass_term_probe
    python -m lab.mass_term_probe --reps 10
"""

import argparse
import time

import pandas as pd
import pm4py_config as pm4py

from lab.discovery import COMBOS
from lab.exp_disco_degrade import CLASSICAL_ALIGNMENT_TIMEOUT
from lab.params import ALL_LOGS
from lab.timing import TimingListener
from process_voids.coveragemass import total_node_count
from process_voids.metric_context import CellContext, ProcessMetric, STAGES as _PRODUCTION_STAGES
from process_voids.voidmass_pn import build_id_net
from skipalignments.alignall import align_sk_all

# (log key into lab.params.ALL_LOGS, combo key into lab.discovery.COMBOS) -
# spans a range of variant counts and tree shapes deliberately, not just
# one pairing: a single comparison can't tell us whether the relationship
# holds generally or is specific to one model's shape.
DEFAULT_PAIRS = [
    ('payment_approval', 'inductive_noise20'),
    ('rtfm', 'inductive_noise20'),
    ('rtfm', 'toothpaste'),
]


def _distinct_variant_lists(log):
    """Distinct activity-sequence variants in `log`, as align_sk_all's
    own List[List[str]] shape (see its docstring). Duplicated from
    process_voids.metric_context._variant_probs rather than importing it
    (that name is underscore-private, and this package's own convention
    is to duplicate a small helper rather than reach into another
    module's private one - see that function's docstring for the prior
    example). Only the distinct variants matter here, not their
    probabilities - align_sk_all aligns each one once regardless of how
    many cases share it."""
    variants = set()
    for _case, group in log.groupby('case:concept:name', sort=False):
        variants.add(tuple(group.sort_values('time:timestamp')['concept:name']))
    return [list(v) for v in variants]


def _dv_mass_term_stage(ctx):
    """The mass-term computation alone: align_sk_all, the direct
    analogue of the production 'classical' stage's align_variant_all -
    skip-alignments' own lumped optimal-skip-alignment search, with no
    skip_prob derivation mixed in (unlike the real 'dv' stage, which
    wraps this same search inside DerivationPipeline.compute() alongside
    the probability-derivation work this probe deliberately excludes).

    align_sk_all returns FUTURES (it uses multiprocessing internally -
    see its own docstring) - .result() is what actually blocks until the
    computation finishes, so it's what must be timed, not the submit
    call, which returns almost immediately regardless of how much work
    it queued.
    """
    variant_strings = _distinct_variant_lists(ctx.log)
    futures = align_sk_all(variant_strings, ctx.tree, timeout=CLASSICAL_ALIGNMENT_TIMEOUT)
    return [f.result() for f in futures]


DV_MASS_TERM_METRIC = ProcessMetric(
    id='dv_mass_term', scope='root', needs=('dv_mass_term',),
    compute=lambda ctx, node: ctx.stage('dv_mass_term'))


class DiagnosticCellContext(CellContext):
    """CellContext with one extra opt-in stage ('dv_mass_term') beyond
    the production STAGES - every real stage still behaves identically,
    since this only ADDS an id nothing else references. Kept as a
    subclass here (lab-only) rather than added to process_voids.
    metric_context.STAGES directly, since it's not a metric any real
    run needs - see module docstring."""
    STAGES = {**_PRODUCTION_STAGES, 'dv_mass_term': _dv_mass_term_stage}


def _cell(log_path, combo_name):
    """(log, tree, classical_net) for one (log, combo) pair - the
    reusable setup shared across every repetition, computed once."""
    log = pm4py.read_xes(log_path)
    tree = COMBOS[combo_name].discover(log).tree
    classical_net = build_id_net(tree)
    return log, tree, classical_net


def run_probe(pairs=DEFAULT_PAIRS, reps=5):
    """
    Runs both methods `reps` times per (log, combo) pair, ALTERNATING
    order each repetition (classical-then-mass-term, then mass-term-
    then-classical, ...) so a systematic first-call/warm-cache bias
    can't masquerade as a real difference between the two methods.

    Returns a DataFrame in the SAME long-form shape as lab.timing.
    TimingListener's own rows (metric_or_stage/seconds/status) plus this
    probe's own (log, combo, rep, n_variants, n_nodes) columns - directly
    analysable with the same pandas idioms already used on every other
    timings_df this project produces, not a bespoke format.
    """
    rows = []
    for log_key, combo_name in pairs:
        log_path = ALL_LOGS[log_key]
        log, tree, classical_net = _cell(log_path, combo_name)
        n_variants = len(_distinct_variant_lists(log))
        n_nodes = total_node_count(tree)

        for rep in range(reps):
            methods = (['classical', 'mass_term'] if rep % 2 == 0
                       else ['mass_term', 'classical'])
            for method in methods:
                listener = TimingListener()
                ctx = DiagnosticCellContext(log=log, tree=tree, listeners=[listener],
                                            classical_net=classical_net,
                                            classical_timeout=CLASSICAL_ALIGNMENT_TIMEOUT)
                started = time.monotonic()
                if method == 'classical':
                    ctx.stage('classical')
                else:
                    ctx.score(DV_MASS_TERM_METRIC)
                elapsed = time.monotonic() - started
                rows.append({
                    'log': log_key, 'combo': combo_name, 'rep': rep,
                    'n_variants': n_variants, 'n_nodes': n_nodes,
                    'metric_or_stage': method, 'seconds': elapsed,
                    'status': 'ok',
                })
                print(f'{log_key} / {combo_name} / rep {rep} / {method}: {elapsed:.3f}s')

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description='Probe: is the skip-alignment mass term (align_sk_all) actually '
                     'faster than the classical one (align_variant_all)? See module '
                     'docstring - a direct, controlled, repeated comparison, not a '
                     'deduction from separately-measured component runtimes.')
    parser.add_argument('--reps', type=int, default=5,
                         help='Repetitions per (log, combo) pair (default: 5)')
    parser.add_argument('--out', default='var/lab/results/mass_term_probe.csv')
    args = parser.parse_args()

    df = run_probe(reps=args.reps)
    df.to_csv(args.out, index=False)
    print(f'\nWrote {args.out} ({len(df)} rows)\n')
    summary = df.groupby(['log', 'combo', 'metric_or_stage']).seconds.agg(
        ['count', 'mean', 'std', 'min', 'max'])
    print(summary.to_string())


if __name__ == '__main__':
    main()
