"""
Degradation sweep on the claims fixture (lab.claims_fixture), against
its known GENERATING tree - unlike exp_disco_degrade, which references
a discovered model, ground truth is known exactly here (see
claims-fixture-spec.md), so this checks metrics against a known-correct
model rather than whatever a discovery algorithm produced.

Computes every metric currently available side by side at each level:
the three in lab.metrics.compute_metrics (weight_coverage,
skipprob, salign_coverage - all skip-alignments-based; duration_coverage
dropped from the roster, see that module's docstring) plus deficit/
voidmass_subprocess/voidmass_process from the classical-alignment path
in process_voids.voidmass_pn (see that module's docstring - still a
prototype, not yet folded into compute_metrics/coveragemass.py, but the
validated direction for voidmass specifically - see session notes on
the E1 lumping fix and the tied-alignment dedup fix).

Targets (see claims-fixture-spec.md's "Pre-identified ablation
targets"): the loop block (mandatory, large, wide receive_docs
duration), assess (mandatory, small), appeal_seq (OPTIONAL - the
specificity negative control: ablating it should show no void, since
its absence from a case is normal model behaviour, not missing data).

degrade_target_subprocess (lab.degradation) removes only the target
subprocess's own events from n chosen cases, not the whole case -
other activities in those cases stay intact.

Usage:
    python -m lab.exp_claims_degrade --target assess loop_block appeal_seq \\
        --n-drops 0 2 4 6 8 [--out ...]
"""

import argparse
import logging
import time
from pathlib import Path

import pandas as pd
import pm4py_config as pm4py

from lab.claims_fixture import CLAIMS_GROUND_TRUTH_CSV, CLAIMS_XES, build_claims_tree
from lab.degradation import degrade_target_subprocess
from lab.logconfig import configure
from lab.metrics import compute_metrics
from lab.exp_disco_degrade import CLASSICAL_METRIC_KEYS
from process_voids.voidmass_pn import build_id_net, voidmass_table_pn, coverage_by_alignment_pn

logger = logging.getLogger(__name__)

CELL_COLS = ['target', 'n_drop_cases']
SLPN_DIR = Path('var/lab/claims_degrade')

TARGETS = {
    'assess': {'assess'},
    'loop_block': {'request_docs', 'receive_docs'},
    'appeal_seq': {'lodge_appeal', 'decide_appeal'},
}


def _merge_write(df, path, cell_cols=CELL_COLS):
    """Upsert by cell - see lab.exp_surprise._merge_write for the full rationale."""
    path = Path(path)
    if df.empty:
        return pd.read_csv(path) if path.exists() else df
    cells = set(df[cell_cols].fillna('').apply(tuple, axis=1))
    if path.exists():
        existing = pd.read_csv(path)
        existing_cells = existing[cell_cols].fillna('').apply(tuple, axis=1)
        existing = existing[~existing_cells.isin(cells)]
        combined = pd.concat([existing, df], ignore_index=True)
    else:
        combined = df
    path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(path, index=False)
    return combined


def _variant_probs(log):
    """{trace variant (activity tuple): probability} from a log's case frequencies."""
    n_cases = log['case:concept:name'].nunique()
    variants = {}
    for _case, group in log.groupby('case:concept:name', sort=False):
        variant = tuple(group.sort_values('time:timestamp')['concept:name'])
        variants[variant] = variants.get(variant, 0) + 1
    return {v: c / n_cases for v, c in variants.items()}


def _target_node(tree, target_activities):
    """
    The (topmost, i.e. coarsest) tree node whose leaf-label set exactly
    equals target_activities - eg {'lodge_appeal', 'decide_appeal'}
    finds appeal_seq itself, not its Xor parent (which also has the
    Tau leaf 'skip-appeal') or either leaf individually.
    """
    target_activities = set(target_activities)

    def _walk(node):
        if set(node.get_leaf_labels()) == target_activities:
            return node
        for child in node.children:
            found = _walk(child)
            if found is not None:
                return found
        return None

    node = _walk(tree)
    if node is None:
        raise ValueError(f'no tree node with leaf set {target_activities}')
    return node


def _deviated_cases(ground_truth_csv):
    """
    Case ids the fixture deliberately deviated from the model (lab.
    claims_fixture's 'deviation' column) - excluded from ablation
    eligibility in every target below. Otherwise an ablation can
    accidentally sweep away the log's only genuine deviation as a side
    effect of dropping cases, confounding "voidmass correctly ignores
    an optional subprocess" with "the deviation happened to get
    deleted" - see session notes, the original appeal_seq run.
    """
    ground_truth = pd.read_csv(ground_truth_csv)
    deviated = ground_truth.loc[ground_truth['deviation'].notna()
                                 & (ground_truth['deviation'] != ''), 'case']
    return set(deviated)


def run_claims_degradation(target_names=None, n_drops=None,
                            out_csv='var/lab/results/exp_claims_degrade.csv',
                            ground_truth_csv=CLAIMS_GROUND_TRUTH_CSV):
    tree = build_claims_tree()
    base_log = pm4py.read_xes(CLAIMS_XES)
    net, im, fm, activity_to_id, tau_ids, id_loop_list = build_id_net(tree)
    exclude_cases = _deviated_cases(ground_truth_csv)

    target_names = target_names or list(TARGETS)
    logger.info('Experiment: claims degrade | targets=%s | n_drops=%s | excluding=%s',
                target_names, n_drops, exclude_cases)

    rows = []
    for target_name in target_names:
        target_activities = TARGETS[target_name]
        target_node = _target_node(tree, target_activities)

        eligible_cases = set(base_log.loc[
            base_log['concept:name'].isin(target_activities), 'case:concept:name'
        ].unique()) - exclude_cases
        max_drops = len(eligible_cases)
        levels = n_drops if n_drops is not None else list(range(max_drops + 1))

        for n in levels:
            cell = f'{target_name} / n_drop_cases={n}'
            started = time.monotonic()
            try:
                if n == 0:
                    degraded_log, dropped = base_log, set()
                else:
                    degraded_log, dropped = degrade_target_subprocess(
                        base_log, target_activities, n, exclude_cases=exclude_cases)

                slpn_path = SLPN_DIR / f'claims_{target_name}_n{n}.slpn'
                metrics, dv = compute_metrics(degraded_log, tree, str(slpn_path), return_dv=True)

                variant_probs = _variant_probs(degraded_log)
                vm_table = voidmass_table_pn(tree, variant_probs, net, im, fm,
                                              activity_to_id, tau_ids,
                                              id_loop_list=id_loop_list, timeout=60)
                target_row = vm_table[target_node]
                classical_values = (
                    target_row['deficit'],
                    target_row['movecount'],
                    target_row['voidmass_subprocess'],
                    target_row['voidmass_process'],
                    coverage_by_alignment_pn(target_node, dv.skip_probs[target_node], vm_table),
                )

                rows.append({
                    'target': target_name, 'n_drop_cases': n,
                    'dropped_case_count': len(dropped), 'status': 'ok',
                    **metrics,
                    **dict(zip(CLASSICAL_METRIC_KEYS, classical_values)),
                    'elapsed_s': time.monotonic() - started,
                })
                logger.info('%s - done in %.1fs', cell, time.monotonic() - started)
            except Exception as e:
                rows.append({'target': target_name, 'n_drop_cases': n,
                              'status': f'error: {e}'})
                logger.warning('%s - error: %s', cell, e)

    return _merge_write(pd.DataFrame(rows), out_csv)


def main():
    configure()
    logger.info('Starting exp_claims_degrade')
    parser = argparse.ArgumentParser(
        description='Degradation sweep on the claims fixture against its known '
                     'generating tree, across every available metric.')
    parser.add_argument('--target', nargs='+', choices=list(TARGETS), default=None,
                         help=f'Ablation targets (default: all of {list(TARGETS)})')
    parser.add_argument('--n-drops', type=int, nargs='+', default=None,
                         help='Explicit case-drop counts (default: every level from '
                              '0 to the number of eligible cases, per target)')
    parser.add_argument('--out', default='var/lab/results/exp_claims_degrade.csv')
    args = parser.parse_args()

    df = run_claims_degradation(target_names=args.target, n_drops=args.n_drops,
                                 out_csv=args.out)
    print(df.to_string(index=False))
    logger.info('Wrote %s', args.out)


if __name__ == '__main__':
    main()
