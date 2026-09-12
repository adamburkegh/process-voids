'''
The extremes criterion (docs/DESIGN.md, section 4) for every live metric
that declares a scale in lab.metric_registry: each is scored at submodel
b of seq(a, b), against logs in which b is never missing, always
missing, and missing from half the traces. The registry's module
docstring says what each scale promises in each case.

A metric that breaks its promise today is pinned in KNOWN_FAILURES at
the value it actually reads. The suite stays green while the failure
stays visible, and turns red when that value changes, so the list can't
quietly go stale.

Values come from the producing scripts' own declarations rather than
being re-derived here: exp_disco_degrade's ALL_METRICS scored through a
CellContext, and exp_voidmass's node and summary rows (voidmass_table
over that same context's dv stage, with b as the target). Skip
probabilities come from the real ebi-backed pipeline.
'''

import unittest

from skipalignments import Activity, Sequence

from lab.exp_disco_degrade import ALL_METRICS, CLASSICAL_ALIGNMENT_TIMEOUT
from lab.metric_registry import METRICS
from process_voids import dtlog
from process_voids.coveragemass import voidmass_table
from process_voids.metric_context import CellContext, METRIC_ERROR, score_all
from process_voids.voidmass_pn import build_id_net

ACT_COST = 100000
TOLERANCE = 1e-6

# Traces in dtlog.convert_timed's activity:hour notation.
CASES = {
    'nothing_missing': ('a:0 b:1', 'a:0 b:1'),
    'always_missing': ('a:0', 'a:0'),
    'half_missing': ('a:0 b:1', 'a:0'),
}

# (metric id, case) -> the value it reads today, where that breaks its promise.
KNOWN_FAILURES = {
    # (1 - skip_prob) and the match-ratio mass both count b's absence, so
    # coverage falls as (1 - p)^2 rather than 1 - p.
    ('salign_coverage', 'half_missing'): 0.25,
    # skip_prob * voidmass: the voidmass has already counted the absence.
    ('voidage_subprocess', 'half_missing'): 0.25,
    ('target_voidage_subprocess', 'half_missing'): 0.25,
    ('voidage_process', 'half_missing'): 0.125,
    ('target_voidage_process', 'half_missing'): 0.125,
    # b is the last activity: a skip move with no following event gets no
    # duration, so an always-missing b has no aligned duration at all.
    ('voidsat', 'always_missing'): 0.0,
    ('voidsat', 'half_missing'): 0.5,
}

_READINGS = {
    'coverage': {'nothing_missing': 1.0, 'always_missing': 0.0, 'half_missing': 0.5},
    'void': {'nothing_missing': 0.0, 'always_missing': 1.0, 'half_missing': 0.5},
}


def _scaled_live_ids():
    return {mid for mid, m in METRICS.items() if m.status == 'live' and m.scale is not None}


def _promise(scale, case, readings):
    '''
    (kept, promised): whether a metric on `scale` keeps its promise for
    `case`, given its readings in every case, and what it promised.
    '''
    value = readings[case]
    if scale == 'void_share':
        share = readings['always_missing']
        if case == 'always_missing':
            return value > TOLERANCE, 'above 0'
        promised = 0.0 if case == 'nothing_missing' else share / 2
    else:
        promised = _READINGS[scale][case]
    return abs(value - promised) < TOLERANCE, f'{promised:g}'


def _seq_ab():
    a = Activity(None, 'a', ACT_COST)
    a.id = '1'
    b = Activity(None, 'b', ACT_COST)
    b.id = '2'
    tree = Sequence(None, [a, b])
    tree.id = '3'
    a.set_parent(tree)
    b.set_parent(tree)
    return tree, b


def _score(case, traces):
    '''{script: {metric id: value at b}} for one case.'''
    tree, b = _seq_ab()
    log = dtlog.convert_timed(*traces, names=[f'c{i}' for i in range(len(traces))])
    ctx = CellContext(log=log, tree=tree,
                      slpn_path=f'var/lab/test_metric_extremes_{case}.slpn',
                      classical_net=build_id_net(tree),
                      classical_timeout=CLASSICAL_ALIGNMENT_TIMEOUT)
    disco = score_all(ctx, ALL_METRICS, node=b)
    errored = sorted(mid for mid, value in disco.items() if value is METRIC_ERROR)
    if errored:
        raise AssertionError(f'{case}: these metrics raised while scoring: {errored}')

    # exp_voidmass's node row, and its summary row for target b
    dv = ctx.stage('dv')
    table = voidmass_table(tree, dv.skip_dict_backup, dv.pl, skip_probs=dv.skip_probs)
    voidmass = {'skip_prob': dv.skip_probs[b], **table[b]}
    voidmass.update({f'target_{key}': table[b][key]
                     for key in ('voidmass_subprocess', 'voidmass_process',
                                 'voidage_subprocess', 'voidage_process')})
    return {'exp_disco_degrade': disco, 'exp_voidmass': voidmass}


class ExtremesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.values = {case: _score(case, traces) for case, traces in CASES.items()}

    def test_every_scaled_live_metric_is_scored(self):
        for metric_id in sorted(_scaled_live_ids()):
            for script in METRICS[metric_id].scripts:
                with self.subTest(metric=metric_id, script=script):
                    self.assertIn(metric_id, self.values['half_missing'][script],
                                  'no reading for this metric - add it to _score')

    def test_scaled_live_metrics_keep_their_promise(self):
        for metric_id in sorted(_scaled_live_ids()):
            metric = METRICS[metric_id]
            for script in metric.scripts:
                if metric_id not in self.values['half_missing'][script]:
                    continue  # reported by test_every_scaled_live_metric_is_scored
                readings = {case: self.values[case][script][metric_id] for case in CASES}
                for case in CASES:
                    kept, promised = _promise(metric.scale, case, readings)
                    with self.subTest(metric=metric_id, case=case):
                        if (metric_id, case) in KNOWN_FAILURES:
                            self.assertFalse(kept, f'now reads {promised} as promised - '
                                                   'remove it from KNOWN_FAILURES')
                            self.assertAlmostEqual(readings[case],
                                                   KNOWN_FAILURES[metric_id, case], places=6)
                        else:
                            self.assertTrue(kept, f'reads {readings[case]:g}, '
                                                  f'promised {promised}')


class KnownFailuresTest(unittest.TestCase):

    def test_known_failures_name_scaled_live_metrics_and_cases(self):
        for metric_id, case in KNOWN_FAILURES:
            with self.subTest(metric=metric_id, case=case):
                self.assertIn(metric_id, _scaled_live_ids())
                self.assertIn(case, CASES)

    def test_known_failures_say_so_in_the_registry(self):
        # docs/DESIGN.md section 4: a metric that fails a criterion says so
        # in its registry description.
        for metric_id in sorted({metric_id for metric_id, _case in KNOWN_FAILURES}):
            with self.subTest(metric=metric_id):
                self.assertIn('extremes criterion', METRICS[metric_id].description)


if __name__ == '__main__':
    unittest.main()
