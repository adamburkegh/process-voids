"""
Tests for process_voids.util.registry_diff, the squash-review report on
what a change did to lab.metric_registry.

Builds real temp git repos (same style as test_release_check.py) rather
than mocking git: the thing under test is largely "read this file at
that ref", and the index-vs-HEAD default is exactly the behaviour a
mock would assume rather than check.
"""

import subprocess
import tempfile
import unittest
from pathlib import Path

from process_voids.util.registry_diff import (
    diff_registries, format_report, load_metrics,
)

REGISTRY_PATH = 'src/lab/metric_registry.py'

_REGISTRY_TEMPLATE = '''\
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Metric:
    id: str
    description: str
    source: str
    scripts: tuple
    status: str = 'live'
    superseded_by: str = None
    history: dict = field(default_factory=dict)
    scale: str = None


METRICS = {{
{entries}
}}
'''


def _registry_source(metrics):
    """metrics: list of dicts of Metric kwargs."""
    entries = []
    for m in metrics:
        args = ', '.join(f'{k}={v!r}' for k, v in m.items())
        entries.append(f"    {m['id']!r}: Metric({args}),")
    return _REGISTRY_TEMPLATE.format(entries='\n'.join(entries))


def _metric(id, **overrides):
    base = dict(id=id, description='d', source='s', scripts=('exp',))
    base.update(overrides)
    return base


class _Repo:
    def __init__(self, tmp):
        self.path = Path(tmp)
        self._git('init', '-q')
        self._git('config', 'user.email', 'test@example.com')
        self._git('config', 'user.name', 'Test')

    def _git(self, *args):
        return subprocess.run(['git', *args], cwd=self.path, check=True,
                              capture_output=True, text=True)

    def write_registry(self, metrics):
        path = self.path / REGISTRY_PATH
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_registry_source(metrics), encoding='utf-8')

    def commit_registry(self, metrics, message='x'):
        self.write_registry(metrics)
        self._git('add', REGISTRY_PATH)
        self._git('commit', '-q', '-m', message)

    def stage_registry(self, metrics):
        self.write_registry(metrics)
        self._git('add', REGISTRY_PATH)


class LoadMetricsTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = _Repo(self.tmp.name)

    def test_loads_from_a_named_ref(self):
        self.repo.commit_registry([_metric('a'), _metric('b')])
        loaded = load_metrics('HEAD', cwd=self.repo.path)
        self.assertEqual(sorted(loaded), ['a', 'b'])

    def test_loads_staged_content_from_the_index(self):
        self.repo.commit_registry([_metric('a')])
        self.repo.stage_registry([_metric('a'), _metric('b')])
        self.assertEqual(sorted(load_metrics('INDEX', cwd=self.repo.path)), ['a', 'b'])
        self.assertEqual(sorted(load_metrics('HEAD', cwd=self.repo.path)), ['a'])


def _as_objects(metrics):
    namespace = {}
    exec(compile(_registry_source(metrics), 'registry', 'exec'), namespace)
    return namespace['METRICS']


def _diff(old_metrics, new_metrics):
    return diff_registries(_as_objects(old_metrics), _as_objects(new_metrics))


class DiffRegistriesTest(unittest.TestCase):

    def test_reports_added_and_removed_ids(self):
        result = _diff([_metric('a'), _metric('gone')],
                            [_metric('a'), _metric('new')])
        self.assertEqual(result.added, ['new'])
        self.assertEqual(result.removed, ['gone'])

    def test_reports_field_changes_on_surviving_ids(self):
        result = _diff([_metric('a', status='live', scale='coverage')],
                            [_metric('a', status='retired', scale='coverage',
                                     superseded_by='b')])
        changes = dict((field, (old, new)) for field, old, new in result.changed['a'])
        self.assertEqual(changes['status'], ('live', 'retired'))
        self.assertEqual(changes['superseded_by'], (None, 'b'))
        self.assertNotIn('scale', changes)

    def test_an_extended_history_entry_is_not_a_violation(self):
        result = _diff([_metric('a', history={'v0.5.0': 'was X.'})],
                            [_metric('a', history={'v0.5.0': 'was X. And Y.'})])
        self.assertEqual(result.history_notes['a'], [('v0.5.0', 'extended')])
        self.assertEqual(result.violations, [])

    def test_an_altered_history_entry_is_a_violation(self):
        result = _diff([_metric('a', history={'v0.5.0': 'was X.'})],
                            [_metric('a', history={'v0.5.0': 'was Z.'})])
        self.assertEqual(result.history_notes['a'], [('v0.5.0', 'altered')])
        self.assertIn(('a', 'v0.5.0', 'altered'), result.violations)

    def test_a_dropped_history_entry_is_a_violation(self):
        """The failure this tool was written for: a retirement that
        silently dropped an id's pre-p4 history."""
        result = _diff([_metric('a', history={'v0.5.0': 'was X.'})],
                            [_metric('a', status='retired', history={})])
        self.assertEqual(result.history_notes['a'], [('v0.5.0', 'dropped')])
        self.assertIn(('a', 'v0.5.0', 'dropped'), result.violations)

    def test_a_removed_id_taking_its_history_with_it_is_a_violation(self):
        result = _diff([_metric('a', history={'v0.5.0': 'was X.'})], [])
        self.assertEqual(result.removed, ['a'])
        self.assertIn(('a', 'v0.5.0', 'dropped'), result.violations)

    def test_a_new_history_entry_is_not_a_violation(self):
        result = _diff([_metric('a', history={'v0.5.0': 'was X.'})],
                            [_metric('a', history={'v0.5.0': 'was X.',
                                                   'v0.6.0': 'was Y.'})])
        self.assertEqual(result.history_notes['a'], [('v0.6.0', 'added')])
        self.assertEqual(result.violations, [])

    def test_no_changes_reports_nothing(self):
        result = _diff([_metric('a')], [_metric('a')])
        self.assertFalse(result.added or result.removed or result.changed
                         or result.history_notes or result.violations)
        self.assertIn('no registry changes', format_report(result).lower())


class FormatReportTest(unittest.TestCase):

    def test_names_every_finding(self):
        result = _diff(
            [_metric('gone', history={'v0.5.0': 'was X.'}), _metric('changed')],
            [_metric('changed', status='retired'), _metric('fresh')])
        report = format_report(result)
        for expected in ('gone', 'changed', 'fresh', 'retired', 'dropped'):
            with self.subTest(expected=expected):
                self.assertIn(expected, report)


if __name__ == '__main__':
    unittest.main()
