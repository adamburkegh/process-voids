"""
Tests for process_voids.util.squash_review, the single-command squash
review that combines the staged diff stat, registry_diff, test_name_diff,
a public-repo hygiene scan of added lines, and a CHANGELOG check.

Each section is pinned on both its findings and its clean case. The
parsing and scanning are pure functions tested on synthetic input; the
report as a whole runs against a throwaway git repository, the same way
the sibling review tools' tests do.
"""

import subprocess
import tempfile
import unittest
from pathlib import Path

from process_voids.util.squash_review import (
    ALLOWANCES, FORBIDDEN_PATTERNS, build_report, changelog_note, main,
    parse_added_lines, scan_added_lines,
)

PATTERN_NAMES = {name for name, _regex, _why in FORBIDDEN_PATTERNS}


class _Repo:
    def __init__(self, tmp):
        self.path = Path(tmp)
        self._git('init', '-q')
        self._git('config', 'user.email', 'test@example.com')
        self._git('config', 'user.name', 'Test')

    def _git(self, *args):
        return subprocess.run(['git', *args], cwd=self.path, check=True,
                              capture_output=True, text=True)

    def write(self, rel_path, content):
        path = self.path / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding='utf-8')
        self._git('add', rel_path)

    def commit(self, message='x'):
        self._git('commit', '-q', '-m', message)


def _registry(*ids):
    entries = ''.join(
        f"    {i!r}: Metric(id={i!r}, description='d', source='s', scripts=()),\n"
        for i in ids)
    return ('from dataclasses import dataclass, field\n\n'
            '@dataclass(frozen=True)\n'
            'class Metric:\n'
            '    id: str\n    description: str\n    source: str\n    scripts: tuple\n'
            "    status: str = 'live'\n    superseded_by: str = None\n"
            '    history: dict = field(default_factory=dict)\n'
            '    scale: str = None\n\n'
            f'METRICS = {{\n{entries}}}\n')


def _hits(added):
    return {(path, lineno, name) for path, lineno, name, _text in scan_added_lines(added)}


class ParseAddedLinesTest(unittest.TestCase):

    def test_reads_path_new_side_line_number_and_text(self):
        diff = ('diff --git a/src/m.py b/src/m.py\n'
                '--- a/src/m.py\n'
                '+++ b/src/m.py\n'
                '@@ -3,0 +4,2 @@\n'
                '+first\n'
                '+second\n')
        self.assertEqual(parse_added_lines(diff),
                         [('src/m.py', 4, 'first'), ('src/m.py', 5, 'second')])

    def test_removed_lines_do_not_advance_the_new_side_count(self):
        diff = ('diff --git a/f b/f\n--- a/f\n+++ b/f\n'
                '@@ -10,2 +10,1 @@\n'
                '-gone\n'
                '-also gone\n'
                '+kept\n')
        self.assertEqual(parse_added_lines(diff), [('f', 10, 'kept')])

    def test_spans_several_files_and_hunks(self):
        diff = ('diff --git a/a b/a\n--- a/a\n+++ b/a\n'
                '@@ -1,0 +2,1 @@\n+one\n'
                '@@ -8,0 +20,1 @@\n+two\n'
                'diff --git a/b b/b\nnew file mode 100644\n--- /dev/null\n+++ b/b\n'
                '@@ -0,0 +1,1 @@\n+three\n')
        self.assertEqual(parse_added_lines(diff),
                         [('a', 2, 'one'), ('a', 20, 'two'), ('b', 1, 'three')])

    def test_an_added_line_that_itself_begins_with_plus_signs_is_content(self):
        """Inside a hunk every '+' line is an addition - one whose own
        text starts '++ ' must not be mistaken for a file header."""
        diff = ('diff --git a/f b/f\n--- a/f\n+++ b/f\n'
                '@@ -0,0 +1,1 @@\n'
                '+++ not a header\n')
        self.assertEqual(parse_added_lines(diff), [('f', 1, '++ not a header')])

    def test_ignores_the_no_newline_marker(self):
        diff = ('diff --git a/f b/f\n--- a/f\n+++ b/f\n'
                '@@ -0,0 +1,1 @@\n+x\n\\ No newline at end of file\n')
        self.assertEqual(parse_added_lines(diff), [('f', 1, 'x')])


class ScanAddedLinesTest(unittest.TestCase):

    def test_each_forbidden_pattern_is_found_with_its_file_and_line(self):
        added = [
            ('src/a.py', 1, '# see labnotes.md for the derivation'),
            ('src/a.py', 2, "PAPER = 'var/papers/void.pdf'"),
            ('src/a.py', 3, '# drafted in var/reports/2026-09-12-x.md'),
            ('src/a.py', 4, '# filed in pvoid-lab'),
            ('src/a.py', 5, "CONFIG = '.claude/settings.json'"),
            ('src/a.py', 6, "LOG = 'C:/working/data/rtfm.xes'"),
            ('src/a.py', 7, r"TOOL = 'C:\working\tools\ebi\ebi.exe'"),
            ('src/a.py', 8, '<' * 7 + ' HEAD'),
            ('src/a.py', 9, '>' * 7 + ' main'),
        ]
        self.assertEqual(_hits(added), {
            ('src/a.py', 1, 'labnotes'),
            ('src/a.py', 2, 'var_papers'),
            ('src/a.py', 3, 'var_reports'),
            ('src/a.py', 4, 'pvoid_lab'),
            ('src/a.py', 5, 'claude_dir'),
            ('src/a.py', 6, 'absolute_path'),
            ('src/a.py', 7, 'absolute_path'),
            ('src/a.py', 8, 'conflict_marker'),
            ('src/a.py', 9, 'conflict_marker'),
        })

    def test_every_pattern_is_exercised_above(self):
        """So a pattern added to the constant without a test is noticed."""
        self.assertEqual(PATTERN_NAMES, {'labnotes', 'var_papers', 'var_reports', 'pvoid_lab',
                                         'claude_dir', 'absolute_path', 'conflict_marker'})

    def test_clean_lines_produce_no_hits(self):
        added = [('src/a.py', 1, "LOG = 'data/payment_approval.xes'"),
                 ('src/a.py', 2, '# see lab.metric_registry'),
                 ('src/a.py', 3, "URL = 'https://example.com/path'")]
        self.assertEqual(scan_added_lines(added), [])

    def test_a_url_scheme_is_not_an_absolute_windows_path(self):
        self.assertEqual(scan_added_lines([('f', 1, "'https://github.com/x/y'"),
                                           ('f', 2, "'http://localhost:8080/'")]), [])

    def test_a_marker_not_at_the_start_of_the_line_is_not_a_conflict(self):
        """Real conflict markers start the line; one indented inside a
        string or docstring is text about markers, not a leftover."""
        self.assertEqual(scan_added_lines([('f', 1, "    '" + '<' * 7 + " HEAD'")]), [])

    def test_params_py_is_held_to_the_same_rule_as_every_other_file(self):
        """Machine paths live in pvoid.toml now, so an absolute path
        reappearing in lab/params.py is a regression, not a convention."""
        added = [('src/lab/params.py', 14, "    'rtfm': 'C:/working/data/rtfm.xes',")]
        self.assertEqual(_hits(added), {('src/lab/params.py', 14, 'absolute_path')})

    def test_an_absolute_path_elsewhere_is_flagged(self):
        added = [('src/lab/run.py', 9, "LOG = 'C:/working/data/rtfm.xes'")]
        self.assertEqual(_hits(added), {('src/lab/run.py', 9, 'absolute_path')})

    def test_only_the_tool_and_its_tests_are_allowed_anything(self):
        """They have to spell out each pattern; nothing else has a reason to."""
        self.assertEqual(set(ALLOWANCES), {'src/process_voids/util/squash_review.py',
                                           'tests/process_voids/test_squash_review.py'})

    def test_allowances_name_only_real_patterns(self):
        for path, allowed in ALLOWANCES.items():
            with self.subTest(path=path):
                self.assertTrue(allowed <= PATTERN_NAMES)


class ChangelogNoteTest(unittest.TestCase):

    def test_src_changed_without_a_changelog_entry_is_reported(self):
        note = changelog_note(['src/lab/run.py', 'tests/lab/test_run.py'])
        self.assertIsNotNone(note)
        self.assertIn('src/lab/run.py', note)

    def test_src_changed_with_a_changelog_entry_says_nothing(self):
        self.assertIsNone(changelog_note(['src/lab/run.py', 'CHANGELOG.md']))

    def test_a_change_outside_src_needs_no_entry(self):
        self.assertIsNone(changelog_note(['tests/lab/test_run.py', 'AGENTS.md']))

    def test_the_note_reports_rather_than_judges(self):
        """A docstring-only change legitimately needs no entry, so the
        note must not read as a failure."""
        note = changelog_note(['src/lab/run.py'])
        self.assertNotIn('FAIL', note.upper())
        self.assertIn('docstring', note)


class BuildReportTest(unittest.TestCase):
    """The whole report against a real staged change."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = _Repo(self.tmp.name)
        self.repo.write('src/lab/metric_registry.py', _registry('a'))
        self.repo.write('tests/test_x.py',
                        'def test_kept():\n    pass\n\ndef test_dropped():\n    pass\n')
        self.repo.write('CHANGELOG.md', '# Changelog\n')
        self.repo.commit()

    def _report(self):
        return build_report('HEAD', 'INDEX', cwd=self.repo.path)

    def test_every_section_is_present_even_when_clean(self):
        self.repo.write('README.md', 'hello\n')
        report = self._report()
        for heading in ('Diff stat', 'Registry', 'Test names',
                        'Public-repo hygiene', 'CHANGELOG'):
            with self.subTest(heading=heading):
                self.assertIn(heading, report)

    def test_a_clean_change_reports_no_findings_in_each_section(self):
        self.repo.write('README.md', 'hello\n')
        report = self._report()
        self.assertIn('No registry changes.', report)
        self.assertIn('No test names dropped.', report)
        self.assertIn('No forbidden references in added lines.', report)
        self.assertIn('No src/ change without a CHANGELOG entry.', report)

    def test_the_diff_stat_names_the_staged_files(self):
        self.repo.write('README.md', 'hello\n')
        self.assertIn('README.md', self._report())

    def test_registry_changes_come_through(self):
        self.repo.write('src/lab/metric_registry.py', _registry('a', 'b'))
        self.repo.write('CHANGELOG.md', '# Changelog\n\n* added b\n')
        self.assertIn('+ b', self._report())

    def test_a_dropped_test_name_comes_through(self):
        self.repo.write('tests/test_x.py', 'def test_kept():\n    pass\n')
        self.assertIn('test_dropped', self._report())

    def test_a_forbidden_reference_is_reported_with_file_and_line(self):
        self.repo.write('src/lab/new.py', 'X = 1\n# see labnotes.md\n')
        self.repo.write('CHANGELOG.md', '# Changelog\n\n* new\n')
        report = self._report()
        self.assertIn('src/lab/new.py:2', report)
        self.assertIn('labnotes', report)

    def test_src_without_changelog_is_reported(self):
        self.repo.write('src/lab/new.py', 'X = 1\n')
        self.assertIn('src/lab/new.py', self._report().split('CHANGELOG')[-1])

    def test_explicit_refs_compare_two_commits(self):
        self.repo.write('src/lab/new.py', 'X = 1  # see labnotes.md\n')
        self.repo.commit()
        report = build_report('HEAD~1', 'HEAD', cwd=self.repo.path)
        self.assertIn('src/lab/new.py:1', report)

    def test_no_merge_base_note_where_the_base_has_not_moved(self):
        """The base is already the merge base here, so the note would
        only suggest a difference that is not there."""
        self.repo.write('src/lab/new.py', 'X = 1\n')
        self.repo.commit()
        report = build_report('HEAD~1', 'HEAD', cwd=self.repo.path)
        self.assertNotIn('merge base', report.splitlines()[0])


class BranchBehindTrunkTest(unittest.TestCase):
    """A branch cut before trunk moved on. Reviewing it against trunk
    must report the branch's own changes - not trunk's newer work, which
    a plain two-tree comparison shows as the branch removing it."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = _Repo(self.tmp.name)
        git = self.repo._git

        self.repo.write('src/lab/metric_registry.py', _registry('a'))
        self.repo.write('tests/test_x.py', 'def test_kept():\n    pass\n')
        self.repo.write('CHANGELOG.md', '# Changelog\n')
        self.repo.commit('base')
        git('branch', '-M', 'trunk')

        git('checkout', '-q', '-b', 'feature')
        self.repo.write('src/lab/feature.py', 'X = 1\n')
        self.repo.write('CHANGELOG.md', '# Changelog\n\n* feature\n')
        self.repo.commit('feature work')

        git('checkout', '-q', 'trunk')
        self.repo.write('src/lab/metric_registry.py', _registry('a', 'b'))
        self.repo.write('tests/test_x.py',
                        'def test_kept():\n    pass\n\ndef test_on_trunk():\n    pass\n')
        self.repo.write('CHANGELOG.md', '# Changelog\n\n* trunk\n')
        self.repo.commit('trunk moves on')

        self.report = build_report('trunk', 'feature', cwd=self.repo.path)

    def test_trunks_newer_registry_id_is_not_reported_as_removed(self):
        self.assertNotIn('- b', self.report)
        self.assertIn('No registry changes.', self.report)

    def test_trunks_newer_test_is_not_reported_as_dropped(self):
        self.assertNotIn('test_on_trunk', self.report)
        self.assertIn('No test names dropped.', self.report)

    def test_the_diff_stat_holds_only_the_branchs_own_files(self):
        stat = self.report.split('== Diff stat ==')[1].split('==')[0]
        self.assertIn('src/lab/feature.py', stat)
        self.assertNotIn('metric_registry.py', stat)
        self.assertNotIn('test_x.py', stat)

    def test_the_header_names_the_merge_base_it_compared_from(self):
        """So a reader can see the comparison was not against trunk's tip."""
        merge_base = self.repo._git('merge-base', 'trunk', 'feature').stdout.strip()
        self.assertIn(merge_base[:7], self.report.splitlines()[0])


class MainTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = _Repo(self.tmp.name)
        self.repo.write('src/lab/metric_registry.py', _registry('a'))
        self.repo.write('CHANGELOG.md', '# Changelog\n')
        self.repo.commit()

    def test_findings_do_not_make_it_exit_non_zero(self):
        """It reports; release_check is the gate."""
        self.repo.write('src/lab/new.py', '# see labnotes.md\n')
        self.assertEqual(main(['HEAD', 'INDEX'], cwd=self.repo.path), 0)

    def test_a_bad_ref_is_an_internal_error(self):
        self.assertEqual(main(['no-such-ref', 'INDEX'], cwd=self.repo.path), 1)


if __name__ == '__main__':
    unittest.main()
