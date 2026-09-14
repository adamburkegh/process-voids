"""
Tests for process_voids.util.branch_survey.

The point of the tool is that "N commits ahead" is meaningless in a
squash workflow: a branch can be ahead and content-identical to main,
and re-squashing such a branch re-adds what already landed. So these
build real repos where a branch is ahead but landed, ahead and genuinely
different, and conflicting, and pin what the survey says about each.
"""

import subprocess
import tempfile
import unittest
from pathlib import Path

from process_voids.util.branch_survey import (
    format_report, squash_state, survey, worktrees_by_branch,
)


class _Repo:
    def __init__(self, tmp):
        self.path = Path(tmp)
        self.git('init', '-q', '-b', 'main')
        self.git('config', 'user.email', 'test@example.com')
        self.git('config', 'user.name', 'Test')
        self.write('base.txt', 'base\n')
        self.commit('initial')

    def git(self, *args):
        return subprocess.run(['git', *args], cwd=self.path, check=True,
                              capture_output=True, text=True)

    def write(self, rel_path, content):
        path = self.path / rel_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding='utf-8')
        self.git('add', rel_path)

    def commit(self, message='x'):
        self.git('commit', '-q', '-m', message)

    def branch(self, name, start='main'):
        self.git('checkout', '-q', '-b', name, start)

    def checkout(self, name):
        self.git('checkout', '-q', name)


class SquashStateTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = _Repo(self.tmp.name)

    def test_a_branch_whose_content_is_already_on_main_is_landed(self):
        """
        The case a commit count gets wrong: the branch has its own
        commit, but main was squash-merged from it, so the merged tree
        is main's tree and re-squashing would re-add landed content.
        """
        self.repo.branch('feature')
        self.repo.write('feature.txt', 'work\n')
        self.repo.commit('feature work')
        self.repo.checkout('main')
        self.repo.git('merge', '--squash', 'feature')
        self.repo.commit('squashed feature')

        state = squash_state('feature', base='main', cwd=self.repo.path)
        self.assertTrue(state.landed)
        self.assertTrue(state.clean)
        self.assertEqual(state.conflicts, [])

    def test_an_unlanded_branch_that_merges_cleanly(self):
        self.repo.branch('feature')
        self.repo.write('feature.txt', 'work\n')
        self.repo.commit('feature work')
        self.repo.checkout('main')

        state = squash_state('feature', base='main', cwd=self.repo.path)
        self.assertFalse(state.landed)
        self.assertTrue(state.clean)
        self.assertEqual(state.conflicts, [])

    def test_a_conflicting_branch_names_the_conflicted_files(self):
        self.repo.branch('feature')
        self.repo.write('shared.txt', 'branch version\n')
        self.repo.commit('branch edit')
        self.repo.checkout('main')
        self.repo.write('shared.txt', 'main version\n')
        self.repo.commit('main edit')

        state = squash_state('feature', base='main', cwd=self.repo.path)
        self.assertFalse(state.landed)
        self.assertFalse(state.clean)
        self.assertEqual(state.conflicts, ['shared.txt'])

    def test_main_against_itself_is_landed(self):
        state = squash_state('main', base='main', cwd=self.repo.path)
        self.assertTrue(state.landed)


class SurveyTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = _Repo(self.tmp.name)

    def test_covers_every_local_branch(self):
        self.repo.branch('feature')
        self.repo.write('feature.txt', 'work\n')
        self.repo.commit('feature work')
        self.repo.checkout('main')

        rows = survey(base='main', cwd=self.repo.path)
        self.assertEqual(sorted(row.branch for row in rows), ['feature', 'main'])
        feature = next(row for row in rows if row.branch == 'feature')
        self.assertFalse(feature.landed)
        self.assertIsNone(feature.worktree)

    def test_reports_a_branch_checked_out_in_a_worktree_with_its_dirt(self):
        worktree_dir = Path(self.tmp.name).parent / 'survey_wt'
        self.repo.branch('feature')
        self.repo.checkout('main')
        self.repo.git('worktree', 'add', '-q', str(worktree_dir), 'feature')
        self.addCleanup(lambda: subprocess.run(
            ['git', 'worktree', 'remove', '--force', str(worktree_dir)],
            cwd=self.repo.path, capture_output=True, text=True))
        (worktree_dir / 'scratch.txt').write_text('untracked\n', encoding='utf-8')

        by_branch = worktrees_by_branch(cwd=self.repo.path)
        self.assertIn('feature', by_branch)
        self.assertEqual(by_branch['feature'].uncommitted, 1)

        rows = survey(base='main', cwd=self.repo.path)
        feature = next(row for row in rows if row.branch == 'feature')
        self.assertIsNotNone(feature.worktree)
        self.assertEqual(feature.worktree.uncommitted, 1)

    def test_report_marks_landed_branches(self):
        self.repo.branch('feature')
        self.repo.write('feature.txt', 'work\n')
        self.repo.commit('feature work')
        self.repo.checkout('main')
        self.repo.git('merge', '--squash', 'feature')
        self.repo.commit('squashed feature')

        report = format_report(survey(base='main', cwd=self.repo.path))
        self.assertIn('feature', report)
        self.assertIn('landed', report.lower())


if __name__ == '__main__':
    unittest.main()
