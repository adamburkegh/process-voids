import subprocess
import tempfile
import unittest
from pathlib import Path

from lab.logconfig import (
    dependency_version_line, editable_source_dir, git_version,
    installed_package_version, project_version,
)


class ProjectVersionTest(unittest.TestCase):
    """project_version('.') against this real checkout's pyproject.toml."""

    def test_real_checkout_returns_a_dotted_version_string(self):
        version = project_version('.')
        self.assertIsInstance(version, str)
        self.assertRegex(version, r'^\d+\.\d+(\.\d+)?')

    def test_missing_pyproject_returns_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(project_version(tmp))

    def test_pyproject_without_project_version_returns_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'pyproject.toml').write_text('[tool.other]\nx = 1\n')
            self.assertIsNone(project_version(tmp))


class GitVersionTest(unittest.TestCase):
    """git_version('.') against this real checkout - not mocked, so a
    change to its subprocess invocations gets caught against real git
    output shape, not an assumption about it."""

    def test_real_checkout_returns_a_short_hash_and_bool_dirty(self):
        commit, dirty = git_version('.')
        self.assertIsInstance(commit, str)
        self.assertEqual(len(commit), 12)
        self.assertIsInstance(dirty, bool)

    def test_non_git_directory_returns_none_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(git_version(tmp), (None, None))

    def test_dirty_true_when_working_tree_has_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run(['git', 'init'], cwd=tmp, capture_output=True, check=True)
            subprocess.run(['git', 'config', 'user.email', 'a@b.c'], cwd=tmp, check=True)
            subprocess.run(['git', 'config', 'user.name', 'Test'], cwd=tmp, check=True)
            (Path(tmp) / 'f.txt').write_text('one')
            subprocess.run(['git', 'add', 'f.txt'], cwd=tmp, check=True)
            subprocess.run(['git', 'commit', '-m', 'init'], cwd=tmp, capture_output=True,
                            check=True)
            commit, dirty = git_version(tmp)
            self.assertIsNotNone(commit)
            self.assertFalse(dirty)

            (Path(tmp) / 'f.txt').write_text('two')
            _, dirty = git_version(tmp)
            self.assertTrue(dirty)


class SkipAlignmentsDependencyVersionTest(unittest.TestCase):
    """skip-alignments is installed non-editable in this dev env (a
    published release or a git tag, per pyproject.toml), so
    installed_package_version resolves it but editable_source_dir has no
    source tree to inspect, and its git state reports as 'unknown'."""

    def test_installed_package_version_is_a_dotted_string(self):
        version = installed_package_version('skipalignments')
        self.assertIsInstance(version, str)
        self.assertRegex(version, r'^\d+\.\d+(\.\d+)?')

    def test_uninstalled_package_returns_none(self):
        self.assertIsNone(installed_package_version('not-a-real-package-xyz'))

    def test_pinned_release_has_no_editable_source_dir(self):
        self.assertIsNone(editable_source_dir('skipalignments'))

    def test_uninstalled_package_has_no_editable_source_dir(self):
        self.assertIsNone(editable_source_dir('not-a-real-package-xyz'))

    def test_dependency_version_line_has_version_but_no_git_state(self):
        line = dependency_version_line('skipalignments')
        self.assertEqual(line, f"{installed_package_version('skipalignments')} (git unknown)")

    def test_uninstalled_package_line_says_unknown(self):
        line = dependency_version_line('not-a-real-package-xyz')
        self.assertEqual(line, 'unknown (git unknown)')


class EditableDependencyVersionTest(unittest.TestCase):
    """process-voids is itself editable-installed in this dev env (`pip
    install -e .` - see pip show process-voids) - the stand-in for
    exercising editable_source_dir/dependency_version_line's editable/
    direct_url.json path, since skip-alignments isn't installed editable
    (see SkipAlignmentsDependencyVersionTest)."""

    def test_editable_source_dir_is_a_real_directory_with_pyproject(self):
        source_dir = editable_source_dir('process-voids')
        self.assertIsNotNone(source_dir)
        self.assertTrue((source_dir / 'pyproject.toml').is_file())

    def test_dependency_version_line_has_version_and_git_parts(self):
        line = dependency_version_line('process-voids')
        self.assertRegex(line, r'^\d+\.\d+(\.\d+)? \(git [0-9a-f]{12}( dirty)?\)$')


if __name__ == '__main__':
    unittest.main()
