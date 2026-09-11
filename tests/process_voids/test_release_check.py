import subprocess
import tempfile
import unittest
from pathlib import Path

from process_voids.util.release_check import (
    _next_steps, _pyproject_field, check_no_conflict_markers, check_no_direct_dependencies,
)


def _git_repo_with_file(tmp, filename, content):
    subprocess.run(['git', 'init', '-q'], cwd=tmp, check=True)
    # a fresh `git init` has no identity configured on CI-style machines;
    # commits below would fail without one, so set a throwaway local identity.
    subprocess.run(['git', 'config', 'user.email', 'test@example.com'], cwd=tmp, check=True)
    subprocess.run(['git', 'config', 'user.name', 'Test'], cwd=tmp, check=True)
    (Path(tmp) / filename).write_text(content, encoding='utf-8')
    subprocess.run(['git', 'add', filename], cwd=tmp, check=True)
    subprocess.run(['git', 'commit', '-q', '-m', 'x'], cwd=tmp, check=True)


def _pyproject(tmp, dependencies):
    path = Path(tmp) / 'pyproject.toml'
    lines = ',\n'.join(f'    "{dep}"' for dep in dependencies)
    path.write_text(f'[project]\nname = "x"\nversion = "0.1.0"\n'
                    f'dependencies = [\n{lines}\n]\n', encoding='utf-8')
    return path


class DirectDependencyCheckTest(unittest.TestCase):
    """
    A git-tag or file-path pin (a PEP 508 direct reference, 'name @ url')
    is fine while developing but must not reach a release - it has to
    depend on a published version. The check fails and names each such
    line.
    """

    def test_git_tag_dependency_is_flagged(self):
        dep = 'skipalignments @ git+https://github.com/adamburkegh/skip-alignments@v0.2.3+p1'
        with tempfile.TemporaryDirectory() as tmp:
            passed, message = check_no_direct_dependencies(_pyproject(tmp, [dep, 'pm4py>=2.7']))
        self.assertFalse(passed)
        self.assertIn(dep, message)
        self.assertNotIn('pm4py', message)

    def test_file_path_dependency_is_flagged(self):
        dep = 'skipalignments @ file:../skip-alignments'
        with tempfile.TemporaryDirectory() as tmp:
            passed, message = check_no_direct_dependencies(_pyproject(tmp, [dep]))
        self.assertFalse(passed)
        self.assertIn(dep, message)

    def test_published_versions_pass(self):
        deps = ['pm4py>=2.7,<=2.8', 'skipalignments==0.2.2', "rustxes>=0.2.10; python_version >= '3.10'"]
        with tempfile.TemporaryDirectory() as tmp:
            passed, _message = check_no_direct_dependencies(_pyproject(tmp, deps))
        self.assertTrue(passed)


class ConflictMarkerCheckTest(unittest.TestCase):
    """
    A tracked file left with unresolved '<<<<<<<'/'>>>>>>>' merge markers
    must block a release rather than slip through unnoticed.
    """

    def test_conflict_markers_are_flagged(self):
        content = '<<<<<<< HEAD\nours\n=======\ntheirs\n>>>>>>> branch\n'
        with tempfile.TemporaryDirectory() as tmp:
            _git_repo_with_file(tmp, 'a.txt', content)
            passed, message = check_no_conflict_markers(cwd=tmp)
        self.assertFalse(passed)
        self.assertIn('a.txt', message)

    def test_a_bare_separator_line_is_not_flagged(self):
        # '=======' alone is common outside conflicts (e.g. a markdown
        # heading underline) and isn't itself checked.
        content = 'Title\n=======\nbody text\n'
        with tempfile.TemporaryDirectory() as tmp:
            _git_repo_with_file(tmp, 'a.txt', content)
            passed, _message = check_no_conflict_markers(cwd=tmp)
        self.assertTrue(passed)

    def test_clean_tree_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            _git_repo_with_file(tmp, 'a.txt', 'nothing to see here\n')
            passed, _message = check_no_conflict_markers(cwd=tmp)
        self.assertTrue(passed)


class NextStepsTest(unittest.TestCase):
    """
    The printed next steps match how releases are actually done here:
    commit and push locally, then create the release and its tag in the
    GitHub web interface - never git tag / git push --tags.
    """

    def setUp(self):
        self.steps = _next_steps()
        self.version = _pyproject_field("version")

    def test_ends_with_a_github_release_step(self):
        self.assertIn(f"create release v{self.version} on GitHub", self.steps)

    def test_never_tags_locally(self):
        self.assertNotIn("git tag", self.steps)
        self.assertNotIn("--tags", self.steps)

    def test_still_commits_and_pushes(self):
        self.assertIn(f'git commit -m "Release v{self.version}"', self.steps)
        self.assertIn("git push", self.steps)


if __name__ == '__main__':
    unittest.main()
