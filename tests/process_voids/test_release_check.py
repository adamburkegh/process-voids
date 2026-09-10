import tempfile
import unittest
from pathlib import Path

from process_voids.util.release_check import (
    _next_steps, _pyproject_field, check_no_direct_dependencies,
)


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
