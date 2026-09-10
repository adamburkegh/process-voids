import unittest

from process_voids.util.release_check import _next_steps, _pyproject_field


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
