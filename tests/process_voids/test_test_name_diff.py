"""
Tests for process_voids.util.test_name_diff, the squash-review report on
which test names a change drops.

Counts are not the question: a refactor that moves code between modules
can drop behaviour tests while the total rises, so these pin that a name
present at the base ref and absent at the target is reported however the
totals move.
"""

import subprocess
import tempfile
import unittest
from pathlib import Path

from process_voids.util.test_name_diff import (
    collect_test_names, diff_test_names, format_report, test_names_in_source,
)


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

    def remove(self, rel_path):
        self._git('rm', '-q', rel_path)

    def commit(self, message='x'):
        self._git('commit', '-q', '-m', message)


def _test_file(*names):
    body = '\n\n'.join(f'    def {name}(self):\n        pass' for name in names)
    return f'import unittest\n\n\nclass T(unittest.TestCase):\n{body}\n'


class TestNamesInSourceTest(unittest.TestCase):

    def test_finds_test_methods_and_functions(self):
        source = ('def test_module_level():\n    pass\n\n'
                  'class T:\n    def test_method(self):\n        pass\n')
        self.assertEqual(test_names_in_source(source), {'test_module_level', 'test_method'})

    def test_ignores_helpers_and_mere_mentions(self):
        source = ('def helper_test_thing():\n    pass\n\n'
                  '# def test_commented_out(self):\n'
                  'NAMES = ["test_in_a_string"]\n')
        self.assertEqual(test_names_in_source(source), set())


class CollectAndDiffTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = _Repo(self.tmp.name)

    def test_collects_across_a_directory_at_a_ref(self):
        self.repo.write('tests/test_a.py', _test_file('test_one'))
        self.repo.write('tests/sub/test_b.py', _test_file('test_two'))
        self.repo.commit()
        names = collect_test_names('HEAD', ['tests'], cwd=self.repo.path)
        self.assertEqual(sorted(names), ['test_one', 'test_two'])
        self.assertEqual(names['test_two'], {'tests/sub/test_b.py'})

    def test_reports_a_name_dropped_even_as_the_total_rises(self):
        """
        The failure this tool was written for: a refactor moved a cell
        loop to a new module, dropping behaviour tests while the total
        count went up.
        """
        self.repo.write('tests/test_a.py', _test_file('test_kept', 'test_dropped'))
        self.repo.commit()
        self.repo.write('tests/test_a.py', _test_file('test_kept'))
        self.repo.write('tests/test_new.py', _test_file('test_added_1', 'test_added_2'))
        self.repo.commit()

        before = collect_test_names('HEAD~1', ['tests'], cwd=self.repo.path)
        after = collect_test_names('HEAD', ['tests'], cwd=self.repo.path)
        self.assertGreater(len(after), len(before))

        result = diff_test_names(before, after)
        self.assertEqual(result.dropped, ['test_dropped'])
        self.assertEqual(result.added, ['test_added_1', 'test_added_2'])

    def test_a_name_moved_between_files_is_not_dropped(self):
        self.repo.write('tests/test_a.py', _test_file('test_moves'))
        self.repo.commit()
        self.repo.remove('tests/test_a.py')
        self.repo.write('tests/test_b.py', _test_file('test_moves'))
        self.repo.commit()

        result = diff_test_names(
            collect_test_names('HEAD~1', ['tests'], cwd=self.repo.path),
            collect_test_names('HEAD', ['tests'], cwd=self.repo.path))
        self.assertEqual(result.dropped, [])
        self.assertEqual(result.moved, [('test_moves', {'tests/test_a.py'},
                                         {'tests/test_b.py'})])

    def test_collects_staged_content_from_the_index(self):
        self.repo.write('tests/test_a.py', _test_file('test_committed'))
        self.repo.commit()
        self.repo.write('tests/test_a.py', _test_file('test_staged'))
        names = collect_test_names('INDEX', ['tests'], cwd=self.repo.path)
        self.assertEqual(sorted(names), ['test_staged'])

    def test_a_missing_path_at_one_ref_is_not_an_error(self):
        self.repo.write('tests/test_a.py', _test_file('test_one'))
        self.repo.commit()
        self.assertEqual(collect_test_names('HEAD', ['no_such_dir'],
                                            cwd=self.repo.path), {})

    def test_report_names_the_dropped_tests(self):
        result = diff_test_names({'test_gone': {'tests/test_a.py'}},
                                 {'test_new': {'tests/test_a.py'}})
        report = format_report(result)
        self.assertIn('test_gone', report)
        self.assertIn('tests/test_a.py', report)

    def test_report_says_so_when_nothing_was_dropped(self):
        result = diff_test_names({'test_a': {'f.py'}}, {'test_a': {'f.py'}})
        self.assertIn('no test names dropped', format_report(result).lower())


if __name__ == '__main__':
    unittest.main()
