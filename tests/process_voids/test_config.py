"""
Tests for process_voids.config, the machine-specific settings read from
pvoid.toml.

Lookup is exercised against real throwaway git repositories, including a
linked worktree, since finding the main checkout's file from inside a
worktree is the case the lookup exists for. Nothing here reads this
machine's own pvoid.toml.
"""

import re
import subprocess
import tempfile
import tomllib
import unittest
from pathlib import Path

from process_voids.config import (
    CONFIG_FILENAME, EXAMPLE_FILENAME, SCHEMA, ConfigError, ebi_executable, find_config,
    load_config, value,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


def _git(cwd, *args):
    subprocess.run(['git', *args], cwd=cwd, check=True, capture_output=True, text=True)


def _repo(path):
    path.mkdir(parents=True, exist_ok=True)
    _git(path, 'init', '-q')
    _git(path, 'config', 'user.email', 'test@example.com')
    _git(path, 'config', 'user.name', 'Test')
    (path / 'README.md').write_text('x\n', encoding='utf-8')
    _git(path, 'add', 'README.md')
    _git(path, 'commit', '-q', '-m', 'x')
    return path


def _write(path, text):
    path.write_text(text, encoding='utf-8')
    return path


class FindConfigTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.main = _repo(self.root / 'main')

    def test_finds_the_file_at_the_checkout_root(self):
        config = _write(self.main / CONFIG_FILENAME, '')
        self.assertEqual(find_config(cwd=self.main), config)

    def test_finds_it_from_a_subdirectory(self):
        config = _write(self.main / CONFIG_FILENAME, '')
        (self.main / 'src' / 'deep').mkdir(parents=True)
        self.assertEqual(find_config(cwd=self.main / 'src' / 'deep'), config)

    def test_a_worktree_without_its_own_file_uses_the_main_checkouts(self):
        """Every session works in a worktree; one file in the main
        checkout has to serve them all, or each new worktree breaks on
        first use the way a missing venv does."""
        config = _write(self.main / CONFIG_FILENAME, '')
        worktree = self.root / 'wt'
        _git(self.main, 'worktree', 'add', '-q', str(worktree))
        self.assertEqual(find_config(cwd=worktree).resolve(), config.resolve())

    def test_a_worktrees_own_file_takes_precedence(self):
        _write(self.main / CONFIG_FILENAME, '')
        worktree = self.root / 'wt'
        _git(self.main, 'worktree', 'add', '-q', str(worktree))
        own = _write(worktree / CONFIG_FILENAME, '')
        self.assertEqual(find_config(cwd=worktree).resolve(), own.resolve())

    def test_none_where_there_is_no_file(self):
        self.assertIsNone(find_config(cwd=self.main))

    def test_none_outside_any_git_repository(self):
        """An installed package has no checkout at all."""
        outside = self.root / 'not_a_repo'
        outside.mkdir()
        self.assertIsNone(find_config(cwd=outside))


class LoadConfigTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def test_no_file_is_an_empty_config(self):
        self.assertEqual(load_config(None), {})

    def test_reads_declared_keys(self):
        path = _write(self.dir / CONFIG_FILENAME,
                      '[paths]\ndata_dir = "/srv/logs"\n')
        self.assertEqual(load_config(path), {'paths': {'data_dir': '/srv/logs'}})

    def test_an_unknown_key_is_rejected_rather_than_ignored(self):
        """A typo would otherwise read as an unset key, and the error
        would point at the wrong problem."""
        path = _write(self.dir / CONFIG_FILENAME, '[paths]\ndata_dri = "/srv/logs"\n')
        with self.assertRaises(ConfigError) as caught:
            load_config(path)
        self.assertIn('data_dri', str(caught.exception))

    def test_agent_keys_load_alongside_the_codes_own(self):
        """Machine facts only agents read - the interpreter to build a
        venv with, where the papers are - share the one file. Declared,
        or the loader would reject every pvoid.toml that set them."""
        path = _write(self.dir / CONFIG_FILENAME,
                      '[paths]\ndata_dir = "/srv/logs"\n\n'
                      '[agent]\npython = "/opt/python314/bin/python"\n'
                      'papers_dir = "/srv/papers"\n')
        self.assertEqual(load_config(path)['agent'],
                         {'python': '/opt/python314/bin/python', 'papers_dir': '/srv/papers'})

    def test_an_unknown_section_is_rejected(self):
        path = _write(self.dir / CONFIG_FILENAME, '[pathz]\ndata_dir = "/srv/logs"\n')
        with self.assertRaises(ConfigError) as caught:
            load_config(path)
        self.assertIn('pathz', str(caught.exception))

    def test_malformed_toml_names_the_file(self):
        path = _write(self.dir / CONFIG_FILENAME, '[paths\n')
        with self.assertRaises(ConfigError) as caught:
            load_config(path)
        self.assertIn(CONFIG_FILENAME, str(caught.exception))


class ValueTest(unittest.TestCase):

    def test_returns_a_configured_value(self):
        self.assertEqual(value('paths', 'data_dir', {'paths': {'data_dir': '/srv/logs'}}),
                         '/srv/logs')

    def test_falls_back_to_a_declared_default(self):
        self.assertEqual(value('tools', 'ebi', {}), 'ebi')

    def test_an_unset_key_without_a_default_says_what_it_is_and_where_to_set_it(self):
        with self.assertRaises(ConfigError) as caught:
            value('paths', 'data_dir', {})
        message = str(caught.exception)
        self.assertIn('data_dir', message)
        self.assertIn(EXAMPLE_FILENAME, message)
        self.assertIn(SCHEMA['paths']['data_dir'].description, message)

    def test_an_undeclared_key_is_a_programming_error_not_a_config_one(self):
        with self.assertRaises(KeyError):
            value('paths', 'no_such_key', {})


class EbiExecutableTest(unittest.TestCase):

    def test_defaults_to_the_bare_name_on_path(self):
        """skip-alignments' own default: resolved through PATH, so an
        installed package with no pvoid.toml still finds ebi."""
        self.assertEqual(ebi_executable({}), 'ebi')

    def test_a_configured_executable_wins(self):
        self.assertEqual(ebi_executable({'tools': {'ebi': '/opt/ebi/bin/ebi'}}),
                         '/opt/ebi/bin/ebi')


class AgentSectionTest(unittest.TestCase):
    """[agent] holds facts an agent or person reads from the file - the
    interpreter needed before any Python runs, the papers to read - so
    no code reads them. That is what keeps the section honest: a key the
    code starts depending on belongs in [paths] or [tools], with the
    fail-loudly handling those get."""

    def test_no_source_module_reads_an_agent_key(self):
        reads = re.compile(r"""value\(\s*['"]agent['"]""")
        offenders = [str(path.relative_to(REPO_ROOT))
                     for path in (REPO_ROOT / 'src').rglob('*.py')
                     if reads.search(path.read_text(encoding='utf-8'))]
        self.assertEqual(offenders, [])

    def test_every_agent_key_says_it_is_not_read_by_code(self):
        for key, declared in SCHEMA['agent'].items():
            with self.subTest(key=key):
                self.assertIn('not read by code', declared.description)


class ExampleFileTest(unittest.TestCase):
    """The tracked example is how a new machine learns what to set, so it
    must hold exactly the schema's keys - no more, no fewer."""

    def setUp(self):
        with open(REPO_ROOT / EXAMPLE_FILENAME, 'rb') as f:
            self.example = tomllib.load(f)

    def test_declares_exactly_the_schemas_sections_and_keys(self):
        self.assertEqual({section: set(keys) for section, keys in self.example.items()},
                         {section: set(keys) for section, keys in SCHEMA.items()})

    def test_the_example_itself_loads(self):
        self.assertIsInstance(load_config(REPO_ROOT / EXAMPLE_FILENAME), dict)

    def test_the_real_file_is_gitignored(self):
        """It holds this machine's absolute paths, which must not reach
        a public repository."""
        result = subprocess.run(['git', 'check-ignore', '-q', CONFIG_FILENAME],
                                cwd=REPO_ROOT, capture_output=True)
        self.assertEqual(result.returncode, 0, f'{CONFIG_FILENAME} is not gitignored')


if __name__ == '__main__':
    unittest.main()
