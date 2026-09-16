"""
Machine-specific settings, read from pvoid.toml.

Where a large log or an external tool lives is a fact about one machine,
not about the code, so it belongs neither in a module constant nor in an
environment variable. The file is gitignored, since it holds exactly the
absolute paths that must not reach a public repository; the tracked
pvoid.example.toml documents every key and is held to SCHEMA by test.

Lookup, first match wins:
  1. the root of the checkout the process is running in;
  2. the root of the main checkout, when running in a linked worktree.
So one pvoid.toml in the main checkout serves every worktree, and a
worktree can still override it with its own. Outside any git repository
- an installed package - there is no file, and only keys with a default
can be read.

Reading a key that is unset and has no default raises ConfigError naming
the key, what it is for, and the example file to copy. Nothing reads a
key until it is needed, so a machine missing one it never uses is not
affected.
"""

import subprocess
import tomllib
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

CONFIG_FILENAME = 'pvoid.toml'
EXAMPLE_FILENAME = 'pvoid.example.toml'


class ConfigError(RuntimeError):
    """pvoid.toml is malformed, declares an unknown key, or lacks one a
    caller needs."""


@dataclass(frozen=True)
class Key:
    description: str
    default: object = None   # None: no fallback, so reading it unset is an error


SCHEMA = {
    'paths': {
        'data_dir': Key('the directory holding the large event logs kept outside '
                        'the repository'),
        'toothpaste_dir': Key('the toothpaste checkout that lab.toothpaste_bridge '
                              'mines with'),
    },
    'tools': {
        'ebi': Key('the ebi executable skip-alignments calls to estimate skip '
                   'probabilities', default='ebi'),
    },
    # Facts an agent or person reads from this file, never code: the
    # interpreter is needed before any Python runs, and the papers are
    # for reading. Declared here only because load_config rejects
    # undeclared sections.
    'agent': {
        'python': Key('the Python 3.14 interpreter to create a worktree venv '
                      'with - not read by code'),
        'papers_dir': Key('the directory holding the relevant papers - not read '
                          'by code'),
    },
}


def _git_path(cwd, *args):
    """A path git reports, or None outside a repository or without git."""
    try:
        result = subprocess.run(['git', *args], cwd=cwd, capture_output=True, text=True)
    except OSError:
        return None
    if result.returncode != 0:
        return None
    return Path(result.stdout.strip())


def search_paths(cwd=None) -> list:
    """Candidate pvoid.toml locations in lookup order - see the module
    docstring. Empty outside any git repository."""
    checkout = _git_path(cwd, 'rev-parse', '--show-toplevel')
    if checkout is None:
        return []
    candidates = [checkout / CONFIG_FILENAME]
    common = _git_path(cwd, 'rev-parse', '--path-format=absolute', '--git-common-dir')
    if common is not None:
        main_checkout = common.parent
        if main_checkout.resolve() != checkout.resolve():
            candidates.append(main_checkout / CONFIG_FILENAME)
    return candidates


def find_config(cwd=None):
    """The pvoid.toml that applies from `cwd`, or None."""
    for candidate in search_paths(cwd):
        if candidate.is_file():
            return candidate
    return None


def load_config(path) -> dict:
    """
    The parsed contents of `path`, or {} where it is None. Any section or
    key SCHEMA does not declare is rejected: a misspelt key would
    otherwise read as unset, and the resulting error would point at the
    wrong problem.
    """
    if path is None:
        return {}
    try:
        with open(path, 'rb') as f:
            config = tomllib.load(f)
    except tomllib.TOMLDecodeError as e:
        raise ConfigError(f'{path} is not valid TOML: {e}') from e

    for section, keys in config.items():
        if section not in SCHEMA:
            raise ConfigError(f'{path}: unknown section [{section}]; '
                              f'known sections: {sorted(SCHEMA)}')
        if not isinstance(keys, dict):
            raise ConfigError(f'{path}: [{section}] must be a table')
        for key in keys:
            if key not in SCHEMA[section]:
                raise ConfigError(f'{path}: unknown key {key!r} in [{section}]; '
                                  f'known keys: {sorted(SCHEMA[section])}')
    return config


@lru_cache(maxsize=None)
def _default_config():
    return load_config(find_config())


def value(section, key, config=None):
    """
    The configured value of [section] key, else its declared default.

    `config` defaults to this process's pvoid.toml. An undeclared key is
    a KeyError - a mistake in the calling code, not in anyone's
    configuration.
    """
    declared = SCHEMA[section][key]
    if config is None:
        config = _default_config()
    configured = config.get(section, {}).get(key)
    if configured is not None:
        return configured
    if declared.default is not None:
        return declared.default
    raise ConfigError(
        f'[{section}] {key} is not set - {declared.description}. Copy '
        f'{EXAMPLE_FILENAME} to {CONFIG_FILENAME} in the main checkout (or this '
        f'worktree) and set it. Searched: '
        f'{[str(p) for p in search_paths()] or "no git checkout found"}')


def ebi_executable(config=None) -> str:
    """The ebi executable to call - bare 'ebi', found on PATH, unless
    pvoid.toml names one."""
    return value('tools', 'ebi', config)
