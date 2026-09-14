"""
Reading tracked file content at a git ref, for the squash-review tools
in this package (registry_diff, test_name_diff, branch_survey).

A "ref" here is anything git resolves - a branch, a tag, a commit - or
the literal INDEX, meaning the staging area rather than any commit.
That pair (HEAD against the index) is the squash-review case these
tools default to: what the staged change does, before it is committed.

None of this belongs in the package's analysis API; it lives beside
release_check for the same reason that does - maintainer tooling about
the project itself, not usage of it.
"""

import subprocess
from pathlib import Path

# src/process_voids/util/gitread.py -> repo root is three levels up
PROJECT_ROOT = Path(__file__).resolve().parents[3]

INDEX = 'INDEX'


class GitError(RuntimeError):
    """A git command failed - a broken ref, or not a repository."""


def is_index(ref: str) -> bool:
    return ref.upper() == INDEX


def run_git(args, cwd=None) -> str:
    """stdout of `git <args>`, raising GitError on a non-zero exit."""
    result = subprocess.run(['git', *args], cwd=cwd or PROJECT_ROOT,
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise GitError(f"git {' '.join(args)} failed: {result.stderr.strip()}")
    return result.stdout


def read_file(ref: str, path: str, cwd=None) -> str:
    """
    Content of `path` at `ref`. `git show :<path>` reads the index, so
    the index is just a ref whose prefix is empty - no separate code
    path, and no need to consult the working tree (which may hold
    unstaged edits that are not part of what is being reviewed).
    """
    spec = f':{path}' if is_index(ref) else f'{ref}:{path}'
    return run_git(['show', spec], cwd=cwd)


def list_files(ref: str, path: str, cwd=None) -> list:
    """
    Tracked files under `path` at `ref`, as repo-relative paths. `path`
    may name a directory or a single file; a path that doesn't exist at
    that ref yields an empty list rather than an error, since a
    directory legitimately may not exist on one side of a comparison.
    """
    try:
        if is_index(ref):
            out = run_git(['ls-files', '--cached', '--', path], cwd=cwd)
        else:
            out = run_git(['ls-tree', '-r', '--name-only', ref, '--', path], cwd=cwd)
    except GitError:
        return []
    return [line.strip() for line in out.splitlines() if line.strip()]
