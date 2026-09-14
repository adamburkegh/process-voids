"""
Every local branch: where it is checked out, whether its content has
already landed, and whether it would squash cleanly.

"N commits ahead" answers none of that in a squash workflow. A branch
squash-merged into main keeps its own commits and still reads as ahead,
while being content-identical to main - and re-squashing such a branch
re-adds what already landed. So "landed" here means the merged tree
equals the base's tree (`git merge-tree --write-tree`), not any
relationship between commit counts.

Reports, does not gate: it exits non-zero only if it cannot do its job.

    bash run.sh python -m process_voids.util.branch_survey
    bash run.sh python -m process_voids.util.branch_survey --base main
"""

import argparse
import subprocess
import sys
from dataclasses import dataclass, field

from process_voids.util.gitread import PROJECT_ROOT, GitError, run_git


def _git_allowing_conflict(args, cwd=None):
    """
    git, without raising on exit 1 - merge-tree uses it to mean
    "conflicts", which is a finding here rather than an error.
    """
    return subprocess.run(['git', *args], cwd=cwd or PROJECT_ROOT,
                          capture_output=True, text=True)


@dataclass
class Worktree:
    path: str
    branch: str
    uncommitted: int


@dataclass
class SquashState:
    landed: bool          # merging into base yields base's own tree
    clean: bool           # merges without conflict
    conflicts: list = field(default_factory=list)


@dataclass
class BranchRow:
    branch: str
    landed: bool
    clean: bool
    conflicts: list
    worktree: Worktree = None


def local_branches(cwd=None) -> list:
    out = run_git(['for-each-ref', '--format=%(refname:short)', 'refs/heads/'], cwd=cwd)
    return [line.strip() for line in out.splitlines() if line.strip()]


def _uncommitted_count(path) -> int:
    """Changed-or-untracked files in that worktree, as `git status` counts them."""
    try:
        out = run_git(['-C', str(path), 'status', '--porcelain'], cwd=path)
    except GitError:
        return 0
    return len([line for line in out.splitlines() if line.strip()])


def worktrees_by_branch(cwd=None) -> dict:
    """
    {branch -> Worktree} for every checked-out worktree. A detached
    worktree has no branch and is skipped: nothing in this survey is
    keyed by it.
    """
    out = run_git(['worktree', 'list', '--porcelain'], cwd=cwd)
    result = {}
    path = None
    for line in out.splitlines():
        if line.startswith('worktree '):
            path = line[len('worktree '):].strip()
        elif line.startswith('branch ') and path is not None:
            branch = line[len('branch '):].strip()
            branch = branch.removeprefix('refs/heads/')
            result[branch] = Worktree(path=path, branch=branch,
                                      uncommitted=_uncommitted_count(path))
    return result


def squash_state(branch: str, base: str = 'main', cwd=None) -> SquashState:
    """
    What squashing `branch` onto `base` would do now.

    `git merge-tree --write-tree` performs the merge in memory: exit 0
    prints the merged tree's oid, exit 1 means conflicts. Landed is that
    merged tree being byte-identical to the base's own tree - ie the
    merge would add nothing.

    On conflict the output is the oid, then (with --name-only) the
    conflicted paths, then a BLANK LINE and git's own informational
    messages ('Auto-merging x', 'CONFLICT (add/add): ...'). Only the
    section before that blank line is a list of paths.
    """
    result = _git_allowing_conflict(
        ['merge-tree', '--write-tree', '--name-only', base, branch], cwd=cwd)
    if result.returncode > 1:
        raise GitError(f'git merge-tree {base} {branch} failed: {result.stderr.strip()}')

    lines = result.stdout.splitlines()
    merged_tree = lines[0].strip() if lines else ''
    if result.returncode == 1:
        conflicts = []
        for line in lines[1:]:
            if not line.strip():
                break
            conflicts.append(line.strip())
        return SquashState(landed=False, clean=False, conflicts=conflicts)

    base_tree = run_git(['rev-parse', f'{base}^{{tree}}'], cwd=cwd).strip()
    return SquashState(landed=merged_tree == base_tree, clean=True)


def survey(base: str = 'main', cwd=None) -> list:
    worktrees = worktrees_by_branch(cwd=cwd)
    rows = []
    for branch in local_branches(cwd=cwd):
        state = squash_state(branch, base=base, cwd=cwd)
        rows.append(BranchRow(branch=branch, landed=state.landed, clean=state.clean,
                              conflicts=state.conflicts,
                              worktree=worktrees.get(branch)))
    return rows


def format_report(rows, base: str = 'main') -> str:
    lines = [f'{len(rows)} local branches against {base}:', '']
    for row in sorted(rows, key=lambda r: (not r.landed, r.branch)):
        if row.landed:
            status = 'landed (content already on base - re-squashing would re-add it)'
        elif row.clean:
            status = 'unlanded, squashes cleanly'
        else:
            status = f'unlanded, CONFLICTS in {len(row.conflicts)} file(s)'
        lines.append(f'  {row.branch}: {status}')
        for path in row.conflicts:
            lines.append(f'      conflict: {path}')
        if row.worktree:
            dirt = (f'{row.worktree.uncommitted} uncommitted file(s)'
                    if row.worktree.uncommitted else 'clean')
            lines.append(f'      worktree: {row.worktree.path} ({dirt})')
    return '\n'.join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description='Local branches: worktree, whether landed, whether they squash '
                    'cleanly. A report, not a gate - exits non-zero only on error.')
    parser.add_argument('--base', default='main',
                        help='branch to compare against (default: main)')
    args = parser.parse_args(argv)

    try:
        rows = survey(base=args.base)
    except GitError as e:
        print(f'error: {e}', file=sys.stderr)
        return 1

    print(format_report(rows, base=args.base))
    return 0


if __name__ == '__main__':
    sys.exit(main())
