"""
The whole squash review in one command, between two refs.

Combines, in order: the diff stat; registry_diff's report on
lab.metric_registry; test_name_diff's report on dropped test names; a
scan of ADDED lines for anything that must not reach this public
repository; and a note where src/ changed without CHANGELOG.md. The
first two reports are called as functions, not shelled out to, so this
reads refs exactly as they do.

Reports, does not gate: it exits non-zero only if it cannot do its job
(a bad ref, an unparseable registry), never on a finding. release_check
is the gate. It deliberately does not run the test suite, which is slow,
runs in the background, and has its own fixed, allowlisted command.

    bash run.sh python -m process_voids.util.squash_review
    bash run.sh python -m process_voids.util.squash_review main my-branch

Defaults to HEAD against the index: the staged change under review.
"""

import argparse
import re
import sys

from process_voids.util import registry_diff, test_name_diff
from process_voids.util.gitread import INDEX, GitError, is_index, run_git

TEST_PATHS = ['tests']

# Things that must not reach this public repository, looked for in added
# lines only - a line that merely survives a change was already reviewed
# when it arrived. Each entry is (name, regex, what it catches). Add a
# pattern here and a case to tests/process_voids/test_squash_review.py.
FORBIDDEN_PATTERNS = (
    ('labnotes', r'labnotes',
     'the private lab notebook, kept in a separate private repo'),
    ('var_papers', r'var/papers',
     'gitignored local copies of papers'),
    ('var_reports', r'var/reports',
     'gitignored private reports and dispatches'),
    ('pvoid_lab', r'pvoid-lab',
     'the private lab archive repository'),
    ('claude_dir', r'\.claude/',
     'per-machine agent configuration and worktrees'),
    # A drive letter then a separator. The \b stops a URL scheme matching:
    # the 's' of 'https:/' follows a word character, so no boundary.
    ('absolute_path', r'\b[A-Za-z]:[\\/]',
     'an absolute Windows path, meaningful on one machine only'),
    # Anchored: a real leftover marker starts its line, where one inside
    # a string or a docstring is text about markers.
    ('conflict_marker', r'^(<<<<<<<|>>>>>>>)',
     'a leftover merge conflict marker'),
)

_COMPILED = tuple((name, re.compile(regex), why) for name, regex, why in FORBIDDEN_PATTERNS)

_ALL = frozenset(name for name, _regex, _why in FORBIDDEN_PATTERNS)

# path -> pattern names that file is allowed to contain.
#
# Only this module and its tests, which are allowed everything: the module
# has to spell out each pattern it looks for, and the tests have to
# exercise them, so scanning either against itself could only ever report
# its own pattern list. Machine paths belong in pvoid.toml, so no other
# file has a reason to hold one.
ALLOWANCES = {
    'src/process_voids/util/squash_review.py': _ALL,
    'tests/process_voids/test_squash_review.py': _ALL,
}

_HUNK = re.compile(r'^@@ -\d+(?:,\d+)? \+(\d+)(?:,\d+)? @@')


def _diff_args(base, target):
    """
    `git diff` arguments comparing `base` to `target`, where either may
    be the index. The index is reached through --cached, which compares
    a commit to the index; with the index as the base the direction is
    reversed.
    """
    if is_index(target):
        return ['diff', '--cached', base]
    if is_index(base):
        return ['diff', '--cached', '-R', target]
    return ['diff', base, target]


def diff_stat(base, target, cwd=None) -> str:
    return run_git([*_diff_args(base, target), '--stat'], cwd=cwd).rstrip()


def changed_files(base, target, cwd=None) -> list:
    out = run_git([*_diff_args(base, target), '--name-only'], cwd=cwd)
    return [line.strip() for line in out.splitlines() if line.strip()]


def _header_path(line):
    """'+++ b/src/x.py' -> 'src/x.py'; None for '+++ /dev/null'."""
    path = line[4:].strip().strip('"')
    if path == '/dev/null':
        return None
    return path[2:] if path.startswith('b/') else path


def parse_added_lines(diff_text) -> list:
    """
    [(path, new-side line number, text), ...] for every added line in a
    zero-context unified diff.

    A '+++ ' line is a file header only before that file's first hunk;
    inside a hunk every '+' line is an addition, including one whose own
    text begins with '++ '.
    """
    added = []
    path = None
    lineno = None
    in_hunk = False
    for line in diff_text.splitlines():
        if line.startswith('diff --git '):
            path, lineno, in_hunk = None, None, False
            continue
        if not in_hunk and line.startswith('+++ '):
            path = _header_path(line)
            continue
        hunk = _HUNK.match(line)
        if hunk:
            lineno = int(hunk.group(1))
            in_hunk = True
            continue
        if not in_hunk or path is None:
            continue
        if line.startswith('+'):
            added.append((path, lineno, line[1:]))
            lineno += 1
        # '-' removals do not advance the new side; with no context there
        # is nothing else to count, and '\ No newline' is a marker.
    return added


def added_lines(base, target, cwd=None) -> list:
    return parse_added_lines(run_git([*_diff_args(base, target), '-U0'], cwd=cwd))


def scan_added_lines(added, allowances=ALLOWANCES) -> list:
    """
    [(path, line number, pattern name, text), ...] for each added line
    matching a pattern its file is not allowed.
    """
    hits = []
    for path, lineno, text in added:
        allowed = allowances.get(path, frozenset())
        for name, regex, _why in _COMPILED:
            if name not in allowed and regex.search(text):
                hits.append((path, lineno, name, text))
    return hits


def changelog_note(changed):
    """
    A note naming the src/ files changed where CHANGELOG.md was not, or
    None. Worded as a report: a docstring-only change legitimately needs
    no entry, so whether one is owed is the reviewer's call.
    """
    if 'CHANGELOG.md' in changed:
        return None
    src = sorted(path for path in changed if path.startswith('src/'))
    if not src:
        return None
    listed = '\n'.join(f'  {path}' for path in src)
    return (f'src/ changed and CHANGELOG.md did not ({len(src)} file(s)) - '
            'check whether an entry is owed; a docstring-only change needs '
            f'none:\n{listed}')


def _format_hits(hits):
    if not hits:
        return 'No forbidden references in added lines.'
    why = {name: reason for name, _regex, reason in FORBIDDEN_PATTERNS}
    lines = [f'{len(hits)} hit(s):']
    for path, lineno, name, text in hits:
        lines.append(f'  {path}:{lineno}  [{name}] {why[name]}')
        lines.append(f'      {text.strip()}')
    return '\n'.join(lines)


def _section(title, body):
    return f'== {title} ==\n{body}'


def comparison_base(base, target, cwd=None):
    """
    The ref to compare `target` from: the merge base of the two when both
    are commits, else `base` itself.

    Reviewing a branch against a trunk that has moved on since the branch
    was cut would otherwise compare two trees, so the trunk's newer work
    reads as the branch removing it - registry ids reported removed, tests
    reported dropped. From the merge base, only the branch's own changes
    show. The index is built on its base commit, so a comparison involving
    it needs no adjustment.
    """
    if is_index(base) or is_index(target):
        return base
    return run_git(['merge-base', base, target], cwd=cwd).strip()


def build_report(base='HEAD', target=INDEX, cwd=None) -> str:
    """The combined report. Raises GitError or SyntaxError if a ref or
    the registry at it cannot be read."""
    start = comparison_base(base, target, cwd=cwd)
    registry = registry_diff.diff_registries(
        registry_diff.load_metrics(start, cwd=cwd),
        registry_diff.load_metrics(target, cwd=cwd))
    names = test_name_diff.diff_test_names(
        test_name_diff.collect_test_names(start, TEST_PATHS, cwd=cwd),
        test_name_diff.collect_test_names(target, TEST_PATHS, cwd=cwd))

    stat = diff_stat(start, target, cwd=cwd)
    changelog = changelog_note(changed_files(start, target, cwd=cwd))

    header = f'Squash review: {base} -> {target}'
    if start != base and start != run_git(['rev-parse', base], cwd=cwd).strip():
        header += f' (from their merge base {start[:12]})'

    return '\n\n'.join([
        header,
        _section('Diff stat', stat or 'No changes.'),
        _section('Registry', registry_diff.format_report(registry)),
        _section('Test names', test_name_diff.format_report(names)),
        _section('Public-repo hygiene',
                 _format_hits(scan_added_lines(added_lines(start, target, cwd=cwd)))),
        _section('CHANGELOG', changelog or 'No src/ change without a CHANGELOG entry.'),
    ])


def main(argv=None, cwd=None) -> int:
    parser = argparse.ArgumentParser(
        description='The whole squash review between two refs, as one report. '
                    'A report, not a gate - exits non-zero only on error.')
    parser.add_argument('base', nargs='?', default='HEAD',
                        help='ref to compare FROM (default: HEAD)')
    parser.add_argument('target', nargs='?', default=INDEX,
                        help=f'ref to compare TO (default: {INDEX}, the staging area)')
    args = parser.parse_args(argv)

    try:
        report = build_report(args.base, args.target, cwd=cwd)
    except (GitError, SyntaxError) as e:
        print(f'error: {e}', file=sys.stderr)
        return 1
    print(report)
    return 0


if __name__ == '__main__':
    sys.exit(main())
