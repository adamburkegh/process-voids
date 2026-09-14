"""
Which test names a change drops, between two refs.

Test counts don't answer this: a refactor that moves code between
modules can drop behaviour tests while the total rises, so a suite that
grew by ten can still have lost a regression test guarding a shipped
fix. This compares the SET of `def test_*` names, so a name that
disappears is reported however the totals move, and a name that merely
moved file is reported as moved rather than lost.

Reports, does not gate: it exits non-zero only if it cannot do its job.
release_check is the gate, and it runs the suite.

    bash run.sh python -m process_voids.util.test_name_diff
    bash run.sh python -m process_voids.util.test_name_diff main my-branch
    bash run.sh python -m process_voids.util.test_name_diff --paths tests/lab

Defaults to HEAD against the index: the staged change under review.
"""

import argparse
import re
import sys
from dataclasses import dataclass, field

from process_voids.util.gitread import INDEX, GitError, list_files, read_file

# A def whose name starts with test_, method or module-level function.
# Anchored to the line's indentation so a name inside a string or a
# comment is not collected.
_TEST_DEF = re.compile(r'^[ \t]*def[ \t]+(test_\w*)[ \t]*\(', re.MULTILINE)


@dataclass
class TestNameDiff:
    dropped: list = field(default_factory=list)
    added: list = field(default_factory=list)
    # [(name, files_before, files_after)] - same name, different file(s)
    moved: list = field(default_factory=list)
    # name -> files, at each ref
    before: dict = field(default_factory=dict)
    after: dict = field(default_factory=dict)


def test_names_in_source(source: str) -> set:
    """Every `def test_*` name defined in `source`."""
    return set(_TEST_DEF.findall(source))


def collect_test_names(ref: str, paths, cwd=None) -> dict:
    """
    {test name -> {files defining it}} across `paths` at `ref`. A path
    may be a directory or a single file, and one that doesn't exist at
    that ref contributes nothing - a directory legitimately may not
    exist on both sides of a comparison.
    """
    names = {}
    for path in paths:
        for file_path in list_files(ref, path, cwd=cwd):
            if not file_path.endswith('.py'):
                continue
            source = read_file(ref, file_path, cwd=cwd)
            for name in test_names_in_source(source):
                names.setdefault(name, set()).add(file_path)
    return names


def diff_test_names(before: dict, after: dict) -> TestNameDiff:
    moved = [(name, before[name], after[name])
             for name in sorted(set(before) & set(after))
             if before[name] != after[name]]
    return TestNameDiff(
        dropped=sorted(set(before) - set(after)),
        added=sorted(set(after) - set(before)),
        moved=moved,
        before=before,
        after=after,
    )


def format_report(result: TestNameDiff) -> str:
    lines = [f'{len(result.before)} test names before, {len(result.after)} after.']
    if result.dropped:
        lines.append('')
        lines.append(f'DROPPED ({len(result.dropped)}) - present at the base ref, '
                     'absent at the target:')
        for name in result.dropped:
            where = ', '.join(sorted(result.before[name]))
            lines.append(f'  - {name}  ({where})')
    else:
        lines.append('No test names dropped.')
    if result.moved:
        lines.append('')
        lines.append(f'Moved ({len(result.moved)}) - same name, different file:')
        for name, before_files, after_files in result.moved:
            lines.append(f'  ~ {name}: {", ".join(sorted(before_files))} -> '
                         f'{", ".join(sorted(after_files))}')
    if result.added:
        lines.append('')
        lines.append(f'Added ({len(result.added)}):')
        lines.extend(f'  + {name}' for name in result.added)
    return '\n'.join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description='Which test names a change drops, between two refs. A report, '
                    'not a gate - exits non-zero only on error.')
    parser.add_argument('base', nargs='?', default='HEAD',
                        help='ref to compare FROM (default: HEAD)')
    parser.add_argument('target', nargs='?', default=INDEX,
                        help=f'ref to compare TO (default: {INDEX}, the staging area)')
    parser.add_argument('--paths', nargs='+', default=['tests'],
                        help='test files or directories to compare (default: tests)')
    args = parser.parse_args(argv)

    try:
        before = collect_test_names(args.base, args.paths)
        after = collect_test_names(args.target, args.paths)
    except GitError as e:
        print(f'error: {e}', file=sys.stderr)
        return 1

    print(f'{", ".join(args.paths)}: {args.base} -> {args.target}')
    print(format_report(diff_test_names(before, after)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
