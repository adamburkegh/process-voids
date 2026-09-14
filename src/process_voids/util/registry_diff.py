"""
What a change did to lab.metric_registry, between two refs.

The registry is an append-only dictionary for every result CSV column
this project has ever written, so a review needs to see more than "the
file changed": which ids appeared or disappeared, what moved on the ids
that survived, and - the part no test catches - whether any history
entry was altered or dropped rather than added to. Extending an entry's
text is fine; rewriting or deleting one silently rewrites the meaning of
columns in result files that were written years earlier and are never
regenerated.

Reports, does not gate: it exits non-zero only if it cannot do its job
(a bad ref, an unparseable registry). release_check is the gate.

    bash run.sh python -m process_voids.util.registry_diff
    bash run.sh python -m process_voids.util.registry_diff main my-branch

Defaults to HEAD against the index: the staged change under review.
"""

import argparse
import sys
from dataclasses import dataclass, field

from process_voids.util.gitread import INDEX, GitError, read_file

REGISTRY_PATH = 'src/lab/metric_registry.py'

# Compared per surviving id. description is deliberately absent: it is
# prose that gets reworded constantly, and a diff of it is what `git
# diff` is for.
COMPARED_FIELDS = ('status', 'superseded_by', 'scripts', 'source', 'scale')

# A history entry that gained text is fine (the append-only rule allows
# adding to the record); one that was altered or dropped is not.
VIOLATION_KINDS = ('altered', 'dropped')


@dataclass
class RegistryDiff:
    added: list = field(default_factory=list)
    removed: list = field(default_factory=list)
    # id -> [(field, old_value, new_value)]
    changed: dict = field(default_factory=dict)
    # id -> [(history_key, 'added' | 'extended' | 'altered' | 'dropped')]
    history_notes: dict = field(default_factory=dict)
    # [(id, history_key, kind)] for the kinds that break the rule
    violations: list = field(default_factory=list)


def load_metrics(ref: str, cwd=None) -> dict:
    """
    lab.metric_registry's METRICS as it stands at `ref`.

    Executed in a throwaway namespace rather than imported: the point is
    to hold two versions of the same module in memory at once, which
    importing cannot do, and the registry is a literal dict of
    dataclasses whose only import is dataclasses itself.
    """
    source = read_file(ref, REGISTRY_PATH, cwd=cwd)
    namespace = {}
    exec(compile(source, f'{ref}:{REGISTRY_PATH}', 'exec'), namespace)
    try:
        return namespace['METRICS']
    except KeyError:
        raise GitError(f'{ref}:{REGISTRY_PATH} defines no METRICS')


def _history_note(old_text, new_text):
    if new_text == old_text:
        return None
    if new_text.startswith(old_text):
        return 'extended'
    return 'altered'


def _compare_history(old_metric, new_metric):
    old_history = getattr(old_metric, 'history', None) or {}
    new_history = getattr(new_metric, 'history', None) or {}
    notes = []
    for key, old_text in old_history.items():
        if key not in new_history:
            notes.append((key, 'dropped'))
            continue
        note = _history_note(old_text, new_history[key])
        if note:
            notes.append((key, note))
    notes.extend((key, 'added') for key in new_history if key not in old_history)
    return notes


def diff_registries(old: dict, new: dict) -> RegistryDiff:
    """
    What changed between two METRICS dicts - see RegistryDiff. A removed
    id takes its whole history with it, so every entry it held counts as
    dropped.
    """
    result = RegistryDiff(
        added=sorted(set(new) - set(old)),
        removed=sorted(set(old) - set(new)),
    )

    for metric_id in result.removed:
        history = getattr(old[metric_id], 'history', None) or {}
        notes = [(key, 'dropped') for key in history]
        if notes:
            result.history_notes[metric_id] = notes
            result.violations.extend((metric_id, key, kind) for key, kind in notes)

    for metric_id in sorted(set(old) & set(new)):
        old_metric, new_metric = old[metric_id], new[metric_id]
        changes = [(f, getattr(old_metric, f, None), getattr(new_metric, f, None))
                   for f in COMPARED_FIELDS
                   if getattr(old_metric, f, None) != getattr(new_metric, f, None)]
        if changes:
            result.changed[metric_id] = changes
        notes = _compare_history(old_metric, new_metric)
        if notes:
            result.history_notes[metric_id] = notes
            result.violations.extend((metric_id, key, kind) for key, kind in notes
                                     if kind in VIOLATION_KINDS)
    return result


def format_report(result: RegistryDiff) -> str:
    lines = []
    if result.added:
        lines.append(f'Added ids ({len(result.added)}):')
        lines.extend(f'  + {metric_id}' for metric_id in result.added)
    if result.removed:
        lines.append(f'Removed ids ({len(result.removed)}):')
        lines.extend(f'  - {metric_id}' for metric_id in result.removed)
    if result.changed:
        lines.append(f'Changed ids ({len(result.changed)}):')
        for metric_id, changes in result.changed.items():
            lines.append(f'  {metric_id}:')
            lines.extend(f'    {f}: {old!r} -> {new!r}' for f, old, new in changes)
    if result.history_notes:
        lines.append('History entries:')
        for metric_id, notes in result.history_notes.items():
            for key, kind in notes:
                lines.append(f'  {metric_id} [{key}]: {kind}')
    if result.violations:
        lines.append('')
        lines.append(f'APPEND-ONLY VIOLATIONS ({len(result.violations)}) - a history '
                     'entry is the only record of what an already-written CSV column '
                     'holds:')
        lines.extend(f'  {metric_id} [{key}]: {kind}'
                     for metric_id, key, kind in result.violations)
    if not lines:
        return 'No registry changes.'
    return '\n'.join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="What a change did to lab.metric_registry between two refs. "
                    "A report, not a gate - exits non-zero only on error.")
    parser.add_argument('base', nargs='?', default='HEAD',
                        help='ref to compare FROM (default: HEAD)')
    parser.add_argument('target', nargs='?', default=INDEX,
                        help=f'ref to compare TO (default: {INDEX}, the staging area)')
    args = parser.parse_args(argv)

    try:
        old = load_metrics(args.base)
        new = load_metrics(args.target)
    except (GitError, SyntaxError) as e:
        print(f'error: {e}', file=sys.stderr)
        return 1

    print(f'lab.metric_registry: {args.base} -> {args.target}')
    print(format_report(diff_registries(old, new)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
