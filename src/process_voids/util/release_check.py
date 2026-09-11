"""
Release hygiene script -- a handful of checks/actions that would
otherwise be hand-cranked (and hand-remembered) before every release:
regenerating requirements.txt, confirming pyproject.toml's version and
CHANGELOG.md's top entry agree, flagging git/file dependency pins,
running the full test suite, flagging untracked files sitting in the
repo that might get missed, and scanning tracked files for leftover
merge conflict markers.

Adapted from skip-alignments' own release_check.py, cut down for a
project that isn't published to PyPI: no `python -m build`/`twine
check` step (there's no dist artifact to validate here).

This performs those checks/actions directly, but deliberately does NOT
perform the release itself (no git tag/push, no git commit) -- those
are irreversible/external actions that stay a human decision. If every
check passes, it prints the exact commands for the remaining steps
instead of running them.

None of this belongs in the package's own analysis API -- it lives
here (rather than tests/ or a top-level scripts/ dir) because it's
maintainer tooling about the project itself, not usage of it.

Run from the project root with the project's own venv, e.g.:
    pvoid/Scripts/python.exe -m process_voids.util.release_check

Exits 0 if every check passes, 1 otherwise -- suitable for a
pre-release step, not (yet) wired into CI.
"""
import re
import subprocess
import sys
import tomllib
from pathlib import Path

# src/process_voids/util/release_check.py -> repo root is three levels up
PROJECT_ROOT = Path(__file__).resolve().parents[3]


def _run(cmd, cwd=PROJECT_ROOT):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def _pyproject_field(name: str) -> str:
    text = (PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(rf'^{name}\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if match is None:
        raise ValueError(f"pyproject.toml: could not find a top-level {name} = \"...\" line")
    return match.group(1)


def regenerate_requirements() -> str:
    """
    Overwrites requirements.txt's dependency lines with the current
    venv's `pip freeze`, stripping this package's own self-reference
    line (pip freeze emits a `-e git+...#egg=<name>` entry for whatever
    is installed editable in the same venv it's run from -- that's the
    project being released, not a dependency of it). The file's own
    explanatory header comment (everything before the first dependency
    line) is preserved rather than overwritten.

    The egg name in that self-reference is the underscored import name
    (process_voids), not necessarily pyproject.toml's own hyphenated
    `name` field (process-voids) -- matched with either separator
    rather than assuming they agree.
    """
    requirements_path = PROJECT_ROOT / "requirements.txt"
    existing_lines = requirements_path.read_text(encoding="utf-8").splitlines() \
        if requirements_path.exists() else []
    header = [line for line in existing_lines if line.startswith("#")]

    package_name = _pyproject_field("name")
    result = _run([sys.executable, "-m", "pip", "freeze"])
    if result.returncode != 0:
        raise RuntimeError(f"pip freeze failed:\n{result.stderr}")

    # Build the [-_] substitution BEFORE escaping, not after - escaping
    # first turns '-' into the literal two-character sequence '\-',
    # and a naive str.replace('-', '[-_]') on that corrupts it into
    # '\[-_]' (an escaped '[', not a character class) instead of the
    # intended one.
    name_pattern = "[-_]".join(re.escape(part) for part in package_name.split("-"))
    self_reference = re.compile(rf"#egg={name_pattern}\b")
    dep_lines = [line for line in result.stdout.splitlines() if not self_reference.search(line)]

    requirements_path.write_text("\n".join(header + dep_lines) + "\n", encoding="utf-8")
    return f"requirements.txt regenerated ({len(dep_lines)} entries, self-reference stripped)"


def check_version_consistency() -> tuple[bool, str]:
    """
    pyproject.toml's version must match CHANGELOG.md's topmost entry, and
    that entry must be a real dated release, not a lingering [Unreleased]
    (or missing a date -- Keep a Changelog's convention this project
    follows, e.g. `## [0.4.1] - 2026-09-10`).
    """
    pyproject_version = _pyproject_field("version")

    changelog_text = (PROJECT_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    heading = re.search(r"^## \[([^\]]+)\](?:\s*-\s*(\d{4}-\d{2}-\d{2}))?", changelog_text, re.MULTILINE)
    if heading is None:
        return False, "CHANGELOG.md: no '## [version] - date' heading found at all"

    changelog_version, changelog_date = heading.group(1), heading.group(2)

    if changelog_version == "Unreleased":
        return False, (
            f"CHANGELOG.md's top entry is still [Unreleased] -- pyproject.toml is at "
            f"{pyproject_version}; move the release notes under a dated "
            f"'## [{pyproject_version}] - YYYY-MM-DD' heading before releasing"
        )
    if changelog_date is None:
        return False, f"CHANGELOG.md's top entry ([{changelog_version}]) has no date"
    if changelog_version != pyproject_version:
        return False, (
            f"version mismatch: pyproject.toml says {pyproject_version}, "
            f"CHANGELOG.md's top entry says {changelog_version}"
        )
    return True, f"pyproject.toml and CHANGELOG.md agree: {pyproject_version} ({changelog_date})"


def check_no_direct_dependencies(pyproject_path=None) -> tuple[bool, str]:
    """
    A release depends on published versions only. A PEP 508 direct
    reference ('name @ git+https://...', 'name @ file:...') is fine while
    developing against an unreleased dependency, but has to be swapped for
    a published version before releasing.
    """
    path = pyproject_path or PROJECT_ROOT / "pyproject.toml"
    with open(path, "rb") as f:
        dependencies = tomllib.load(f).get("project", {}).get("dependencies", [])
    direct = [dep for dep in dependencies if "@" in dep.split(";")[0]]
    if direct:
        return False, ("direct-reference dependencies in pyproject.toml - pin a published "
                       "release before releasing: " + "; ".join(direct))
    return True, "all dependencies are published versions"


def run_tests() -> tuple[bool, str]:
    """Full test suite -- the release gate that actually matters here,
    with no PyPI packaging step (build/twine) to validate alongside it.

    -t (top-level dir) is passed explicitly as the repo root, not left
    to default to -s's value (tests/) - discover infers each module's
    dotted name and what to put on sys.path relative to -t, and without
    it 'tests' itself gets treated as the top level, breaking any test
    module's `from lab.xxx import ...`/`from process_voids.xxx import
    ...` that (correctly) assumes the repo root's layout.
    """
    result = _run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-t", "."])
    if result.returncode != 0:
        # unittest's own summary is on stderr
        tail = "\n".join(result.stderr.strip().splitlines()[-15:])
        return False, f"test suite failed:\n{tail}"
    return True, "full test suite passed"


def check_stray_files() -> tuple[bool, str]:
    """
    Untracked files are easy to either forget (they never make it into a
    commit and quietly vanish) or accidentally sweep in (a careless
    `git add -A`). Neither is what you want right before a release --
    flag them so it's a decision, not an accident.
    """
    result = _run(["git", "status", "--porcelain"])
    if result.returncode != 0:
        return False, f"git status failed:\n{result.stderr}"

    untracked = [line[3:] for line in result.stdout.splitlines() if line.startswith("??")]
    if untracked:
        return False, "untracked files present: " + ", ".join(untracked)
    return True, "no untracked files"


def check_no_conflict_markers(cwd=None) -> tuple[bool, str]:
    """
    An unresolved merge can leave '<<<<<<<'/'>>>>>>>' markers sitting in a
    tracked file - easy to miss in a large diff, and disastrous if it
    reaches a release. Only those two markers are checked: '=======' alone
    is common outside conflicts (markdown headings, ASCII separators), so
    checking it too would flag plenty of files that were never touched by
    a conflict.

    git grep only searches tracked files, so this can't flag markers
    sitting in an untracked or ignored file. check_stray_files flags
    untracked files separately; ignored files are checked by neither.
    """
    result = _run(["git", "grep", "-n", "-E", r"^(<{7}|>{7})"], cwd=cwd or PROJECT_ROOT)
    if result.returncode == 0:
        return False, "merge conflict markers found:\n" + result.stdout.strip()
    if result.returncode == 1:
        return True, "no merge conflict markers in tracked files"
    return False, f"git grep failed:\n{result.stderr}"


def _next_steps() -> str:
    version = _pyproject_field("version")
    return "\n".join([
        "All checks passed. Remaining steps are yours to run:",
        "",
        "  git add -A",
        f'  git commit -m "Release v{version}"',
        "  git push",
        f"  then create release v{version} on GitHub (tag v{version}, target main)",
    ])


def main() -> int:
    ok = True

    # requirements.txt regeneration and the stray-files check both act on
    # git-visible state, so run stray-files last: regenerating
    # requirements.txt only ever *modifies* a tracked file (never adds an
    # untracked one), so it doesn't change what counts as "stray" -- but
    # checking last still means the report reflects the repo exactly as
    # this run leaves it.
    print(regenerate_requirements())

    for check in (check_version_consistency, check_no_direct_dependencies, run_tests,
                  check_stray_files, check_no_conflict_markers):
        passed, message = check()
        print(("PASS: " if passed else "FAIL: ") + message)
        ok = ok and passed

    print()
    print(_next_steps() if ok else "One or more checks failed -- not printing release commands.")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
