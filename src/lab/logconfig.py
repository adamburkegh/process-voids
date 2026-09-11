"""
Console + file logging setup shared by the experiment entrypoints.

Times/status go through logging (INFO by default, DEBUG for the more
chatty per-cell start markers) rather than print(), so verbosity is a
level away instead of a code change.
"""

import importlib.metadata
import json
import logging
import subprocess
import sys
import time
import tomllib
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import url2pathname

from skipalignments.progress import disable_progress_bars

LOG_DIR = Path('var/lab/logs')


def project_version(repo_dir='.'):
    """
    pyproject.toml's [project].version at repo_dir, or None if missing/
    unreadable/absent - the PRIMARY "what code is this" identifier for
    a run's log (deliberately bumped, meaningful across an experiment
    runner checkout that sits still for a while - see configure). The
    commit hash from git_version is supplementary: always precise, but
    changes on every commit rather than marking a deliberate release
    point.
    """
    try:
        with open(Path(repo_dir) / 'pyproject.toml', 'rb') as f:
            data = tomllib.load(f)
        return data['project']['version']
    except (OSError, KeyError, tomllib.TOMLDecodeError):
        return None


def git_version(repo_dir='.'):
    """
    (short_hash, dirty) for the git repo at repo_dir - dirty is True if
    there are any uncommitted changes (tracked or untracked). Returns
    (None, None) if unavailable (not a git checkout, git not on PATH,
    etc) - a diagnostic nicety, never worth failing a run over.
    """
    try:
        commit = subprocess.run(
            ['git', '-C', str(repo_dir), 'rev-parse', '--short=12', 'HEAD'],
            capture_output=True, text=True, timeout=5, check=True).stdout.strip()
        status = subprocess.run(
            ['git', '-C', str(repo_dir), 'status', '--porcelain'],
            capture_output=True, text=True, timeout=5, check=True).stdout
        return commit, bool(status.strip())
    except (OSError, subprocess.SubprocessError):
        return None, None


def installed_package_version(package_name):
    """
    pip-installed version of `package_name` (from importlib.metadata),
    or None if it isn't installed. For an editable install this is
    whatever version its OWN pyproject.toml declared at install time -
    see editable_source_dir for that package's current git state, which
    can have moved on since without the pip-recorded version changing
    (an editable dependency is exactly the "shifting dev env" case this
    matters for).
    """
    try:
        return importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        return None


def editable_source_dir(package_name):
    """
    Source directory for `package_name` if it's installed editable (PEP
    660's direct_url.json, written by pip -e), else None - lets a
    dependency's OWN git_version be logged alongside process-voids' own,
    without process-voids hardcoding a sibling-repo path that would
    break the moment that checkout moves (a different machine, a
    worktree, a clone under a different name).
    """
    try:
        dist = importlib.metadata.distribution(package_name)
        direct_url = json.loads(dist.read_text('direct_url.json') or '{}')
    except (importlib.metadata.PackageNotFoundError, OSError, ValueError):
        return None
    if not direct_url.get('dir_info', {}).get('editable'):
        return None
    url = direct_url.get('url', '')
    if not url.startswith('file:'):
        return None
    return Path(url2pathname(urlparse(url).path))


def _version_line(version, commit, dirty):
    git_part = f'git {commit}{" dirty" if dirty else ""}' if commit is not None else 'git unknown'
    return f'{version if version is not None else "unknown"} ({git_part})'


def dependency_version_line(package_name):
    """
    '<pip version> (git <hash>[ dirty])' for `package_name` - same
    shape as configure()'s own process-voids line, built from
    installed_package_version/editable_source_dir+git_version instead
    of project_version/git_version('.') directly, since a dependency
    isn't running from repo_dir='.'.
    """
    version = installed_package_version(package_name)
    source_dir = editable_source_dir(package_name)
    commit, dirty = git_version(source_dir) if source_dir is not None else (None, None)
    return _version_line(version, commit, dirty)


def configure(level=logging.INFO, log_dir=LOG_DIR):
    """
    Logs to console AND a timestamped file under log_dir - relying on
    console output alone means the only record of a run is wherever its
    stdout happened to be captured, which is exactly what's fragile
    (a `| tail` swallowing everything until exit, a background task's
    capture getting discarded once it's stopped or superseded). The file
    handler makes a run's log durable and independent of how the
    command was invoked, and the timestamp in its name (script name +
    YYYYMMDD-HHMMSS) means a later run can never silently overwrite an
    earlier one's record the way a fixed path could.

    Also logs the process-voids version this run executed against,
    immediately after the log path - useful any time results might
    outlive the code that produced them (an experiment runner checkout
    that sits still for days while dev moves on elsewhere, a run
    resumed from a background task days later), so the log itself says
    which version ran rather than relying on remembering, or on a
    results CSV's mtime. pyproject.toml's version is primary (a
    deliberate, meaningful release point); the git commit/dirty state
    is supplementary, in parentheses - always precise, but moves on
    every commit rather than marking anything deliberate.

    A second line logs skip-alignments the same way (dependency_version_
    line) - its installed version, plus its OWN git state when it's
    installed editable (which can move independently of process-voids');
    installed from a release or git tag, the git part reads 'unknown'.
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    script_name = Path(sys.argv[0]).stem or 'run'
    stamp = time.strftime('%Y%m%d-%H%M%S')
    log_path = log_dir / f'{script_name}_{stamp}.log'

    logging.basicConfig(
        level=level, format='%(asctime)s %(levelname)s %(message)s', datefmt='%H:%M:%S',
        handlers=[logging.StreamHandler(), logging.FileHandler(log_path, encoding='utf-8')])
    disable_progress_bars()
    log = logging.getLogger(__name__)
    log.info('Logging to %s', log_path)

    log.info('process-voids version: %s', _version_line(project_version(), *git_version()))
    log.info('skip-alignments version: %s', dependency_version_line('skipalignments'))


def enable_skipalignments_debug():
    """
    skip-alignments' execution.py logs waste-ratio/per-variant-timing
    detail at DEBUG under its own logger, silent by default. Call this
    to opt in when investigating a slow run.
    """
    logging.getLogger('skipalignments.execution').setLevel(logging.DEBUG)
