# Working rules for agents

Local conventions for coding agents working in this repo. The user sets
design direction and handles all source control on `main`.

## Environment

- Windows. Python 3.14 at `C:\working\tools\python\python314`.
- The project venv is `pvoid/` at the repo root (gitignored). Run everything
  through `run.sh`, which activates it: `bash run.sh python -m ...`. Run from
  the repo root; don't prefix commands with `cd`.
- **ebi**: skip-alignments calls `ebi` by bare name, so it must be on PATH. It
  is installed at `C:\working\tools\ebi\ebi.exe` and hard-linked into
  `pvoid\Scripts\ebi.exe`, so activating the venv puts it on PATH. Don't set
  environment variables or override `EBI_EXECUTABLE` per command.
- **toothpaste**: `lab/toothpaste_bridge.py` expects it at
  `C:\working\tools\toothpaste`.
- **Logs**: small fixtures are in `data/`; large logs (rtfm, sepsis, BPI) live
  outside the repo in `C:/working/data`.
- **Dependencies**: declared in `pyproject.toml`. After changing it, reinstall
  with `bash run.sh pip install -e .`. `requirements.txt` is a `pip freeze`
  record of the environment - regenerate it, don't hand-edit it. Don't add
  new dependencies without discussing it with the user.

## Worktrees

- Each worktree needs its own venv, or imports silently resolve to the main
  checkout:

  ```
  C:\working\tools\python\python314\python.exe -m venv pvoid
  bash run.sh pip install -e .
  ```

  then hard-link ebi into it (PowerShell):

  ```
  New-Item -ItemType HardLink -Path pvoid\Scripts\ebi.exe -Target C:\working\tools\ebi\ebi.exe
  ```

- Check imports resolve inside the worktree before trusting any result:
  `bash run.sh python -c "import process_voids, lab; print(process_voids.__file__, lab.__file__)"`
- Commit freely on your worktree branch. The user lands work on `main` with
  `git merge --squash`. Never run `git add`/`commit`/`push` on `main` - print
  the commands for the user instead.
- Work is exchanged between sessions as branches (and findings as messages to
  the user), never as patch files applied in another session's tree.
- Don't sign your commits with the Claude Code guff

## Tests

- Full suite, exactly like this with nothing chained after it (this shape is
  allowlisted; redirection, `tail` or `echo` trigger a permission prompt):

  ```
  bash run.sh python -m unittest discover -s tests -t .
  ```

  `-t .` is required. The output is long - tests of error handling log their
  expected tracebacks - so read the final `Ran N tests` / `OK` lines.
- One module: `bash run.sh python -m unittest tests.lab.test_exp_disco_degrade -v`
- Test first. For a bugfix, write the regression test and see it fail before
  writing the fix.
- Run the full suite before proposing anything to land.

## Experiments

- Current entry point: `bash run.sh python -m lab.exp_disco_degrade --run smoke`
  (`--dry-run` prints the resolved configuration). This is being consolidated
  into a single runner; don't add new `exp_*` entry points.
- Each run logs to `var/lab/logs/<experiment>_<timestamp>.log`, opening with
  the process-voids and skip-alignments versions; results go to
  `var/lab/results/`, timestamped. `var/` is gitignored.
- Probe one cell on a new log before committing to a full sweep.

## Code conventions

- `process_voids` is the product package; `process_voids.pvoid` is the product
  entry point - don't modify it without asking. `lab` is the research harness.
- The metric registry (`lab/metric_registry.py`) is append-only: it's the
  dictionary for every result file ever written. Retire ids, never delete them.
- Comments and docstrings describe the code as it is. Change history and
  rationale go in `CHANGELOG.md` or commit messages.
- This repo is public: don't reference files outside source control in code,
  comments or commit messages.
- Don't modify sibling repos (e.g. skip-alignments); report findings to the
  user.
- If an implementation would narrow or change an agreed design, raise it
  with the user - don't document the deviation as though it were the design.
- Australian/British spelling.

## Releases

- `CHANGELOG.md` follows Keep a Changelog; changes accumulate under
  `## Unreleased`.
- `bash run.sh python -m process_voids.util.release_check` checks version and
  CHANGELOG agreement, dependency pins, the full suite and untracked files.
- The user creates releases and tags in the GitHub web interface.
