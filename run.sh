#!/bin/bash
# Activates the project venv, then execs whatever command follows.
# Exists purely so commands needing the venv (ebi.exe on PATH, installed
# packages) don't need a "source ... && ..." compound command, which
# permission allowlisting can't reliably match on.
#
# Usage: bash run.sh [--seed N] <command> [args...]
#
# --seed N sets PYTHONHASHSEED before invoking the command. It has to
# live here rather than as a Python CLI flag: the interpreter reads
# PYTHONHASHSEED at startup, so by the time any main() could parse a
# flag, its own hash seed is already fixed. Pinning it makes pm4py's
# Inductive cut selection reproducible where a noise_threshold cut sits
# near a tie. Recognised in first position only, so it can't swallow an
# argument meant for the command itself.
source pvoid/Scripts/activate
if [ "$1" = "--seed" ]; then
    export PYTHONHASHSEED="$2"
    shift 2
fi
exec "$@"
