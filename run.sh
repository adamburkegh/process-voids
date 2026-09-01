#!/bin/bash
# Activates the project venv, then execs whatever command follows.
# Exists purely so commands needing the venv (ebi.exe on PATH, installed
# packages) don't need a "source ... && ..." compound command, which
# permission allowlisting can't reliably match on.
#
# Usage: bash run.sh <command> [args...]
source pvoid/Scripts/activate
exec "$@"
