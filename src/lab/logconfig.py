"""
Minimal console logging setup shared by the experiment entrypoints.

Times/status go through logging (INFO by default, DEBUG for the more
chatty per-cell start markers) rather than print(), so verbosity is a
level away instead of a code change.
"""

import logging

from skipalignments.progress import disable_progress_bars


def configure(level=logging.INFO):
    logging.basicConfig(level=level, format='%(asctime)s %(levelname)s %(message)s',
                         datefmt='%H:%M:%S')
    disable_progress_bars()


def enable_skipalignments_debug():
    """
    skip-alignments' execution.py logs waste-ratio/per-variant-timing
    detail at DEBUG under its own logger, silent by default. Call this
    to opt in when investigating a slow run.
    """
    logging.getLogger('skipalignments.execution').setLevel(logging.DEBUG)
