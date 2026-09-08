"""
Console + file logging setup shared by the experiment entrypoints.

Times/status go through logging (INFO by default, DEBUG for the more
chatty per-cell start markers) rather than print(), so verbosity is a
level away instead of a code change.
"""

import logging
import sys
import time
from pathlib import Path

from skipalignments.progress import disable_progress_bars

LOG_DIR = Path('var/lab/logs')


def configure(level=logging.INFO, log_dir=LOG_DIR):
    """
    Logs to console AND a timestamped file under log_dir - relying on
    console output alone means the only record of a run is wherever its
    stdout happened to be captured, which is exactly what's fragile
    (a `| tail` swallowing everything until exit, a background task's
    capture getting discarded once it's stopped or superseded - see
    session notes on the rtfm run that left nothing behind). The file
    handler makes a run's log durable and independent of how the
    command was invoked, and the timestamp in its name (script name +
    YYYYMMDD-HHMMSS) means a later run can never silently overwrite an
    earlier one's record the way a fixed path could.
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
    logging.getLogger(__name__).info('Logging to %s', log_path)


def enable_skipalignments_debug():
    """
    skip-alignments' execution.py logs waste-ratio/per-variant-timing
    detail at DEBUG under its own logger, silent by default. Call this
    to opt in when investigating a slow run.
    """
    logging.getLogger('skipalignments.execution').setLevel(logging.DEBUG)
