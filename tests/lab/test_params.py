"""
lab.params' external logs: registered by filename, with the directory
coming from pvoid.toml.

lab.runs builds its named runs from ALL_LOGS at import, and lab.params
is imported almost everywhere, so an external log must not need its
directory until something actually opens or names it - otherwise a
machine without pvoid.toml could not import the harness or run the
suite at all.
"""

import os
import unittest
from pathlib import Path
from unittest.mock import patch

from lab.params import ALL_LOGS, ExternalLog
from process_voids.config import ConfigError

CONFIG = {'paths': {'data_dir': '/srv/logs'}}


def _with_config(config):
    return patch('process_voids.config._default_config', return_value=config)


class ExternalLogTest(unittest.TestCase):

    def test_constructing_one_needs_no_configuration(self):
        with _with_config({}):
            ExternalLog('rtfm.xes')

    def test_resolves_under_the_configured_data_dir(self):
        with _with_config(CONFIG):
            self.assertEqual(os.fspath(ExternalLog('rtfm.xes')), '/srv/logs/rtfm.xes')

    def test_resolves_with_forward_slashes(self):
        """So a path recorded from it stays legible off Windows, as
        run_history's paths already are."""
        with _with_config({'paths': {'data_dir': r'shared\logs'}}):
            self.assertNotIn('\\', os.fspath(ExternalLog('rtfm.xes')))

    def test_names_the_log_by_its_stem_like_a_plain_path(self):
        with _with_config(CONFIG):
            self.assertEqual(Path(ExternalLog('rtfm.xes')).stem, 'rtfm')

    def test_an_unset_data_dir_fails_when_the_log_is_used(self):
        with _with_config({}):
            with self.assertRaises(ConfigError) as caught:
                os.fspath(ExternalLog('rtfm.xes'))
        self.assertIn('data_dir', str(caught.exception))

    def test_repr_does_not_need_configuration(self):
        """Listing registered logs must not fail on a machine without
        pvoid.toml."""
        with _with_config({}):
            self.assertIn('rtfm.xes', repr(ExternalLog('rtfm.xes')))


class AllLogsTest(unittest.TestCase):

    def test_repository_fixtures_stay_plain_relative_paths(self):
        self.assertEqual(ALL_LOGS['payment_approval'], 'data/payment_approval.xes')

    def test_no_registered_log_hardcodes_an_absolute_path(self):
        for name, path in ALL_LOGS.items():
            with self.subTest(log=name):
                if isinstance(path, ExternalLog):
                    continue
                self.assertFalse(Path(path).is_absolute(), path)
                self.assertNotRegex(path, r'^[A-Za-z]:')


if __name__ == '__main__':
    unittest.main()
