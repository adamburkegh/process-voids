"""
Tests for lab.discovery's shared discovered-tree cache.

The discovered tree depends only on (log, combo), so every experiment
script caches it to var/lab/tree_cache/ rather than rediscovering. That
also makes the tree reproducible across processes, which matters because
pm4py's Inductive cut selection is hash-seed dependent where a
noise_threshold cut sits near a tie - rtfm draws a 12- or 13-node tree
from the same log depending on the process.

The cache holds a (tree, ppt_weights) pair: only exp_disco_degrade
consumes ppt_weights (toothpaste's fixed PPT weights, threaded into
pvoid.skipprob), but all three scripts share these files, so the shape
has to serve the richest caller.
"""

import pickle
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from lab.discovery import DiscoveryCombo, DiscoveryResult, discover_cached


class DiscoverCachedTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.cache_dir = Path(self.tmp.name)
        patcher = patch('lab.discovery.TREE_CACHE_DIR', self.cache_dir)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.discover = Mock(return_value=DiscoveryResult('TREE', 'WEIGHTS'))
        self.combo = DiscoveryCombo('combo', self.discover)

    def test_miss_discovers_and_returns_tree_and_weights(self):
        self.assertEqual(discover_cached('log', 'combo', self.combo, 'LOG'),
                         ('TREE', 'WEIGHTS'))
        self.discover.assert_called_once_with('LOG')

    def test_hit_returns_the_cached_pair_without_rediscovering(self):
        first = discover_cached('log', 'combo', self.combo, 'LOG')
        self.discover.return_value = DiscoveryResult('OTHER_TREE', 'OTHER_WEIGHTS')
        second = discover_cached('log', 'combo', self.combo, 'LOG')
        self.assertEqual(second, first)
        self.discover.assert_called_once()

    def test_cache_is_keyed_by_log_and_combo(self):
        discover_cached('log_a', 'combo', self.combo, 'LOG')
        self.discover.return_value = DiscoveryResult('TREE_B', 'WEIGHTS_B')
        self.assertEqual(discover_cached('log_b', 'combo', self.combo, 'LOG'),
                         ('TREE_B', 'WEIGHTS_B'))
        self.assertEqual(discover_cached('log_a', 'combo', self.combo, 'LOG'),
                         ('TREE', 'WEIGHTS'))

    def test_combo_without_ppt_weights_caches_none(self):
        self.discover.return_value = DiscoveryResult('TREE')
        self.assertEqual(discover_cached('log', 'combo', self.combo, 'LOG'),
                         ('TREE', None))

    def test_a_bare_tree_cache_file_is_not_read_as_a_pair(self):
        """
        Files written before the cache held pairs contain a bare tree.
        Unpickling one as a (tree, ppt_weights) pair would silently
        misread it - a string tree would unpack into its first two
        characters - so the pair cache must not read them.
        """
        legacy = self.cache_dir / 'log__combo.pkl'
        with open(legacy, 'wb') as f:
            pickle.dump('LEGACY_TREE', f)
        self.assertEqual(discover_cached('log', 'combo', self.combo, 'LOG'),
                         ('TREE', 'WEIGHTS'))

    def test_discovery_failure_propagates_and_caches_nothing(self):
        """
        exp_disco_degrade turns NotImplementedError into a
        'not_implemented' discover_status at its own call site, so the
        helper must not swallow it - nor leave a file behind that a
        later run would read as a successful discovery.
        """
        self.discover.side_effect = NotImplementedError('no such miner')
        with self.assertRaises(NotImplementedError):
            discover_cached('log', 'combo', self.combo, 'LOG')
        self.assertEqual(list(self.cache_dir.glob('*.pkl')), [])


if __name__ == '__main__':
    unittest.main()
