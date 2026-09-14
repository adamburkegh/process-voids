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

from lab.discovery import (
    CachedDiscovery, DiscoveryCombo, DiscoveryResult, discover_cached)


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
        found = discover_cached('log', 'combo', self.combo, 'LOG')
        self.assertEqual((found.tree, found.ppt_weights), ('TREE', 'WEIGHTS'))
        self.discover.assert_called_once_with('LOG')

    def test_hit_returns_the_cached_pair_without_rediscovering(self):
        first = discover_cached('log', 'combo', self.combo, 'LOG')
        self.discover.return_value = DiscoveryResult('OTHER_TREE', 'OTHER_WEIGHTS')
        second = discover_cached('log', 'combo', self.combo, 'LOG')
        self.assertEqual((second.tree, second.ppt_weights),
                         (first.tree, first.ppt_weights))
        self.discover.assert_called_once()

    def test_cache_is_keyed_by_log_and_combo(self):
        discover_cached('log_a', 'combo', self.combo, 'LOG')
        self.discover.return_value = DiscoveryResult('TREE_B', 'WEIGHTS_B')
        b = discover_cached('log_b', 'combo', self.combo, 'LOG')
        self.assertEqual((b.tree, b.ppt_weights), ('TREE_B', 'WEIGHTS_B'))
        a = discover_cached('log_a', 'combo', self.combo, 'LOG')
        self.assertEqual((a.tree, a.ppt_weights), ('TREE', 'WEIGHTS'))

    def test_combo_without_ppt_weights_caches_none(self):
        self.discover.return_value = DiscoveryResult('TREE')
        found = discover_cached('log', 'combo', self.combo, 'LOG')
        self.assertEqual((found.tree, found.ppt_weights), ('TREE', None))

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
        found = discover_cached('log', 'combo', self.combo, 'LOG')
        self.assertEqual((found.tree, found.ppt_weights), ('TREE', 'WEIGHTS'))

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

    def test_a_miss_reports_it_discovered_now_and_which_file_it_wrote(self):
        """Which tree a result was scored against is the point of the
        cache, so the caller is told rather than left to infer it from
        whether the file happened to exist beforehand."""
        found = discover_cached('log', 'combo', self.combo, 'LOG')
        self.assertEqual(found.source, 'discovered')
        self.assertEqual(found.cache_path, self.cache_dir / 'log__combo__pair.pkl')

    def test_a_hit_reports_it_came_from_the_cache(self):
        discover_cached('log', 'combo', self.combo, 'LOG')
        found = discover_cached('log', 'combo', self.combo, 'LOG')
        self.assertEqual(found.source, 'cached')
        self.assertEqual(found.cache_path, self.cache_dir / 'log__combo__pair.pkl')

    def test_a_legacy_bare_tree_file_counts_as_a_miss(self):
        """The legacy file is ignored, so this run really did discover
        the tree - reporting 'cached' would name a file it did not read."""
        with open(self.cache_dir / 'log__combo.pkl', 'wb') as f:
            pickle.dump('LEGACY_TREE', f)
        self.assertEqual(discover_cached('log', 'combo', self.combo, 'LOG').source,
                         'discovered')


if __name__ == '__main__':
    unittest.main()
