'''
lab.print_tree is a thin wrapper around lab.discovery.discover_cached -
the tree itself already has a readable str() (skip-alignments' own
process-tree notation), so this only needs to test that it resolves a
log/combo name correctly and hands back that string, not the notation
itself.
'''

import unittest
from pathlib import Path
from unittest.mock import patch, MagicMock

from lab.discovery import CachedDiscovery
from lab.print_tree import discovered_tree_text


def _found(tree):
    """What discover_cached really returns. A bare (tree, weights) tuple
    here would let a caller that unpacks it as a pair pass, while the
    real call raises."""
    return CachedDiscovery(tree, None, 'cached', Path('x__y__pair.pkl'))


class DiscoveredTreeTextTest(unittest.TestCase):
    def test_resolves_a_registered_log_name_and_returns_the_tree_str(self):
        fake_tree = MagicMock()
        fake_tree.__str__.return_value = 'FAKE TREE TEXT'
        fake_log = object()

        with patch('lab.print_tree.pm4py.read_xes', return_value=fake_log) as mock_read, \
             patch('lab.print_tree.discover_cached',
                   return_value=_found(fake_tree)) as mock_discover, \
             patch('lab.print_tree.ALL_LOGS', {'payment_approval': 'data/payment_approval.xes'}), \
             patch('lab.print_tree.ALL_COMBOS', {'inductive_noise20': 'fake_combo'}):
            text = discovered_tree_text('payment_approval', 'inductive_noise20')

        mock_read.assert_called_once_with('data/payment_approval.xes')
        mock_discover.assert_called_once_with(
            'payment_approval', 'inductive_noise20', 'fake_combo', fake_log)
        self.assertEqual(text, 'FAKE TREE TEXT')

    def test_an_unregistered_log_is_treated_as_a_path_directly(self):
        fake_tree = MagicMock()
        fake_tree.__str__.return_value = 'FAKE TREE TEXT'

        with patch('lab.print_tree.pm4py.read_xes', return_value=object()) as mock_read, \
             patch('lab.print_tree.discover_cached',
                   return_value=_found(fake_tree)) as mock_discover, \
             patch('lab.print_tree.ALL_LOGS', {}), \
             patch('lab.print_tree.ALL_COMBOS', {'inductive_noise20': 'fake_combo'}):
            discovered_tree_text('some/ad_hoc.xes', 'inductive_noise20')

        mock_read.assert_called_once_with('some/ad_hoc.xes')
        # log_name passed to discover_cached (the cache key) is the
        # path's stem, not the full ad hoc path, matching
        # exp_disco_degrade's own convention
        self.assertEqual(mock_discover.call_args.args[0], 'ad_hoc')


if __name__ == '__main__':
    unittest.main()
