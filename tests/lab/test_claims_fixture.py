'''
Tests for lab.claims_fixture - the synthetic claims-handling fixture
(see claims-fixture-spec.md): a mid-sized (~30 trace) log/tree pair,
bigger than the six-trace payment running example but fast enough for
iteration, with a non-zero baseline deficit and a known-optional
subprocess for a specificity check.
'''

import unittest

from lab.claims_fixture import (
    ACTIVITY_ALPHABET, build_claims_tree, generate_claims_log, summarize_ground_truth,
)


class ClaimsTreeStructureTest(unittest.TestCase):
    '''The tree matches the spec's model exactly - every leaf reachable,
    the right activities mandatory vs optional.'''

    def setUp(self):
        self.tree = build_claims_tree()

    def test_leaf_alphabet_matches_spec(self):
        self.assertEqual(set(self.tree.get_leaf_labels()), ACTIVITY_ALPHABET)

    def test_root_is_a_sequence_of_four(self):
        self.assertEqual(len(self.tree.children), 4)

    def test_register_and_close_are_mandatory_leaves(self):
        labels = [c.name for c in self.tree.children if not c.children]
        self.assertIn('register', labels)
        self.assertIn('close', labels)

    def test_appeal_subprocess_is_optional_via_xor_with_tau(self):
        # the Xor(appeal_seq, tau) node - one child has no 'children'
        # attribute issue; find the Xor by its two children, one of
        # which is a Tau leaf.
        from skipalignments.processtree import Xor, Tau
        xor_nodes = [c for c in self.tree.children if isinstance(c, Xor)]
        self.assertEqual(len(xor_nodes), 1)
        xor = xor_nodes[0]
        self.assertTrue(any(isinstance(c, Tau) for c in xor.children))

    def test_loop_block_is_mandatory_and_parallel_with_assess(self):
        from skipalignments.processtree import And, Loop
        and_nodes = [c for c in self.tree.children if isinstance(c, And)]
        self.assertEqual(len(and_nodes), 1)
        parallel = and_nodes[0]
        loop_nodes = [c for c in parallel.children if isinstance(c, Loop)]
        self.assertEqual(len(loop_nodes), 1)
        self.assertEqual(set(loop_nodes[0].get_leaf_labels()),
                          {'request_docs', 'receive_docs'})


class ClaimsLogGenerationTest(unittest.TestCase):
    '''Structural properties of the generated log/ground-truth pair -
    not pinning exact sampled values (that would be fragile and not
    the point), but the invariants the spec cares about.'''

    def setUp(self):
        self.log, self.ground_truth = generate_claims_log(n_traces=30, seed=42)

    def test_trace_count(self):
        self.assertEqual(self.log['case:concept:name'].nunique(), 30)
        self.assertEqual(len(self.ground_truth), 30)

    def test_every_trace_starts_with_register_and_ends_with_close(self):
        for case, group in self.log.groupby('case:concept:name', sort=False):
            ordered = group.sort_values('time:timestamp')
            self.assertEqual(ordered.iloc[0]['concept:name'], 'register')
            self.assertEqual(ordered.iloc[-1]['concept:name'], 'close')

    def test_only_alphabet_and_injected_extra_activities_appear(self):
        observed = set(self.log['concept:name'].unique())
        deviated_extra = set(self.ground_truth.loc[
            self.ground_truth['deviation'].str.startswith('extra_activity', na=False),
            'deviation'
        ].str.split(':').str[1])
        self.assertTrue(observed.issubset(ACTIVITY_ALPHABET | deviated_extra))

    def test_loop_iteration_counts_are_consistent_with_request_docs_events(self):
        for case, group in self.log.groupby('case:concept:name', sort=False):
            n_iterations = self.ground_truth.loc[
                self.ground_truth['case'] == case, 'loop_iterations'].iloc[0]
            n_request_docs = (group['concept:name'] == 'request_docs').sum()
            self.assertEqual(n_request_docs, n_iterations)

    def test_appeal_flag_matches_presence_of_appeal_activities(self):
        for case, group in self.log.groupby('case:concept:name', sort=False):
            appeal = self.ground_truth.loc[
                self.ground_truth['case'] == case, 'appeal'].iloc[0]
            has_appeal_activity = group['concept:name'].isin(
                {'lodge_appeal', 'decide_appeal'}).any()
            self.assertEqual(bool(appeal), bool(has_appeal_activity))

    def test_a_few_traces_are_deviated(self):
        n_deviated = (self.ground_truth['deviation'] != '').sum()
        self.assertIn(n_deviated, (2, 3))

    def test_deterministic_given_seed(self):
        log2, gt2 = generate_claims_log(n_traces=30, seed=42)
        self.assertTrue(self.log.equals(log2))
        self.assertTrue(self.ground_truth.equals(gt2))

    def test_reordered_deviation_deliberately_swaps_the_appeal_pair(self):
        # Not "some random adjacent pair, wherever the seeded rng
        # stream happens to land" - that was fragile (only became a
        # real voidmass-visible deviation because the swap incidentally
        # hit lodge_appeal/decide_appeal, a strict Sequence; if it had
        # landed in the order-tolerant parallel/loop region instead, it
        # would have cost nothing). The target must be deliberate: the
        # appeal pair is the one place in the model where order is
        # actually enforced, so that's what gets swapped, every time -
        # checked across several seeds, since seed=42 alone could pass
        # by the same kind of accident this is meant to rule out.
        for seed in (42, 1, 7, 123, 9999):
            with self.subTest(seed=seed):
                log, ground_truth = generate_claims_log(n_traces=30, seed=seed)
                reordered = ground_truth[ground_truth['deviation'].str.startswith(
                    'reordered', na=False)]
                self.assertEqual(len(reordered), 1)
                case = reordered.iloc[0]['case']
                self.assertTrue(reordered.iloc[0]['appeal'])

                group = log[log['case:concept:name'] == case].sort_values('time:timestamp')
                ordered_activities = list(group['concept:name'])
                i_decide = ordered_activities.index('decide_appeal')
                i_lodge = ordered_activities.index('lodge_appeal')
                self.assertLess(i_decide, i_lodge)  # swapped: decide now precedes lodge

    def test_summary_reports_deviated_case_and_counts(self):
        summary = summarize_ground_truth(self.ground_truth)
        n_appeal = int(self.ground_truth['appeal'].sum())
        deviated_case = self.ground_truth.loc[
            self.ground_truth['deviation'] != '', 'case'].iloc[0]
        self.assertIn(f'n_traces: {len(self.ground_truth)}', summary)
        self.assertIn(f'appeal: {n_appeal} of {len(self.ground_truth)}', summary)
        self.assertIn(deviated_case, summary)


if __name__ == '__main__':
    unittest.main()
