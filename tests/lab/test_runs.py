import unittest

from lab.runs import Experiment, RUNS


class ExperimentTest(unittest.TestCase):
    def test_cell_count_is_the_cartesian_product(self):
        experiment = Experiment(name='x', log_paths=['a.xes', 'b.xes'],
                                 combos={'c1': None, 'c2': None},
                                 degradations={'d1': None}, levels=[0.0, 0.5, 1.0])
        self.assertEqual(experiment.cell_count(), 2 * 2 * 1 * 3)

    def test_describe_lists_log_stems_not_full_paths(self):
        experiment = Experiment(name='x', log_paths=['data/rtfm.xes'],
                                 combos={'inductive': None}, degradations={'trace': None},
                                 levels=[0.0])
        text = experiment.describe()
        self.assertIn('rtfm', text)
        self.assertNotIn('data/rtfm.xes', text)

    def test_describe_flags_the_level_zero_overcount(self):
        experiment = Experiment(name='x', log_paths=['a.xes'], combos={'c': None},
                                 degradations={'d1': None, 'd2': None}, levels=[0.0, 0.5])
        self.assertIn('overcounts', experiment.describe())

    def test_describe_does_not_flag_overcount_with_one_degradation_dim(self):
        experiment = Experiment(name='x', log_paths=['a.xes'], combos={'c': None},
                                 degradations={'d1': None}, levels=[0.0, 0.5])
        self.assertNotIn('overcounts', experiment.describe())

    def test_describe_does_not_flag_overcount_without_level_zero(self):
        experiment = Experiment(name='x', log_paths=['a.xes'], combos={'c': None},
                                 degradations={'d1': None, 'd2': None}, levels=[0.5])
        self.assertNotIn('overcounts', experiment.describe())


class RunsRegistryTest(unittest.TestCase):
    def test_every_run_is_an_experiment_instance(self):
        for name, run in RUNS.items():
            self.assertIsInstance(run, Experiment, f'{name} is not an Experiment')

    def test_every_run_has_at_least_one_log_combo_degradation_and_level(self):
        for name, run in RUNS.items():
            self.assertTrue(run.log_paths, f'{name} has no logs')
            self.assertTrue(run.combos, f'{name} has no combos')
            self.assertTrue(run.degradations, f'{name} has no degradations')
            self.assertTrue(run.levels, f'{name} has no levels')


if __name__ == '__main__':
    unittest.main()
