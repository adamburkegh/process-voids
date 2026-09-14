"""
run.sh's --seed flag.

PYTHONHASHSEED is read by the interpreter at startup, so nothing inside
a Python main() can set it for its own process. The flag therefore lives
in the wrapper script that already exists to prepare the environment
before invoking python, alongside activating the venv.

Driven through a real `bash run.sh`, since what is being tested is the
script's own argument handling - a Python-level reimplementation of it
would pass while the script itself was broken.
"""

import subprocess
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Prints the seed this interpreter actually started under, then its own
# argv, so one invocation can check both the environment and that --seed
# was consumed rather than forwarded.
PROBE = ('import os, sys; '
         "print(os.environ.get('PYTHONHASHSEED', 'unset')); "
         "print(' '.join(sys.argv[1:]))")


def _run(*args):
    result = subprocess.run(['bash', 'run.sh', *args], cwd=REPO_ROOT,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    # no .strip() on the whole output first: an empty argv line is a real
    # result (the command got no arguments) and strip would swallow it
    seed, argv = result.stdout.splitlines()[:2]
    return seed.strip(), argv.strip()


class SeedFlagTest(unittest.TestCase):

    def test_seed_reaches_the_interpreter_that_runs_the_command(self):
        seed, _argv = _run('--seed', '0', 'python', '-c', PROBE)
        self.assertEqual(seed, '0')

    def test_any_seed_value_is_passed_through_not_just_zero(self):
        seed, _argv = _run('--seed', '12345', 'python', '-c', PROBE)
        self.assertEqual(seed, '12345')

    def test_the_flag_is_consumed_not_forwarded_to_the_command(self):
        """The command must see its own arguments only - a stray
        '--seed 0' reaching an argparse CLI would be a hard error."""
        _seed, argv = _run('--seed', '0', 'python', '-c', PROBE, 'a', 'b')
        self.assertEqual(argv, 'a b')

    def test_without_the_flag_the_seed_is_left_alone(self):
        """Unset is the honest default - run_history records it as
        'unset', and pretending to a seed nobody chose would be worse
        than saying so."""
        seed, _argv = _run('python', '-c', PROBE)
        self.assertEqual(seed, 'unset')

    def test_a_command_named_like_the_flag_is_still_runnable(self):
        """--seed is only recognised in first position, so it can't
        swallow an argument meant for the command."""
        _seed, argv = _run('python', '-c', PROBE, '--seed', '0')
        self.assertEqual(argv, '--seed 0')


if __name__ == '__main__':
    unittest.main()
