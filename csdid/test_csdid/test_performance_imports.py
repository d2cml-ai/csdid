"""Importing estimators must not eagerly request pyplot from csdid modules."""
import os
from pathlib import Path
import subprocess
import sys


def test_pyplot_import_is_deferred_and_plot_helpers_remain_available():
    # A fresh process avoids pytest's own plotting imports. Dependencies may
    # import matplotlib themselves, so check imports initiated by csdid only.
    source = """
import builtins
original_import = builtins.__import__
def check_import(name, globals=None, locals=None, fromlist=(), level=0):
    caller = (globals or {}).get('__name__', '')
    if caller.startswith('csdid') and name == 'matplotlib.pyplot':
        raise AssertionError('csdid imported pyplot before a plotting call')
    return original_import(name, globals, locals, fromlist, level)
builtins.__import__ = check_import
from csdid.att_gt import ATTgt, gplot, splot
assert callable(gplot) and callable(splot)
assert callable(ATTgt.plot_attgt) and callable(ATTgt.plot_aggte)
"""
    repo = Path(__file__).resolve().parents[2]
    env = dict(os.environ, PYTHONPATH=str(repo))
    result = subprocess.run([sys.executable, "-c", source], cwd=str(repo),
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
