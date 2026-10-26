"""Unit tests of the shared modules in bin/ (run: pytest tests/). Synthetic data only."""
import os
import sys

BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, os.path.abspath(BIN))
ASSETS = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'assets'))


def getter(d):
    """get_field callable over a dict (missing field = '')."""
    return lambda k: d.get(k, '')
