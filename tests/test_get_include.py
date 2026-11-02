"""Smoke test of :func:`composites.get_include`, which must run without error
"""
import sys
sys.path.append('..')

import composites

def test_get_include():
    composites.get_include()
