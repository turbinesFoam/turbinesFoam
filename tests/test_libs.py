#!/usr/bin/env python
"""Tests for shared libraries."""

import ctypes
import platform


def test_turbinesfoam():
    """Test that `libturbinesFoam.so` can be loaded."""
    if platform.system() == "Darwin":
        ext = ".dylib"
    else:
        ext = ".so"
    ctypes.cdll.LoadLibrary("libturbinesFoam" + ext)
