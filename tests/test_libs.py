#!/usr/bin/env python
"""Tests for shared libraries."""

import ctypes


def test_turbinesfoam():
    """Test that `libturbinesFoam.so` can be loaded."""
    ctypes.cdll.LoadLibrary("libturbinesFoam.so")


def test_turbulence_models():
    """Test that `libturbinesFoamTurbulenceModels.so` can be loaded."""
    ctypes.cdll.LoadLibrary("libturbinesFoamTurbulenceModels.so")
