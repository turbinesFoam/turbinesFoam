#!/usr/bin/env python
"""Tests for `actuatorLineSource`."""

from __future__ import division, print_function
import subprocess
import pandas as pd
import os
import numpy as np


alpha_deg = 8.0
output_fpath = "postProcessing/actuatorLines/0/foil.csv"


def setup():
    os.chdir("actuatorLineSource")
    output_clean = subprocess.check_output("./Allclean")
    output_run = subprocess.check_output(["./Allrun", str(alpha_deg)])


def load_output():
    """Load `actuatorLineSource` output file from CSV."""
    return pd.read_csv(output_fpath)


def test_created():
    """Test that actuatorLineSource was created."""
    txt = "Selecting finite volume options model type actuatorLineSource"
    subprocess.check_output(["grep", txt, "log.simpleFoam"])


def test_output_file_exists():
    """Test that actuatorLineSource output file exists."""
    assert os.path.isfile(output_fpath)


def test_geometric_alpha():
    """
    Test that actuatorLineSource geometric angle of attack was set properly.
    """
    df = load_output()
    assert np.all(df.alpha_geom_deg == alpha_deg)


def test_alpha_sweep():
    """Test 2-D actuator line angle of attack sweep."""
    out = subprocess.check_output(["python", "paramsweep.py", "-15", "16", "5"])
    df = pd.read_csv("processed/alpha_sweep.csv")
    mse = np.mean((df.alpha_geom_deg - df.alpha_deg)**2)
    print("Mean square error between geometric and detected alpha (deg):", mse)
    assert mse < 0.1


def teardown():
    """Move back into tests directory."""
    os.chdir("../")
