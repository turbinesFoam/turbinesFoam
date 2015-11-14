#!/usr/bin/env python
"""Tests for `actuatorLineSource`."""

from __future__ import division, print_function
import subprocess
import pandas as pd
import os
import numpy as np
from nose.tools import timed


alpha_deg = 10.0
output_fpath = "postProcessing/actuatorLines/0/foil.csv"


def setup():
    os.chdir("actuatorLineSource")


def load_output():
    """Load `actuatorLineSource` output file from CSV."""
    return pd.read_csv(output_fpath)


def check_created():
    """Test that actuatorLineSource was created."""
    txt = "Selecting finite volume options model type actuatorLineSource"
    subprocess.check_output(["grep", txt, "log.simpleFoam"])


def check_output_file_exists():
    """Test that actuatorLineSource output file exists."""
    assert os.path.isfile(output_fpath)


def check_geometric_alpha():
    """
    Test that actuatorLineSource geometric angle of attack was set properly.
    """
    df = load_output()
    assert np.all(df.alpha_geom_deg == alpha_deg)


def test_2d():
    """Test 2-D actuatorLineSource."""
    output_clean = subprocess.check_output("./Allclean")
    output_run = subprocess.check_output(["./Allrun", str(alpha_deg)])
    check_created()
    check_output_file_exists()
    check_geometric_alpha()


def test_alpha_sweep():
    """Test 2-D actuatorLineSource angle of attack sweep."""
    out = subprocess.check_output(["python", "paramsweep.py", "-10", "11", "5"])
    df = pd.read_csv("processed/alpha_sweep.csv")
    mse = np.mean((df.alpha_geom_deg - df.alpha_deg)**2)
    print("Mean square error between geometric and detected alpha (deg):", mse)
    assert mse < 0.1


def test_3d():
    """Test 3-D actuatorLineSource."""
    out = subprocess.check_output("./Allclean")
    out = subprocess.check_output(["./Allrun3D", str(alpha_deg)])
    log_end = subprocess.check_output(["tail", "log.simpleFoam"]).decode()
    print(log_end)
    assert log_end.split()[-1] == "End"


@timed(30) # Test must run faster than 30 seconds
def test_parallel():
    """Test 3-D actuatorLineSource in parallel."""
    out = subprocess.check_output("./Allclean")
    out = subprocess.check_output(["./Allrun3D", "-parallel", str(alpha_deg)])
    log_end = subprocess.check_output(["tail", "log.simpleFoam"]).decode()
    print(log_end)
    assert "Finalising parallel run" in log_end


def teardown():
    """Move back into tests directory."""
    os.chdir("../")
