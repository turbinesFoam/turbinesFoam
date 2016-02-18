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
element_dir = "postProcessing/actuatorLineElements/0/"


def setup():
    os.chdir("actuatorLineSource")


def get_tutorial_files(case="static"):
    cmd = "./getTutorialFiles.sh {}".format(case)
    out = subprocess.check_output(cmd, shell=True)


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


def check_element_file_exists():
    """Test that the element perf file was created."""
    assert os.path.isfile(os.path.join(element_dir, "foil.element0.csv"))


def check_geometric_alpha():
    """Test that actuatorLineSource geometric angle of attack was set properly.
    """
    df = load_output()
    assert np.all(df.alpha_geom_deg == alpha_deg)


def check_re_corrections():
    """Parse and check data from Reynolds number corrections."""
    cmd_temp = "grep '{} {} {} coefficient' log.simpleFoam"
    def make_array(when=None, minmax=None, quantity=None, cmd=None):
        if cmd is None:
            cmd = cmd_temp.format(when, minmax, quantity)
        raw = subprocess.check_output(cmd, shell=True).decode().split("\n")
        arr = []
        for i in raw:
            try:
                arr.append(float(i.split()[-1]))
            except IndexError:
                pass
        return np.array(arr)
    df = pd.DataFrame()
    for when in ["Initial", "Corrected"]:
        for minmax, quantity in zip(["maximum","minimum"], ["lift", "drag"]):
            name = when.lower() + "_" + minmax + "_" + quantity
            df[name] = make_array(when, minmax, quantity)
    Re = make_array(cmd="grep Re: log.simpleFoam")[-1]
    Re_ref = make_array(cmd="grep ReRef: log.simpleFoam")[-1]
    cdi = df.initial_minimum_drag.iloc[-1]
    cdc = df.corrected_minimum_drag.iloc[-1]
    cli = df.initial_maximum_lift.iloc[-1]
    clc = df.corrected_maximum_lift.iloc[-1]
    if Re > Re_ref:
        assert clc > cli
        assert cdc < cdi
    elif Re < Re_ref:
        assert clc < cli
        assert cdc > cdi
    elif Re == Re_ref:
        assert clc == cli
        assert cdc == cdi


def test_2d():
    """Test 2-D actuatorLineSource."""
    get_tutorial_files(case="static")
    output_clean = subprocess.check_output("./Allclean")
    output_run = subprocess.check_output(["./Allrun", "2D", str(alpha_deg)])
    check_created()
    check_output_file_exists()
    check_element_file_exists()
    check_geometric_alpha()
    check_re_corrections()


def test_alpha_sweep():
    """Test 2-D actuatorLineSource angle of attack sweep."""
    get_tutorial_files(case="static")
    out = subprocess.check_output(["python", "paramsweep.py", "-10", "11", "5"])
    df = pd.read_csv("processed/alpha_sweep.csv")
    mse = np.mean((df.alpha_geom_deg - df.alpha_deg)**2)
    print("Mean square error between geometric and detected alpha (deg):", mse)
    assert mse < 0.1


def test_3d():
    """Test 3-D actuatorLineSource."""
    get_tutorial_files(case="static")
    out = subprocess.check_output("./Allclean")
    out = subprocess.check_output(["./Allrun", "3D", str(alpha_deg)])
    log_end = subprocess.check_output(["tail", "log.simpleFoam"]).decode()
    print(log_end)
    assert log_end.split()[-1] == "End"


@timed(30) # Test must run faster than 30 seconds
def test_parallel():
    """Test 3-D actuatorLineSource in parallel."""
    get_tutorial_files(case="static")
    out = subprocess.check_output("./Allclean")
    try:
        out = subprocess.check_output(["./Allrun", "3D", str(alpha_deg),
                                       "-parallel"])
    except subprocess.CalledProcessError:
        print(subprocess.check_output(["tail", "-n", "200",
                                       "log.simpleFoam"]).decode())
        assert False
    log_end = subprocess.check_output(["tail", "log.simpleFoam"]).decode()
    print(log_end)
    assert "Finalising parallel run" in log_end


def check_pitching_geom():
    """Check that the geometric values of the pitch match those prescribed."""
    # Read the pitching parameters from fvOptions
    txt = subprocess.check_output("grep reducedFreq system/fvOptions",
                                  shell=True).decode()
    reduced_freq = float(txt.replace(";", " ").split()[1])
    txt = subprocess.check_output("grep amplitude system/fvOptions",
                                  shell=True).decode()
    amp_deg = float(txt.replace(";", " ").split()[1])
    free_stream = 29.9 # Not detected automatically
    chord = 1.22 # Not detected automatically
    offset = 15.0 # Not detected automatically
    omega = reduced_freq*2.0*free_stream/chord
    # Read values from CSV
    df = pd.read_csv("postProcessing/actuatorLines/0/foil.csv")
    # Calculate the prescribed values for alpha_geom_deg
    amp_rad = np.deg2rad(amp_deg)
    alpha_rad_prescribed = amp_rad*np.sin(omega*df.time) + np.deg2rad(offset)
    alpha_deg_prescribed = np.rad2deg(alpha_rad_prescribed)
    np.testing.assert_allclose(df.alpha_geom_deg, alpha_deg_prescribed)


def test_pitching():
    """Test actuatorLineSource harmonic pitching."""
    get_tutorial_files(case="pitching")
    out = subprocess.check_output("./Allclean")
    out = subprocess.check_output("./Allrun")
    check_pitching_geom()


def teardown():
    """Move back into tests directory."""
    os.chdir("../")
