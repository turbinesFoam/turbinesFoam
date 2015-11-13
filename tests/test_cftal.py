#!/usr/bin/env python
"""Tests for `crossFlowTurbineALSource`."""

from __future__ import division, print_function
import subprocess
import pandas as pd
import os
import numpy as np
# import sys
# sys.path = [os.path.abspath("./crossFlowTurbineALSource")] + sys.path
# from modules import processing


def setup():
    os.chdir("crossFlowTurbineALSource")
    out = subprocess.check_output("./getTutorialFiles.sh", shell=True)


def check_created():
    """Test that crossFlowTurbineALSource was created."""
    txt = "Selecting finite volume options model type crossFlowTurbineALSource"
    subprocess.check_output(["grep", txt, "log.pimpleFoam"])


def check_perf(angle0=540.0):
    """Test CFTAL performance was written and in reasonable range."""
    df = pd.read_csv("postProcessing/turbines/0/turbine.csv")
    df = df.drop_duplicates("time", take_last=True)
    if df.angle_deg.max() < angle0:
        angle0 = 0.0
    mean_tsr = df.tsr[df.angle_deg >= angle0].mean()
    mean_cp = df.cp[df.angle_deg >= angle0].mean()
    mean_cd = df.cd[df.angle_deg >= angle0].mean()
    print("Performance from {:.1f}--{:.1f} degrees:".format(angle0,
          df.angle_deg.max()))
    print("Mean TSR = {:.2f}".format(mean_tsr))
    print("Mean C_P = {:.2f}".format(mean_cp))
    print("Mean C_D = {:.2f}".format(mean_cd))
    assert 1.6 < mean_tsr < 2.1
    assert 0.2 < mean_cp < 1.0
    assert 0.8 < mean_cd < 1.8


def test_serial():
    """Test crossFlowTurbineALSource in serial."""
    output_clean = subprocess.check_output("./Allclean")
    output_run = subprocess.check_output("./Allrun")
    check_created()
    check_perf()


def test_parallel():
    """Test crossFlowTurbineALSource in parallel."""
    output_clean = subprocess.check_output("./Allclean")
    output_run = subprocess.check_output(["./Allrun", "-parallel"])
    check_created()
    check_perf()


def teardown():
    """Move back into tests directory."""
    os.chdir("../")
