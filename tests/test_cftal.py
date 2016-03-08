#!/usr/bin/env python
"""Tests for `crossFlowTurbineALSource`."""

from __future__ import division, print_function
import subprocess
import pandas as pd
import os
import numpy as np
from numpy.testing import assert_array_almost_equal


element_dir = "postProcessing/actuatorLineElements/0/"
al_dir = "postProcessing/actuatorLines/0/"


def setup():
    os.chdir("crossFlowTurbineALSource")
    out = subprocess.check_output("./getTutorialFiles.sh", shell=True)


def check_created():
    """Test that crossFlowTurbineALSource was created."""
    txt = "Selecting finite volume options model type crossFlowTurbineALSource"
    subprocess.check_output(["grep", txt, "log.pimpleFoam"])


def check_al_file_exists():
    """Test that the actuator line perf file was created."""
    assert os.path.isfile(os.path.join(al_dir, "turbine.blade1.csv"))


def check_element_file_exists():
    """Test that the element perf file was created."""
    assert os.path.isfile(os.path.join(element_dir,
                                       "turbine.blade1.element0.csv"))


def check_perf(angle0=540.0):
    """Test CFTAL performance was written and in reasonable range."""
    df = pd.read_csv("postProcessing/turbines/0/turbine.csv")
    df = df.drop_duplicates("time", keep="last")
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


def check_periodic_tsr():
    """Check that the periodic TSR oscillation was set properly."""
    df = pd.read_csv("postProcessing/turbines/0/turbine.csv")
    df = df.drop_duplicates("time", keep="last")
    theta_rad = np.deg2rad(df.angle_deg)
    tsr_amp = 0.19
    tsr_phase = 1.9
    tsr_mean = 1.9
    nblades = 3
    tsr_ref = tsr_amp*np.cos(nblades*(theta_rad - tsr_phase)) + tsr_mean
    assert_array_almost_equal(df.tsr, tsr_ref)


def test_serial():
    """Test crossFlowTurbineALSource in serial."""
    output_clean = subprocess.check_output("./Allclean")
    try:
        output_run = subprocess.check_output("./Allrun")
    except subprocess.CalledProcessError:
        print(subprocess.check_output(["tail", "-n", "200",
                                       "log.pimpleFoam"]).decode())
        assert False
    check_created()
    check_al_file_exists()
    check_element_file_exists()
    check_perf()
    check_periodic_tsr()
    log_end = subprocess.check_output(["tail", "log.pimpleFoam"]).decode()
    print(log_end)
    assert log_end.split()[-1] == "End"


def test_parallel():
    """Test crossFlowTurbineALSource in parallel."""
    output_clean = subprocess.check_output("./Allclean")
    try:
        output_run = subprocess.check_output(["./Allrun", "-parallel"])
    except subprocess.CalledProcessError:
        print(subprocess.check_output(["tail", "-n", "200",
                                       "log.pimpleFoam"]).decode())
        assert False
    check_created()
    check_al_file_exists()
    check_element_file_exists()
    check_perf()
    check_periodic_tsr()
    log_end = subprocess.check_output(["tail", "log.pimpleFoam"]).decode()
    print(log_end)
    assert "Finalising parallel run" in log_end


def test_multi_re():
    """Test crossFlowTurbineALSource with multiRe profileData."""
    # Switch on change profileData tableType to multiRe
    with open("system/fvOptions") as f:
        txt = f.read()
    txt = txt.replace("tableType   singleRe;", "tableType   multiRe;")
    with open("system/fvOptions", "w") as f:
        f.write(txt)
    test_serial()


def teardown():
    """Move back into tests directory."""
    os.chdir("../")
