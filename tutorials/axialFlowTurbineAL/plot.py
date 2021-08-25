#!/usr/bin/env python
"""
This script plots mean power coefficient from the turbinesFoam axial-flow
turbine actuator line tutorial.
"""

from __future__ import division, print_function, absolute_import
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
try:
    import seaborn
except ImportError:
    plt.style.use("ggplot")


def plot_cp(angle0=2160.0):
    df = pd.read_csv("postProcessing/turbines/0/turbine.csv")
    df = df.drop_duplicates("time", keep="last")
    if df.angle_deg.max() < angle0:
        angle0 = 0.0
    print("Performance from {:.1f}--{:.1f} degrees:".format(
            angle0, df.angle_deg.max()))
    print("Mean TSR = {:.2f}".format(df.tsr[df.angle_deg >= angle0].mean()))
    print("Mean C_P = {:.2f}".format(df.cp[df.angle_deg >= angle0].mean()))
    print("Mean C_D = {:.2f}".format(df.cd[df.angle_deg >= angle0].mean()))
    plt.plot(df.angle_deg, df.cp)
    plt.xlabel("Azimuthal angle (degrees)")
    plt.ylabel("$C_P$")
    plt.tight_layout()


def plot_spanwise():
    """Plot spanwise distribution of angle of attack and relative velocity."""
    elements_dir = "postProcessing/actuatorLineElements/0"
    elements = os.listdir(elements_dir)
    dfs = {}
    urel = np.zeros(len(elements))
    alpha_deg = np.zeros(len(elements))
    cl = np.zeros(len(elements))
    f = np.zeros(len(elements))
    root_dist = np.zeros(len(elements))
    for e in elements:
        i = int(e.replace("turbine.blade1.element", "").replace(".csv", ""))
        df = pd.read_csv(os.path.join(elements_dir, e))
        root_dist[i] = df.root_dist.iloc[-1]
        urel[i] = df.rel_vel_mag.iloc[-1]
        alpha_deg[i] = df.alpha_deg.iloc[-1]
        f[i] = df.end_effect_factor.iloc[-1]
        cl[i] = df.cl.iloc[-1]
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7.5, 3.25))
    ax[0].plot(root_dist, cl)
    ax[0].set_ylabel(r"$C_l$")
    ax[1].plot(root_dist, f)
    ax[1].set_ylabel(r"$f$")
    for a in ax:
        a.set_xlabel("$r/R$")
    fig.tight_layout()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "perf":
            plot_cp()
        elif sys.argv[1] == "blade":
            plot_blade_perf()
        elif sys.argv[1] == "spanwise":
            plot_spanwise()
    else:
        plot_cp()
    plt.show()
