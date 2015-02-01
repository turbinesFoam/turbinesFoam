#!/usr/bin/env python
"""
This script plots mean power coefficient from the turbinesFoam cross-flow
turbine actuator line tutorial.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

plt.style.use("ggplot")
savefig = False

def plot_cp(angle0=540.0):
    df = pd.read_csv("postProcessing/turbines/0/turbine.csv")
    if df.angle_deg.max() < angle0:
        angle0 = 0.0
    print("Performance from {} degrees onward:".format(angle0))
    print("Mean TSR =", df.tsr[df.angle_deg >= angle0].mean())
    print("Mean C_P =", df.cp[df.angle_deg >= angle0].mean())
    print("Mean C_D =", df.cd[df.angle_deg >= angle0].mean())
    plt.plot(df.angle_deg, df.cp)
    plt.xlabel("Azimuthal angle (degrees)")
    plt.ylabel("$C_P$")
    if savefig:
        if not os.path.isdir("figures"):
            os.makedirs("figures")
        plt.savefig("cp.pdf")
    plt.show()

plot_cp()
