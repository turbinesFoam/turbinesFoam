#!/usr/bin/env python
"""
This script plots mean power coefficient from the turbinesFoam cross-flow
turbine actuator line tutorial.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
from pxl.timeseries import smooth

plt.style.use("ggplot")
savefig = False

def plot_cp(angle0=540.0):
    df = pd.read_csv("postProcessing/turbines/0/turbine.csv")
    df = df.drop_duplicates("time", take_last=True)
    df.cp = smooth(df.cp, 4)
    if df.angle_deg.max() < angle0:
        angle0 = 0.0
    print("Performance from {:.1f}--{:.1f} degrees:".format(angle0, df.angle_deg.max()))
    print("Mean TSR = {:.2f}".format(df.tsr[df.angle_deg >= angle0].mean()))
    print("Mean C_P = {:.2f}".format(df.cp[df.angle_deg >= angle0].mean()))
    print("Mean C_D = {:.2f}".format(df.cd[df.angle_deg >= angle0].mean()))
    plt.plot(df.angle_deg, df.cp)
    plt.xlabel("Azimuthal angle (degrees)")
    plt.ylabel("$C_P$")
    if savefig:
        if not os.path.isdir("figures"):
            os.makedirs("figures")
        plt.savefig("cp.pdf")
    plt.show()

if __name__ == "__main__":
    plot_cp()
