#!/usr/bin/env python
"""
This script plots mean power coefficient from the turbinesFoam axial-flow
turbine actuator line tutorial.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from modules.processing import *

plt.style.use("fivethirtyeight")

def plot_meancontquiv():
    data = calcwake(t1=0.4)
    y_R = data["y/R"]
    z_R = data["z/R"]
    u = data["meanu"]
    v = data["meanv"]
    w = data["meanw"]
    plt.figure(figsize=(7, 9))
    # Add contours of mean velocity
    cs = plt.contourf(y_R, z_R, u/U, 20, cmap=plt.cm.coolwarm)
    cb = plt.colorbar(cs, shrink=1, extend="both", 
                      orientation="horizontal", pad=0.1)
                      #ticks=np.round(np.linspace(0.44, 1.12, 10), decimals=2))
    cb.set_label(r"$U/U_{\infty}$")
    # Make quiver plot of v and w velocities
    Q = plt.quiver(y_R, z_R, v/U, w/U, angles="xy", width=0.0022,
                   edgecolor="none", scale=3.0)
    plt.xlabel(r"$y/R$")
    plt.ylabel(r"$z/R$")
    plt.quiverkey(Q, 0.8, 0.21, 0.1, r"$0.1 U_\infty$",
               labelpos="E",
               coordinates="figure",
               fontproperties={"size": "small"})
    ax = plt.axes()
    ax.set_aspect(1)
    # Plot circle to represent turbine frontal area
    circ = plt.Circle((0, 0), radius=1, facecolor="none", edgecolor="gray",
                      linewidth=3.0)
    ax.add_patch(circ)
    plt.tight_layout()

def plot_cp(angle0=2160.0):
    df = pd.read_csv("postProcessing/turbines/0/turbine.csv")
    df = df.drop_duplicates("time", take_last=True)
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

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "wake":
            plot_meancontquiv()
        elif sys.argv[1] == "perf":
            plot_cp()
        elif sys.argv[1] == "blade":
            plot_blade_perf()
    else:
        plot_cp()
    plt.show()
