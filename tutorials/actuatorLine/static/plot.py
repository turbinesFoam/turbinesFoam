#!/usr/bin/env python
"""This script plots results from the turbinesFoam actuator line tutorial."""

from __future__ import division, print_function
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns


H = 1.0
U_infty = 1.0


def plot_spanwise():
    """Plot spanwise distribution of angle of attack and relative velocity."""
    elements_dir = "postProcessing/actuatorLineElements/0"
    elements = os.listdir(elements_dir)
    dfs = {}
    z_H = np.zeros(len(elements))
    urel = np.zeros(len(elements))
    alpha_deg = np.zeros(len(elements))
    for e in elements:
        i = int(e.replace("foil.element", "").replace(".csv", ""))
        df = pd.read_csv(os.path.join(elements_dir, e))
        z_H[i] = df.z.iloc[-1]/H
        urel[i] = df.rel_vel_mag.iloc[-1]/U_infty
        alpha_deg[i] = df.alpha_deg.iloc[-1]
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7.5, 3.25))
    ax[0].plot(z_H, alpha_deg)
    ax[0].set_ylabel(r"$\alpha$ (deg)")
    ax[1].plot(z_H, urel)
    ax[1].set_ylabel(r"$ | U_{\mathrm{rel}} | / U_\infty $")
    for a in ax:
        a.set_xlabel("$z/H$")
        a.grid(True)
    fig.tight_layout()


def plot_sweep():
    """Plot results from the angle of attack sweep."""
    df = pd.read_csv("processed/alpha_sweep.csv")
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(7.5, 3))
    ax1.plot(df.alpha_geom_deg, df.alpha_deg, "o", label="Detected")
    ax1.plot(df.alpha_geom_deg, df.alpha_geom_deg, "--", label="Geometric")
    ax1.set_xlabel(r"$\alpha$ (geometric, degrees)")
    ax1.set_ylabel(r"$\alpha$ (detected, degrees)")
    ax1.legend(loc="lower right")
    ax2.plot(df.alpha_deg, df.rel_vel_mag, "o", label="Detected")
    ax2.plot(df.alpha_geom_deg, np.ones(len(df)), "--", label="Geometric",
             lw=2)
    ax2.set_xlabel(r"$\alpha$ (detected, degrees)")
    ax2.set_ylabel(r"$|U_\mathrm{rel}|$")
    fig.tight_layout()


if __name__ == "__main__":
    plot_spanwise()
    plt.show()
