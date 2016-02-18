#!/usr/bin/env python
"""Visualization for OpenFOAM actuatorLine simulation."""

from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def loadperf():
    df = pd.read_csv("postProcessing/actuatorLines/0/foil.csv")
    df = df.drop_duplicates("time", keep="last")
    df["alpha_rad"] = df.alpha_deg/180.0*np.pi
    df["cn"] =  df.cl*np.cos(df.alpha_rad) - df.cd*np.sin(df.alpha_rad)
    df["cc"] = df.cl*np.sin(df.alpha_rad) - df.cd*np.cos(df.alpha_rad)
    return df


def plot_alpha():
    df = loadperf()
    plt.figure()
    plt.plot(df.time, df.alpha_geom_deg, label="Geometric")
    plt.plot(df.time, df.alpha_deg, label="Actual")
    plt.xlabel("Time (s)")
    plt.ylabel("Angle of attack (deg)")
    plt.legend(loc="best")
    plt.tight_layout()


def plot_cn(t0=0.5):
    df = loadperf()
    plt.figure()
    ind = df.time >= t0
    plt.plot(df.alpha_geom_deg[ind], df.cn[ind])
    plt.xlabel(r"$\alpha$ (geometric, degrees)")
    plt.ylabel(r"$C_N$")
    plt.xlim((0, None))
    plt.ylim((0, None))
    plt.tight_layout()


def plot_cc(t0=0.5):
    df = loadperf()
    plt.figure()
    ind = df.time >= t0
    plt.plot(df.alpha_geom_deg[ind], df.cc[ind])
    plt.xlabel(r"$\alpha$ (geometric, degrees)")
    plt.ylabel(r"$C_C$")
    plt.tight_layout()


if __name__ == "__main__":
    plot_alpha()
    plot_cn()
    plot_cc()
    plt.show()
