#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Processing functions for crossFlowTurbineAL tutorial.
"""
from __future__ import division, print_function
import matplotlib.pyplot as plt
import re
import numpy as np
import os
import sys
import foampy
from pxl import fdiff
import pandas as pd

# Some constants
R = 0.45
U = 10.0
U_infty = U
D = R*2
A = np.pi*R**2
rho = 1000.0

ylabels = {"meanu" : r"$U/U_\infty$",
           "stdu" : r"$\sigma_u/U_\infty$",
           "meanv" : r"$V/U_\infty$",
           "meanw" : r"$W/U_\infty$",
           "meanuv" : r"$\overline{u'v'}/U_\infty^2$"}
           

class WakeMap(object):
    """
    Object that represents a wake map or statistics.
    """
    def __init__(self):
        self.load()
        
    def load_single_time(self, time):
        """
        Loads data from a single time step.
        """
        timedir = "postProcessing/sets/{}".format(time)


def loadwake(time):
    """Loads wake data and returns y/R and statistics."""
    # Figure out if time is an int or float
    if not isinstance(time, str):
        if time % 1 == 0:
            folder = str(int(time))
        else:
            folder = str(time)
    else:
        folder = time
    flist = os.listdir("postProcessing/sets/"+folder)
    data = {}
    for fname in flist:
        fpath = "postProcessing/sets/"+folder+"/"+fname
        z_R = float(fname.split("_")[1])
        data[z_R] = np.loadtxt(fpath, unpack=True)
    return data
    
def calcwake(t1=0.0):
    times = os.listdir("postProcessing/sets")
    times = [float(time) for time in times]
    times.sort()
    times = np.asarray(times)
    data = loadwake(times[0])
    z_R = np.asarray(sorted(data.keys()))
    y_R = data[z_R[0]][0]/R
    # Find first timestep from which to average over
    t = times[times >= t1]
    # Assemble 3-D arrays, with time as first index
    u = np.zeros((len(t), len(z_R), len(y_R)))
    v = np.zeros((len(t), len(z_R), len(y_R)))
    w = np.zeros((len(t), len(z_R), len(y_R)))
    xvorticity = np.zeros((len(t), len(z_R), len(y_R)))
    # Loop through all times
    for n in range(len(t)):
        data = loadwake(t[n])
        for m in range(len(z_R)):
            u[n,m,:] = data[z_R[m]][1]
            v[n,m,:] = data[z_R[m]][2]
            w[n,m,:] = data[z_R[m]][3]
            try:
                xvorticity[n,m,:] = data[z_R[m]][4]
            except IndexError:
                pass
    meanu = u.mean(axis=0)
    meanv = v.mean(axis=0)
    meanw = w.mean(axis=0)
    xvorticity = xvorticity.mean(axis=0)
    return {"meanu" : meanu,
            "meanv" : meanv,
            "meanw" : meanw,
            "xvorticity" : xvorticity,
            "y/R" : y_R, 
            "z/R" : z_R}
    
def plot_al_perf(name="blade1"):
    df_turb = pd.read_csv("postProcessing/turbines/0/turbine.csv")
    df_turb = df_turb.drop_duplicates("time", keep="last")
    df = pd.read_csv("postProcessing/actuatorLines/0/{}.csv".format(name))
    df = df.drop_duplicates("time", keep="last")
    df["angle_deg"] = df_turb.angle_deg
    plt.figure()
    plt.plot(df.angle_deg, df.alpha_deg)
    plt.xlabel("Azimuthal angle (degrees)")
    plt.ylabel("Angle of attack (degrees)")
    plt.tight_layout()
    plt.figure()
    plt.plot(df.angle_deg, df.rel_vel_mag)
    plt.xlabel("Azimuthal angle (degrees)")
    plt.ylabel("Relative velocity (m/s)")
    plt.tight_layout()
    
def plot_blade_perf():
    plot_al_perf("blade1")
    
def plot_spanwise():
    elements_dir = "postProcessing/actuatorLineElements/0"
    elements = os.listdir(elements_dir)
    dfs = {}
    r_R = np.zeros(len(elements))
    fx = np.zeros(len(elements))
    ft = np.zeros(len(elements))
    for e in elements:
        i = int(e.replace("blade1Element", "").replace(".csv", ""))
        df = pd.read_csv(os.path.join(elements_dir, e))
        r_R[i] = np.sqrt(df.y**2 + df.z**2).iloc[-1]/R
        fx[i] = df.fx.iloc[-1]
        ft[i] = np.sqrt(df.fy**2 + df.fz**2).iloc[-1]
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7.5, 3.25))
    ax[0].plot(r_R, 2*ft/(U_infty**2*R))
    ax[0].set_ylabel(r"$F_\theta / (\rho R U_\infty^2)$")
    ax[1].plot(r_R, 2*fx/(R*U_infty**2))
    ax[1].set_ylabel(r"$F_x / (\rho R U_\infty^2 )$")
    for a in ax:
        a.set_xlabel("$r/R$")
    fig.tight_layout()
    
