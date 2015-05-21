#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Processing for OpenFOAM actuatorLine simulation.
"""
from __future__ import division, print_function
import matplotlib.pyplot as plt
import re
import numpy as np
import os
from styleplot import styleplot
import sys
import foampy
import fdiff
    
# Some constants
R = 0.5
U = 1.0
H = 0.05
D = 1.0
A = H*D
rho = 1000.0

ylabels = {"meanu" : r"$U/U_\infty$",
           "stdu" : r"$\sigma_u/U_\infty$",
           "meanv" : r"$V/U_\infty$",
           "meanw" : r"$W/U_\infty$",
           "meanuv" : r"$\overline{u'v'}/U_\infty^2$"}
    
def loadwake():
    """Loads wake data and returns y/R and statistics."""
    folder = os.listdir("postProcessing/sets")[0]
    flist = os.listdir("postProcessing/sets/"+folder)
    data = {}
    for fname in flist:
        fpath = "postProcessing/sets/"+folder+"/"+fname
        z_H = float(fname.split("_")[1])
        data_s = np.loadtxt(fpath, unpack=True)
        data[z_H] = data_s
    return data
    
def plotwake(plotlist=["meanu"], save=False, savepath="", savetype=".pdf"):
    data = loadwake()
    y_R = data[0][0]/R
    z_H = np.asarray(sorted(data.keys()))
    # Assemble 2-D arrays
    u = np.zeros((len(z_H), len(y_R)))
    v = np.zeros((len(z_H), len(y_R)))
    w = np.zeros((len(z_H), len(y_R)))
    xvorticity = np.zeros((len(z_H), len(y_R)))
    for n in range(len(z_H)):
        u[n,:] = data[z_H[n]][1]
        v[n,:] = data[z_H[n]][2]
        w[n,:] = data[z_H[n]][3]
        xvorticity[n,:] = data[z_H[n]][4]
    def turb_lines():
        plt.hlines(0.5, -1, 1, linestyles='solid', linewidth=2)
        plt.vlines(-1, 0, 0.5, linestyles='solid', linewidth=2)
        plt.vlines(1, 0, 0.5, linestyles='solid', linewidth=2)
    if "meanu" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, u, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.2)
        cb.set_label(r'$U/U_{\infty}$')
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.grid(True)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
    if "meanv" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y/0.5, z, v, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.3)
        cb.set_label(r'$V/U_{\infty}$')
        #turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.grid(True)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
    if "v-wquiver" in plotlist or "all" in plotlist:
        # Make quiver plot of v and w velocities
        plt.figure(figsize=(10,5))
        Q = plt.quiver(y_R, z_H, v, w, angles='xy')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(Q, 0.75, 0.2, 0.1, r'$0.1$ m/s',
                   labelpos='E',
                   coordinates='figure',
                   fontproperties={'size': 'small'})
        plt.tight_layout()
        plt.hlines(0.5, -1, 1, linestyles='solid', colors='r',
                   linewidth=2)
        plt.vlines(-1, -0.2, 0.5, linestyles='solid', colors='r',
                   linewidth=2)
        plt.vlines(1, -0.2, 0.5, linestyles='solid', colors='r',
                   linewidth=2)
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'v-wquiver'+savetype)
    if "xvorticity" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, xvorticity, 10, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.26)
        cb.set_label(r"$\Omega_x$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+'/xvorticity_AD'+savetype)
    if "meancomboquiv" in plotlist or "all" in plotlist:
        plt.figure(figsize=(9, 8))
        # Add contours of mean velocity
        cs = plt.contourf(y_R, z_H, u, 20, cmap=plt.cm.coolwarm)
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.12)
                          #ticks=np.round(np.linspace(0.44, 1.12, 10), decimals=2))
        cb.set_label(r'$U/U_{\infty}$')
        plt.hold(True)
        # Make quiver plot of v and w velocities
        Q = plt.quiver(y_R, z_H, v, w, angles='xy', width=0.0022)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        #plt.ylim(-0.2, 0.78)
        #plt.xlim(-3.2, 3.2)
        plt.xlim(-3.66, 3.66)
        plt.ylim(-1.22, 1.22)
        plt.quiverkey(Q, 0.8, 0.22, 0.1, r'$0.1 U_\infty$',
                   labelpos='E',
                   coordinates='figure',
                   fontproperties={'size': 'small'})
        plt.hlines(0.5, -1, 1, linestyles='solid', colors='gray',
                   linewidth=3)
        plt.hlines(-0.5, -1, 1, linestyles='solid', colors='gray',
                   linewidth=3)
        plt.vlines(-1, -0.5, 0.5, linestyles='solid', colors='gray',
                   linewidth=3)
        plt.vlines(1, -0.5, 0.5, linestyles='solid', colors='gray',
                   linewidth=3)
        ax = plt.axes()
        ax.set_aspect(2.0)
        styleplot()
        if save:
            plt.savefig(savepath+"\\meancomboquiv_AD"+savetype)
    plt.show()
        
def plotexpwake(Re_D, quantity, z_H=0.0, save=False, savepath="", 
                savetype=".pdf", newfig=True, marker="--ok",
                fill="none", figsize=(10, 5)):
    """Plots the transverse wake profile of some quantity. These can be
      * meanu
      * meanv
      * meanw
      * stdu
    """
    U = Re_D/1e6
    label = "Exp."
    folder = exp_path + "/Wake/U_" + str(U) + "/Processed/"
    z_H_arr = np.load(folder + "z_H.npy")
    i = np.where(z_H_arr==z_H)
    q = np.load(folder + quantity + ".npy")[i]
    y_R = np.load(folder + "y_R.npy")[i]
    if newfig:
        plt.figure(figsize=figsize)
    plt.plot(y_R, q/U, marker, markerfacecolor=fill, label=label)
    plt.xlabel(r"$y/R$")
    plt.ylabel(ylabels[quantity])
    plt.grid(True)
    styleplot()

def main():
    p = "Google Drive/Research/Papers/JOT CFT near-wake/Figures"
    if "linux" in sys.platform:
        p = "/home/pete/" + p
    elif "win" in sys.platform:
        p = "C:/Users/Pete/" + p
    plt.close("all")
    
    plotwake(plotlist=["meancomboquiv"], save=True, savepath=p)

if __name__ == "__main__":
    main()
