#!/usr/bin/env python
"""
This script plots results from `paramsweep.py`.
"""
from __future__ import division, print_function
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns

U_infty = 1.0


if __name__ == "__main__":
    sns.set(style="white", context="paper", font_scale=1.5,
            rc={"axes.grid": True, "legend.frameon": True})
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
    plt.show()
