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

def plot_cp():
    df = pd.read_csv("postProcessing/turbines/0/turbine.csv")
    print("Mean C_P =", df.cp.mean())
    plt.plot(df.time, df.cp)
    plt.xlabel("Time (s)")
    plt.ylabel("$C_P$")
    if savefig:
        if not os.path.isdir("figures"):
            os.makedirs("figures")
        plt.savefig("cp.pdf")
    plt.show()

plot_cp()
