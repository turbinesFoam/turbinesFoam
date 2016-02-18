#!/usr/bin/env python
"""Run multiple 2-D simulations varying the angle of attack."""

import numpy as np
from subprocess import call
import os
import pandas as pd
import argparse


def read_force_coeffs():
    """Read force coefficients from `postProcessing/actuatorLines`."""
    df = pd.read_csv("postProcessing/actuatorLines/0/foil.csv")
    df = df.iloc[-1]
    return df[["time", "rel_vel_mag", "alpha_geom_deg", "alpha_deg", "cl",
               "cd", "cm"]]


def read_turbulence_fields():
    """Dummy function for now."""
    return {"k": np.nan, "omega": np.nan, "epsilon": np.nan, "nut": np.nan,
            "z_turbulence": np.nan}


def alpha_sweep(start, stop, step, append=False):
    """Vary the foil angle of attack and log results."""
    alpha_list = np.arange(start, stop, step)
    df_fname = "processed/alpha_sweep.csv"
    if append:
        df = pd.read_csv(df_fname)
    else:
        df = pd.DataFrame(columns=["time", "rel_vel_mag", "alpha_geom_deg",
                                   "alpha_deg", "cl", "cd", "cm", "k", "omega",
                                   "epsilon", "nut", "z_turbulence"])
    for alpha in alpha_list:
        call("./Allclean")
        call(["./Allrun", "2D", str(alpha)])
        d = dict(read_force_coeffs())
        d.update(read_turbulence_fields())
        df = df.append(d, ignore_index=True)
        df.to_csv(df_fname, index=False)


if __name__ == "__main__":
    if not os.path.isdir("processed"):
        os.mkdir("processed")

    parser = argparse.ArgumentParser(description="Vary the foil angle of \
                                     attack and log results.")
    parser.add_argument("start", type=float, help="Start angle of sweep.",
                        nargs="?", default=-15.0)
    parser.add_argument("stop", type=float, help="End angle of sweep. The sweep\
                        does not include this value.", nargs="?", default=15.0)
    parser.add_argument("step", type=float, default=1.0, nargs="?",
                        help="Spacing between values.")
    parser.add_argument("--append", "-a", action="store_true", default=False,
                        help="Append to previous results")
    args = parser.parse_args()

    alpha_sweep(args.start, args.stop, args.step, append=args.append)
