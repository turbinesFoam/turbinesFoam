#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alpha_deg", nargs="?", type=float, default=10.0)
    parser.add_argument("--three-dim", "-d", action="store_true", default=False,
                        help="3-D actuator line")

    args = parser.parse_args()
    alpha_deg = args.alpha_deg
    print("Setting angle of attack to {} degrees".format(alpha_deg))

    if args.three_dim:
        semispan = 0.5
        n_elements = 12
    else:
        semispan = 0.05
        n_elements = 1

    with open("system/fvOptions", "w") as f:
        with open("system/fvOptions.template") as template:
            txt = template.read()
        f.write(txt.format(n_elements=n_elements, semispan=semispan,
                           alpha_deg=alpha_deg))
