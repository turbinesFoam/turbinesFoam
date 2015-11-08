#!/usr/bin/env python

from __future__ import print_function
import sys

if len(sys.argv) > 1:
    alpha_deg = sys.argv[1]
else:
    alpha_deg = 4.0

print("Setting angle of attack to {} degrees".format(alpha_deg))

with open("system/fvOptions", "w") as f:
    with open("system/fvOptions.template") as template:
        txt = template.read()
    f.write(txt.format(alpha_deg=alpha_deg))
