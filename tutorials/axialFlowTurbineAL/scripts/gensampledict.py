#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate sampleDict for multiple cross-stream profiles
"""
from __future__ import division, print_function
import numpy as np
import os
import sys

# Input parameters
setformat = "raw"
interpscheme = "cellPoint"
fields = ["U", "vorticity"]
R = 0.45
D = R*2
x_D = 1.0
y_R_max = 1.6
y_R_min = -1.6
ny = 41
z_R_max = 1.6
z_R_min = -1.6
nz = 21

header = r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.4.x                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      sampleDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def main():
    zmin = z_R_min*R
    zmax = z_R_max*R
    z_array = np.linspace(zmin, zmax, nz)
    
    ymax = y_R_max*R
    ymin = y_R_min*R
    
    x = x_D*D
    
    txt = header + "\n" 
    txt += "setFormat " + setformat + "; \n\n"
    txt += "interpolationScheme " + interpscheme + "; \n\n"
    txt += "sets \n ( \n"
    
    for z in z_array:
        txt += "    " + "profile_" + str(z/R) + "\n"
        txt += "    { \n"
        txt += "        type        uniform; \n"
        txt += "        axis        y; \n"
        txt += "        start       (" + str(x) + " " + str(ymin) + " " + str(z) + ");\n"
        txt += "        end         (" + str(x) + " " + str(ymax) + " " + str(z) + ");\n"
        txt += "        nPoints     " + str(ny) + ";\n    }\n\n"
        
    txt += ");\n\n"
    txt += "fields \n(\n"
    
    for field in fields:
        txt += "    " + field + "\n"
        
    txt += "); \n\n"
    txt += "// *********************************************************************** // \n"
    
    with open("system/sampleDict", "w") as f:
        f.write(txt)

if __name__ == "__main__":
    main()
