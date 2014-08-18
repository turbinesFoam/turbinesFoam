#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate sampleDict for multiple cross-stream profiles

by Pete Bachant (petebachant@gmail.com)

"""
from __future__ import division, print_function
import numpy as np
import os
import sys
import foampy

# Input parameters
setformat = "raw"
interpscheme = "cellPoint"
fields = ["U", "vorticity"]
x = 1.0
ymax = 1.83
ymin = -1.83
ny = 41
zmax = 1.22
zmin = -1.22
nz = 21

header = r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.x                                 |
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
    z_array = np.linspace(zmin, zmax, nz)
    
    txt = header + "\n" 
    txt += "setFormat " + setformat + "; \n\n"
    txt += "interpolationScheme " + interpscheme + "; \n\n"
    txt += "sets \n ( \n"
    
    for z in z_array:
        txt += "    " + "profile_" + str(z) + "\n"
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
