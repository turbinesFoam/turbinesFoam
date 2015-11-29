#!/usr/bin/env sh
# This script copies necessary tutorial files to this directory

cd ${0%/*} || exit 1    # run from this directory

TUTORIAL_DIR="../../tutorials/axialFlowTurbineAL"

cp -rfT $TUTORIAL_DIR/0.org 0.org
cp -rfT $TUTORIAL_DIR/system system
cp -rfT $TUTORIAL_DIR/constant constant

# Get run and plotting scripts
cp $TUTORIAL_DIR/Allrun Allrun
cp $TUTORIAL_DIR/Allclean Allclean
cp $TUTORIAL_DIR/plot.py plot.py

# Fix line in fvOptions
sed -i 's/..\/..\/resources\/foilData/..\/..\/..\/tutorials\/resources\/foilData/g' system/fvOptions

# Reduce endTime
sed -i '/endTime /c\endTime         0.003;' system/controlDict

# Add debug switches to controlDict
cat debugSwitches >> system/controlDict

# Delete reconstructPar from Allrun to save time
sed -i '/reconstructPar/d' Allrun
