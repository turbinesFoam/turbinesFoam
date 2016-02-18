#!/usr/bin/env sh
# This script copies necessary tutorial files to this directory
# The first arg should be either `static` or `pitching`

cd ${0%/*} || exit 1    # run from this directory

TUTORIAL_DIR="../../tutorials/actuatorLine/$1"

cp -rfT $TUTORIAL_DIR/0.org 0.org
cp -rfT $TUTORIAL_DIR/system system
cp -rfT $TUTORIAL_DIR/constant constant

# Get run and plotting scripts
cp $TUTORIAL_DIR/Allrun Allrun
cp $TUTORIAL_DIR/Allclean Allclean
cp $TUTORIAL_DIR/plot.py plot.py

# Get paramsweep script if static
if [ "$1" = "static" ]
then
    cp $TUTORIAL_DIR/paramsweep.py paramsweep.py
    cp -rfT $TUTORIAL_DIR/scripts scripts
fi

# If running pitching, reduce endTime and correct foilData dir
if [ "$1" = "pitching" ]
then
    sed -i '/endTime /c\endTime         0.1;' system/controlDict
    sed -i 's/..\/..\/..\/resources\/foilData/..\/..\/..\/tutorials\/resources\/foilData/g' system/fvOptions
fi

# Add debug switches to controlDict
cat debugSwitches >> system/controlDict

# Delete reconstructPar from Allrun to save time
sed -i '/reconstructPar/d' Allrun
