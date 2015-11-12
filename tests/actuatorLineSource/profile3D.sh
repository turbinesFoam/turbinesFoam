#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

TUTORIAL_DIR="../../tutorials/actuatorLine/simpleFoam"

# Copy files from tutorial
cp -rf $TUTORIAL_DIR/0.org 0
cp $TUTORIAL_DIR/constant/polyMesh/blockMeshDict constant/polyMesh/blockMeshDict-3D

runApplication blockMesh -dict constant/polyMesh/blockMeshDict-3D
runApplication snappyHexMesh -overwrite
# runApplication refineMesh -overwrite -all
runApplication topoSet

if [ "$1" = "-parallel" ]
    then
    python scripts/set_alpha.py $2 -d
    nProc=$(getNumberOfProcessors)
    runApplication decomposePar
    ls -d processor* | xargs -I {} rm -rf ./{}/0
    ls -d processor* | xargs -I {} cp -r $TUTORIAL_DIR/0.org ./{}/0
    mpirun -np $nProc valgrind --tool=callgrind simpleFoam -parallel > log.simpleFoam
    runApplication reconstructPar
else
    python scripts/set_alpha.py $1 -d
    valgrind --tool=callgrind simpleFoam > log.simpleFoam
fi
