#!/usr/bin/env bash
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Copy files from tutorial
./getTutorialFiles.sh static

cp -rf 0.org 0

runApplication blockMesh -dict system/blockMeshDict-3D
runApplication snappyHexMesh -overwrite
# runApplication refineMesh -overwrite -all
runApplication topoSet

if [ "$1" = "-parallel" ]
    then
    python scripts/set_alpha.py $2 -d
    if [[ $WM_PROJECT_VERSION == "3."* ]]
        then
        nProc=$(getNumberOfProcessors)
    else
        nProc=""
    fi
    runApplication decomposePar
    ls -d processor* | xargs -I {} rm -rf ./{}/0
    ls -d processor* | xargs -I {} cp -r 0.org ./{}/0
    mpirun -np $nProc valgrind --tool=callgrind simpleFoam -parallel > log.simpleFoam
    runApplication reconstructPar
else
    python scripts/set_alpha.py $1 -d
    valgrind --tool=callgrind simpleFoam > log.simpleFoam
fi
