#!/usr/bin/env bash
# Install OpenFOAM for Ubuntu--handy for CI

OF_VERS="$1"

if [ "$OF_VERS" = "1706" ]; then
    wget https://sourceforge.net/projects/openfoamplus/files/v1706/OpenFOAM-v1706-windows10.tgz
    sudo tar -xzf  OpenFOAM-v1706-windows10.tgz -C /opt/
    sudo mv /opt/OpenFOAM/OpenFOAM-v1706 /opt/openfoam1706
    sudo apt-get install -qq bison flex m4
else
    sudo add-apt-repository http://dl.openfoam.org/ubuntu
    sudo sh -c "wget -O - http://dl.openfoam.org/gpg.key | apt-key add -"
    sudo apt-get update -qq
    sudo apt-get install -qq openfoam${OF_VERS}
fi
