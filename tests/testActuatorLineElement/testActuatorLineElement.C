/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam

Application
    testActuatorLineElement

Description
    Unit tests for actuatorLineElement

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "actuatorLineElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    word name = "testElement";
    dictionary dict;
    dict.add("profileName", "testProfile");
    dictionary profileData;
    dict.add("profileData", profileData); 
    
    Foam::fv::actuatorLineElement element(name, dict, mesh);
    
    Info<< element.name() << endl;

    return 0;
}


// ************************************************************************* //
