    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//Read the turbine Options
    Info << "Reading Turbine Options On or Off..." << endl;
    IOdictionary TurbineOption
    (
        IOobject
        (
            "TurbineOption",
            runTime.time().constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Are turbines to be activated?
    bool turbineArrayOn(TurbineOption.lookupOrDefault<bool>("turbineArrayOn", false));
   if(turbineArrayOn==1) Info <<"All The Turbine are Switch On    "<<turbineArrayOn<<endl;
   if(turbineArrayOn==0) Info <<"All The Turbine are Switch Off   "<<turbineArrayOn<<endl;

  
   // Create an object of the horizontalWindTurbineArray class if there  is to be a turbine array
    turbineModels::horizontalAxisWindTurbinesALM turbines(U);


#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
