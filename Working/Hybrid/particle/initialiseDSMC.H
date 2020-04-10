#include "fvCFD.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


IOdictionary dsmcInitialiseDict
(
    IOobject
    (
        "dsmcInitialiseDict",
        particleRegions[0].time().system(),
        particleRegions[0],
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Initialising dsmc for Time = " << runTime.timeName() << nl << endl;

dsmcCloud dsmcInit("dsmc", particleRegions[0], dsmcInitialiseDict);

label totalMolecules = dsmcInit.size();

if (Pstream::parRun())
{
    reduce(totalMolecules, sumOp<label>());
}

Info<< nl << "Total number of molecules added: " << totalMolecules
    << nl << endl;

IOstream::defaultPrecision(15);

if (!particleRegions[0].write())
{
    FatalErrorInFunction
        << "Failed writing dsmcCloud."
        << nl << exit(FatalError);
}

Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;



