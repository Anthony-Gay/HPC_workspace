const wordList fluidNames(rp["fluid"]);

    Info<< "Create mesh for time = "
    << runTime.timeName() << nl << endl;

autoPtr<fvMesh> fldMeshPtr
(
    new fvMesh
    (
        IOobject
        (
            fluidNames[0],
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    )
);

fvMesh& fluidMesh = fldMeshPtr();

    