
    const wordList particleNames(rp["particle"]);

    Info<< "Create mesh for time = "
    << runTime.timeName() << nl << endl;

autoPtr<fvMesh> partMeshPtr
(
    new fvMesh
    (
        IOobject
        (
            particleNames[0],
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    )
);

fvMesh& particleMesh = partMeshPtr();
