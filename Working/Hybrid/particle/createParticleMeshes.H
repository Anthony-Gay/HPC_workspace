    const wordList particleNames(rp["particle"]);

    PtrList<fvMesh> particleRegions(particleNames.size());

    forAll(particleNames, i)
    {
        Info<< "Create particle mesh for region " << particleNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        particleRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    particleNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
