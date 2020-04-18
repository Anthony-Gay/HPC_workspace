bool LTS = fv::localEulerDdt::enabled(fluidMesh);

tmp<volScalarField> trDeltaT;

if (LTS)
{
    Info<< "Using LTS" << endl;

    trDeltaT = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTName,
                runTime.timeName(),
                fluidMesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidMesh,
            dimensionedScalar(dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
}