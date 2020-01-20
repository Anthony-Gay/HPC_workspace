/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::word& Foam::DSMCCloudV1<ParcelType>::cloudName() const
{
    return cloudName_;
}


template<class ParcelType>
inline const Foam::fvMesh& Foam::DSMCCloudV1<ParcelType>::mesh() const
{
    return mesh_;
}


template<class ParcelType>
inline const Foam::IOdictionary&
Foam::DSMCCloudV1<ParcelType>::particleProperties() const
{
    return particleProperties_;
}


template<class ParcelType>
inline const Foam::List<Foam::word>&
Foam::DSMCCloudV1<ParcelType>::typeIdList() const
{
    return typeIdList_;
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloudV1<ParcelType>::nParticle() const
{
    return nParticle_;
}


template<class ParcelType>
inline const Foam::List<Foam::DynamicList<ParcelType*>>&
Foam::DSMCCloudV1<ParcelType>::cellOccupancy() const
{
    return cellOccupancy_;
}


template<class ParcelType>
inline Foam::volScalarField& Foam::DSMCCloudV1<ParcelType>::sigmaTcRMax()
{
    return sigmaTcRMax_;
}


template<class ParcelType>
inline Foam::scalarField&
Foam::DSMCCloudV1<ParcelType>::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}


template<class ParcelType>
inline const Foam::List<typename ParcelType::constantProperties>&
Foam::DSMCCloudV1<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::DSMCCloudV1<ParcelType>::constProps
(
    label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorInFunction
            << "constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}


template<class ParcelType>
inline Foam::Random& Foam::DSMCCloudV1<ParcelType>::rndGen()
{
    return rndGen_;
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudV1<ParcelType>::qBF()
{
    return q_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volVectorField::Boundary&
Foam::DSMCCloudV1<ParcelType>::fDBF()
{
    return fD_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudV1<ParcelType>::rhoNBF()
{
    return rhoN_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudV1<ParcelType>::rhoMBF()
{
    return rhoM_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudV1<ParcelType>::linearKEBF()
{
    return linearKE_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudV1<ParcelType>::internalEBF()
{
    return internalE_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudV1<ParcelType>::iDofBF()
{
    return iDof_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volVectorField::Boundary&
Foam::DSMCCloudV1<ParcelType>::momentumBF()
{
    return momentum_.boundaryFieldRef();
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudV1<ParcelType>::boundaryT() const
{
    return boundaryT_;
}


template<class ParcelType>
inline const Foam::volVectorField&
Foam::DSMCCloudV1<ParcelType>::boundaryU() const
{
    return boundaryU_;
}


template<class ParcelType>
inline const Foam::BinaryCollisionModel<Foam::DSMCCloudV1<ParcelType>>&
Foam::DSMCCloudV1<ParcelType>::binaryCollision() const
{
    return binaryCollisionModel_;
}


template<class ParcelType>
inline Foam::BinaryCollisionModel<Foam::DSMCCloudV1<ParcelType>>&
Foam::DSMCCloudV1<ParcelType>::binaryCollision()
{
    return binaryCollisionModel_();
}


template<class ParcelType>
inline const Foam::WallInteractionModel<Foam::DSMCCloudV1<ParcelType>>&
Foam::DSMCCloudV1<ParcelType>::wallInteraction() const
{
    return wallInteractionModel_;
}


template<class ParcelType>
inline Foam::WallInteractionModel<Foam::DSMCCloudV1<ParcelType>>&
Foam::DSMCCloudV1<ParcelType>::wallInteraction()
{
    return wallInteractionModel_();
}


template<class ParcelType>
inline const Foam::InflowBoundaryModel<Foam::DSMCCloudV1<ParcelType>>&
Foam::DSMCCloudV1<ParcelType>::inflowBoundary() const
{
    return inflowBoundaryModel_;
}


template<class ParcelType>
inline Foam::InflowBoundaryModel<Foam::DSMCCloudV1<ParcelType>>&
Foam::DSMCCloudV1<ParcelType>::inflowBoundary()
{
    return inflowBoundaryModel_();
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloudV1<ParcelType>::massInSystem() const
{
    scalar sysMass = 0.0;

    forAllConstIter(typename DSMCCloudV1<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        sysMass += cP.mass();
    }

    return nParticle_*sysMass;
}


template<class ParcelType>
inline Foam::vector Foam::DSMCCloudV1<ParcelType>::linearMomentumOfSystem() const
{
    vector linearMomentum(Zero);

    forAllConstIter(typename DSMCCloudV1<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearMomentum += cP.mass()*p.U();
    }

    return nParticle_*linearMomentum;
}


template<class ParcelType>
inline Foam::scalar
Foam::DSMCCloudV1<ParcelType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename DSMCCloudV1<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        const typename ParcelType::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearKineticEnergy += 0.5*cP.mass()*(p.U() & p.U());
    }

    return nParticle_*linearKineticEnergy;
}


template<class ParcelType>
inline Foam::scalar
Foam::DSMCCloudV1<ParcelType>::internalEnergyOfSystem() const
{
    scalar internalEnergy = 0.0;

    forAllConstIter(typename DSMCCloudV1<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        internalEnergy += p.Ei();
    }

    return nParticle_*internalEnergy;
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloudV1<ParcelType>::maxwellianAverageSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


template<class ParcelType>
inline Foam::scalarField Foam::DSMCCloudV1<ParcelType>::maxwellianAverageSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
    return tfld();
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloudV1<ParcelType>::maxwellianRMSSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField Foam::DSMCCloudV1<ParcelType>::maxwellianRMSSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        sqrt(3.0*physicoChemical::k.value()*temperature/mass);
    return tfld();
}


template<class ParcelType>
inline Foam::scalar
Foam::DSMCCloudV1<ParcelType>::maxwellianMostProbableSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField
Foam::DSMCCloudV1<ParcelType>::maxwellianMostProbableSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    tmp<scalarField> tfld =
        sqrt(2.0*physicoChemical::k.value()*temperature/mass);
    return tfld();
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::DSMCCloudV1<ParcelType>::q() const
{
    return q_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::DSMCCloudV1<ParcelType>::fD() const
{
    return fD_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudV1<ParcelType>::rhoN() const
{
    return rhoN_;
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::DSMCCloudV1<ParcelType>::rhoM() const
{
    return rhoM_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudV1<ParcelType>::dsmcRhoN() const
{
    return dsmcRhoN_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudV1<ParcelType>::linearKE() const
{
    return linearKE_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudV1<ParcelType>::internalE() const
{
    return internalE_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudV1<ParcelType>::iDof() const
{
    return iDof_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::DSMCCloudV1<ParcelType>::momentum() const
{
    return momentum_;
}


template<class ParcelType>
inline void Foam::DSMCCloudV1<ParcelType>::clear()
{
    return IDLList<ParcelType>::clear();
}


// ************************************************************************* //
