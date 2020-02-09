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
inline const Foam::word& Foam::DSMCCloudMod<ParcelType>::cloudName() const
{
    return cloudName_;
}


template<class ParcelType>
inline const Foam::fvMesh& Foam::DSMCCloudMod<ParcelType>::mesh() const
{
    return mesh_;
}


template<class ParcelType>
inline const Foam::IOdictionary&
Foam::DSMCCloudMod<ParcelType>::particleProperties() const
{
    return particleProperties_;
}


template<class ParcelType>
inline const Foam::List<Foam::word>&
Foam::DSMCCloudMod<ParcelType>::typeIdList() const
{
    return typeIdList_;
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloudMod<ParcelType>::nParticle() const
{
    return nParticle_;
}


template<class ParcelType>
inline const Foam::List<Foam::DynamicList<ParcelType*>>&
Foam::DSMCCloudMod<ParcelType>::cellOccupancy() const
{
    return cellOccupancy_;
}


template<class ParcelType>
inline Foam::volScalarField& Foam::DSMCCloudMod<ParcelType>::sigmaTcRMax()
{
    return sigmaTcRMax_;
}


template<class ParcelType>
inline Foam::scalarField&
Foam::DSMCCloudMod<ParcelType>::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}


template<class ParcelType>
inline const Foam::List<typename ParcelType::constantProperties>&
Foam::DSMCCloudMod<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::DSMCCloudMod<ParcelType>::constProps
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
inline Foam::Random& Foam::DSMCCloudMod<ParcelType>::rndGen()
{
    return rndGen_;
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudMod<ParcelType>::qBF()
{
    return q_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volVectorField::Boundary&
Foam::DSMCCloudMod<ParcelType>::fDBF()
{
    return fD_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudMod<ParcelType>::rhoNBF()
{
    return rhoN_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudMod<ParcelType>::rhoMBF()
{
    return rhoM_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudMod<ParcelType>::linearKEBF()
{
    return linearKE_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudMod<ParcelType>::internalEBF()
{
    return internalE_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volScalarField::Boundary&
Foam::DSMCCloudMod<ParcelType>::iDofBF()
{
    return iDof_.boundaryFieldRef();
}


template<class ParcelType>
inline Foam::volVectorField::Boundary&
Foam::DSMCCloudMod<ParcelType>::momentumBF()
{
    return momentum_.boundaryFieldRef();
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudMod<ParcelType>::boundaryT() const
{
    return boundaryT_;
}


template<class ParcelType>
inline const Foam::volVectorField&
Foam::DSMCCloudMod<ParcelType>::boundaryU() const
{
    return boundaryU_;
}


template<class ParcelType>
inline const Foam::BinaryCollisionModel<Foam::DSMCCloudMod<ParcelType>>&
Foam::DSMCCloudMod<ParcelType>::binaryCollision() const
{
    return binaryCollisionModel_;
}


template<class ParcelType>
inline Foam::BinaryCollisionModel<Foam::DSMCCloudMod<ParcelType>>&
Foam::DSMCCloudMod<ParcelType>::binaryCollision()
{
    return binaryCollisionModel_();
}


template<class ParcelType>
inline const Foam::WallInteractionModel<Foam::DSMCCloudMod<ParcelType>>&
Foam::DSMCCloudMod<ParcelType>::wallInteraction() const
{
    return wallInteractionModel_;
}


template<class ParcelType>
inline Foam::WallInteractionModel<Foam::DSMCCloudMod<ParcelType>>&
Foam::DSMCCloudMod<ParcelType>::wallInteraction()
{
    return wallInteractionModel_();
}


template<class ParcelType>
inline const Foam::InflowBoundaryModel<Foam::DSMCCloudMod<ParcelType>>&
Foam::DSMCCloudMod<ParcelType>::inflowBoundary() const
{
    return inflowBoundaryModel_;
}


template<class ParcelType>
inline Foam::InflowBoundaryModel<Foam::DSMCCloudMod<ParcelType>>&
Foam::DSMCCloudMod<ParcelType>::inflowBoundary()
{
    return inflowBoundaryModel_();
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloudMod<ParcelType>::massInSystem() const
{
    scalar sysMass = 0.0;

    forAllConstIter(typename DSMCCloudMod<ParcelType>, *this, iter)
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
inline Foam::vector Foam::DSMCCloudMod<ParcelType>::linearMomentumOfSystem() const
{
    vector linearMomentum(Zero);

    forAllConstIter(typename DSMCCloudMod<ParcelType>, *this, iter)
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
Foam::DSMCCloudMod<ParcelType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename DSMCCloudMod<ParcelType>, *this, iter)
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
Foam::DSMCCloudMod<ParcelType>::internalEnergyOfSystem() const
{
    scalar internalEnergy = 0.0;

    forAllConstIter(typename DSMCCloudMod<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        internalEnergy += p.Ei();
    }

    return nParticle_*internalEnergy;
}


template<class ParcelType>
inline Foam::scalar Foam::DSMCCloudMod<ParcelType>::maxwellianAverageSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


template<class ParcelType>
inline Foam::scalarField Foam::DSMCCloudMod<ParcelType>::maxwellianAverageSpeed
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
inline Foam::scalar Foam::DSMCCloudMod<ParcelType>::maxwellianRMSSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField Foam::DSMCCloudMod<ParcelType>::maxwellianRMSSpeed
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
Foam::DSMCCloudMod<ParcelType>::maxwellianMostProbableSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


template<class ParcelType>
inline Foam::scalarField
Foam::DSMCCloudMod<ParcelType>::maxwellianMostProbableSpeed
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
inline const Foam::volScalarField& Foam::DSMCCloudMod<ParcelType>::q() const
{
    return q_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::DSMCCloudMod<ParcelType>::fD() const
{
    return fD_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudMod<ParcelType>::rhoN() const
{
    return rhoN_;
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::DSMCCloudMod<ParcelType>::rhoM() const
{
    return rhoM_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudMod<ParcelType>::dsmcRhoN() const
{
    return dsmcRhoN_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudMod<ParcelType>::linearKE() const
{
    return linearKE_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudMod<ParcelType>::internalE() const
{
    return internalE_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::DSMCCloudMod<ParcelType>::iDof() const
{
    return iDof_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::DSMCCloudMod<ParcelType>::momentum() const
{
    return momentum_;
}


template<class ParcelType>
inline void Foam::DSMCCloudMod<ParcelType>::clear()
{
    return IDLList<ParcelType>::clear();
}


// ************************************************************************* //
