/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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

Class
    Foam::RASModels::SSGMod

Description
    Speziale, Sarkar and Gatski Reynolds-stress turbulence model for
    incompressible and compressible flows.

    Reference:
    \verbatim
        Speziale, C. G., Sarkar, S., & Gatski, T. B. (1991).
        Modelling the pressure–strain correlation of turbulence:
        an invariant dynamical systems approach.
        Journal of Fluid Mechanics, 227, 245-272.
    \endverbatim

    Including the generalized gradient diffusion model of
    Daly and Harlow:
    \verbatim
        Daly, B. J., & Harlow, F. H. (1970).
        Transport equations in turbulence.
        Physics of Fluids (1958-1988), 13(11), 2634-2649.
    \endverbatim
	
	Optional Shir (1973) and Gibson-Launder wall-reflection is also provided for the free surface damping:
    \verbatim
        Gibson, M. M., & Launder, B. E. (1978).
        Ground effects on pressure fluctuations in the
        atmospheric boundary layer.
        Journal of Fluid Mechanics, 86(03), 491-511.
    \endverbatim

    The default model coefficients are:// ref. Cokljat (1993), Cokljat & Younis (1995)
    \verbatim
        SSGModCoeffs
        {
            Cmu             0.09;

            C1              3.4;
            C1s             1.8;
            C2              4.2;
            C3              0.8;
            C3s             1.3;
            C4              1.25;
            C5              0.4;

            Ceps1           1.45;
            Ceps2           1.9;
            Cs              0.22;
            Ceps            0.18;
			
			freeSurfaceReflection  yes; //added
            kappa           0.41;
            Cref1           0.5;
            Cref2           0.1;

            couplingFactor  1.0;
        }
    \endverbatim

SourceFiles
    SSGMod.C

\*---------------------------------------------------------------------------*/

#ifndef SSGMod_H
#define SSGMod_H

#include "RASModel.H"
#include "ReynoldsStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class SSGMod Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class SSGMod
:
    public ReynoldsStress<RASModel<BasicMomentumTransportModel>>
{
protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;

            dimensionedScalar C1_;
            dimensionedScalar C1s_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar C3s_;
            dimensionedScalar C4_;
            dimensionedScalar C5_;

            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Cs_;
            dimensionedScalar Ceps_;
			
		// freeSurface-refection coefficients // added

            Switch freeSurfaceReflection_;
            dimensionedScalar kappa_;
            dimensionedScalar Cref1_;
            dimensionedScalar Cref2_;

        // Fields

            volScalarField k_;
            volScalarField epsilon_;
			volVectorField CC_;//added


    // Protected Member Functions

        //- Update the eddy-viscosity
        virtual void correctNut();
		
		//- Source term for the epsilon equation
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("SSGMod");


    // Constructors

        //- Construct from components
        SSGMod
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        SSGMod(const SSGMod&) = delete;


    //- Destructor
    virtual ~SSGMod()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }
		
		//- Return the cell centes//added
        virtual tmp<volVectorField> cellCentres() const
        {
            return CC_;
        }
		
		//- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return volScalarField::New
            (
                "omega",
                epsilon_/(Cmu_*k_),
                epsilon_.boundaryField().types()
            );
        }

        //- Return the effective diffusivity for R
        tmp<volSymmTensorField> DREff() const;

        //- Return the effective diffusivity for epsilon
        tmp<volSymmTensorField> DepsilonEff() const;

        //- Solve the turbulence equations and correct eddy-Viscosity and
        //  related properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SSGMod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SSGMod.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
