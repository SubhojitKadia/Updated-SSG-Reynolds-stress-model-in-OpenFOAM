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

\*---------------------------------------------------------------------------*/

#include "SSGMod.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "wallFvPatch.H"

//added
#include "wallDist.H"

#include "bound.H"
#include "vectorIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void SSGMod<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = this->Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> SSGMod<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
SSGMod<BasicMomentumTransportModel>::SSGMod
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    ReynoldsStress<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            3.4
        )
    ),
    C1s_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1s",
            this->coeffDict_,
            1.8
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            4.2
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0.8
        )
    ),
    C3s_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3s",
            this->coeffDict_,
            1.3
        )
    ),
    C4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4",
            this->coeffDict_,
            1.25
        )
    ),
    C5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C5",
            this->coeffDict_,
            0.4
        )
    ),

    Ceps1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1",
            this->coeffDict_,
            //1.44
			1.45
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            //1.92
			1.9
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            //0.25
			0.22
        )
    ),
    Ceps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps",
            this->coeffDict_,
            //0.15
			0.18
        )
    ),
	
	freeSurfaceReflection_// added
    (
        Switch::lookupOrAddToDict
        (
            "freeSurfaceReflection",
            this->coeffDict_,
            true
        )
    ),
    kappa_// added
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    Cref1_// added
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cref1",
            this->coeffDict_,
            0.5
        )
    ),
    Cref2_// added
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cref2",
            this->coeffDict_,
            0.1
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*tr(this->R_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
	CC_//added
    (
        IOobject
        (
            "cellCentres",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1.0*this->mesh_.C()
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);

        this->boundNormalStress(this->R_);
        bound(epsilon_, this->epsilonMin_);
        k_ = 0.5*tr(this->R_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool SSGMod<BasicMomentumTransportModel>::read()
{
    if (ReynoldsStress<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C1s_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        C3s_.readIfPresent(this->coeffDict());
        C4_.readIfPresent(this->coeffDict());
        C5_.readIfPresent(this->coeffDict());

        Ceps1_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
        Ceps_.readIfPresent(this->coeffDict());
		
		freeSurfaceReflection_.readIfPresent("freeSurfaceReflection", this->coeffDict());//added
        kappa_.readIfPresent(this->coeffDict());//added
        Cref1_.readIfPresent(this->coeffDict());//added
        Cref2_.readIfPresent(this->coeffDict());//added

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> SSGMod<BasicMomentumTransportModel>::DREff() const
{
    return volSymmTensorField::New
    (
        "DREff",
        (Cs_*(this->k_/this->epsilon_))*this->R_ + I*this->nu()
    );
}


template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> SSGMod<BasicMomentumTransportModel>::DepsilonEff() const
{
    return volSymmTensorField::New
    (
        "DepsilonEff",
        (Ceps_*(this->k_/this->epsilon_))*this->R_ + I*this->nu()
    );
}


template<class BasicMomentumTransportModel>
void SSGMod<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volSymmTensorField& R = this->R_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    ReynoldsStress<RASModel<BasicMomentumTransportModel>>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField P(-twoSymm(R & gradU));
    volScalarField G(this->GName(), 0.5*mag(tr(P)));

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        Ceps1_*alpha*rho*G*epsilon_/k_
      - fvm::Sp(Ceps2_*alpha*rho*epsilon_/k_, epsilon_)
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Correct the trace of the tensorial production to be consistent
    // with the near-wall generation from the wall-functions
    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label celli = curPatch.faceCells()[facei];
                P[celli] *= min
                (
                    G[celli]/(0.5*mag(tr(P[celli])) + small),
                    1.0
                );
            }
        }
    }

    volSymmTensorField b(dev(R)/(2*k_));
    volSymmTensorField S(symm(gradU));
    volTensorField Omega(skew(gradU));

    // Reynolds stress equation
    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(alpha, rho, R)
      + fvm::div(alphaRhoPhi, R)
      - fvm::laplacian(alpha*rho*DREff(), R)
      + fvm::Sp(((C1_/2)*epsilon_ + (C1s_/2)*G)*alpha*rho/k_, R)// G part is rapid pressure strain
     ==
        alpha*rho*P
      - ((1.0/3.0)*I)*(((2.0 - C1_)*epsilon_ - C1s_*G)*alpha*rho)// G part is rapid pressure strain
      + (C2_*(alpha*rho*epsilon_))*dev(innerSqr(b))//nonlinear part of the slow pressure strain
      + alpha*rho*k_
       *(
           (C3_ - C3s_*mag(b))*S//was (C3_ - C3s_*mag(b))*dev(S)
          + C4_*dev(twoSymm(b&S))
          + C5_*twoSymm(b&Omega)
        )//rapid pressure strain
      + this->RSource()
      + fvModels.source(alpha, rho, R)
    );
	
	if (freeSurfaceReflection_)//added//
    {
        const volVectorField& n_(wallDist::New(this->mesh_).n());// modified later
        const volScalarField& y_(wallDist::New(this->mesh_).y());// modified later
		//// added
		// inputs
		const dimensionedScalar h = max(this->mesh_.Cf().component(1));//flow depth in m
		const dimensionedScalar b1 = -2.0*min(this->mesh_.Cf().component(2)); // channel width in m //-2 as half channel width is simulated along the -z direction
		
		// Top width of the channel
		const dimensionedScalar bt = b1;// top width of the flow area, if the channel filling is lower than the maximum width
		// or
		//const dimensionedScalar bt = 0.866*b1;// top width of the flow area, if the channel filling is at the 3/4ths of the circle
		//
		
		volScalarField yc("yc", y_);
		volScalarField zc("zc", y_);
		// initialization of the average distance from the free surface
		volScalarField y1("y1", y_);
		// initialization of the unit normal vector 
		volVectorField n1("n1", n_);
		forAll(n_, i)
		{
			n1[i] = vector(0, 1, 0);// for free surface
			yc[i] = CC_[i].component(1);
			zc[i] = CC_[i].component(2);
			scalar y_1 = h.value() - yc[i]; // cell centre from the top boundary
			scalar z_1 = 0.5*bt.value() - mag(zc[i]); // horizontal distance b/w the cell centre and the mixed corner
			// to obtain the average distance from the free surface
			scalar nn = 500; // top width is divided into 500 segments (501 nodes)
			scalar ii = 0;
			scalar sum = 0.0; // initialising the integral value deltaTheta/(distance)^2
			for(ii=0; ii<= nn-1; ii++)
			{
				scalar z_2 = bt.value()*ii/nn; // free surface points from the mixed corner
				scalar z_2_next = bt.value()*(ii + 1)/nn; // next free surface points from the mixed corner
				scalar deltaTheta = 0.0;
				if(z_2 < z_1 && z_2_next > z_1)
				{
					deltaTheta = atan(mag(z_1 - z_2)/y_1) + atan(mag(z_1 - z_2_next)/y_1);
				}
				else
				{
					deltaTheta = mag(atan(mag(z_1 - z_2)/y_1) - atan(mag(z_1 - z_2_next)/y_1));
				}
				scalar sqrDist = sqr(y_1) + sqr(z_1 - z_2);
				sum += deltaTheta/sqrDist;
			}
			y1[i] = sqrt(1.0/(2.0*sum/constant::mathematical::pi)); // average distance from the free surface [see Naot and Rodi (1982))
		}	
		const volScalarField Lt1 = (pow(Cmu_, 0.75)/kappa_)*(pow(k_, 1.5)/epsilon_);// lurbulent length scale Lt1 added	
		const volScalarField fw1 = sqr(Lt1/(y1 + 0.16*Lt1));
        const volSymmTensorField reflect
        (
			(Cref1_*(epsilon_/k_))*R //+ Cref2_*Rapid pressure-strain part
			//- (Cref2_*C1s_*G)*b // addition of this part causes instability near the free surface
			+ (Cref2_*k_)
			*(
				(C3_ - C3s_*mag(b))*S//in place of (C3_ - C3s_*mag(b))*dev(S)
				+ C4_*dev(twoSymm(b&S))
				+ C5_*twoSymm(b&Omega)
			)
        );
		Info<< reflect.size() << endl;//n1.write();//y1.write();
		REqn.ref() +=
            (3.0*alpha*rho*fw1)
           *(dev(symm((n1 & reflect)*n1)));
    }

    REqn.ref().relax();
    fvConstraints.constrain(REqn.ref());
    solve(REqn);
    fvConstraints.constrain(R);

    this->boundNormalStress(R);

    k_ = 0.5*tr(R);

    correctNut();

    // Correct wall shear-stresses when applying wall-functions//
    this->correctWallShearStress(R);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
