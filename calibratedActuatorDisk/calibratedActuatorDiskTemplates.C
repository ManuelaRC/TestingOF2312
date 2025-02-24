/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd
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

#include "calibratedActuatorDisk.H"
#include "myLookupTable4Turbines.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "volFields.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::calibratedActuatorDisk::calc
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    switch (forceMethod_)
    {
        case forceMethodType::constantParameters_:
        {
            calcConstantParameters(alpha, rho, eqn);
            break;
        }

        case forceMethodType::constantCt_:
        {
            calcConstantCt(alpha, rho, eqn);
            break;
        }
        
        case forceMethodType::calafMethod_:
        {
            calcCalafMethod(alpha, rho, eqn);
            break;
        }

        case forceMethodType::interpolatedFromTable:
        {
            calcInterpolatedFromTable(alpha, rho, eqn);
            break;
        }

        default:
            break;
    }
}

// * * * * * * * * * * * * *                    * * * * * * * * * * * //



template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::calibratedActuatorDisk::monitorData
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    const vectorField& U,
    const scalarField& cellsV
)
{
    vector Udisk(Zero);
    rhoRef_ = 0.0; 

    if (monitorMethod_ == "actuatorDisk")
    {
        for (const label celli : cells_)
        {
            Udisk += U[celli]*cellsV[celli];
            rhoRef_ += rho[celli]*cellsV[celli];
        }
        reduce(Udisk, sumOp<vector>());
        reduce(rhoRef_, sumOp<scalar>());
        Udisk /= V(); 
        rhoRef_ /= V(); 
        Umonitor_ = Udisk & diskDir_;
    }

    else if (monitorMethod_ == "upstreamPoint")
    {
        label localCelli = mesh_.findCell(monitorPoint_);
        label globalCelli = returnReduce(localCelli, maxOp<label>());
        
        if (globalCelli == -1)
        {
            FatalErrorInFunction
                << "No owner cell found for point - Verify input coordinates: "
                << monitorPoint_ << endl;
        }
        else
        {
            Umonitor_ = U[globalCelli] & diskDir_;
            rhoRef_ = rho[globalCelli];
        }

        // Optionally, further reduce these values if needed (depending on your algorithm)
        reduce(Umonitor_, minOp<scalar>());
        reduce(rhoRef_, minOp<scalar>());
    }


    else if (monitorMethod_ == "upstreamCellSet")
    {
        labelList monitorCellSetIDs_ = cellSet(mesh(), monitorCellSetName_).sortedToc();
        returnReduce(monitorCellSetIDs_ , sumOp<labelList>());

        label szMonitorCells = monitorCellSetIDs_.size();
        reduce(szMonitorCells, sumOp<label>());

        if (szMonitorCells == 0)
        {
            FatalErrorInFunction
                << "No cell is available for incoming velocity monitoring in the cellSet: "
                << monitorCellSetName_
                << exit(FatalError);
        }

        scalar totalVofMonitoredCellSet_ = 0.0; 

        for (const label i : monitorCellSetIDs_)
        {
            Udisk += U[i]*cellsV[i];
            rhoRef_ = rhoRef_ + rho[i]*cellsV[i];
            totalVofMonitoredCellSet_ += cellsV[i];
        }
        reduce(Udisk, sumOp<vector>());
        reduce(rhoRef_, sumOp<scalar>());
        reduce(totalVofMonitoredCellSet_, sumOp<scalar>());
        
        Udisk /= totalVofMonitoredCellSet_;
        rhoRef_ /= totalVofMonitoredCellSet_;
        
        Umonitor_ = Udisk & diskDir_;
    }
    else
    {
        FatalErrorInFunction
        << "Unavailable monitor method requested"
        << monitorMethod_
        << exit(FatalError);
    }
}

void Foam::fv::calibratedActuatorDisk::getHubCenter
(
    const labelList& cells_,
    const scalarField& cellsV
)
{
    //static bool computed = false;
    if (!computed_)
    {  
        const Field<vector> zoneCellCentres(mesh().cellCentres(), cells_);
        const Field<scalar> zoneCellVolumes(mesh().cellVolumes(), cells_);

        avgCentre = gSum(zoneCellVolumes * zoneCellCentres) / V();
        maxR = gMax(mag(zoneCellCentres - avgCentre));

        totalV = 0.0; 

        for (const label celli : cells_)
        {
            totalV += cellsV[celli];
        }
        reduce(totalV, sumOp<scalar>());

        Info << "The actuator disk: " << this->name() << endl; 
        Info << "Total cellSet volume is : " << totalV << endl;
        Info << "The cellSet center is : " << avgCentre << endl;
        Info << "The max radius is : " << maxR << endl; 
        Info << "Monitored method applied : " << monitorMethod_ << endl;
        
        if (monitorMethod_ == "upstreamPoint")
        {
            Info << "Point is : " << monitorPoint_ << endl;
        }
        else if (monitorMethod_ == "upstreamCellSet")
        {
            Info << "Cellzone is : " << monitorCellSetName_ << endl;
        }
        computed_ = true; // Prevent recomputation
    }
}


void Foam::fv::calibratedActuatorDisk::writeDataToOutputFile()
{
    if
    (
        mesh_.time().timeOutputValue() >= writeFileStart_
     && mesh_.time().timeOutputValue() <= writeFileEnd_
    )
    {
        Ostream& os = file();
        writeCurrentTime(os);
        os.precision(3);
        os << Uref_ << tab << tab << Ct_ref_ << tab << Cp_ref_ << tab
         << Thrust_ << tab << tab << magSource_X << tab  << errorRatioT << tab;   
        if (rotSpeed_ != 0.0)
        {
        os << Power_ << tab << magSource_YZ << tab << errorRatioP << tab;
        }        
        os << endl;
    }
}


void Foam::fv::calibratedActuatorDisk::applyingSourceTerm
(
    const labelList& cells_,
    const scalarField& cellsV,
    vectorField& Usource
)
{
    switch (distributionMethod_)
    {
        case distributionMethodType::uniformlyDistributed_:
        {
            magSource_X = 0.0;
            magSource_YZ = 0.0;
            errorRatioT = 0.0;
            errorRatioP = 0.0;
            
            Thrust_ = 0.5*rhoRef_*diskArea_*pow(Uref_,2)*Ct_ref_;

            if (rotSpeed_ != 0.0)
            {
                Power_ = 0.5*rhoRef_*diskArea_*pow3(mag(Uref_))*Cp_ref_;
                
                for (const label celli : cells_)
                {
                    //-Computing local position of the cell
                    const scalar rPos = mag(mesh().cellCentres()[celli] - avgCentre);
                    //Calc of distribution for tangential (rotational) force
                    scalar tRatio = (rPos*(rPos-maxR))/pow(maxR,2);
                    //Calc of positional vector and its direction
                    vector vecPos = mesh().cellCentres()[celli] - avgCentre;
                    scalar dY = vecPos.y()/rPos;
                    scalar dZ = vecPos.z()/rPos;

                    Usource[celli].x() += (Thrust_) *(cellsV[celli]/totalV);
                    Usource[celli].y() += ((Power_/(rotSpeed_)) * tRatio *(cellsV[celli]/totalV)) * dZ;
                    Usource[celli].z() += - ((Power_/(rotSpeed_)) * tRatio *(cellsV[celli]/totalV)) * dY;
                    magSource_X += mag(Usource[celli].x());
                    magSource_YZ += mag(Usource[celli].y())  + mag(Usource[celli].z());
                }    
                reduce(magSource_X, sumOp<scalar>());
                reduce(magSource_YZ, sumOp<scalar>()); 
                errorRatioT = (1.0 - (Thrust_  / (magSource_X + 0.0001)));
                errorRatioP = (1.0 - (Power_  / (magSource_YZ + 0.0001)));
                
                writeDataToOutputFile();
            }
            else if (rotSpeed_ == 0.0)
            {
                for (const label celli : cells_)
                {
                    Usource[celli] += (Thrust_) *(cellsV[celli]/totalV)*diskDir_;
                    magSource_X += mag(Usource[celli] & diskDir_);
                } 
                reduce(magSource_X, sumOp<scalar>());
                errorRatioT = (1.0 - (Thrust_  / (magSource_X + 0.0001)))*100;
                
                writeDataToOutputFile();
            }  
            break;
        }

        case distributionMethodType::NRELdistribution_:
        {
            magSource_X = 0.0;
            magSource_YZ = 0.0;
            errorRatioT = 0.0;
            errorRatioP = 0.0;

            Thrust_ = 0.5*rhoRef_*diskArea_*pow(Uref_,2)*Ct_ref_;
            
            if (rotSpeed_ != 0.0)
            {
                Power_ = 0.5*rhoRef_*diskArea_*pow3(mag(Uref_))*Cp_ref_;

                for (const label celli : cells_)
                {
                    //-Computing local position of the cell
                    const scalar rPos = mag(mesh().cellCentres()[celli] - avgCentre);
                    //Calc of distribution for tangential (rotational) force
                    scalar tRatio = (rPos*(rPos-maxR))/pow(maxR,2);
                    //Calc of positional vector and its direction
                    vector vecPos = mesh().cellCentres()[celli] - avgCentre;
                    scalar dY = vecPos.y()/rPos;
                    scalar dZ = vecPos.z()/rPos;

                    //-Calc of local distribution for thrust - ORIGINAL NREL distribution
                    const scalar rRatio3 = -(27.0/4.0) * (pow(rPos,2.0)*(rPos-maxR))/(pow(maxR,3.0));
                    const scalar rRatio9 = -(387420489.0/16777216.0)*(pow(rPos,8.0)*(rPos-maxR))/(pow(maxR,9.0));
                    const scalar rRatio = 1.52*1.19*(0.5*rRatio3 + 0.5*rRatio9);

                    Usource[celli].x() += Thrust_*rRatio*(cellsV[celli]/totalV);
                    Usource[celli].y() += ((Power_/(rotSpeed_)) * tRatio *(cellsV[celli]/totalV)) * dZ;
                    Usource[celli].z() += - ((Power_/(rotSpeed_)) * tRatio *(cellsV[celli]/totalV)) * dY;

                    magSource_X += mag(Usource[celli].x());
                    magSource_YZ += mag(Usource[celli].y())  + mag(Usource[celli].z());
                }
                
                reduce(magSource_X, sumOp<scalar>());
                reduce(magSource_YZ, sumOp<scalar>());

                errorRatioT = (1.0 - (Thrust_  / (magSource_X + 0.0001)))*100;
                errorRatioP = (1.0 - (Power_  / (magSource_YZ + 0.0001)))*100;
                
                writeDataToOutputFile();
            }
            else if (rotSpeed_ == 0.0)
            {                
                for (const label celli : cells_)
                {
                    //-Computing local position of the cell
                    const scalar rPos = mag(mesh().cellCentres()[celli] - avgCentre);
                    //-Calc of local distribution for thrust - ORIGINAL NREL distribution
                    const scalar rRatio3 = -(27.0/4.0) * (pow(rPos,2.0)*(rPos-maxR))/(pow(maxR,3.0));
                    const scalar rRatio9 = -(387420489.0/16777216.0)*(pow(rPos,8.0)*(rPos-maxR))/(pow(maxR,9.0));
                    const scalar rRatio = 1.52*1.19*(0.5*rRatio3 + 0.5*rRatio9);
                    
                    Usource[celli] += (Thrust_)*rRatio *(cellsV[celli]/totalV)*diskDir_;
                    magSource_X += mag(Usource[celli] & diskDir_);
                }
                reduce(magSource_X, sumOp<scalar>());
                errorRatioT = (1.0 - (Thrust_  / (magSource_X + 0.0001)))*100;

                writeDataToOutputFile();
            }  
            break;
        }
        default:
            break;
    }
} 


// * * * * * * * * * * * * *                    * * * * * * * * * * * //


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::calibratedActuatorDisk::calcConstantParameters
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();

    getHubCenter(cells_, cellsV);

    if (V() < SMALL)
    {
        FatalErrorInFunction
            << "No cells in the cellSet/cellZone of the actuator disk."
            << exit(FatalError);
    }

    // This case, from reading dictionary is updated the values of
    //      Uref_, Ct_ref_ and Cp_ref which are constants thorought the simulation

    rhoRef_ = 0.0;
    for (const label celli : cells_)
    {
        rhoRef_ += rho[celli]*cellsV[celli];
    }
    reduce(rhoRef_, sumOp<scalar>());
    rhoRef_ /= totalV; 

    applyingSourceTerm(cells_, cellsV, Usource);
    
}



template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::calibratedActuatorDisk::calcConstantCt
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    const vectorField& U = eqn.psi();
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();

    getHubCenter(cells_, cellsV);

    if (V() < SMALL)
    {
        FatalErrorInFunction
            << "No cells in the cellSet/cellZone of the actuator disk."
            << exit(FatalError);
    }
    if (monitorMethod_ == "actuatorDisk")
    {
        FatalErrorInFunction
        << "This monitoring method is not suitable for the current thrust calculation"
        << exit(FatalError);
    }
    else if (monitorMethod_ == "upstreamPoint" || monitorMethod_ == "upstreamCellSet")
    {
        monitorData(alpha, rho, U, cellsV);
    }
    else
    {
        FatalErrorInFunction
        << "Wrong name entered of forceMethod in dictionary."
        << exit(FatalError);
    }    
    Uref_ = Umonitor_;


    // This case, from reading dictionary is updated the values of
    //      Ct_ref_ and Cp_ref which are constants thorought the simulation

    applyingSourceTerm(cells_, cellsV, Usource);

}



template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::calibratedActuatorDisk::calcCalafMethod 
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    const vectorField& U = eqn.psi();
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();

    getHubCenter(cells_, cellsV);

    if (V() < SMALL)
    {
        FatalErrorInFunction
            << "No cells in the cellSet/cellZone of the actuator disk."
            << exit(FatalError);
    }
    
    // Forcing readings on actuator disk cellset only
    monitorMethod_ = "actuatorDisk";
    monitorData(alpha, rho, U, cellsV);

    //Info << "Velocity at monitored region: " << Umonitor_ << endl;
    
    Uref_ = Umonitor_;
    if (Ct_star_ < VSMALL)
    {
        FatalErrorInFunction
            << "Scaled thrust coefficient is equal or below zero"
            << "Please check dictionary "
        << exit(FatalError);
    }

    Ct_ref_ = Ct_star_;

    if (rotSpeed_ != 0.0)
    {
        if (Cp_star_ < VSMALL)
        {
            FatalErrorInFunction
                << "Scaled power coefficient is equal or below zero"
                << "Please check dictionary "
            << exit(FatalError);
        }    
        Cp_ref_ = Cp_star_;
    }

    // This case, from reading dictionary Ct_star, it's stored as Ct_ref_
    
    applyingSourceTerm(cells_, cellsV, Usource);
}



template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::calibratedActuatorDisk::calcInterpolatedFromTable
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn
)
{
    const vectorField& U = eqn.psi();
    vectorField& Usource = eqn.source();
    const scalarField& cellsV = mesh_.V();

    getHubCenter(cells_, cellsV);

    if (V() < SMALL)
    {
        FatalErrorInFunction
            << "No cells in the cellSet/cellZone of the actuator disk."
            << exit(FatalError);
    }

    if (monitorMethod_ == "actuatorDisk")
    {
        monitorData(alpha, rho, U, cellsV);

        if (!csvFile_.empty())
        {
            Foam::myLookupTable4Turbines lookup(csvFile_);
            if(!lookup.findByUad(Umonitor_, Uref_, Ct_ref_, Cp_ref_))
            {
                Foam::Info << "No valid interpolation found for Umonitor_ = " << Umonitor_ << Foam::nl;
                Foam::Info << "Actuatordisk: " << this->name() << " turned off." << Foam::nl;
                Uref_ = 0.0;
                Ct_ref_ = 0.0;
                Cp_ref_ = 0.0;
            }

        }
        else
        {
            FatalErrorInFunction
            << "-   pathToCsvTable entered not valid (check path, format and content)."
            << exit(FatalError);
        }
    }
    else if (monitorMethod_ == "upstreamPoint")
    {
        monitorData(alpha, rho, U, cellsV);

        if (!csvFile_.empty())
        {
            Foam::myLookupTable4Turbines lookup(csvFile_);
            if(!lookup.findByU(Umonitor_, Ct_ref_, Cp_ref_))
            {
                Foam::Info << "No valid interpolation found for Umonitor_ = " << Umonitor_ << Foam::nl;
                Foam::Info << "Actuatordisk: " << this->name() << " turned off." << Foam::nl;
                Uref_ = 0.0;
                Ct_ref_ = 0.0;
                Cp_ref_ = 0.0;
            }
            Uref_ = Umonitor_;

        }
        else
        {
            FatalErrorInFunction
            << "-   pathToCsvTable entered not valid (check path, format and content)."
            << exit(FatalError);
        }
    }
    else if (monitorMethod_ == "upstreamCellSet")
    {
        monitorData(alpha, rho, U, cellsV);

        if (!csvFile_.empty())
        {
            Foam::myLookupTable4Turbines lookup(csvFile_);
            if(!lookup.findByU(Umonitor_, Ct_ref_, Cp_ref_))
            {
                Foam::Info << "No valid interpolation found for Umonitor_ = " << Umonitor_ << Foam::nl;
                Foam::Info << "Actuatordisk: " << this->name() << " turned off." << Foam::nl;
                Uref_ = 0.0;
                Ct_ref_ = 0.0;
                Cp_ref_ = 0.0;
            }
            Uref_ = Umonitor_;

        }
        else
        {
            FatalErrorInFunction
            << "-   pathToCsvTable entered not valid (check path, format and content)."
            << exit(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
        << "Wrong name entered of forceMethod in dictionary."
        << exit(FatalError);
    }  

    applyingSourceTerm(cells_, cellsV, Usource);
}