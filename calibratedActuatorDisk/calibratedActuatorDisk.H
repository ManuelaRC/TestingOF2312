/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Built-by:
    María Manuela Rosales (https://www.linkedin.com/in/mmrosales/)
    February 2025
    version: 0.0

    Class
    Foam::fv::calibratedActuatorDisk

Group
    grpFvOptionsSources

Description
    
    User-built explicit source term to model the effects of wind turbines inspired 
        in the openfoam class actuationDiskSource. This class looks to provide a single
        tool to run     single-wake wind turbine and miltorow wind farm simulations with
        known or unknown working conditions.

    To better understand the purposes the several configurations, it's recommended
        to review the following literature:

        [Calaf] Marc Calaf, Charles Meneveau, Johan Meyers; Large eddy simulation study
            of fully developed wind-turbine array boundary layers. Physics of Fluids 1 
            22 (1): 015110. https://doi.org/10.1063/1.3291077
        
        [Réthoré] Réthoré, P.-E., van der Laan, P., Troldborg, N., Zahle, F. and Sørensen,
            N.N. (2014), Verification and validation of an actuator disc model. Wind Energ.,
            17: 919-937. https://doi.org/10.1002/we.1607
        
        [PVDL-1] van der Laan, M. P., Sørensen, N. N., Réthoré, P.-E., Mann, J., 
            Kelly, M. C., Troldborg, N., Schepers, J. G., and Machefaux, E. (2015),
            An improved k- ϵ model applied to a wind turbine wake in atmospheric
            turbulence. Wind Energ., 18: 889–907. doi: 10.1002/we.1736.
        
        [PVDL-2] van der Laan, M. P., Sørensen, N. N., Réthoré, P.-E., Mann, J.,
            Kelly, M. C., and Troldborg, N. (2015) The k-ε-fP model applied to double
            wind turbine wakes using different actuator disk force methods. 
            Wind Energ., 18: 2223–2240. doi: 10.1002/we.1816. 
        
        [Navarro] Gonzalo P. Navarro Diaz, Matias Avila and Arnau Folch (2018),
            An annual energy production estimation methodology for onshore wind farms
            over complex terrain using a RANS model with actuator discs. 
            Journal of Physics: Conference Series, Volume 1037, Issue 7. 
            doi: 10.1088/1742-6596/1037/7/072018


    The present class provides 4 methods to value the magnitude of the thrust 
        force based on the equation:  
        Thrust_ = 0.5 *rhoRef_ * diskArea * Uref_ ^2 * Ct_ref_

    1) constantParameters: [PVDL-1]user inputs in dictionary values of diskArea, Uref, Ct_ref

    2) constantCt: user inputs in dictionary values of diskArea, Ct_ref and chooses 
        to extract Uref_ either from an upstreamPoint or from a upstreamCellSet

    3) calafMethod: [Calaf][PVDL-2] user inputs in dictionary values of diskArea and Ct_star
        and thrust force magnitude bases on the average-volume disk streamwise velocity

    4) interpolatedFromTable: [PVDL-2][Navarro] user inputs in dictionary values of diskArea, the 
        monitoring method of choice (the current actuatordisk, upstreamPoint or 
        upstreamCellSet) and provides the path to an external '.csv' file which carries
        information of performance data of an isolated single-wake turbine: 
        Uref, Ctref, Cpref, Uad. The actuator disk will calibrate itself based on readings
        on the mesh and interpolation within data of the .csv file.
    
    The present class provides 2 methods to distribute the thrust force across the disk: 

    1) uniformlyDistributed: [Calaf][Navarro] uniform thrust force distribution
        across the radius
    2) NRELdistribution: [Réthoré][PVDL-1] predefined radial-loading based on the
        loads found in NREL 5 MW wind turbine at 8 m/s. 
    At those 2 distributions, it can be added the effects of a radial-rotating force
        if rotSpeed is differente from 0.00 (default value).

    Important:  
        1) This class is built to handle actuatordisk which disk-surface is orthogonal
        the streamwise direction in x-axis ONLY. This class does not handle another domain
        axis configuration. 
        
        2) The present class requires an auxiliar class user-built named
        'myLookupTable4Turbines', that should be found inside this class directory.
        If either upsteamPoit or upstreamCellSet are used monitoring methods, 
        the search in table data will only based on the Uref values listed.
        If actuatorDisk method is used, the search will be based on the column 
        Uad and extract data or interpolate: Uref, Ct and Cp.
        If wind turbines Uref is lower than 4.0 m/s or higher than 25.0 m/s,
        Uref will set to 0.00 since at these ranges wind turbines are not 
        operating. A message will be prompted to screen or to solver's log.  
    
    Due to the extent of possible configurations and amount of input data, an example file 
    of how to set the fvOption is provided inside the current's class directory.
    

\*---------------------------------------------------------------------------*/

#ifndef calibratedActuatorDisk_H
#define calibratedActuatorDisk_H

#include "cellSetOption.H"
#include "writeFile.H"
#include "cellSet.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class calibratedActuatorDisk Declaration
\*---------------------------------------------------------------------------*/

class calibratedActuatorDisk
:
    public fv::cellSetOption,
    public functionObjects::writeFile
{
protected:

    // Protected Enumerations

        //- Options for the force computation method types
        enum forceMethodType
        {
            constantParameters_,
            constantCt_,
            calafMethod_, 
            interpolatedFromTable
        };
        
        //- Names for forceMethodType
        static const Enum<forceMethodType> forceMethodTypeNames;
        
        //- Options for the force computation method types
        enum distributionMethodType
        {
            uniformlyDistributed_,
            NRELdistribution_
        };
        //- Names for distributionMethodType
        static const Enum<distributionMethodType> distributionMethodTypeNames;

    // Protected Data

        //- The type of the force computation method
        const enum forceMethodType forceMethod_;
        
        //- The type of incoming velocity monitoring method
        const enum distributionMethodType distributionMethod_;

        //- Start time for file output
        scalar writeFileStart_;

        //- End time for file output
        scalar writeFileEnd_;

        //- Actuator disk planar surface area [m2]
        scalar diskArea_;

        //- Surface-normal vector of the actuator disk pointing downstream
        vector diskDir_;

        //- Wind speed [m/s] predetermined or monitored used to compute thrust force 
        scalar  Uref_;

        //- Coefficients [-] predetermined or monitored used to compute thrust force 
        scalar  Ct_ref_;
        scalar  Cp_ref_;
        
        //- Rotary wind speed in [rad/s] used to compute rotary force 
        scalar  rotSpeed_;

        //- (if required) Method name to be assigned by the user and gather data of Uref_
        word    monitorMethod_; 
        
        //- sumplementary info for each method
        point   monitorPoint_;
        word    monitorCellSetName_; 

        //- Input data for forceMethod_: calafMethod_
        scalar  Ct_star_;
        scalar  Cp_star_;

        //- Path of the .csv file where to sort in table monitored values to compute thrust method
        string  csvFile_;

        //- Flag use to compute once geometrical features of the cellSet selected
        bool    computed_;
        vector  avgCentre;
        scalar  maxR;
        scalar  totalV;

        //- Variables to store monitored values
        scalar  Umonitor_; 
        scalar  rhoRef_;

        //- Variables used to store forces, sourceterm magnitudes and error based on the difference between these.
        scalar  Thrust_;
        scalar  Power_;
        scalar  magSource_X;
        scalar  magSource_YZ;
        scalar  errorRatioT;
        scalar  errorRatioP;
       

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

private:

    // Private Member Functions

        //- Locate the set of cells whereat the incoming velocity is monitored
        template<class AlphaFieldType, class RhoFieldType>
        void monitorData
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            const vectorField& U,
            const scalarField& cellsV
        );
       
        //- Get rotor's center, max. radius and total volume of the cellSet
        void getHubCenter
        (
            const labelList& cells_,
            const scalarField& cellsV
        );
        
        //- Writes data to the output file of the source term
        void writeDataToOutputFile();      
        
        //- Computes and apply the drag forces (axial and rotational) as sourceterm
        void applyingSourceTerm
        (
            const labelList& cells_,
            const scalarField& cellsV,
            vectorField& Usource
        );

        //- Compute the selected force computation method type
        template<class AlphaFieldType, class RhoFieldType>
        void calc
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn
        );

        //- Compute Froude's ideal actuator disk method
        template<class AlphaFieldType, class RhoFieldType>
        void calcConstantParameters
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn
        );

        //- Compute variable-scaling actuator disk method
        template<class AlphaFieldType, class RhoFieldType>
        void calcConstantCt
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn
        );
        //- Compute variable-scaling actuator disk method
        template<class AlphaFieldType, class RhoFieldType>
        void calcCalafMethod 
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn
        );

        //- Compute variable-scaling actuator disk method
        template<class AlphaFieldType, class RhoFieldType>
        void calcInterpolatedFromTable
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn
        );

public:

    //- Runtime type information
    TypeName("calibratedActuatorDisk");


    // Constructors

        //- No default construct
        calibratedActuatorDisk() = delete;

        //- Construct from components
        calibratedActuatorDisk
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        calibratedActuatorDisk(const calibratedActuatorDisk&) = delete;

        //- No copy assignment
        void operator=(const calibratedActuatorDisk&) = delete;


    //- Destructor
    virtual ~calibratedActuatorDisk() = default;


    // Member Functions

        // Access

            //- Normal disk direction
            const vector& diskDir() const
            {
                return diskDir_;
            }

            //- Disk area
            scalar diskArea() const
            {
                return diskArea_;
            }


         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to phase-compressible momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "calibratedActuatorDiskTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
