/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
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

Application
    myLookupTable4Turbines

Description
    Tool used to read '.csv' files or any that uses commas as value separators
    skipping the first line (header). This tool is invoked in the fvOptions
    source 'calibratedActuatorDisk' to read from the characteristics of the 
    wind turbine performances curve and provide input data for the thrust
    calculation according to local readings.
    Inner structure of the ".csv" file is n-rows and 4 columns:
            0      1      2     3
        0   Uinf, Ctinf, Cpinf, Uad
        1   float, float, float, float
            .
            .
            .
        n   float, float, float, float
    Where:
        Uinf: reference velocity in [m/s]
        Ctinf: reference thrust coefficient [-] @ Uinf
        Cpinf: reference power coefficient [-] @ Uinf
        Uad: average wind velocity at cellset of the actuator disk [m/s]
    IMPORTANT:
    When searching Uad
    1) If the lowest value from .csv of either Uinf or Uad are below the conditions
     for 3 m/s, Uad and Uinf are set to 0 m/s 
    2) If the lowest value from .csv of either Uinf or Uad are below the conditions
     for 3 m/s, Uad and Uinf are set to 0 m/s 

\*---------------------------------------------------------------------------*/
#ifndef MYLOOKUPTABLE4TURBINES_H
#define MYLOOKUPTABLE4TURBINES_H

#include "scalar.H"
#include "vector.H"
#include "HashTable.H"
#include "IFstream.H"
#include "IOdictionary.H"
#include "word.H"
#include <vector>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myLookupTable4Turbines
{
private:
    std::vector<scalar> U, ct, cp, Uad;

    // Helper function for linear interpolation
    scalar interpolate(scalar x1, scalar x2, scalar y1, scalar y2, scalar x) const;

public:
    //myLookupTable4Turbines(const word& fileName);  // Constructor loads CSV
    myLookupTable4Turbines(const std::string& stdFileName);  // Constructor loads CSV

    // Search functions
    bool findByUad(scalar UadValue, scalar& Uout, scalar& ctout, scalar& cpout) const;
    bool findByU(scalar UValue, scalar& ctout, scalar& cpout) const;
};

} // namespace Foam

#endif



// ************************************************************************* //
