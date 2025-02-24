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

\*---------------------------------------------------------------------------*/

#include "myLookupTable4Turbines.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "ListOps.H"
#include <algorithm>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

// Utility function to split a string by a given delimiter (comma for CSV)
std::vector<std::string> split(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;

    while (std::getline(ss, token, delimiter))
    {
        tokens.push_back(token);
    }

    return tokens;
}

// Constructor to read the CSV file while skipping the first row (header)
//myLookupTable4Turbines::myLookupTable4Turbines(const word& fileName)
myLookupTable4Turbines::myLookupTable4Turbines(const std::string& stdFileName)
{
    // Convert std::string to Foam::word for OpenFOAM compatibility
    Foam::word fileName(stdFileName);

    IFstream file(fileName);
    if (!file.good())
    {
        FatalErrorIn("myLookupTable4Turbines::myLookupTable4Turbines()")
            << "Cannot open file: " << fileName << exit(FatalError);
    }

    // Skip the first line (header row)
    string headerLine;
    file.getLine(headerLine);

    string line;
    while (file.getLine(line)) // Read each line
    {
        std::vector<std::string> values = split(line, ','); // Split by comma

        if (values.size() != 4) // Ensure correct number of columns
        {
            FatalErrorIn("myLookupTable4Turbines::myLookupTable4Turbines()")
                << "Invalid row in CSV file: " << line << exit(FatalError);
        }

        // Convert string values to scalars
        scalar u = Foam::readScalar(IStringStream(values[0])());
        scalar ctVal = Foam::readScalar(IStringStream(values[1])());
        scalar cpVal = Foam::readScalar(IStringStream(values[2])());
        scalar uad = Foam::readScalar(IStringStream(values[3])());

        // Store values
        U.push_back(u);
        ct.push_back(ctVal);
        cp.push_back(cpVal);
        Uad.push_back(uad);
    }
}

// Helper function for linear interpolation 
scalar myLookupTable4Turbines::interpolate(scalar x1, scalar x2, scalar y1, scalar y2, scalar x) const
{
    if (x1 == x2) return y1; // Avoid division by zero
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

// Search by Uad and interpolate
bool myLookupTable4Turbines::findByUad(scalar UadValue, scalar& Uout, scalar& ctout, scalar& cpout) const
{
    if (Uad.empty()) return false;

    // Check first if the monitored values is within the limits of the content of the .csv file:
    label sizeOfData = Uad.size();
    if (Uad[0] >= UadValue)
    {
        Uout = interpolate(Uad[0], Uad[1], U[0], U[1], UadValue);
        
        if (Uout < 4.0)
        {
            // Non-working conditions for a turbine
            Uout = 0.0;
            ctout = 0.0;
            cpout = 0.0;
            Info << "Non-operable conditions detected on the turbine" << endl; 
        }
        else
        {
            ctout = interpolate(Uad[0], Uad[1], ct[0], ct[1], UadValue);
            cpout = interpolate(Uad[0], Uad[1], cp[0], cp[1], UadValue);
        }

        return true;
    }
    else if (Uad[sizeOfData -1] <= UadValue)
    {
        Uout = interpolate(Uad[sizeOfData -2], Uad[sizeOfData -1], U[sizeOfData -2], U[sizeOfData -1], UadValue);
        if (Uout > 25.1)
        {
            // Non-working conditions for a turbine
            Uout = 0.0;
            ctout = 0.0;
            cpout = 0.0;
            Info << "Non-operable conditions detected on the turbine" << endl; 
        }
        else
        {
            ctout = interpolate(Uad[sizeOfData -2], Uad[sizeOfData -1], ct[sizeOfData -2], ct[sizeOfData -1], UadValue);
            cpout = interpolate(Uad[sizeOfData -2], Uad[sizeOfData -1], cp[sizeOfData -2], cp[sizeOfData -1], UadValue);
        }
        return true;
    }
    else if (Uad[0] <= UadValue && Uad[sizeOfData -1] >= UadValue)
    {
        for (size_t i = 0; i < Uad.size() - 1; ++i)
        {
            Uout = interpolate(Uad[i], Uad[i + 1], U[i], U[i + 1], UadValue);
            ctout = interpolate(Uad[i], Uad[i + 1], ct[i], ct[i + 1], UadValue);
            cpout = interpolate(Uad[i], Uad[i + 1], cp[i], cp[i + 1], UadValue);
            return true;
        }
    }
    return false; // No match found
}

// Search by U and interpolate
bool myLookupTable4Turbines::findByU(scalar UValue, scalar& ctout, scalar& cpout) const
{
    if (U.empty()) return false;
    // Check first if the monitored values is within the limits of the content of the .csv file:
    label sizeOfData = Uad.size();

    if (UValue <= U[0])
    {
        if (UValue < 4.0)
        {
            // Non-working conditions for a turbine
            ctout = 0.0;
            cpout = 0.0;
            Info << "Non-operable conditions detected on the turbine" << endl; 
        }
        else
        {
            ctout = interpolate(U[0], U[1], ct[0], ct[1], UValue);
            cpout = interpolate(U[0], U[1], cp[0], cp[1], UValue);
        }
        return true;
    }
    else if (U[sizeOfData -1] <= UValue)
    {
        if (UValue > 25.1)
        {
            // Non-working conditions for a turbine
            ctout = 0.0;
            cpout = 0.0;
            Info << "Non-operable conditions detected on the turbine" << endl; 
        }
        else
        {
            ctout = interpolate(U[sizeOfData -2], U[sizeOfData -1], ct[sizeOfData -2], ct[sizeOfData -1], UValue);
            cpout = interpolate(U[sizeOfData -2], U[sizeOfData -1], cp[sizeOfData -2], cp[sizeOfData -1], UValue);
        }
        return true;
    }

    
    else if (U[0] <= UValue && U[sizeOfData -1] >= UValue)
    {
        for (size_t i = 0; i < U.size() - 1; ++i)
        {
            ctout = interpolate(U[i], U[i + 1], ct[i], ct[i + 1], UValue);
            cpout = interpolate(U[i], U[i + 1], cp[i], cp[i + 1], UValue);
            return true;
        }
    }

    return false; // No match found
}

} // namespace Foam
