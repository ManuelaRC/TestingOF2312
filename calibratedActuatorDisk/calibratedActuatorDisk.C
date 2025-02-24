/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd
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
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(calibratedActuatorDisk, 0);
    addToRunTimeSelectionTable(option, calibratedActuatorDisk, dictionary);
}
}


const Foam::Enum
<
    Foam::fv::calibratedActuatorDisk::forceMethodType
>
Foam::fv::calibratedActuatorDisk::forceMethodTypeNames
({
    { forceMethodType::constantParameters_, "constantParameters" },
    { forceMethodType::constantCt_, "constantCt" },
    { forceMethodType::calafMethod_, "calafMethod" },
    { forceMethodType::interpolatedFromTable, "interpolatedFromTable" }
});

const Foam::Enum
<
    Foam::fv::calibratedActuatorDisk::distributionMethodType
>
Foam::fv::calibratedActuatorDisk::distributionMethodTypeNames
({
    { distributionMethodType::uniformlyDistributed_, "uniformlyDistributed" },
    { distributionMethodType::NRELdistribution_, "NRELdistribution" }
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::calibratedActuatorDisk::writeFileHeader(Ostream& os)
{
    writeFile::writeHeader(os, "Calibrated actuator disk source");
    writeFile::writeHeaderValue(os, "Turbine ID:                ", this->name());
    writeFile::writeHeaderValue(os, "Thrust method:             ", forceMethodTypeNames[forceMethod_]);
    writeFile::writeHeaderValue(os, "Monitor reading method:    ", monitorMethod_);
    writeFile::writeHeaderValue(os, "Rotation applied:          ", rotSpeed_);
    writeFile::writeHeaderValue(os, "Load distribution:         ", distributionMethodTypeNames[distributionMethod_]);
    if (monitorMethod_ == "upstreamPoint")
    {
        writeFile::writeHeaderValue(os, "Monitor point:         ", monitorPoint_);
    }
    if (monitorMethod_ == "upstreamCellSet")
    {
        writeFile::writeHeaderValue(os, "Monitor cellSet:         ", monitorCellSetName_);
    }
    writeFile::writeHeader(os, " ");
    writeFile::writeCommented(os, "Time");
    writeFile::writeCommented(os, "Uref");
    writeFile::writeCommented(os, "Ctref");
    writeFile::writeCommented(os, "Cpref");
    writeFile::writeCommented(os, "Thr_Calc");
    writeFile::writeCommented(os, "Thr_Applied");
    writeFile::writeCommented(os, "errorT (%)");
    if (rotSpeed_ != 0.0)
    {
        writeFile::writeCommented(os, "RotF_Calc");
        writeFile::writeCommented(os, "RotF_Applied");
        writeFile::writeCommented(os, "errorP (%)");
    }
    writeFile::writeHeader(os, " ");
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::calibratedActuatorDisk::calibratedActuatorDisk
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    writeFile(mesh, name, modelType, coeffs_),
    forceMethod_
    (
        forceMethodTypeNames.getOrDefault
        (
            "forceMethod",
            coeffs_,
            forceMethodType::constantParameters_
        )
    ),
    distributionMethod_
    (
        distributionMethodTypeNames.getOrDefault
        (
            "loadDistribution",
            coeffs_,
            distributionMethodType::uniformlyDistributed_
        )
    ),
    writeFileStart_(coeffs_.getOrDefault<scalar>("writeFileStart", 0)),
    writeFileEnd_(coeffs_.getOrDefault<scalar>("writeFileEnd", VGREAT)),
    diskArea_
    (
        coeffs_.getCheck<scalar>
        (
            "diskArea",
            scalarMinMax::ge(VSMALL)
        )
    ),
    diskDir_
    (
        coeffs_.getCheck<vector>
        (
            "diskDir",
            [&](const vector& vec){ return mag(vec) > VSMALL; }
        ).normalise()
    ),
    Uref_
    (
        coeffs_.getOrDefault<scalar>
        (
            "U_inf",
            0.0
        )
    ),
    Ct_ref_
    (
        coeffs_.getOrDefault<scalar>
        (
            "Ct_inf",
            0.0
        )
    ),
    Cp_ref_
    (
        coeffs_.getOrDefault<scalar>
        (
            "Cp_inf",
            0.0
        )
    ),
    rotSpeed_
    (
        coeffs_.getOrDefault<scalar>
        (
            "rotSpeed",
            0.0
        )
    ),
    monitorMethod_
    (
        coeffs_.getOrDefault<word>
        (
            "monitorMethod",
            "actuatorDisk"
        )
    ),
    monitorPoint_
    (
        coeffs_.getOrDefault<point>
        (
            "upstreamPoint",
            vector::uniform(NAN)
        )
    ),
    monitorCellSetName_
    (
        coeffs_.getOrDefault<word>
        (
            "upstreamCellSet",
            "None"
        )
    ),
    Ct_star_
    (
        coeffs_.getOrDefault<scalar>
        (
            "Ct_star",
            0.0
        )
    ),
    Cp_star_
    (
        coeffs_.getOrDefault<scalar>
        (
            "Cp_star",
            0.0
        )
    ),
    csvFile_
    (
        coeffs_.getOrDefault<string>
        (
            "pathToCsvTable",
            "None"
        )
    ),
    computed_(false)
{
    fieldNames_.resize(1, "U");

    fv::option::resetApplied();

    Info<< "    - creating actuation disk zone: " << this->name() << endl;

    Info<< "    - force computation method: "
        << forceMethodTypeNames[forceMethod_]  << " with load distribution: " << distributionMethodTypeNames[distributionMethod_] << endl;

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::calibratedActuatorDisk::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (V() > VSMALL)
    {
        calc(geometricOneField(), geometricOneField(), eqn);
    }
}


void Foam::fv::calibratedActuatorDisk::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (V() > VSMALL)
    {
        calc(geometricOneField(), rho, eqn);
    }
}


void Foam::fv::calibratedActuatorDisk::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (V() > VSMALL)
    {
        calc(alpha, rho, eqn);
    }
}


bool Foam::fv::calibratedActuatorDisk::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict) && writeFile::read(dict))
    {
        dict.readIfPresent("pathToCsvTable", csvFile_);
        dict.readIfPresent("monitorMethod", monitorMethod_);
        dict.readIfPresent("upstreamPoint", monitorPoint_);
        dict.readIfPresent("upstreamCellSet", monitorCellSetName_);
        dict.readIfPresent("writeFileStart", writeFileStart_);
        dict.readIfPresent("writeFileEnd", writeFileEnd_);
        dict.readIfPresent("diskArea", diskArea_);
        if (diskArea_ < VSMALL)
        {
            FatalIOErrorInFunction(dict)
                << "Actuator disk has zero area: "
                << "diskArea = " << diskArea_
                << exit(FatalIOError);
        }
        dict.readIfPresent("diskDir", diskDir_);
        diskDir_.normalise();
        if (mag(diskDir_) < VSMALL)
        {
            FatalIOErrorInFunction(dict)
                << "Actuator disk surface-normal vector is zero: "
                << "diskDir = " << diskDir_
                << exit(FatalIOError);
        }
        return true;
    }

    return false;
}


// ************************************************************************* //
