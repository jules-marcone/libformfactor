//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      lib/Vectors3D.h
//! @brief     Defines basic vectors in Z^3, R^3, C^3.
//!
//! @homepage  http://www.bornagainproject.org
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2018
//! @authors   Scientific Computing Group at MLZ (see CITATION, AUTHORS)
//
//  ************************************************************************************************

#ifndef LIBFORMFACTOR_LIB_VECTORS3D_H
#define LIBFORMFACTOR_LIB_VECTORS3D_H

#include "Base/Vector/BasicVector3D.h"

using I3 = BasicVector3D<int>;
using R3 = BasicVector3D<double>;
using C3 = BasicVector3D<complex_t>;

#endif // LIBFORMFACTOR_LIB_VECTORS3D_H
