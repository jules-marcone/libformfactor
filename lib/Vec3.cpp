//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      lib/Vec3.cpp
//! @brief     Implements type-specific functions from template class Vec3.
//!
//! @homepage  http://www.bornagainproject.org
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2021
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

#include "Vec3.h"
#include <stdexcept>

//! Returns complex conjugate vector
template <> Vec3<double> Vec3<double>::conj() const
{
    return *this;
}

template <> Vec3<complex_t> Vec3<complex_t>::conj() const
{
    return {std::conj(x()), std::conj(y()), std::conj(z())};
}

//! Returns this, trivially converted to complex type.
template <> Vec3<complex_t> Vec3<double>::complex() const
{
    return {x(), y(), z()};
}

//! Returns real parts.
template <> Vec3<double> Vec3<double>::real() const
{
    return *this;
}

template <> Vec3<double> Vec3<complex_t>::real() const
{
    return {x().real(), y().real(), z().real()};
}

//! Returns unit vector in direction of this. Throws for null vector.
template <> Vec3<double> Vec3<double>::unit() const
{
    double len = mag();
    if (len == 0.0)
        throw std::runtime_error("Cannot normalize zero vector");
    return {x() / len, y() / len, z() / len};
}

template <> Vec3<complex_t> Vec3<complex_t>::unit() const
{
    double len = mag();
    if (len == 0.0)
        throw std::runtime_error("Cannot normalize zero vector");
    return {x() / len, y() / len, z() / len};
}
