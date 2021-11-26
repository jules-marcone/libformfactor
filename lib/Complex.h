//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      lib/Complex.h
//! @brief     Defines complex_t, and a few elementary functions
//!
//! @homepage  http://www.bornagainproject.org
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2021
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

#ifndef LIBFORMFACTOR_LIB_COMPLEX_H
#define LIBFORMFACTOR_LIB_COMPLEX_H

#include <complex>

using complex_t = std::complex<double>;
constexpr complex_t I = complex_t(0.0, 1.0);

//! Returns product I*z, where I is the imaginary unit.
inline complex_t mul_I(complex_t z)
{
    return complex_t(-z.imag(), z.real());
}

//! Returns exp(I*z), where I is the imaginary unit.
inline complex_t exp_I(complex_t z)
{
    return std::exp(complex_t(-z.imag(), z.real()));
}

#endif // LIBFORMFACTOR_LIB_COMPLEX_H