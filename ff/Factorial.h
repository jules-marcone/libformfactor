//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      lib/Factorial.h
//! @brief     Precomputes 1/n!
//!
//! @homepage  http://www.bornagainproject.org
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2021
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

#ifdef SWIG
#error no need to expose this header to Swig
#endif

#ifndef USER_API
#ifndef LIBFORMFACTOR_LIB_FACTORIAL_H
#define LIBFORMFACTOR_LIB_FACTORIAL_H

#include <array>
#include <cstddef>
#include <utility>

namespace {

template <size_t N> struct ReciprocalFactorial {
    static constexpr double value = ReciprocalFactorial<N - 1>::value / N;
};

template <> struct ReciprocalFactorial<0> {
    static constexpr double value = 1.0;
};

template <template <size_t> class F, size_t... I>
constexpr std::array<double, sizeof...(I)> generateArrayHelper(std::index_sequence<I...>)
{
    return {F<I>::value...};
}

} // namespace


namespace ff_aux {

//! Returns a compile-time generated std::array of reciprocal factorials.

template <size_t N, typename Indices = std::make_index_sequence<N>>
constexpr std::array<double, N> generateReciprocalFactorialArray()
{
    return generateArrayHelper<ReciprocalFactorial>(Indices{});
}

} // namespace ff_aux

#endif // LIBFORMFACTOR_LIB_FACTORIAL_H
#endif // USER_API
