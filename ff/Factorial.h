//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      ff/Factorial.h
//! @brief     Precomputes 1/n!
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

#ifndef FORMFACTOR_FF_FACTORIAL_H
#define FORMFACTOR_FF_FACTORIAL_H

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

#endif // FORMFACTOR_FF_FACTORIAL_H
