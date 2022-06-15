//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      demo/tribipy.cpp
//! @brief     Computes the formfactor of a triangular bipyramid
//!
//! @homepage  https://github.com/jules-marcone/libformfactor
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @author    Jules Marcone, June 2022
//
//  ************************************************************************************************

#include <ff/Tri.h>
#include <iostream>

constexpr double twopi = 6.28318530718;

//! Prints list t vs |F(q(t))| for a logarithmic range of t values

int main()
{
    std::cout << "# Tribipy form factor, for different q scans\n";
    ff::tri::TriangularBipyramid TriangularBipyramid(1.);

    std::cout << "# q vs |F(q)| for q in direction 111, perpendicular to two faces\n";
    for (double t = 0.2; t < 200; t *= 1.002) {
        C3 q(t/sqrt(3), t/sqrt(3), t/sqrt(3));
        std::cout << t << " " << std::abs(TriangularBipyramid.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "# q vs |F(q)| for q in direction 110, perpendicular to two edges\n";
    for (double t = 0.2; t < 200; t *= 1.002) {
        C3 q(t/sqrt(2), t/sqrt(2), 0);
        std::cout << t << " " << std::abs(TriangularBipyramid.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;

}
