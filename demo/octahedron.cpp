//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      demo/octahedron.cpp
//! @brief     Computes the formfactor of a platonic octahedron
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

#include <ff/Platonic.h>
#include <iostream>

//! Prints list t vs |F(q(t))| for a logarithmic range of t values

int main()
{
    std::cout << "# Octahedral form factor, for different q scans\n";
    ff::platonic::Octahedron octahedron(1.);

    std::cout << "# q vs |F(q)| for q in direction 111, perpendicular to two faces\n";
    for (double t = 0.2; t < 200; t *= 1.002) {
        C3 q(t/sqrt(3), t/sqrt(3), t/sqrt(3));
        std::cout << t << " " << std::abs(octahedron.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "# q vs |F(q)| for q in direction 110, perpendicular to two edges\n";
    for (double t = 0.2; t < 200; t *= 1.002) {
        C3 q(t/sqrt(2), t/sqrt(2), 0);
        std::cout << t << " " << std::abs(octahedron.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "# q vs |F(q)| for q in direction 345, no special symmetry\n";
    for (double t = 0.2; t < 200; t *= 1.002) {
        C3 q(3*t/sqrt(50), 4*t/sqrt(50), 5*t/sqrt(50));
        std::cout << t << " " << std::abs(octahedron.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;
}
