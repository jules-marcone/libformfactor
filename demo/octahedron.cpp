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
    ff::platonic::Octahedron octahedron(1.);
    for (double t = 0.2; t < 200; t *= 1.002) {
        // choose q perpendicular to two opposite faces
        C3 q(0, sqrt(2. / 3) * t, sqrt(1. / 3) * t);
        std::cout << t << " " << std::abs(octahedron.formfactor(q)) << std::endl;
    }
}
