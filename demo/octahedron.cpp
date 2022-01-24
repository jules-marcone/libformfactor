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

#include <ff/Polyhedron.h>
#include <iostream>

const ff::PolyhedralTopology octahedron_topology = {
    {{{0, 2, 1}, false},
     {{0, 3, 2}, false},
     {{0, 4, 3}, false},
     {{0, 1, 4}, false},
     {{2, 3, 5}, false},
     {{1, 2, 5}, false},
     {{4, 1, 5}, false},
     {{3, 4, 5}, false}},
    true};

std::vector<R3> octahedron_vertices(const double edge)
{

    const double a = edge / 2;
    const double h = a * sqrt(2);

    return {
        {0, 0, -h},
        {-a, -a, 0},
        {a, -a, 0},
        {a, a, 0},
        {-a, a, 0},
        {0, 0, h}};
}

//! Prints list t vs |F(q(t))| for a logarithmic range of t values

int main() {
    ff::Polyhedron octahedron(octahedron_topology, octahedron_vertices(1.));
    for (double t=0.2; t<200;  t *= 1.002) {
        // choose q perpendicular to two opposite faces
        C3 q(0, sqrt(2./3)*t, sqrt(1./3)*t);
        std::cout << t << " " << std::abs(octahedron.formfactor(q)) << std::endl;
    }
}
