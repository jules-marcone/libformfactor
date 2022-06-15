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

#include <ff/Cuboid.h>
#include <iostream>

constexpr double twopi = 6.28318530718;

//! Prints list t vs |F(q(t))| for a logarithmic range of t values

int main()
{
    std::cout << "# Cube form factor, for different q scans\n";
    ff::cuboid::Cube cube(1.);

    std::cout << "# q vs |F(q)| for q in direction 111, perpendicular to two faces\n";
    for (double t = 0.2; t < 200; t *= 1.002) {
        C3 q(t/sqrt(3), t/sqrt(3), t/sqrt(3));
        std::cout << t << " " << std::abs(cube.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "# q vs |F(q)| for q in direction 110, perpendicular to two edges\n";
    for (double t = 0.2; t < 200; t *= 1.002) {
        C3 q(t/sqrt(2), t/sqrt(2), 0);
        std::cout << t << " " << std::abs(cube.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "# q vs |F(q)| for q in direction 345, no special symmetry\n";
    for (double t = 0.2; t < 200; t *= 1.002) {
        C3 q(3*t/sqrt(50), 4*t/sqrt(50), 5*t/sqrt(50));
        std::cout << t << " " << std::abs(cube.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "# q vs |F(q)| for |q|=50, q on grand cercle through 111 and -1,-1,1 directions\n";
    const C3 a1(1/sqrt(3), 1/sqrt(3), 1/sqrt(3));
    const C3 a2(-1/sqrt(6), -1/sqrt(6), 2/sqrt(6));
    for (double t = 0; t < twopi; t += twopi/500) {
        C3 q = 50.*(cos(t)*a1+sin(t)*a2);
        std::cout << t << " " << std::abs(cube.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "# q vs |F(q)| for |q|=50, q on grand cercle through 111 and -2,3,-5 directions\n";
    const C3 b1(1/sqrt(3), 1/sqrt(3), 1/sqrt(3));
    const C3 b2(-8/sqrt(98), 3/sqrt(98), 5/sqrt(98));
    for (double t = 0; t < twopi; t += twopi/500) {
        C3 q = 50.*(cos(t)*b1+sin(t)*b2);
        std::cout << t << " " << std::abs(cube.formfactor(q)) << std::endl;
    }
    std::cout << std::endl;
}
