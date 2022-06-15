//  ************************************************************************************************
//
//  BornAgain: simulate and fit reflection and scattering
//
//! @file      ff/Tri.cpp
//! @brief     trigonal shapes class, bipyramids, bifrustum, etc ...
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @authors  Marianne Imperor-Clerc, Jules Marcone June 2022
//
//  ************************************************************************************************

#include "ff/Tri.h"

namespace ff::tri {


//  ************************************************************************************************
//  Triangular Bipyramid
//  ************************************************************************************************


ff::PolyhedralTopology TriangularBipyramid::topology()

{
    return {{{{0, 1, 3}, false},
             {{1, 2, 3}, false},
             {{2, 0, 3}, false},
             {{1, 0, 4}, false},
             {{2, 1, 4}, false},
             {{0, 2, 4}, false}},
            false};
}

std::vector<R3> TriangularBipyramid::vertices(const double edge)
{
    const double a = edge/sqrt(3.);
    const double x = a / 2.;
    const double y = sqrt(3.) * a / 2.;
    const double h = sqrt(2.) * a;

    return {{-x, y, 0.}, {-x, -y, 0.}, {a, 0., 0.}, {0., 0., h}, {0., 0., -h}};
}

TriangularBipyramid::TriangularBipyramid(const double edge) : ff::Polyhedron(topology(), vertices(edge)) {}

//  ************************************************************************************************
//  Elongated Triangular Bipyramid (Triangular Bypramid with an added parameter for anisotropy)
//  ************************************************************************************************


ff::PolyhedralTopology ElongatedTriangularBipyramid::topology()

{
    return {{{{0, 1, 3}, false},
             {{1, 2, 3}, false},
             {{2, 0, 3}, false},
             {{1, 0, 4}, false},
             {{2, 1, 4}, false},
             {{0, 2, 4}, false}},
            false};
}

std::vector<R3> ElongatedTriangularBipyramid::vertices2(const double edge, const double height)
{
    const double a = edge/sqrt(3.);
    const double x = a / 2.;
    const double y = sqrt(3.) * a / 2.;
    const double h = height;

    return {{-x, y, 0.}, {-x, -y, 0.}, {a, 0., 0.}, {0., 0., h}, {0., 0., -h}};
}

ElongatedTriangularBipyramid::ElongatedTriangularBipyramid(const double edge, const double height) : ff::Polyhedron(topology(), vertices2(edge, height)) {}
//  ************************************************************************************************
//  Triangular Bifrustum (parameters are edges of base triangle, total theoritcal height of bipyramid, and height where truncature was operated as ratio of theoretical height)
//  ************************************************************************************************


ff::PolyhedralTopology TriangularBifrustum::topology()

{
    return {{{{0, 1, 4, 3}, false},
             {{1, 2, 5, 4}, false},
             {{2, 0, 3, 5}, false},
             {{1, 0, 6, 7}, false},
             {{2, 1, 7, 8}, false},
             {{0, 2, 8, 6}, false},
             {{3, 4, 5}, false},
             {{7, 6, 8}, false}},
            false};
}

std::vector<R3> TriangularBifrustum::vertices3(const double edge, const double height, const double trunc)
{
    const double a = edge / sqrt(3.);
    const double z = trunc; //z is the new height of bifrustum compared to the original bipyramid (z is between 0 and 1)
    const double x = 0.5 * a;
    const double y = sqrt(3.) * a / 2.;
    const double h = height;

    return {{-x, y, 0.}, {-x, -y, 0.}, {a, 0., 0.}, //middle plane
            {-x*(1-z), y*(1-z), z*h}, {-x*(1-z), -y*(1-z), z*h}, {(1-z)*a, 0., z*h}, //top plane
            {-x*(1-z), y*(1-z), -z*h}, {-x*(1-z), -y*(1-z), -z*h}, {(1-z)*a, 0., -z*h}};//bottom plane
}

TriangularBifrustum::TriangularBifrustum(const double edge, const double height, const double trunc) : ff::Polyhedron(topology(), vertices3(edge, height, trunc)) {}

} // namespace ff::tri
