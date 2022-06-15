//  ************************************************************************************************
//
//  BornAgain: simulate and fit reflection and scattering
//
//! @file      ff/Cuboid.cpp
//! @brief     cuboid shapes, classes ...
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @authors   Marianne Imperor-Clerc, Jules Marcone June 2022
//
//  ************************************************************************************************
#include <ff/Cuboid.h>

namespace ff::cuboid {

//  ************************************************************************************************
//  class Cube
//  ************************************************************************************************

ff::PolyhedralTopology Cube::topology()
{
    return {{{{3, 2, 1, 0}, true},
             {{1, 2, 6, 5}, true},
             {{0, 1, 5, 4}, true},
             {{3, 0, 4, 7}, true},
             {{2, 3, 7, 6}, true},
             {{4, 5, 6, 7}, true}},
            false};
}

std::vector<R3> Cube::vertices(const double edge)
{
    
    const double a = edge / 2.;

    
    return {{a, -a, -a}, {a, a, -a}, {-a, a, -a}, {-a, -a, -a}, {a, -a, a}, {a, a, a}, {-a, a, a}, {-a, -a, a}};
}

Cube::Cube(const double edge) : ff::Polyhedron(topology(), vertices(edge)) {}


//  ************************************************************************************************
//  class Pave
//  ************************************************************************************************

ff::PolyhedralTopology Pave::topology()
{
    return {{{{3, 2, 1, 0}, true},
             {{1, 2, 6, 5}, true},
             {{0, 1, 5, 4}, true},
             {{3, 0, 4, 7}, true},
             {{2, 3, 7, 6}, true},
             {{4, 5, 6, 7}, true}},
            false};
}

std::vector<R3> Pave::vertices3(const double edge_a, const double edge_b, const double edge_c)
{
    
    const double a = edge_a/2.;
    const double b = edge_b/2.;
    const double c = edge_c/2.;

    
    return {{a, -b, -c}, {a, b, -c}, {-a, b, -c}, {-a, -b, -c}, {a, -b, c}, {a, b, c}, {-a, b, c}, {-a, -b, c}};
}

Pave::Pave(const double edge_a, const double edge_b, const double edge_c) : ff::Polyhedron(topology(), vertices3(edge_a, edge_b, edge_c)) {}

} // namespace ff::platonic
