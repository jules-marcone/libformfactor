//  ************************************************************************************************
//
//  BornAgain: simulate and fit reflection and scattering
//
//! @file      ff/Platonic.cpp
//! @brief     Implements platonic solid classes.
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @authors   Scientific Computing Group at MLZ (see CITATION, AUTHORS)
//
//  ************************************************************************************************

#include "ff/Platonic.h"

namespace ff::platonic {

//  ************************************************************************************************
//  class Tetrahedron
//  ************************************************************************************************

const ff::PolyhedralTopology Tetrahedron::topology = {
    {{{2, 1, 0}, false}, {{0, 1, 3}, false}, {{1, 2, 3}, false}, {{2, 0, 3}, false}}, false};

std::vector<R3> Tetrahedron::vertices(const double edge)
{
    const double a = edge;
    const double as = a / 2;
    const double ac = a / sqrt(3) / 2;
    const double ah = a / sqrt(3);
    const double height = sqrt(2. / 3) * edge;
    const double zcom = height / 4; // center of mass

    return{
        {-ac, as, -zcom},
        {-ac, -as, -zcom},
        {ah, 0., -zcom},
        {0, 0., height - zcom}};
}

Tetrahedron::Tetrahedron(const double edge)
    : ff::Polyhedron(topology, vertices(edge))
{}

//  ************************************************************************************************
//  class Octahedron
//  ************************************************************************************************

const ff::PolyhedralTopology Octahedron::topology = {
    {{{0, 2, 1}, false},
     {{0, 3, 2}, false},
     {{0, 4, 3}, false},
     {{0, 1, 4}, false},
     {{2, 3, 5}, false},
     {{1, 2, 5}, false},
     {{4, 1, 5}, false},
     {{3, 4, 5}, false}},
    true};

std::vector<R3> Octahedron::vertices(const double edge)
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

Octahedron::Octahedron(const double edge)
    : ff::Polyhedron(topology, vertices(edge))
{}

} // namespace ff::platonic
