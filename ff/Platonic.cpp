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

//  ************************************************************************************************
//  class Dodecahedron
//  ************************************************************************************************

const ff::PolyhedralTopology Dodecahedron::topology = {
    {// bottom:
        {{0, 4, 3, 2, 1}, false},
        // lower ring:
        {{0, 5, 12, 9, 4}, false},
        {{4, 9, 11, 8, 3}, false},
        {{3, 8, 10, 7, 2}, false},
        {{2, 7, 14, 6, 1}, false},
        {{1, 6, 13, 5, 0}, false},
        // upper ring:
        {{8, 11, 16, 15, 10}, false},
        {{9, 12, 17, 16, 11}, false},
        {{5, 13, 18, 17, 12}, false},
        {{6, 14, 19, 18, 13}, false},
        {{7, 10, 15, 19, 14}, false},
        // top:
        {{15, 16, 17, 18, 19}, false}},
    true};

std::vector<R3> Dodecahedron::vertices(const double a)
{
    return {
        {0.8506508083520399 * a, 0 * a, -1.113516364411607 * a},
        {0.2628655560595668 * a, 0.8090169943749473 * a, -1.113516364411607 * a},
        {-0.6881909602355868 * a, 0.5 * a, -1.113516364411607 * a},
        {-0.6881909602355868 * a, -0.5 * a, -1.113516364411607 * a},
        {0.2628655560595668 * a, -0.8090169943749473 * a, -1.113516364411607 * a},
        {1.376381920471174 * a, 0 * a, -0.2628655560595667 * a},
        {0.42532540417602 * a, 1.309016994374947 * a, -0.2628655560595667 * a},
        {-1.113516364411607 * a, 0.8090169943749475 * a, -0.2628655560595667 * a},
        {-1.113516364411607 * a, -0.8090169943749475 * a, -0.2628655560595667 * a},
        {0.42532540417602 * a, -1.309016994374947 * a, -0.2628655560595667 * a},
        {-1.376381920471174 * a, 0 * a, 0.2628655560595667 * a},
        {-0.42532540417602 * a, -1.309016994374947 * a, 0.2628655560595667 * a},
        {1.113516364411607 * a, -0.8090169943749475 * a, 0.2628655560595667 * a},
        {1.113516364411607 * a, 0.8090169943749475 * a, 0.2628655560595667 * a},
        {-0.42532540417602 * a, 1.309016994374947 * a, 0.2628655560595667 * a},
        {-0.8506508083520399 * a, 0 * a, 1.113516364411607 * a},
        {-0.2628655560595668 * a, -0.8090169943749473 * a, 1.113516364411607 * a},
        {0.6881909602355868 * a, -0.5 * a, 1.113516364411607 * a},
        {0.6881909602355868 * a, 0.5 * a, 1.113516364411607 * a},
        {-0.2628655560595668 * a, 0.8090169943749473 * a, 1.113516364411607 * a}};
}

Dodecahedron::Dodecahedron(const double edge)
    : ff::Polyhedron(topology, vertices(edge))
{}

//  ************************************************************************************************
//  class Icosahedron
//  ************************************************************************************************

const ff::PolyhedralTopology Icosahedron::topology = {
    {// bottom:
        {{0, 2, 1}, false},
        // 1st row:
        {{0, 5, 2}, false},
        {{2, 3, 1}, false},
        {{1, 4, 0}, false},
        // 2nd row:
        {{0, 6, 5}, false},
        {{2, 5, 8}, false},
        {{2, 8, 3}, false},
        {{1, 3, 7}, false},
        {{1, 7, 4}, false},
        {{0, 4, 6}, false},
        // 3rd row:
        {{3, 8, 9}, false},
        {{5, 11, 8}, false},
        {{5, 6, 11}, false},
        {{4, 10, 6}, false},
        {{4, 7, 10}, false},
        {{3, 9, 7}, false},
        // 4th row:
        {{8, 11, 9}, false},
        {{6, 10, 11}, false},
        {{7, 9, 10}, false},
        // top:
        {{9, 11, 10}, false}},
    true};

std::vector<R3> Icosahedron::vertices(const double a)
{
    return {
        {0.5773502691896258 * a, 0 * a, -0.7557613140761708 * a},
        {-0.288675134594813 * a, 0.5 * a, -0.7557613140761708 * a},
        {-0.288675134594813 * a, -0.5 * a, -0.7557613140761708 * a},
        {-0.9341723589627158 * a, 0 * a, -0.1784110448865449 * a},
        {0.467086179481358 * a, 0.8090169943749475 * a, -0.1784110448865449 * a},
        {0.467086179481358 * a, -0.8090169943749475 * a, -0.1784110448865449 * a},
        {0.9341723589627158 * a, 0 * a, 0.1784110448865449 * a},
        {-0.467086179481358 * a, 0.8090169943749475 * a, 0.1784110448865449 * a},
        {-0.467086179481358 * a, -0.8090169943749475 * a, 0.1784110448865449 * a},
        {-0.5773502691896258 * a, 0 * a, 0.7557613140761708 * a},
        {0.288675134594813 * a, 0.5 * a, 0.7557613140761708 * a},
        {0.288675134594813 * a, -0.5 * a, 0.7557613140761708 * a}};
}

Icosahedron::Icosahedron(const double edge)
    : ff::Polyhedron(topology, vertices(edge))
{}

} // namespace ff::platonic
