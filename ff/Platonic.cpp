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

ff::PolyhedralTopology Tetrahedron::topology()
{
    return {{{{2, 1, 0}, false},
             {{0, 1, 3}, false},
             {{1, 2, 3}, false},
             {{2, 0, 3}, false}},
            false};
}

std::vector<R3> Tetrahedron::vertices(const double edge)
{
    const double a = edge;
    const double as = a / 2;
    const double ac = a / sqrt(3) / 2;
    const double ah = a / sqrt(3);
    const double height = sqrt(2. / 3) * edge;
    const double zcom = height / 4; // center of mass

    return {{-ac, as, -zcom}, {-ac, -as, -zcom}, {ah, 0., -zcom}, {0, 0., height - zcom}};
}

Tetrahedron::Tetrahedron(const double edge) : ff::Polyhedron(topology(), vertices(edge)) {}

//  ************************************************************************************************
//  class Octahedron
//  ************************************************************************************************

ff::PolyhedralTopology Octahedron::topology()
{
    return {{{{0, 2, 1}, false},
             {{0, 3, 2}, false},
             {{0, 4, 3}, false},
             {{0, 1, 4}, false},
             {{2, 3, 5}, false},
             {{1, 2, 5}, false},
             {{4, 1, 5}, false},
             {{3, 4, 5}, false}},
            true};
}

/// diagonal position: axes on vertices
//std::vector<R3> Octahedron::vertices(const double edge)
//{
//    const double h = edge / sqrt(2);

//    return {{0, 0, -h}, {h, 0, 0}, {0, h, 0}, {-h, 0, 0}, {0, -h, 0}, {0, 0, h}};
//}

/// rotated position: x and y axes are perpendicular to the edges
std::vector<R3> Octahedron::vertices(const double edge)
{
    const double h = edge / sqrt(2);
    const double a = edge / 2;

    return {{0, 0, -h}, {a, -a, 0}, {a, a, 0}, {-a, a, 0}, {-a, -a, 0}, {0, 0, h}};
}

Octahedron::Octahedron(const double edge) : ff::Polyhedron(topology(), vertices(edge)) {}

//  ************************************************************************************************
//  class Dodecahedron
//  ************************************************************************************************

ff::PolyhedralTopology Dodecahedron::topology()
{
    return {{// bottom:
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
}

std::vector<R3> Dodecahedron::vertices(const double a)
{
    const double r1 = 0.2628655560595668 * a; // sqrt((5-sqrt(5)/40)
    const double r2 = 0.42532540417602 * a; // r1*phi
    const double r3 = 0.5 * a;
    const double r4 = 0.6881909602355868 * a; // r2*phi
    const double r5 = 0.8090169943749473 * a; // r3*phi
    const double r6 = 0.8506508083520399 * a; // r1 * 2 * phi
    const double r7 = 1.113516364411607 * a; // r4*phi
    const double r8 = 1.309016994374947 * a; // r5*phi
    const double r9 = 1.376381920471174 * a; // r6*phi
    return {{r6, 0, -r7},
            {r1, r5, -r7},
            {-r4, r3, -r7},
            {-r4, -r3, -r7},
            {r1, -r5, -r7},
            {r9, 0, -r1},
            {r2, r8, -r1},
            {-r7, r5, -r1},
            {-r7, -r5, -r1},
            {r2, -r8, -r1},
            {-r9, 0, r1},
            {-r2, -r8, r1},
            {r7, -r5, r1},
            {r7, r5, r1},
            {-r2, r8, r1},
            {-r6, 0, r7},
            {-r1, -r5, r7},
            {r4, -r3, r7},
            {r4, r3, r7},
            {-r1, r5, r7}};
}

Dodecahedron::Dodecahedron(const double edge) : ff::Polyhedron(topology(), vertices(edge)) {}

//  ************************************************************************************************
//  class Icosahedron
//  ************************************************************************************************

ff::PolyhedralTopology Icosahedron::topology()
{
    return {{// bottom:
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
}

std::vector<R3> Icosahedron::vertices(const double a)
{
    const double s1 = 0.1784110448865449 * a; // 1/sqrt(6)/sqrt(3+sqrt(5))
    const double s2 = 0.288675134594813 * a; // s1 * phi
    const double s3 = 0.467086179481358 * a; // s2 * phi
    const double s4 = 0.5 * a;
    const double s5 = 0.5773502691896258 * a; // 2 * s2
    const double s6 = 0.7557613140761708 * a; // s3 * phi
    const double s7 = 0.8090169943749473 * a; // phi/2
    const double s8 = 0.9341723589627158 * a; // s5 * phi
    return {{s5, 0, -s6},
            {-s2, s4, -s6},
            {-s2, -s4, -s6},
            {-s8, 0, -s1},
            {s3, s7, -s1},
            {s3, -s7, -s1},
            {s8, 0, s1},
            {-s3, s7, s1},
            {-s3, -s7, s1},
            {-s5, 0, s6},
            {s2, s4, s6},
            {s2, -s4, s6}};
}

Icosahedron::Icosahedron(const double edge) : ff::Polyhedron(topology(), vertices(edge)) {}

} // namespace ff::platonic
