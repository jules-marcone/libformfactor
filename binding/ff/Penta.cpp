//  ************************************************************************************************
//
//  BornAgain: simulate and fit reflection and scattering
//
//! @file      ff/Penta.cpp
//! @brief     pentagonal shapes class, bipyramids, decahedron, etc ...
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @authors   Marianne Imperor-Clerc, Jules Marcone June 2022
//
//  ************************************************************************************************

#include "ff/Penta.h"

namespace ff::penta {


//  ************************************************************************************************
//  regular Decahedron
//  ************************************************************************************************


ff::PolyhedralTopology Decahedron::topology()

{
    return {{{{0, 1, 5}, false}, 
             {{1, 2, 5}, false}, 
             {{2, 3, 5}, false}, 
             {{3, 4, 5}, false}, 
             {{4, 0, 5}, false}, 
             {{1, 0, 6}, false}, 
             {{2, 1, 6}, false}, 
             {{3, 2, 6}, false}, 
             {{4, 3, 6}, false}, 
             {{0, 4, 6}, false}},
            false};
}

std::vector<R3> Decahedron::vertices(const double edge)
{
    const double coeff = 0.8506508083520399;
    const double a = edge*coeff;
    const double ac5 = a*0.30901699437494745;
    const double as5 = a*0.9510565162951535;
    const double a2c5 = -a*0.8090169943749475;
    const double a2s5 = a*0.5877852522924731;
    const double height = edge * sqrt(1-coeff);

    return {{a, 0., 0.},
            {ac5, as5, 0.}, 
            {a2c5, a2s5, 0.},
            {a2c5, -a2s5, 0.},
            {ac5,-as5, 0.},
            {0., 0., height}, 
            {0., 0., -height}};
}

Decahedron::Decahedron(const double edge) : ff::Polyhedron(topology(), vertices(edge)) {}

//  ************************************************************************************************
//  Elongated Decahedron (decahedron with added parameter for anisotropy)
//  ************************************************************************************************


ff::PolyhedralTopology ElongatedDecahedron::topology()

{
    return {{{{0, 1, 5}, false}, 
             {{1, 2, 5}, false}, 
             {{2, 3, 5}, false}, 
             {{3, 4, 5}, false}, 
             {{4, 0, 5}, false}, 
             {{1, 0, 6}, false}, 
             {{2, 1, 6}, false}, 
             {{3, 2, 6}, false}, 
             {{4, 3, 6}, false}, 
             {{0, 4, 6}, false}},
            false};
}

std::vector<R3> ElongatedDecahedron::vertices2(const double edge, const double height)
{
    const double coeff = 0.8506508083520399;
    const double a = edge*coeff;
    const double ac5 = a*0.30901699437494745;
    const double as5 = a*0.9510565162951535;
    const double a2c5 = -a*0.8090169943749475;
    const double a2s5 = a*0.5877852522924731;
    const double h = height;

    return {{a, 0., 0.},
            {ac5, as5, 0.}, 
            {a2c5, a2s5, 0.},
            {a2c5, -a2s5, 0.},
            {ac5,-as5, 0.},
            {0., 0., h}, 
            {0., 0., -h}};
}

ElongatedDecahedron::ElongatedDecahedron(const double edge, const double height) : ff::Polyhedron(topology(), vertices2(edge, height)) {}

//  ************************************************************************************************
//  Pentagonal Bifrustum
//  ************************************************************************************************


ff::PolyhedralTopology PentagonalBifrustum::topology()

{
    return {{//Top face
             {{5, 6, 7, 8, 9}, false},
             //First row of faces
             {{0, 1, 6, 5}, false}, 
             {{1, 2, 7, 6}, false}, 
             {{2, 3, 8, 7}, false}, 
             {{3, 4, 9, 8}, false}, 
             {{4, 0, 5, 9}, false},
             //Second row of faces 
             {{1, 0, 10, 11}, false}, 
             {{2, 1, 11, 12}, false}, 
             {{3, 2, 12, 13}, false}, 
             {{4, 3, 13, 14}, false}, 
             {{0, 4, 14, 10}, false},
             //Bottom face
             {{14, 13, 12, 11, 10}, false}},
            false};
}

std::vector<R3> PentagonalBifrustum::vertices3(const double edge, const double height, const double trunc)
{
    const double coeff = 0.8506508083520399;
    const double a = edge*coeff;
    const double z = trunc;
    const double ac5 = a*0.30901699437494745;
    const double as5 = a*0.9510565162951535;
    const double a2c5 = -a*0.8090169943749475;
    const double a2s5 = a*0.5877852522924731;
    const double h = height;

    return {//middle plane
            {a, 0., 0.},
            {ac5, as5, 0.}, 
            {a2c5, a2s5, 0.},
            {a2c5, -a2s5, 0.},
            {ac5,-as5, 0.},
            //top plane
            {a*(1.-z), 0., z*h},
            {ac5*(1.-z), as5*(1.-z), z*h}, 
            {a2c5*(1.-z), a2s5*(1.-z), z*h},
            {a2c5*(1.-z), -a2s5*(1.-z), z*h},
            {ac5*(1.-z),-as5*(1.-z), z*h},
            //bottom plane 
            {a*(1.-z), 0., -z*h},
            {ac5*(1.-z), as5*(1.-z), -z*h}, 
            {a2c5*(1.-z), a2s5*(1.-z), -z*h},
            {a2c5*(1.-z), -a2s5*(1.-z), -z*h},
            {ac5*(1.-z),-as5*(1.-z), -z*h},};
}

PentagonalBifrustum::PentagonalBifrustum(const double edge, const double height, const double trunc) : ff::Polyhedron(topology(), vertices3(edge, height, trunc)) {}

//  ************************************************************************************************
//  Capped Pentagonal Prism (nanorods)
//  Height is length of prism, capsize is height of pyramids on each side
//  ************************************************************************************************


ff::PolyhedralTopology CappedPentagonalPrism::topology()

{
    return {{//Top Pyramid
             {{0, 1, 10}, false},
             {{1, 2, 10}, false},
             {{2, 3, 10}, false},
             {{3, 4, 10}, false},
             {{4, 0, 10}, false},
             //Central Prism
             {{5, 6, 1, 0}, true},
             {{6, 7, 2, 1}, true},
             {{7, 8, 3, 2}, true},
             {{8, 9, 4, 3}, true},
             {{9, 5, 0, 4}, true},
             //Bottom Pyramid}
             {{6, 5, 11}, false},
             {{7, 6, 11}, false},
             {{8, 7, 11}, false},
             {{9, 8, 11}, false},
             {{5, 9, 11}, false}},
            false};
}

std::vector<R3> CappedPentagonalPrism::vertices3(const double edge, const double height, const double capsize)
{
    const double coeff = 0.8506508083520399;
    const double a = edge*coeff;
    const double ac5 = a*0.30901699437494745;
    const double as5 = a*0.9510565162951535;
    const double a2c5 = -a*0.8090169943749475;
    const double a2s5 = a*0.5877852522924731;
    const double z = capsize;
    const double h = height/2.;

    return {//Top face of prism
            {a, 0., h},
            {ac5, as5, h}, 
            {a2c5, a2s5, h},
            {a2c5, -a2s5, h},
            {ac5,-as5, h},
            //Bottom face of prism
            {a, 0., -h},
            {ac5, as5, -h}, 
            {a2c5, a2s5, -h},
            {a2c5, -a2s5, -h},
            {ac5,-as5, -h},
            //Top vertex of upper pyramid
            {0., 0., h+z},
            //Bottom vertex of bottom pyramid
            {0., 0., -h-z}};
}

CappedPentagonalPrism::CappedPentagonalPrism(const double edge, const double height, const double capsize) : ff::Polyhedron(topology(), vertices3(edge, height, capsize)) {}

} // namespace ff::penta
