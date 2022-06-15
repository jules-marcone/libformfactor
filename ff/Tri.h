//  ************************************************************************************************
//
//  BornAgain: simulate and fit reflection and scattering
//
//! @file      ff/Tri.h
//! @brief     defines trigonal shapes class
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @authors   Marianne Imperor-Clerc, Jules Marcone June 2022
//
//  ************************************************************************************************

#include <ff/Polyhedron.h>

namespace ff::tri {

class TriangularBipyramid : public ff::Polyhedron {
public:
    static ff::PolyhedralTopology topology();
    static std::vector<R3> vertices(const double edge);
    TriangularBipyramid(const double edge);
};

class ElongatedTriangularBipyramid : public ff::Polyhedron {
public:
    static ff::PolyhedralTopology topology();
    static std::vector<R3> vertices2(const double edge, const double height);
    ElongatedTriangularBipyramid(const double edge, const double height);
};

class TriangularBifrustum : public ff::Polyhedron {
public:
    static ff::PolyhedralTopology topology();
    static std::vector<R3> vertices3(const double edge, const double height, const double trunc);
    TriangularBifrustum(const double edge, const double height, const double trunc);
};
} // namespace ff::platonic
