//  ************************************************************************************************
//
//  BornAgain: simulate and fit reflection and scattering
//
//! @file      ff/Penta.h
//! @brief     defines pentagonal shapes class
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @authors   Marianne Imperor-Clerc, Jules Marcone June 2022
//
//  ************************************************************************************************

#include <ff/Polyhedron.h>

namespace ff::penta {

// regular decahedron
class Decahedron : public ff::Polyhedron {
public:
    static ff::PolyhedralTopology topology();
    static std::vector<R3> vertices(const double edge);
    Decahedron(const double edge);

};

class ElongatedDecahedron : public ff::Polyhedron {
public:
    static ff::PolyhedralTopology topology();
    static std::vector<R3> vertices2(const double edge, const double height);
    ElongatedDecahedron(const double edge, const double height);

};

class PentagonalBifrustum : public ff::Polyhedron {
public:
    static ff::PolyhedralTopology topology();
    static std::vector<R3> vertices3(const double edge, const double height, const double trunc);
    PentagonalBifrustum(const double edge, const double height, const double trunc);

};

class CappedPentagonalPrism : public ff::Polyhedron {
public:
    static ff::PolyhedralTopology topology();
    static std::vector<R3> vertices3(const double edge, const double height, const double capsize);
    CappedPentagonalPrism(const double edge, const double height, const double capsize);

};
} // namespace ff::penta
