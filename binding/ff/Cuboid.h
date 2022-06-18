//  ************************************************************************************************
//
//  BornAgain: simulate and fit reflection and scattering
//
//! @file      ff/Cuboid.h
//! @brief     defines cuboid shapes class
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @authors   Marianne Imperor-Clerc, Jules Marcone June 2022
//
//  ************************************************************************************************

#include <ff/Polyhedron.h>

namespace ff::cuboid {

class Cube : public ff::Polyhedron {
public:
    static ff::PolyhedralTopology topology();
    static std::vector<R3> vertices(const double edge);
    Cube(const double edge);
};

class Pave : public ff::Polyhedron {
public:
    static ff::PolyhedralTopology topology();
    static std::vector<R3> vertices3(const double edge_a, const double edge_b, const double edge_c);
    Pave(const double edge_a, const double edge_b, const double edge_c);
};
} // namespace ff::platonic
