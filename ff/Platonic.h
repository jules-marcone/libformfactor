//  ************************************************************************************************
//
//  BornAgain: simulate and fit reflection and scattering
//
//! @file      ff/Platonic.h
//! @brief     Defines platonic solid classes.
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @authors   Scientific Computing Group at MLZ (see CITATION, AUTHORS)
//
//  ************************************************************************************************

#include <ff/Polyhedron.h>

namespace ff::platonic {

class Tetrahedron : public ff::Polyhedron {
 public:
    static const ff::PolyhedralTopology topology;
    static std::vector<R3> vertices(const double edge);
    Tetrahedron(const double edge);
};

class Octahedron : public ff::Polyhedron {
 public:
    static const ff::PolyhedralTopology topology;
    static std::vector<R3> vertices(const double edge);
    Octahedron(const double edge);
};

class Dodecahedron : public ff::Polyhedron {
 public:
    static const ff::PolyhedralTopology topology;
    static std::vector<R3> vertices(const double edge);
    Dodecahedron(const double edge);
};

class Icosahedron : public ff::Polyhedron {
 public:
    static const ff::PolyhedralTopology topology;
    static std::vector<R3> vertices(const double edge);
    Icosahedron(const double edge);
};

} // namespace ff::platonic
