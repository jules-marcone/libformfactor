//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      ff/Polyhedron.h
//! @brief     Defines class Polyhedron.
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

#ifndef FORMFACTOR_FF_POLYHEDRON_H
#define FORMFACTOR_FF_POLYHEDRON_H

#include <array>
#include <memory>
#include <vector>

#include <ff/PolyhedralComponents.h>
#include <ff/PolyhedralTopology.h>
#include <heinz/Complex.h>
#include <heinz/Vectors3D.h>

namespace ff {

//! A polyhedron, implementation class for use in IFormFactorPolyhedron

class Polyhedron {
public:
    Polyhedron(const PolyhedralTopology& topology, const std::vector<R3>& vertices);
    Polyhedron(const Polyhedron&) = delete;

    void assert_platonic() const;
    double volume() const;
    double radius() const;

    complex_t formfactor(const C3& q) const;

private:
    bool m_sym_Ci; //!< if true, then faces obtainable by inversion are not provided

    std::vector<PolyhedralFace> m_faces;
    double m_radius;
    double m_volume;
};

} // namespace ff

#endif // FORMFACTOR_FF_POLYHEDRON_H
