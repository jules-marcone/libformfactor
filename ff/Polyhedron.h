//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      ff/Polyhedron.h
//! @brief     Defines class Polyhedron.
//!
//! @homepage  http://www.bornagainproject.org
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2021
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

#ifdef SWIG
#error no need to expose this header to Swig
#endif

#ifndef USER_API
#ifndef FORMFACTOR_FF_POLYHEDRON_H
#define FORMFACTOR_FF_POLYHEDRON_H

#include <array>
#include <memory>
#include <vector>

#include <heinz/Complex.h>
#include <heinz/Vectors3D.h>

namespace ff {

class PolyhedralFace;
class PolyhedralTopology;

//! A polyhedron, implementation class for use in IFormFactorPolyhedron

class Polyhedron {
public:
    Polyhedron(const PolyhedralTopology& topology, double z_bottom,
               const std::vector<R3>& vertices);
    Polyhedron(const Polyhedron&) = delete;

    void assert_platonic() const;
    double volume() const;
    double radius() const;

    std::vector<R3> vertices() const; //! needed for topZ, bottomZ computation
    complex_t formfactor_at_center(const C3& q) const;
    complex_t formfactor_at_bottom(const C3& q) const;

private:

    double m_z_bottom;
    bool m_sym_Ci; //!< if true, then faces obtainable by inversion are not provided

    std::vector<PolyhedralFace> m_faces;
    double m_radius;
    double m_volume;
    std::vector<R3> m_vertices; //! for topZ, bottomZ computation only
};

} // namespace ff

#endif // FORMFACTOR_FF_POLYHEDRON_H
#endif // USER_API
