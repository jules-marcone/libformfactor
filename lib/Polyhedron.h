//  ************************************************************************************************
//
//  BornAgain: simulate and fit reflection and scattering
//
//! @file      Sample/LibFF/Polyhedron.h
//! @brief     Defines class Polyhedron.
//!
//! @homepage  http://www.bornagainproject.org
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2018
//! @authors   Scientific Computing Group at MLZ (see CITATION, AUTHORS)
//
//  ************************************************************************************************

#ifdef SWIG
#error no need to expose this header to Swig
#endif

#ifndef USER_API
#ifndef BORNAGAIN_SAMPLE_LIBFF_POLYHEDRON_H
#define BORNAGAIN_SAMPLE_LIBFF_POLYHEDRON_H

#include "Sample/LibFF/PolyhedralComponents.h"
#include <memory>

class PolyhedralTopology;

//! A polyhedron, implementation class for use in IFormFactorPolyhedron

class Polyhedron {
public:
    Polyhedron(const PolyhedralTopology& topology, double z_bottom,
               const std::vector<R3>& vertices);
    Polyhedron(const Polyhedron&) = delete;
    ~Polyhedron();

    void assert_platonic() const;
    double volume() const;
    double radius() const;

    const std::vector<R3>& vertices() const; //! needed for topZ, bottomZ computation
    complex_t evaluate_for_q(const C3& q) const;
    complex_t evaluate_centered(const C3& q) const;

private:
    double m_z_bottom;
    bool m_sym_Ci; //!< if true, then faces obtainable by inversion are not provided

    std::vector<PolyhedralFace> m_faces;
    double m_radius;
    double m_volume;
    std::vector<R3> m_vertices; //! for topZ, bottomZ computation only
};

#endif // BORNAGAIN_SAMPLE_LIBFF_POLYHEDRON_H
#endif // USER_API