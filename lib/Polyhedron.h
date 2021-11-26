//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      lib/Polyhedron.h
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
#ifndef LIBFORMFACTOR_LIB_POLYHEDRON_H
#define LIBFORMFACTOR_LIB_POLYHEDRON_H

#include <array>
#include <complex>
#include <memory>
#include <vector>

using complex_t = std::complex<double>;

template<class T> class Vec3;
using R3 = Vec3<double>;
using C3 = Vec3<std::complex<double>>;

using R3api = std::array<double,3>;
using C3api = std::array<complex_t,3>;

class PolyhedralFace;
class PolyhedralTopology;

//! A polyhedron, implementation class for use in IFormFactorPolyhedron

class Polyhedron {
public:
    Polyhedron(const PolyhedralTopology& topology, double z_bottom,
               const std::vector<R3api>& vertices);
    Polyhedron(const Polyhedron&) = delete;

    void assert_platonic() const;
    double volume() const;
    double radius() const;

    const std::vector<R3api> vertices() const; //! needed for topZ, bottomZ computation
    complex_t evaluate_for_q(const C3api& q) const;
    complex_t evaluate_centered(const C3api& q) const;

private:
    double m_z_bottom;
    bool m_sym_Ci; //!< if true, then faces obtainable by inversion are not provided

    std::vector<PolyhedralFace> m_faces;
    double m_radius;
    double m_volume;
    std::vector<R3> m_vertices; //! for topZ, bottomZ computation only
};

#endif // LIBFORMFACTOR_LIB_POLYHEDRON_H
#endif // USER_API
