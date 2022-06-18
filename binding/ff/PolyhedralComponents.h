//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      ff/PolyhedralComponents.h
//! @brief     Defines classes PolyhedralEdge, PolyhedralFace
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

#ifndef FORMFACTOR_FF_POLYHEDRALCOMPONENTS_H
#define FORMFACTOR_FF_POLYHEDRALCOMPONENTS_H

#include <heinz/Complex.h>
#include <heinz/Vectors3D.h>
#include <vector>

namespace ff {

#ifdef ALGORITHM_DIAGNOSTIC
#include <string>

struct PolyhedralDiagnosis {
    int algo;
    int order;
    std::string msg;
    void reset();
    std::string message() const;
    bool operator==(const PolyhedralDiagnosis&) const;
    bool operator!=(const PolyhedralDiagnosis&) const;
};
inline PolyhedralDiagnosis polyhedralDiagnosis;
#endif

//! One edge of a polygon, for form factor computation.

class PolyhedralEdge {
public:
    PolyhedralEdge(R3 Vlow, R3 Vhig);

    R3 E() const { return m_E; }
    R3 R() const { return m_R; }
    complex_t qE(C3 q) const { return m_E.dot(q); }
    complex_t qR(C3 q) const { return m_R.dot(q); }

    complex_t contrib(int m, C3 qpa, complex_t qrperp) const;

private:
    R3 m_E; //!< vector pointing from mid of edge to upper vertex
    R3 m_R; //!< position vector of edge midpoint
};

//! A polygon, for form factor computation.

class PolyhedralFace {
public:
    static double diameter(const std::vector<R3>& V);

    PolyhedralFace(const std::vector<R3>& _V = std::vector<R3>(), bool _sym_S2 = false);

    double area() const { return m_area; }
    double pyramidalVolume() const { return m_rperp * m_area / 3; }
    double radius3d() const { return m_radius_3d; }
    //! Returns conj(q)*normal [BasicVector3D::dot is antilinear in 'this' argument]
    complex_t normalProjectionConj(C3 q) const { return q.dot(m_normal); }
    complex_t ff_n(int n, C3 q) const;
    complex_t ff(C3 q, bool sym_Ci) const;
    complex_t ff_2D(C3 qpa) const;
    complex_t ff_2D_direct(C3 qpa) const;   // for TestTriangle
    complex_t ff_2D_expanded(C3 qpa) const; // for TestTriangle
    void assert_Ci(const PolyhedralFace& other) const;

private:
    static double qpa_limit_series; //!< determines when use power series
    static int n_limit_series;

    bool sym_S2; //!< if true, then edges obtainable by inversion are not provided
    std::vector<PolyhedralEdge> edges;
    double m_area;
    R3 m_normal;        //!< normal vector of this polygon's plane
    double m_rperp;     //!< distance of this polygon's plane from the origin, along 'm_normal'
    double m_radius_2d; //!< radius of enclosing cylinder
    double m_radius_3d; //!< radius of enclosing sphere

    void decompose_q(C3 q, complex_t& qperp, C3& qpa) const;
    complex_t ff_n_core(int m, C3 qpa, complex_t qperp) const;
    complex_t edge_sum_ff(C3 q, C3 qpa, bool sym_Ci) const;
    complex_t expansion(complex_t fac_even, complex_t fac_odd, C3 qpa, double abslevel) const;
};

} // namespace ff

#endif // FORMFACTOR_FF_POLYHEDRALCOMPONENTS_H
