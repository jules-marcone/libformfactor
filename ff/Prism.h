//  ************************************************************************************************
//
//  BornAgain: simulate and fit reflection and scattering
//
//! @file      ff/Prism.h
//! @brief     Defines class Prism.
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
#ifndef FORMFACTOR_FF_PRISM_H
#define FORMFACTOR_FF_PRISM_H

#include "ff/PolyhedralComponents.h"
#include "ff/PolyhedralTopology.h"
#include <memory>

namespace ff {

class Prism {
public:
    Prism(bool symmetry_Ci, double height, const std::vector<R3>& vertices);
    Prism(const Prism&) = delete;

    double area() const;
    const std::vector<R3>& vertices() const; //! needed for topZ, bottomZ computation
    complex_t formfactor_at_center(const C3& q) const;

private:
    std::unique_ptr<ff::PolyhedralFace> m_base;
    double m_height;
    std::vector<R3> m_vertices; //! for topZ, bottomZ computation only
};

} // namespace ff

#endif // FORMFACTOR_FF_PRISM_H
#endif // USER_API
