//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      ff/Prism.h
//! @brief     Defines class Prism.
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

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
    complex_t formfactor(const C3& q) const;

private:
    std::unique_ptr<ff::PolyhedralFace> m_base;
    double m_height;
};

} // namespace ff

#endif // FORMFACTOR_FF_PRISM_H
