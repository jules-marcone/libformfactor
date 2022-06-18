//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      ff/PolyhedralTopology.h
//! @brief     Defines classes PolygonalTopology, PolyhedralTopology
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

#ifndef FORMFACTOR_FF_POLYHEDRALTOPOLOGY_H
#define FORMFACTOR_FF_POLYHEDRALTOPOLOGY_H

#include <vector>

namespace ff {

//! For internal use in PolyhedralFace.
class PolygonalTopology {
public:
    std::vector<int> vertexIndices;
    bool symmetry_S2;
};

//! For internal use in IFormFactorPolyhedron.
class PolyhedralTopology {
public:
    std::vector<PolygonalTopology> faces;
    bool symmetry_Ci;
};

} // namespace ff

#endif // FORMFACTOR_FF_POLYHEDRALTOPOLOGY_H
