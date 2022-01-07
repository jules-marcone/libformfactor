//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      ff/PolyhedralTopology.h
//! @brief     Defines classes PolygonalTopology, PolyhedralTopology
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
#endif // USER_API
