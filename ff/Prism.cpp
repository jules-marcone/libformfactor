//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      ff/Prism.cpp
//! @brief     Implements class Prism.
//!
//! @homepage  https://jugit.fz-juelich.de/mlz/libformfactor
//! @license   GNU General Public License v3 or higher (see COPYING)
//! @copyright Forschungszentrum Jülich GmbH 2022
//! @authors   Scientific Computing Group at MLZ (see CITATION, AUTHORS)
//
//  ************************************************************************************************

//! The mathematics implemented here is described in full detail in a paper
//! by Joachim Wuttke, entitled
//! "Form factor (Fourier shape transform) of polygon and polyhedron."

#include "ff/Prism.h"
#include <stdexcept>

namespace {

complex_t sinc(const complex_t z) // cardinal sine function, sin(x)/x
{
    if (z == complex_t(0., 0.))
        return 1.0;
    return std::sin(z) / z;
}

} // namespace


ff::Prism::Prism(bool symmetry_Ci, double height, const std::vector<R3>& vertices)
{
    m_height = height;

    try {
        m_base = std::make_unique<ff::PolyhedralFace>(vertices, symmetry_Ci);
    } catch (std::invalid_argument& e) {
        throw std::invalid_argument(std::string("Invalid parameterization of Prism: ") + e.what());
    } catch (std::logic_error& e) {
        throw std::logic_error(std::string("Bug in Prism: ") + e.what()
                               + " [please report to the maintainers]");
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("Unexpected exception in Prism: ") + e.what()
                                 + " [please report to the maintainers]");
    }
}

double ff::Prism::area() const
{
    return m_base->area();
}

complex_t ff::Prism::formfactor(const C3& q) const
{
    try {
#ifdef ALGORITHM_DIAGNOSTIC
        polyhedralDiagnosis.reset();
        polyhedralDiagnosis.algo = 500;
#endif
        C3 qxy(q.x(), q.y(), 0.);
        return m_height * sinc(m_height / 2 * q.z()) * m_base->ff_2D(qxy);
    } catch (std::logic_error& e) {
        throw std::logic_error(std::string("Bug in Prism: ") + e.what()
                               + " [please report to the maintainers]");
    } catch (std::runtime_error& e) {
        throw std::runtime_error(std::string("Numeric computation failed in Prism: ") + e.what()
                                 + " [please report to the maintainers]");
    } catch (std::exception& e) {
        throw std::runtime_error(std::string("Unexpected exception in Prism: ") + e.what()
                                 + " [please report to the maintainers]");
    }
}
