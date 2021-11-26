//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      lib/Vec3.h
//! @brief     Defines vectors in R^3, C^3.
//!
//! @homepage  http://www.bornagainproject.org
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2021
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//
//  ************************************************************************************************

#ifndef LIBFORMFACTOR_LIB_VEC3_H
#define LIBFORMFACTOR_LIB_VEC3_H

#include "Complex.h"
#include <array>
#include <ostream>

//! Three-dimensional vector template, for use with integer, double, or complex components.

template <class T> class Vec3 : public std::array<T, 3> {
private:
    using super = std::array<T, 3>;

public:
    // -------------------------------------------------------------------------
    // Constructors and other set functions
    // -------------------------------------------------------------------------

    //! Constructs the null vector.
    Vec3() : super{0., 0., 0.} {}

    //! Constructs a vector from cartesian components.
    Vec3(const T x, const T y, const T z) : super{x, y, z} {}

    // -------------------------------------------------------------------------
    // Component access
    // -------------------------------------------------------------------------

    //! Returns x-component in cartesian coordinate system.
    inline T x() const { return (*this)[0]; }
    //! Returns y-component in cartesian coordinate system.
    inline T y() const { return (*this)[1]; }
    //! Returns z-component in cartesian coordinate system.
    inline T z() const { return (*this)[2]; }

    // -------------------------------------------------------------------------
    // In-place operations
    // -------------------------------------------------------------------------

    //! Adds other vector to this, and returns result.
    Vec3<T>& operator+=(const Vec3<T>& v)
    {
        (*this)[0] += v[0];
        (*this)[1] += v[1];
        (*this)[2] += v[2];
        return *this;
    }

    //! Subtracts other vector from this, and returns result.
    Vec3<T>& operator-=(const Vec3<T>& v)
    {
        (*this)[0] -= v[0];
        (*this)[1] -= v[1];
        (*this)[2] -= v[2];
        return *this;
    }

    //! Multiplies this with a scalar, and returns result.
#ifndef SWIG
    template <class U> auto operator*=(U a)
    {
        (*this)[0] *= a;
        (*this)[1] *= a;
        (*this)[2] *= a;
        return *this;
    }
#endif // USER_API

    //! Divides this by a scalar, and returns result.
#ifndef SWIG
    template <class U> auto operator/=(U a)
    {
        (*this)[0] /= a;
        (*this)[1] /= a;
        (*this)[2] /= a;
        return *this;
    }
#endif // USER_API

    // -------------------------------------------------------------------------
    // Functions of this (with no further argument)
    // -------------------------------------------------------------------------

    //! Returns complex conjugate vector
    Vec3<T> conj() const;

    //! Returns magnitude squared of the vector.
    double mag2() const { return std::norm(x()) + std::norm(y()) + std::norm(z()); }

    //! Returns magnitude of the vector.
    double mag() const { return sqrt(mag2()); }

    //! Returns unit vector in direction of this. Throws for null vector.
    Vec3<T> unit() const;

    //! Returns this, trivially converted to complex type.
    Vec3<complex_t> complex() const;

    //! Returns real parts.
    Vec3<double> real() const;

    // -------------------------------------------------------------------------
    // Functions of this and another vector
    // -------------------------------------------------------------------------

    //! Returns dot product of vectors (antilinear in the first [=self] argument).
#ifndef SWIG
    template <class U> auto dot(const Vec3<U>& v) const;
#endif // USER_API

    //! Returns cross product of vectors (linear in both arguments).
#ifndef SWIG
    template <class U> auto cross(const Vec3<U>& v) const;
#endif // USER_API
};

// =============================================================================
// Non-member functions
// =============================================================================

//! Output to stream.
//! @relates Vec3
template <class T> std::ostream& operator<<(std::ostream& os, const Vec3<T>& a)
{
    return os << "(" << a.x() << "," << a.y() << "," << a.z() << ")";
}

// -----------------------------------------------------------------------------
// Unary operators
// -----------------------------------------------------------------------------

//! Unary plus.
//! @relates Vec3
template <class T> inline Vec3<T> operator+(const Vec3<T>& v)
{
    return v;
}

//! Unary minus.
//! @relates Vec3
template <class T> inline Vec3<T> operator-(const Vec3<T>& v)
{
    return {-v.x(), -v.y(), -v.z()};
}

// -----------------------------------------------------------------------------
// Binary operators
// -----------------------------------------------------------------------------

//! Addition of two vectors.
//! @relates Vec3
template <class T> inline Vec3<T> operator+(const Vec3<T>& a, const Vec3<T>& b)
{
    return {a.x() + b.x(), a.y() + b.y(), a.z() + b.z()};
}

//! Subtraction of two vectors.
//! @relates Vec3
template <class T> inline Vec3<T> operator-(const Vec3<T>& a, const Vec3<T>& b)
{
    return {a.x() - b.x(), a.y() - b.y(), a.z() - b.z()};
}

//! Multiplication vector by scalar.
//! @relates Vec3
#ifndef SWIG
template <class T, class U> inline auto operator*(const Vec3<T>& v, const U a)
{
    return Vec3<decltype(v.x() * v.x())>(v.x() * a, v.y() * a, v.z() * a);
}
#endif // USER_API

//! Multiplication scalar by vector.
//! @relates Vec3
#ifndef SWIG
template <class T, class U> inline auto operator*(const U a, const Vec3<T>& v)
{
    return Vec3<decltype(a * v.x())>(a * v.x(), a * v.y(), a * v.z());
}
#endif // USER_API

//! Division vector by scalar.
//! @relates Vec3
template <class T, class U> inline Vec3<T> operator/(const Vec3<T>& v, U a)
{
    return Vec3<T>(v.x() / a, v.y() / a, v.z() / a);
}

// =============================================================================
// ?? for API generation ??
// =============================================================================

//! Returns dot product of (complex) vectors (antilinear in the first [=self] argument).
#ifndef SWIG
template <class T> template <class U> inline auto Vec3<T>::dot(const Vec3<U>& v) const
{
    Vec3<T> left_star = this->conj();
    return left_star.x() * v.x() + left_star.y() * v.y() + left_star.z() * v.z();
}
#endif // USER_API

//! Returns cross product of (complex) vectors.
#ifndef SWIG
template <class T> template <class U> inline auto Vec3<T>::cross(const Vec3<U>& v) const
{
    return Vec3<decltype(this->x() * v.x())>(y() * v.z() - v.y() * z(), z() * v.x() - v.z() * x(),
                                             x() * v.y() - v.x() * y());
}
#endif // USER_API


using R3 = Vec3<double>;
using C3 = Vec3<complex_t>;

#endif // LIBFORMFACTOR_LIB_VEC3_H
