//  ************************************************************************************************
//
//  libformfactor: efficient and accurate computation of scattering form factors
//
//! @file      lib/BasicVector3D.h
//! @brief     Declares and partly implements template class BasicVector3D.
//!
//! @homepage  http://www.bornagainproject.org
//! @license   GNU General Public License v3 or higher (see LICENSE)
//! @copyright Forschungszentrum Jülich GmbH 2021
//! @author    Joachim Wuttke, Scientific Computing Group at MLZ (see CITATION)
//!
//! Forked from CLHEP/Geometry by E. Chernyaev <Evgueni.Tcherniaev@cern.ch>,
//! then reworked beyond recognition. Removed split of point and vector semantics.
//! Transforms are relegated to a separate class Transform3D.
//
//  ************************************************************************************************

#ifndef LIBFORMFACTOR_LIB_BASICVECTOR3D_H
#define LIBFORMFACTOR_LIB_BASICVECTOR3D_H

#include "Base/Types/Complex.h"
#include <array>

//! Three-dimensional vector template, for use with integer, double, or complex components.

template <class T> class BasicVector3D : public std::array<T,3> {
private:
    using super = std::array<T,3>;

public:
    // -------------------------------------------------------------------------
    // Constructors and other set functions
    // -------------------------------------------------------------------------

    //! Constructs the null vector.
    BasicVector3D() : super{0., 0., 0.} {}

    //! Constructs a vector from cartesian components.
    BasicVector3D(const T x, const T y, const T z) : super{x,y,z} {}

    // -------------------------------------------------------------------------
    // Component access
    // -------------------------------------------------------------------------

    //! Returns x-component in cartesian coordinate system.
    inline T x() const { return (*this)[0]; }
    //! Returns y-component in cartesian coordinate system.
    inline T y() const { return (*this)[1]; }
    //! Returns z-component in cartesian coordinate system.
    inline T z() const { return (*this)[2]; }

    //! Sets x-component in cartesian coordinate system.
    void setX(const T& a) { (*this)[0] = a; }
    //! Sets y-component in cartesian coordinate system.
    void setY(const T& a) { (*this)[1] = a; }
    //! Sets z-component in cartesian coordinate system.
    void setZ(const T& a) { (*this)[2] = a; }

    // -------------------------------------------------------------------------
    // In-place operations
    // -------------------------------------------------------------------------

    //! Adds other vector to this, and returns result.
    BasicVector3D<T>& operator+=(const BasicVector3D<T>& v)
    {
        (*this)[0] += v[0];
        (*this)[1] += v[1];
        (*this)[2] += v[2];
        return *this;
    }

    //! Subtracts other vector from this, and returns result.
    BasicVector3D<T>& operator-=(const BasicVector3D<T>& v)
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
    BasicVector3D<T> conj() const;

    //! Returns magnitude squared of the vector.
    double mag2() const { return std::norm(x()) + std::norm(y()) + std::norm(z()); }

    //! Returns magnitude of the vector.
    double mag() const { return sqrt(mag2()); }

    //! Returns squared distance from z axis.
    double magxy2() const { return std::norm(x()) + std::norm(y()); }

    //! Returns distance from z axis.
    double magxy() const { return sqrt(magxy2()); }

    //! Returns azimuth angle.
    double phi() const;

    //! Returns polar angle.
    double theta() const;

    //! Returns cosine of polar angle.
    double cosTheta() const;

    //! Returns squared sine of polar angle.
    double sin2Theta() const;

    //! Returns unit vector in direction of this. Throws for null vector.
    BasicVector3D<T> unit() const;

    //! Returns this, trivially converted to complex type.
    BasicVector3D<complex_t> complex() const;

    //! Returns real parts.
    BasicVector3D<double> real() const;

    // -------------------------------------------------------------------------
    // Functions of this and another vector
    // -------------------------------------------------------------------------

    //! Returns dot product of vectors (antilinear in the first [=self] argument).
#ifndef SWIG
    template <class U> auto dot(const BasicVector3D<U>& v) const;
#endif // USER_API

    //! Returns cross product of vectors (linear in both arguments).
#ifndef SWIG
    template <class U> auto cross(const BasicVector3D<U>& v) const;
#endif // USER_API

    //! Returns angle with respect to another vector.
    double angle(const BasicVector3D<T>& v) const;

    //! Returns projection of this onto other vector: (this*v)*v/|v|^2.
    BasicVector3D<T> project(const BasicVector3D<T>& v) const
    {
        return dot(v) * v / v.mag2();
    }

    // -------------------------------------------------------------------------
    // Rotations
    // -------------------------------------------------------------------------

    //! Returns result of rotation around x-axis.
    BasicVector3D<T> rotatedX(double a) const;
    //! Returns result of rotation around y-axis.
    BasicVector3D<T> rotatedY(double a) const;
    //! Returns result of rotation around z-axis.
    BasicVector3D<T> rotatedZ(double a) const
    {
        return BasicVector3D<T>(cos(a) * x() + sin(a) * y(), -sin(a) * x() + cos(a) * y(), z());
    }
    //! Returns result of rotation around the axis specified by another vector.
    BasicVector3D<T> rotated(double a, const BasicVector3D<T>& v) const;
};

// =============================================================================
// Non-member functions
// =============================================================================

//! Output to stream.
//! @relates BasicVector3D
template <class T> std::ostream& operator<<(std::ostream& os, const BasicVector3D<T>& a)
{
    return os << "(" << a.x() << "," << a.y() << "," << a.z() << ")";
}

// -----------------------------------------------------------------------------
// Unary operators
// -----------------------------------------------------------------------------

//! Unary plus.
//! @relates BasicVector3D
template <class T> inline BasicVector3D<T> operator+(const BasicVector3D<T>& v)
{
    return v;
}

//! Unary minus.
//! @relates BasicVector3D
template <class T> inline BasicVector3D<T> operator-(const BasicVector3D<T>& v)
{
    return {-v.x(), -v.y(), -v.z()};
}

// -----------------------------------------------------------------------------
// Binary operators
// -----------------------------------------------------------------------------

//! Addition of two vectors.
//! @relates BasicVector3D
template <class T>
inline BasicVector3D<T> operator+(const BasicVector3D<T>& a, const BasicVector3D<T>& b)
{
    return {a.x() + b.x(), a.y() + b.y(), a.z() + b.z()};
}

//! Subtraction of two vectors.
//! @relates BasicVector3D
template <class T>
inline BasicVector3D<T> operator-(const BasicVector3D<T>& a, const BasicVector3D<T>& b)
{
    return {a.x() - b.x(), a.y() - b.y(), a.z() - b.z()};
}

//! Multiplication vector by scalar.
//! @relates BasicVector3D
#ifndef SWIG
template <class T, class U> inline auto operator*(const BasicVector3D<T>& v, const U a)
{
    return BasicVector3D<decltype(v.x() * v.x())>(v.x() * a, v.y() * a, v.z() * a);
}
#endif // USER_API

//! Multiplication scalar by vector.
//! @relates BasicVector3D
#ifndef SWIG
template <class T, class U> inline auto operator*(const U a, const BasicVector3D<T>& v)
{
    return BasicVector3D<decltype(a * v.x())>(a * v.x(), a * v.y(), a * v.z());
}
#endif // USER_API

// vector*vector not supported
//    (We do not provide the operator form a*b of the dot product:
//     Though nice to write, and in some cases perfectly justified,
//     in general it tends to make expressions more difficult to read.)

//! Division vector by scalar.
//! @relates BasicVector3D
template <class T, class U> inline BasicVector3D<T> operator/(const BasicVector3D<T>& v, U a)
{
    return BasicVector3D<T>(v.x() / a, v.y() / a, v.z() / a);
}

// =============================================================================
// ?? for API generation ??
// =============================================================================

//! Returns dot product of (complex) vectors (antilinear in the first [=self] argument).
#ifndef SWIG
template <class T>
template <class U>
inline auto BasicVector3D<T>::dot(const BasicVector3D<U>& v) const
{
    BasicVector3D<T> left_star = this->conj();
    return left_star.x() * v.x() + left_star.y() * v.y() + left_star.z() * v.z();
}
#endif // USER_API

//! Returns cross product of (complex) vectors.
#ifndef SWIG
template <class T>
template <class U>
inline auto BasicVector3D<T>::cross(const BasicVector3D<U>& v) const
{
    return BasicVector3D<decltype(this->x() * v.x())>(
        y() * v.z() - v.y() * z(), z() * v.x() - v.z() * x(), x() * v.y() - v.x() * y());
}
#endif // USER_API

#endif // LIBFORMFACTOR_LIB_BASICVECTOR3D_H