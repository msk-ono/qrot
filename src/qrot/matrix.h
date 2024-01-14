#ifndef QROT_MATRIX_H
#define QROT_MATRIX_H

#include "qrot/boost.h"
#include "qrot/number.h"

namespace qrot {
#pragma region Matrix
template <RingConcept Ring>
class Matrix {
public:
    Matrix() : m_{0, 0, 0, 0} {}
    Matrix(std::int32_t s) : m_{s, 0, 0, s} {}
    Matrix(const Ring& s) : m_{s, 0, 0, s} {}
    Matrix(const Ring& a, const Ring& b, const Ring& c, const Ring& d) : m_{a, b, c, d} {}

    static Matrix Identity() { return {1, 0, 0, 1}; }

    const Ring& Get(std::size_t row, std::size_t col) const { return m_[2 * row + col]; }
    Ring& GetMut(std::size_t row, std::size_t col) { return m_[2 * row + col]; }
    Ring Det() const { return Get(0, 0) * Get(1, 1) - Get(0, 1) * Get(1, 0); }
    Matrix Transpose() const { return {Get(0, 0), Get(1, 0), Get(0, 1), Get(1, 1)}; }
    Matrix Inv() const {
        if constexpr (std::is_same_v<Float, Ring>) {
            const auto d = Det();
            return {Get(1, 1) / d, -Get(0, 1) / d, -Get(1, 0) / d, Get(0, 0) / d};
        } else {
            if (Det() == Ring{1}) {
                return {Get(1, 1), -Get(0, 1), -Get(1, 0), Get(0, 0)};
            } else if (Det() == Ring{-1}) {
                return {-Get(1, 1), Get(0, 1), Get(1, 0), -Get(0, 0)};
            }
#ifdef QROT_VERBOSE
            throw std::runtime_error("Cannot calculate inverse of non-special matrix");
#endif
            return {Get(1, 1), -Get(0, 1), -Get(1, 0), Get(0, 0)};
        }
    }

    Matrix operator+() const { return *this; }
    Matrix operator-() const { return {-Get(0, 0), -Get(0, 1), -Get(1, 0), -Get(1, 1)}; }
    Matrix& operator+=(const Matrix& rhs) {
        for (auto i = std::size_t{0}; i < 4; ++i) { m_[i] += rhs.m_[i]; }
        return *this;
    }
    Matrix& operator-=(const Matrix& rhs) {
        for (auto i = std::size_t{0}; i < 4; ++i) { m_[i] -= rhs.m_[i]; }
        return *this;
    }
    Matrix& MulFromLeft(const Matrix& lhs);
    Matrix& operator*=(const Matrix& rhs);
    Matrix& operator*=(const Ring& s) {
        for (auto i = std::size_t{0}; i < 4; ++i) { m_[i] *= s; }
        return *this;
    }

private:
    //        col 0   col 1
    // row 0  m_[0]   m_[1]
    // row 1  m_[2]   m_[3]
    std::array<Ring, 4> m_;
};
template <RingConcept Ring>
Matrix<Ring> operator+(const Matrix<Ring>& lhs, const Matrix<Ring>& rhs) {
    auto ret = lhs;
    ret += rhs;
    return ret;
}
template <RingConcept Ring>
Matrix<Ring> operator-(const Matrix<Ring>& lhs, const Matrix<Ring>& rhs) {
    auto ret = lhs;
    ret -= rhs;
    return ret;
}
template <RingConcept Ring>
Matrix<Ring> operator*(const Matrix<Ring>& lhs, const Matrix<Ring>& rhs) {
    const Ring a = lhs.Get(0, 0) * rhs.Get(0, 0) + lhs.Get(0, 1) * rhs.Get(1, 0);
    const Ring b = lhs.Get(0, 0) * rhs.Get(0, 1) + lhs.Get(0, 1) * rhs.Get(1, 1);
    const Ring c = lhs.Get(1, 0) * rhs.Get(0, 0) + lhs.Get(1, 1) * rhs.Get(1, 0);
    const Ring d = lhs.Get(1, 0) * rhs.Get(0, 1) + lhs.Get(1, 1) * rhs.Get(1, 1);
    return Matrix<Ring>(a, b, c, d);
}
template <RingConcept Ring>
Matrix<Ring> operator*(const Matrix<Ring>& m, const Ring& s) {
    return Matrix<Ring>(m.Get(0, 0) * s, m.Get(0, 1) * s, m.Get(1, 0) * s, m.Get(1, 1) * s);
}
template <RingConcept Ring>
Matrix<Ring> operator*(const Ring& s, const Matrix<Ring>& m) {
    return Matrix<Ring>(m.Get(0, 0) * s, m.Get(0, 1) * s, m.Get(1, 0) * s, m.Get(1, 1) * s);
}
template <RingConcept Ring>
Matrix<Ring>& Matrix<Ring>::MulFromLeft(const Matrix<Ring>& lhs) {
    *this = lhs * (*this);
    return *this;
}
template <RingConcept Ring>
Matrix<Ring>& Matrix<Ring>::operator*=(const Matrix<Ring>& rhs) {
    *this = (*this) * rhs;
    return *this;
}
template <RingConcept Ring>
bool operator==(const Matrix<Ring>& lhs, const Matrix<Ring>& rhs) {
    return lhs.Get(0, 0) == rhs.Get(0, 0) && lhs.Get(0, 1) == rhs.Get(0, 1) &&
           lhs.Get(1, 0) == rhs.Get(1, 0) && lhs.Get(1, 1) == rhs.Get(1, 1);
}
template <RingConcept Ring>
bool operator!=(const Matrix<Ring>& lhs, const Matrix<Ring>& rhs) {
    return lhs.Get(0, 0) != rhs.Get(0, 0) || lhs.Get(0, 1) != rhs.Get(0, 1) ||
           lhs.Get(1, 0) != rhs.Get(1, 0) || lhs.Get(1, 1) != rhs.Get(1, 1);
}
template <RingConcept Ring>
inline std::ostream& operator<<(std::ostream& out, const Matrix<Ring>& x) {
    return out << "matrix [" << x.Get(0, 0) << ',' << x.Get(0, 1) << ',' << x.Get(1, 0) << ','
               << x.Get(1, 1) << ']';
}
using Mat = Matrix<Float>;
using MatC = Matrix<Complex>;
using MD2 = Matrix<D2>;
using MCD2 = Matrix<CD2>;
namespace constant::mcd2 {
static inline const MCD2 I = MCD2(1, 0, 0, 1);
static inline const MCD2 H = MCD2(cd2::InvSqrt, cd2::InvSqrt, cd2::InvSqrt, -cd2::InvSqrt);
static inline const MCD2 S = MCD2(1, 0, 0, cd2::Imag);
static inline const MCD2 T = MCD2(1, 0, 0, CD2(D2(0, d::Half), D2(0, d::Half)));
static inline const MCD2 X = MCD2(0, 1, 1, 0);
static inline const MCD2 Y = MCD2(0, -cd2::Imag, cd2::Imag, 0);
static inline const MCD2 Z = MCD2(1, 0, 0, -1);
static inline const MCD2 W = MCD2(cd2::Omega, 0, 0, cd2::Omega);
static inline const MCD2 TDag = MCD2(1, 0, 0, CD2(D2(0, d::Half), D2(0, -d::Half)));
}  // namespace constant::mcd2
#pragma endregion Matrix
#pragma region Vector
template <RingConcept Ring>
class Vector {
public:
    Vector() : v_{Ring{0}, Ring{0}} {}
    Vector(const Ring& a, const Ring& b) : v_{a, b} {}

    const Ring& X() const { return v_[0]; }
    const Ring& Y() const { return v_[1]; }
    const Ring& Get(std::int32_t i) const { return v_[i]; }
    Ring& GetMut(std::int32_t i) { return v_[i]; }

    Vector operator+() const { return *this; }
    Vector operator-() const { return {-Get(0), -Get(1)}; }
    Vector& operator+=(const Vector& rhs) {
        for (auto i = std::size_t{0}; i < 2; ++i) { v_[i] += rhs.v_[i]; }
        return *this;
    }
    Vector& operator-=(const Vector& rhs) {
        for (auto i = std::size_t{0}; i < 2; ++i) { v_[i] -= rhs.v_[i]; }
        return *this;
    }
    Vector& operator*=(const Ring& s) {
        for (auto i = std::size_t{0}; i < 2; ++i) { v_[i] *= s; }
        return *this;
    }
    Vector& operator*=(const Matrix<Ring>& m) { return *this = (*this) * m; }
    Vector& MulFromLeft(const Matrix<Ring>& m) { return *this = m * (*this); }

private:
    std::array<Ring, 2> v_;
};
template <RingConcept Ring>
Vector<Ring> operator+(const Vector<Ring>& lhs, const Vector<Ring>& rhs) {
    auto ret = lhs;
    ret += rhs;
    return ret;
}
template <RingConcept Ring>
Vector<Ring> operator-(const Vector<Ring>& lhs, const Vector<Ring>& rhs) {
    auto ret = lhs;
    ret -= rhs;
    return ret;
}
template <RingConcept Ring>
Vector<Ring> operator*(const Vector<Ring>& v, const Matrix<Ring>& m) {
    return Vector<Ring>(v.Get(0) * m.Get(0, 0) + v.Get(1) * m.Get(1, 0),
                        v.Get(0) * m.Get(0, 1) + v.Get(1) * m.Get(1, 1));
}
template <RingConcept Ring>
Vector<Ring> operator*(const Matrix<Ring>& m, const Vector<Ring>& v) {
    return Vector<Ring>(m.Get(0, 0) * v.Get(0) + m.Get(0, 1) * v.Get(1),
                        m.Get(1, 0) * v.Get(0) + m.Get(1, 1) * v.Get(1));
}
template <RingConcept Ring>
Vector<Ring> operator*(const Vector<Ring>& v, const Ring& s) {
    return Vector<Ring>(v.Get(0) * s, v.Get(1) * s);
}
template <RingConcept Ring>
Vector<Ring> operator*(const Ring& s, const Vector<Ring>& v) {
    return Vector<Ring>(v.Get(0) * s, v.Get(1) * s);
}
template <RingConcept Ring>
inline std::ostream& operator<<(std::ostream& out, const Vector<Ring>& v) {
    return out << "vector [" << v.X() << ',' << v.Y() << ']';
}
using Vec = Vector<Float>;
using VecC = Vector<Complex>;
using VD2 = Vector<D2>;
using VCD2 = Vector<CD2>;
#pragma endregion Vector
#pragma region Algorithm
MD2 Adj2(const MD2& m);
Mat ToMat(const MD2& m);
MatC ToMatC(const MCD2& m);
#pragma endregion Algorithm
}  // namespace qrot

#endif  // QROT_MATRIX_H
