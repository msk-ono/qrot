#ifndef QROT_NUMBER_H
#define QROT_NUMBER_H

#include <concepts>
#include <cstdint>
#include <iostream>
#include <type_traits>

#include "qrot/boost.h"

namespace qrot {
#pragma region RingConcept
// clang-format off
// Concepts
template <typename T>
concept RingConcept = requires(T a, T b) {
    T();
    T{0};
    +a;
    -a;
    a + b;
    a - b;
    a * b;
    a += b;
    a -= b;
    a *= b;
    a == b;
    a != b;
};
template <typename T>
concept RealRingConcept = requires(T a, T b) {
    RingConcept<T>;
    a > b;
    a >= b;
    a < b;
    a <= b;
};
// clang-format on
#pragma endregion RingConcept
#pragma region DyadicFraction
/**
 * @brief Dyadic fraction: num / 2^{den_exp} (den_exp >= 0)
 */
class DyadicFraction {
public:
    DyadicFraction() : num_{0}, den_exp_{0} {}
    DyadicFraction(std::int32_t num) : num_{num}, den_exp_{0} {}
    DyadicFraction(const Integer& num) : num_{num}, den_exp_{0} {}
    DyadicFraction(const Integer& num, std::int32_t den_exp) : num_{num}, den_exp_{den_exp} {
        Normalize();
    }

    bool IsInteger() const { return den_exp_ == 0; }
    const Integer& Num() const { return num_; }
    std::int32_t DenExp() const { return den_exp_; }
    Float ToFloat() const { return Float{num_} / Float{Integer{1} << den_exp_}; }

    DyadicFraction operator+() const { return *this; }
    DyadicFraction operator-() const { return {-num_, den_exp_}; }
    DyadicFraction& operator+=(const DyadicFraction& rhs) {
        const auto max_den = std::max(den_exp_, rhs.den_exp_);
        num_ = (num_ << (max_den - den_exp_)) + (rhs.num_ << (max_den - rhs.den_exp_));
        den_exp_ = max_den;
        return Normalize();
    }
    DyadicFraction& operator-=(const DyadicFraction& rhs) {
        const auto exp_max = std::max(den_exp_, rhs.den_exp_);
        num_ = (num_ << (exp_max - den_exp_)) - (rhs.num_ << (exp_max - rhs.den_exp_));
        den_exp_ = exp_max;
        return Normalize();
    }
    DyadicFraction& operator*=(const DyadicFraction& rhs) {
        num_ *= rhs.num_;
        den_exp_ += rhs.den_exp_;
        return Normalize();
    }
    DyadicFraction& operator>>=(const std::size_t n) {
        den_exp_ += n;
        return Normalize();
    }
    DyadicFraction& operator<<=(const std::size_t n) {
        num_ <<= n;
        return Normalize();
    }

private:
    DyadicFraction& Normalize() {
        if (num_ == 0) {
            den_exp_ = 0;
            return *this;
        }
        while (den_exp_ > 0) {
            if ((num_ & 1) == 0) {
                num_ >>= 1;
                den_exp_--;
            } else {
                break;
            }
        }
        return *this;
    }

    Integer num_;
    std::int32_t den_exp_;
};
inline DyadicFraction operator+(const DyadicFraction& lhs, const DyadicFraction& rhs) {
    auto ret = lhs;
    ret += rhs;
    return ret;
}
inline DyadicFraction operator-(const DyadicFraction& lhs, const DyadicFraction& rhs) {
    auto ret = lhs;
    ret -= rhs;
    return ret;
}
inline DyadicFraction operator*(const DyadicFraction& lhs, const DyadicFraction& rhs) {
    auto ret = lhs;
    ret *= rhs;
    return ret;
}
inline DyadicFraction operator>>(const DyadicFraction& lhs, const std::size_t n) {
    auto ret = lhs;
    ret >>= n;
    return ret;
}
inline DyadicFraction operator<<(const DyadicFraction& lhs, const std::size_t n) {
    auto ret = lhs;
    ret <<= n;
    return ret;
}
inline bool operator==(const DyadicFraction& lhs, const DyadicFraction& rhs) {
    return lhs.Num() == rhs.Num() && lhs.DenExp() == rhs.DenExp();
}
inline bool operator!=(const DyadicFraction& lhs, const DyadicFraction& rhs) {
    return lhs.Num() != rhs.Num() || lhs.DenExp() != rhs.DenExp();
}
inline bool operator<(const DyadicFraction& lhs, const DyadicFraction& rhs) {
    const auto exp_max = std::max(lhs.DenExp(), rhs.DenExp());
    return (lhs.Num() << (exp_max - lhs.DenExp())) < (rhs.Num() << (exp_max - rhs.DenExp()));
}
inline bool operator<=(const DyadicFraction& lhs, const DyadicFraction& rhs) {
    const auto exp_max = std::max(lhs.DenExp(), rhs.DenExp());
    return (lhs.Num() << (exp_max - lhs.DenExp())) <= (rhs.Num() << (exp_max - rhs.DenExp()));
}
inline bool operator>(const DyadicFraction& lhs, const DyadicFraction& rhs) {
    const auto exp_max = std::max(lhs.DenExp(), rhs.DenExp());
    return (lhs.Num() << (exp_max - lhs.DenExp())) > (rhs.Num() << (exp_max - rhs.DenExp()));
}
inline bool operator>=(const DyadicFraction& lhs, const DyadicFraction& rhs) {
    const auto exp_max = std::max(lhs.DenExp(), rhs.DenExp());
    return (lhs.Num() << (exp_max - lhs.DenExp())) >= (rhs.Num() << (exp_max - rhs.DenExp()));
}
inline std::ostream& operator<<(std::ostream& out, const DyadicFraction& x) {
    if (x.IsInteger()) {
        out << x.Num();
    } else {
        out << x.Num() << "/2^" << x.DenExp();
    }
    return out;
}
namespace constant::d {
static inline const DyadicFraction Half = DyadicFraction{1, 1};
}  // namespace constant::d
#pragma endregion DyadicFraction
#pragma region SqrtRing
template <RealRingConcept Ring>
std::int32_t GetSign(const Ring& x) {
    auto sign = std::int32_t{0};
    if constexpr (std::is_same_v<Ring, Integer>) {
        sign = mp::sign(x);
    } else if constexpr (std::is_same_v<Ring, DyadicFraction>) {
        sign = mp::sign(x.Num());
    } else {
        if (x < 0) {
            sign = -1;
        } else if (x > 0) {
            sign = 1;
        } else {
            sign = 0;
        }
    }
    return sign;
}
/**
 * @brief RealRing + \sqrt{2} * RealRing
 * @details a + \sqrt{2} * b (a, b \in RealRing)
 */
template <RealRingConcept Ring>
class SqrtRing {
public:
    SqrtRing() : a_{0}, b_{0} {}
    SqrtRing(std::int32_t a) : a_{a} {}
    SqrtRing(const Ring& a) : a_{a}, b_{0} {}
    SqrtRing(const Ring& a, const Ring& b) : a_{a}, b_{b} {}

    const Ring& Int() const { return a_; }
    const Ring& Sqrt() const { return b_; }
    Ring& IntMut() { return a_; }
    Ring& SqrtMut() { return b_; }
    const Ring Norm() const { return a_ * a_ - 2 * b_ * b_; }
    SqrtRing Adj2() const { return {a_, -b_}; }
    Float ToFloat() const {
        using constant::f::Sqrt;
        if constexpr (std::is_same_v<Ring, Integer>) {
            return static_cast<Float>(a_) + static_cast<Float>(b_) * Sqrt;
        } else if constexpr (std::is_same_v<Ring, DyadicFraction>) {
            return a_.ToFloat() + b_.ToFloat() * Sqrt;
        } else {
            return static_cast<Float>(a_) + static_cast<Float>(b_) * Sqrt;
        }
    }

    void Adj2Inplace() { b_ = -b_; }
    void DivSqrt() {
        std::swap(a_, b_);
        b_ >>= 1;
    }

    SqrtRing operator+() const { return *this; }
    SqrtRing operator-() const { return {-a_, -b_}; }
    SqrtRing& operator+=(const SqrtRing& rhs) {
        a_ += rhs.a_;
        b_ += rhs.b_;
        return *this;
    }
    SqrtRing& operator-=(const SqrtRing& rhs) {
        a_ -= rhs.a_;
        b_ -= rhs.b_;
        return *this;
    }
    SqrtRing& operator*=(const SqrtRing& rhs) {
        const Ring a = a_ * rhs.a_ + 2 * b_ * rhs.b_;
        const Ring b = b_ * rhs.a_ + a_ * rhs.b_;
        a_ = a;
        b_ = b;
        return *this;
    }

private:
    Ring a_;
    Ring b_;
};
template <RealRingConcept Ring>
inline SqrtRing<Ring> operator+(const SqrtRing<Ring>& lhs, const SqrtRing<Ring>& rhs) {
    auto ret = lhs;
    ret += rhs;
    return ret;
}
template <RealRingConcept Ring>
inline SqrtRing<Ring> operator-(const SqrtRing<Ring>& lhs, const SqrtRing<Ring>& rhs) {
    auto ret = lhs;
    ret -= rhs;
    return ret;
}
template <RealRingConcept Ring>
inline SqrtRing<Ring> operator*(const SqrtRing<Ring>& lhs, const SqrtRing<Ring>& rhs) {
    auto ret = lhs;
    ret *= rhs;
    return ret;
}
template <RealRingConcept Ring>
inline bool operator==(const SqrtRing<Ring>& lhs, const SqrtRing<Ring>& rhs) {
    return lhs.Int() == rhs.Int() && lhs.Sqrt() == rhs.Sqrt();
}
template <RealRingConcept Ring>
inline bool operator!=(const SqrtRing<Ring>& lhs, const SqrtRing<Ring>& rhs) {
    return lhs.Int() != rhs.Int() || lhs.Sqrt() != rhs.Sqrt();
}
namespace impl {
/**
 * @brief Calculate: a + \sqrt{2} b > 0 (or >= 0 if include_zero is true)
 */
template <RealRingConcept Ring>
inline bool IsPositive(const Ring& a, const Ring& b, bool include_zero) {
    auto a_sign = GetSign(a);
    auto b_sign = GetSign(b);
    if (a_sign < 0) {
        if (b_sign < 0) {
            return false;
        } else if (b_sign == 0) {
            return false;
        } else {
            return Ring{2} * b * b > a * a;
        }
    } else if (a_sign == 0) {
        if (b_sign < 0) {
            return false;
        } else if (b_sign == 0) {
            return include_zero;
        } else {
            return true;
        }
    } else {
        if (b_sign < 0) {
            return a * a > Ring{2} * b * b;
        } else if (b_sign == 0) {
            return true;
        } else {
            return true;
        }
    }
    assert(0 && "Unreachable");
    return false;
}
}  // namespace impl
template <RealRingConcept Ring>
inline bool operator<(const SqrtRing<Ring>& lhs, const SqrtRing<Ring>& rhs) {
    return impl::IsPositive(rhs.Int() - lhs.Int(), rhs.Sqrt() - lhs.Sqrt(), false);
}
template <RealRingConcept Ring>
inline bool operator<=(const SqrtRing<Ring>& lhs, const SqrtRing<Ring>& rhs) {
    return impl::IsPositive(rhs.Int() - lhs.Int(), rhs.Sqrt() - lhs.Sqrt(), true);
}
template <RealRingConcept Ring>
inline bool operator>(const SqrtRing<Ring>& lhs, const SqrtRing<Ring>& rhs) {
    return impl::IsPositive(lhs.Int() - rhs.Int(), lhs.Sqrt() - rhs.Sqrt(), false);
}
template <RealRingConcept Ring>
inline bool operator>=(const SqrtRing<Ring>& lhs, const SqrtRing<Ring>& rhs) {
    return impl::IsPositive(lhs.Int() - rhs.Int(), lhs.Sqrt() - rhs.Sqrt(), true);
}
template <RealRingConcept Ring>
inline std::ostream& operator<<(std::ostream& out, const SqrtRing<Ring>& x) {
    if (x.Int() == Ring{0}) {
        const auto sign = GetSign(x.Sqrt());
        if (sign == 0) {
            out << "0";
        } else {
            out << x.Sqrt() << " √2";
        }
    } else {
        out << x.Int();
        const auto sign = GetSign(x.Sqrt());
        if (sign > 0) {
            out << " + " << x.Sqrt() << " √2";
        } else if (sign < 0) {
            out << " - " << Ring{sign} * x.Sqrt() << " √2";
        }
    }
    return out;
}
using Z2 = SqrtRing<Integer>;
using D2 = SqrtRing<DyadicFraction>;
namespace constant::z2 {
static inline const Z2 Sqrt = Z2{0, 1};
static inline const Z2 Lambda = Z2{1} + Sqrt;
static inline const Z2 InvLambda = Z2{-1} + Sqrt;
}  // namespace constant::z2
namespace constant::d2 {
static inline const D2 Sqrt = D2{0, 1};
static inline const D2 InvSqrt = D2{0, d::Half};
static inline const D2 Lambda = D2{1} + Sqrt;
static inline const D2 InvLambda = D2{-1} + Sqrt;
}  // namespace constant::d2
#pragma endregion SqrtRing
#pragma region ComplexRing
/**
 * @brief RealRing + i * RealRing
 * @details r + imaginary_unit * i (r, i \in RealRing)
 */
template <RealRingConcept Ring>
class ComplexRing {
public:
    ComplexRing() : r_{0}, i_{0} {}
    ComplexRing(std::int32_t r) : r_{r}, i_{0} {}
    ComplexRing(const Ring& r) : r_{r}, i_{0} {}
    ComplexRing(const Ring& r, const Ring& i) : r_{r}, i_{i} {}

    const Ring& Real() const { return r_; }
    const Ring& Imag() const { return i_; }
    Ring& RealMut() { return r_; }
    Ring& ImagMut() { return i_; }

    bool IsReal() const { return i_ == Ring{0}; }
    bool IsImag() const { return r_ == Ring{0}; }
    bool IsComplex() const { return r_ != Ring{0} && i_ != Ring{0}; }

    Ring Norm() const { return r_ * r_ + i_ * i_; }
    ComplexRing Adj() const { return {r_, -i_}; }

    void AdjInplace() { i_ = -i_; }

    ComplexRing operator+() const { return *this; }
    ComplexRing operator-() const { return {-r_, -i_}; }
    ComplexRing& operator+=(const ComplexRing& rhs) {
        r_ += rhs.r_;
        i_ += rhs.i_;
        return *this;
    }
    ComplexRing& operator-=(const ComplexRing& rhs) {
        r_ -= rhs.r_;
        i_ -= rhs.i_;
        return *this;
    }
    ComplexRing& operator*=(const ComplexRing& rhs) {
        const Ring r = r_ * rhs.r_ - i_ * rhs.i_;
        const Ring i = r_ * rhs.i_ + i_ * rhs.r_;
        r_ = r;
        i_ = i;
        return *this;
    }

private:
    Ring r_, i_;
};
template <RealRingConcept Ring>
ComplexRing<Ring> operator+(const ComplexRing<Ring>& lhs, const ComplexRing<Ring>& rhs) {
    auto ret = lhs;
    ret += rhs;
    return ret;
}
template <RealRingConcept Ring>
ComplexRing<Ring> operator-(const ComplexRing<Ring>& lhs, const ComplexRing<Ring>& rhs) {
    auto ret = lhs;
    ret -= rhs;
    return ret;
}
template <RealRingConcept Ring>
ComplexRing<Ring> operator*(const ComplexRing<Ring>& lhs, const ComplexRing<Ring>& rhs) {
    auto ret = lhs;
    ret *= rhs;
    return ret;
}
template <RealRingConcept Ring>
bool operator==(const ComplexRing<Ring>& lhs, const ComplexRing<Ring>& rhs) {
    return lhs.Real() == rhs.Real() && lhs.Imag() == rhs.Imag();
}
template <RealRingConcept Ring>
bool operator!=(const ComplexRing<Ring>& lhs, const ComplexRing<Ring>& rhs) {
    return lhs.Real() != rhs.Real() || lhs.Imag() != rhs.Imag();
}
template <RealRingConcept Ring>
inline std::ostream& operator<<(std::ostream& out, const ComplexRing<Ring>& x) {
    if (x.Real() == Ring{0}) {
        if (x.Imag() == Ring{0}) {
            out << "0";
        } else {
            out << x.Imag() << " i";
        }
    } else {
        if (x.Imag() == Ring{0}) {
            out << x.Real();
        } else {
            out << x.Real() << " + (" << x.Imag() << ") i";
        }
    }
    return out;
}
using CZ = ComplexRing<Integer>;
using CD = ComplexRing<DyadicFraction>;
using CZ2 = ComplexRing<Z2>;
using CD2 = ComplexRing<D2>;
namespace constant::cz2 {
static inline const CZ2 Imag = CZ2{0, 1};
static inline const CZ2 Sqrt = CZ2{Z2{0, 1}};
static inline const CZ2 Lambda = CZ2{1} + Sqrt;
static inline const CZ2 InvLambda = CZ2{-1} + Sqrt;
}  // namespace constant::cz2
namespace constant::cd2 {
static inline const CD2 Imag = CD2{0, 1};
static inline const CD2 Sqrt = CD2{D2{0, 1}};
static inline const CD2 InvSqrt = CD2{D2{0, d::Half}};
static inline const CD2 Lambda = CD2{1} + Sqrt;
static inline const CD2 InvLambda = CD2{-1} + Sqrt;
static inline const CD2 Omega = CD2{D2{0, d::Half}, D2{0, d::Half}};
static inline const CD2 Omega3 = CD2{D2{0, -d::Half}, D2{0, d::Half}};
static inline const CD2 Delta = CD2{1} + Omega;
static inline const CD2 InvDelta = CD2{1} + Omega3;
}  // namespace constant::cd2
#pragma endregion
#pragma region Omega
/**
 * @brief RealRing + omega * RealRing + omega^2 * RealRing + omega^3 * RealRing
 * @details omega = \sqrt{1/2} + \sqrt{1/2} i (omega^4=-1)
 * x_[0] + x_[1] * omega + x_[2] * omega^2 + x_[3] * omega^3
 */
template <RealRingConcept Ring>
class OmegaRing {
public:
    OmegaRing() : x_{Ring{0}, Ring{0}, Ring{0}, Ring{0}} {}
    OmegaRing(std::int32_t x) : x_{x, Ring{0}, Ring{0}, Ring{0}} {}
    OmegaRing(const Ring& x) : x_{x, Ring{0}, Ring{0}, Ring{0}} {}
    OmegaRing(const Ring& a, const Ring& b, const Ring& c, const Ring& d) : x_{a, b, c, d} {}

    const Ring& Get(std::size_t idx) const { return x_[idx]; }

    Ring Norm() const {
        const auto tmp = *this * Adj();
        return (tmp * tmp.Adj2()).x_[0];
    }
    OmegaRing Adj() const { return {x_[0], -x_[3], -x_[2], -x_[1]}; }
    OmegaRing Adj2() const { return {x_[0], -x_[1], x_[2], -x_[3]}; }

    void AdjInplace() { *this = Adj(); }
    void Adj2Inplace() { *this = Adj2(); }

    OmegaRing operator+() const { return *this; }
    OmegaRing operator-() const { return {-x_[0], -x_[1], -x_[2], -x_[3]}; }
    OmegaRing& operator+=(const OmegaRing& rhs) {
        for (auto i = std::size_t{0}; i < 4; ++i) { x_[i] += rhs.x_[i]; }
        return *this;
    }
    OmegaRing& operator-=(const OmegaRing& rhs) {
        for (auto i = std::size_t{0}; i < 4; ++i) { x_[i] -= rhs.x_[i]; }
        return *this;
    }
    OmegaRing& operator*=(const OmegaRing& rhs) {
        const Ring a =
            x_[0] * rhs.x_[0] - x_[1] * rhs.x_[3] - x_[2] * rhs.x_[2] - x_[3] * rhs.x_[1];
        const Ring b =
            x_[0] * rhs.x_[1] + x_[1] * rhs.x_[0] - x_[2] * rhs.x_[3] - x_[3] * rhs.x_[2];
        const Ring c =
            x_[0] * rhs.x_[2] + x_[1] * rhs.x_[1] + x_[2] * rhs.x_[0] - x_[3] * rhs.x_[3];
        const Ring d =
            x_[0] * rhs.x_[3] + x_[1] * rhs.x_[2] + x_[2] * rhs.x_[1] + x_[3] * rhs.x_[0];
        x_ = {a, b, c, d};
        return *this;
    }

private:
    std::array<Ring, 4> x_;
};
template <RealRingConcept Ring>
OmegaRing<Ring> operator+(const OmegaRing<Ring>& lhs, const OmegaRing<Ring>& rhs) {
    auto ret = lhs;
    ret += rhs;
    return ret;
}
template <RealRingConcept Ring>
OmegaRing<Ring> operator-(const OmegaRing<Ring>& lhs, const OmegaRing<Ring>& rhs) {
    auto ret = lhs;
    ret -= rhs;
    return ret;
}
template <RealRingConcept Ring>
OmegaRing<Ring> operator*(const OmegaRing<Ring>& lhs, const OmegaRing<Ring>& rhs) {
    auto ret = lhs;
    ret *= rhs;
    return ret;
}
template <RealRingConcept Ring>
bool operator==(const OmegaRing<Ring>& lhs, const OmegaRing<Ring>& rhs) {
    return lhs.Get(0) == rhs.Get(0) && lhs.Get(1) == rhs.Get(1) && lhs.Get(2) == rhs.Get(2) &&
           lhs.Get(3) == rhs.Get(3);
}
template <RealRingConcept Ring>
bool operator!=(const OmegaRing<Ring>& lhs, const OmegaRing<Ring>& rhs) {
    return lhs.Get(0) != rhs.Get(0) || lhs.Get(1) != rhs.Get(1) || lhs.Get(2) != rhs.Get(2) ||
           lhs.Get(3) != rhs.Get(3);
}
template <RealRingConcept Ring>
inline std::ostream& operator<<(std::ostream& out, const OmegaRing<Ring>& x) {
    return out << "omega [" << x.Get(0) << ',' << x.Get(1) << ',' << x.Get(2) << ',' << x.Get(3)
               << ']';
}
using ZOmega = OmegaRing<Integer>;
using DOmega = OmegaRing<DyadicFraction>;
namespace constant::zom {
static inline const ZOmega Imag = ZOmega{0, 0, 1, 0};
static inline const ZOmega Sqrt = ZOmega{0, 1, 0, -1};
static inline const ZOmega Lambda = ZOmega{1} + Sqrt;
static inline const ZOmega InvLambda = ZOmega{-1} + Sqrt;
static inline const ZOmega Omega = ZOmega{0, 1, 0, 0};
static inline const ZOmega Omega3 = ZOmega{0, 0, 0, 1};
static inline const ZOmega Delta = ZOmega{1} + Omega;
static inline const ZOmega InvDelta = ZOmega{1} + Omega3;
}  // namespace constant::zom
namespace constant::dom {
static inline const DOmega Imag = DOmega{0, 0, 1, 0};
static inline const DOmega Sqrt = DOmega{0, 1, 0, -1};
static inline const DOmega Lambda = DOmega{1} + Sqrt;
static inline const DOmega InvLambda = DOmega{-1} + Sqrt;
static inline const DOmega Omega = DOmega{0, 1, 0, 0};
static inline const DOmega Omega3 = DOmega{0, 0, 0, 1};
static inline const DOmega Delta = DOmega{1} + Omega;
static inline const DOmega InvDelta = DOmega{1} + Omega3;
}  // namespace constant::dom
#pragma endregion Omega
#pragma region Algorithm
D2 ToD2(const Z2& x);
CD2 ToCD2(const ZOmega& x);
CD2 ToCD2(const DOmega& x);
template <RingConcept Ring>
Ring Pow(Ring x, Integer e) {
#ifdef QROT_VERBOSE
    if (e < 0) { throw std::runtime_error("Exponent must be non-negative"); }
#endif

    auto ret = Ring{1};
    if (e == 0) { return ret; }
    while (e > 0) {
        if ((e & 1) == 1) { ret *= x; }
        x *= x;
        e >>= 1;
    }
    return ret;
}
Integer ModPow(Integer x, Integer exp, const Integer& mod);
/**
 * @brief Z2 is Euclidean domain.
 */
Z2 EuclidGCD(const Z2& lhs, const Z2& rhs);
/**
 * @brief ZOmega is Euclidean domain.
 */
ZOmega EuclidGCD(const ZOmega& lhs, const ZOmega& rhs);
/**
 * @brief Solve modular equation x^2 = a mod p using Cipolla algorithm
 * @details TODO: Use more efficient algorithm.
 *
 * @param a 0 <= a < p
 * @param p prime number
 * @return Integer x 0 <= x < p
 */
Integer SqrtMod(const Integer& a, const Integer& p);

#pragma endregion Algorithm
}  // namespace qrot

#endif  // QROT_NUMBER_H
