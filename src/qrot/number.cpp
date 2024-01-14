#include "qrot/number.h"

namespace qrot {
#pragma region Algorithm
namespace {
Integer RoundDiv(const Integer& num, const Integer& den) {
    const auto tmp = num + den / 2;
    Integer q, r;
    mp::divide_qr(tmp, den, q, r);
    return r < 0 ? q - 1 : q;
}
Z2 EuclidGCDImpl(const Z2& lhs, const Z2& rhs) {
    // Assert: Norm of lhs >= Norm of rhs
    if (rhs == Z2{0}) { return lhs; }

    // Calculate lhs / rhs in Q[\sqrt 2]
    const auto den = rhs.Norm();
    const auto num = lhs * rhs.Adj2();
    // Calculate the nearest integer of num.IntPart() / den
    Integer x = RoundDiv(num.Int(), den);
    // Calculate the nearest integer of num.SqrtPart() / den
    Integer y = RoundDiv(num.Sqrt(), den);

    return EuclidGCDImpl(rhs, lhs - Z2(x, y) * rhs);
}
ZOmega EuclidGCDImpl(const ZOmega& lhs, const ZOmega& rhs) {
    // Assert: Norm of lhs >= Norm of rhs
    if (rhs == ZOmega{0}) { return lhs; }

    const auto den = rhs.Norm();
    const auto num = lhs * rhs.Adj() * (rhs * rhs.Adj()).Adj2();
    const Integer a = RoundDiv(num.Get(0), den);
    const Integer b = RoundDiv(num.Get(1), den);
    const Integer c = RoundDiv(num.Get(2), den);
    const Integer d = RoundDiv(num.Get(3), den);
    return EuclidGCDImpl(rhs, lhs - ZOmega(a, b, c, d) * rhs);
}
std::pair<Integer, Integer> ModPow(const Integer& x, const Integer& y, const Integer& sqrt,
                                   Integer exp, const Integer& mod) {
    if (exp == 0) { return {1, 0}; }

    const auto mod_mul = [&sqrt, &mod](const Integer& lhs_int, const Integer& lhs_sqrt,
                                       const Integer& rhs_int,
                                       const Integer& rhs_sqrt) -> std::pair<Integer, Integer> {
        return {(lhs_int * rhs_int + sqrt * lhs_sqrt * rhs_sqrt) % mod,
                (lhs_int * rhs_sqrt + lhs_sqrt * rhs_int) % mod};
    };

    auto ret_int = Integer{1};
    auto ret_sqrt = Integer{0};
    auto pow_int = x;
    auto pow_sqrt = y;
    while (exp > 0) {
        if ((exp & 1) == 1) {
            std::tie(ret_int, ret_sqrt) = mod_mul(ret_int, ret_sqrt, pow_int, pow_sqrt);
        }
        std::tie(pow_int, pow_sqrt) = mod_mul(pow_int, pow_sqrt, pow_int, pow_sqrt);
        exp >>= 1;
    }
    return {ret_int, ret_sqrt};
}
}  // namespace
D2 ToD2(const Z2& x) { return D2(x.Int(), x.Sqrt()); }
CD2 ToCD2(const ZOmega& x) {
    using constant::cd2::Omega, constant::cd2::Imag, constant::cd2::Omega3;
    return CD2(D2(x.Get(0))) + Omega * CD2(D2(x.Get(1))) + Imag * CD2(D2(x.Get(2))) +
           Omega3 * CD2(D2(x.Get(3)));
}
CD2 ToCD2(const DOmega& x) {
    using constant::cd2::Omega, constant::cd2::Imag, constant::cd2::Omega3;
    return CD2(x.Get(0)) + Omega * CD2(x.Get(1)) + Imag * CD2(x.Get(2)) + Omega3 * CD2(x.Get(3));
}
Integer ModPow(Integer x, Integer exp, const Integer& mod) {
#ifdef QROT_VERBOSE
    if (exp < 0) { throw std::runtime_error("Exponent must be non-negative"); }
#endif

    auto ret = Integer{1};
    if (exp == 0) { return ret; }
    while (exp > 0) {
        if ((exp & 1) == 1) { ret = (ret * x) % mod; }
        x = (x * x) % mod;
        exp >>= 1;
    }
    return ret;
}
Z2 EuclidGCD(const Z2& lhs, const Z2& rhs) {
    const auto l_norm = mp::abs(lhs.Norm());
    const auto r_norm = mp::abs(rhs.Norm());
    return l_norm >= r_norm ? EuclidGCDImpl(lhs, rhs) : EuclidGCDImpl(rhs, lhs);
}
ZOmega EuclidGCD(const ZOmega& lhs, const ZOmega& rhs) {
    const auto l_norm = mp::abs(lhs.Norm());
    const auto r_norm = mp::abs(rhs.Norm());
    return l_norm >= r_norm ? EuclidGCDImpl(lhs, rhs) : EuclidGCDImpl(rhs, lhs);
}
Integer SqrtMod(const Integer& a, const Integer& p) {
#ifdef QROT_VERBOSE
    if (p <= 1) { throw std::runtime_error("`p` must be prime number"); }
    if (a < 0) { throw std::runtime_error("`a` must be non-negative"); }
    if (p <= a) { throw std::runtime_error("`a` must be smaller than prime number `p`"); }
#endif

    if (p == 2) { return a; }
    if (a == 0) { return 0; }
    if (ModPow(a, (p - 1) / 2, p) != 1) { return -1; }
    auto b = Integer{0};
    while (ModPow((b * b + p - a) % p, (p - 1) / 2, p) == 1) { b++; }
    return ModPow(b, 1, (b * b + p - a) % p, (p + 1) / 2, p).first;
}
#pragma endregion Algorithm
#pragma region Explicit Instantiation
template class SqrtRing<Integer>;            // Z2
template class SqrtRing<DyadicFraction>;     // D2
template class ComplexRing<Integer>;         // CZ
template class ComplexRing<DyadicFraction>;  // CD
template class ComplexRing<Z2>;              // CZ2
template class ComplexRing<D2>;              // CD2
template class OmegaRing<Integer>;           // ZOmega
template class OmegaRing<DyadicFraction>;    // DOmega
#pragma endregion Explicit Instantiation
}  // namespace qrot
