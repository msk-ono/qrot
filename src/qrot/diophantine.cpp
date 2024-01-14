#include "qrot/diophantine.h"

#include <algorithm>
#include <queue>
#include <unordered_map>

namespace qrot {
namespace {
Z2 CalcUnit(const D2& x, const D2& y) {
    // Assert x = o * y (o is unit in Z2)
    // Calculate x / y
    const auto num_den_exp = std::max(x.Int().DenExp(), x.Sqrt().DenExp());
    const auto den_den_exp = std::max(y.Int().DenExp(), y.Sqrt().DenExp());
    const auto den_exp = std::max(num_den_exp, den_den_exp);
    auto num = Z2((x.Int() << den_exp).Num(), (x.Sqrt() << den_exp).Num());
    auto den = Z2((y.Int() << den_exp).Num(), (y.Sqrt() << den_exp).Num());
    const auto norm = den.Norm();
    den.Adj2Inplace();
    num *= den;
    // std::cout << "Numerator = " << num.Int() << " " << num.Sqrt() << ", Denominator = " << norm
    //           << std::endl;
    assert(num.Int() % norm == 0 && "x/y must be in Z2");
    assert(num.Sqrt() % norm == 0 && "x/y must be in Z2");
    num.IntMut() /= norm;
    num.SqrtMut() /= norm;
    return num;
}
Z2 SqrtOfUnit(const Z2& x) {
#ifdef QROT_VERBOSE
    std::cout << "Calculate sqrt of " << x << ", norm = " << x.Norm() << std::endl;
#endif
    const auto& a = x.Int();
    const auto i1 = static_cast<Integer>(mp::sqrt((a + 1) / 2));  // FIXME: use int sqrt
    const auto i2 = static_cast<Integer>(mp::sqrt((a - 1) / 2));  // FIXME: use int sqrt
    const auto s1 = static_cast<Integer>(mp::sqrt((a - 1) / 4));  // FIXME: use int sqrt
    const auto s2 = static_cast<Integer>(mp::sqrt((a + 1) / 4));  // FIXME: use int sqrt
    const auto y1 = Z2(i1, s1);
    const auto y2 = Z2(i2, s2);
    const auto y3 = Z2(i1, -s1);
    const auto y4 = Z2(i2, -s2);
    if (y1 * y1 == x) {
        return y1;
    } else if (y2 * y2 == x) {
        return y2;
    } else if (y3 * y3 == x) {
        return y3;
    } else if (y4 * y4 == x) {
        return y4;
    }
    assert(0 && "Unreachable if x is unit");
    return Z2(0, 0);
}
bool Dividable(const Z2& x, const Z2& y) {
    // Is x / y is Z2 ?
    const auto norm = y.Norm();
    const auto num = x * y.Adj2();
    return num.Int() % norm == 0 && num.Sqrt() % norm == 0;
}
}  // namespace
Diophantine::Diophantine() {
    constexpr auto SearchLimit = std::size_t{10'000'000};
    auto is_prime = std::vector<bool>(SearchLimit, true);
    is_prime[0] = is_prime[1] = false;
    for (auto i = std::size_t{2}; i < SearchLimit; ++i) {
        if (is_prime[i]) {
            primes_.emplace_back(Integer(i));
            for (auto j = i * i; j < SearchLimit; j += i) { is_prime[j] = false; }
        }
    }
}
bool Diophantine::Solve(const D2& g, CD2& t) const {
    using constant::cd2::Delta;
    if (g < D2(0)) { return false; }
    if (g.Adj2() < D2(0)) { return false; }

    const auto tmp_den_exp = std::max(g.Int().DenExp(), g.Sqrt().DenExp());
    const auto den_exp = tmp_den_exp % 2 == 0 ? tmp_den_exp : tmp_den_exp + 1;
    assert(den_exp % 2 == 0);
    const auto num = Z2(g.Int().Num() << (den_exp - g.Int().DenExp()),
                        g.Sqrt().Num() << (den_exp - g.Sqrt().DenExp()));
    const auto norm = num.Norm();

#ifdef QROT_VERBOSE
    std::cout << "Factorize " << num << " (norm = " << norm << ")" << std::endl;
#endif

    auto fac = std::unordered_map<Integer, std::uint32_t>();
    FactorizeIntoPrime(norm, fac);

#ifdef QROT_VERBOSE
    std::cout << "Prime factorization = ";
    for (const auto& [p, n] : fac) { std::cout << "{" << p << ", " << n << "}"; }
    std::cout << std::endl;
#endif

    // NOTE: Prime factorization of large integers can fail.
    // If the prime factorization is successful, we can always find a solution for the case that
    // goes through the rough check.
    // If the prime factorization fails, the calculation proceeds with errors.
    // We return `g == t.Norm()` so that the program works well even if prime factorization fails.

    // Rough check
    for (const auto& [p, n] : fac) {
        if (p % 8 == 3 && n % 2 != 0) { return false; }
        if (p % 8 == 5 && n % 2 != 0) { return false; }
        if (p % 8 == 7 && n % 2 != 0) { return false; }
    }

    // If the rough check passes, a solution always exists (if prime factorization is successful)
    t = CD2(1, 0);
    try {
        for (const auto& [p, n] : fac) {
            const auto r = p % 8;

            if (n % 2 == 0) {
                if (r == 1) {
                    // r^2 = 2 mod p
                    // xi = gcd(p, r + \sqrt{2})
                    const auto r = SqrtMod(2, p);
                    auto xi = EuclidGCD(p, Z2(r, 1));  // xi or xi^adj
                    if (!Dividable(num, xi)) { xi.Adj2Inplace(); }
                    assert(Dividable(num, xi));
                    for (auto j = 0u; j < n / 2; ++j) { t *= CD2(D2(xi.Int(), xi.Sqrt())); }
                } else if (r == 3) {
                    // u^2 + 2 = 0 mod p
                    // x = gcd (xi, u + i \sqrt{2})
                    const auto u = SqrtMod(p - 2, p);
                    const auto x = ToCD2(EuclidGCD(p, ZOmega(u, 1, 0, 1)));
                    for (auto j = 0u; j < n / 2; ++j) { t *= x; }
                } else if (r == 5) {
                    // u^2 + 1 = 0 mod p
                    // x = gcd (xi, u + i)
                    const auto u = SqrtMod(p - 1, p);
                    const auto x = ToCD2(EuclidGCD(p, ZOmega(u, 0, 1, 0)));
                    for (auto j = 0u; j < n / 2; ++j) { t *= x; }
                } else if (r == 7) {
                    // r^2 = 2 mod p
                    // xi = gcd(p, r + \sqrt{2})
                    const auto r = SqrtMod(2, p);
                    auto xi = EuclidGCD(p, Z2(r, 1));  // xi or xi^adj
                    if (!Dividable(num, xi)) { xi.Adj2Inplace(); }
                    assert(Dividable(num, xi));
                    for (auto j = 0u; j < n / 2; ++j) { t *= CD2(D2(xi.Int(), xi.Sqrt())); }
                } else {
                    // p = 2
                    for (auto j = 0u; j < n; ++j) { t *= Delta; }
                }
                continue;
            }

            if (r == 1) {
                // r^2 = 2 mod p
                // xi = gcd(p, r + \sqrt{2})
                // u^2 + 1 = 0 mod p
                // x = gcd (xi, u + i)
                const auto r = SqrtMod(2, p);
                auto xi = EuclidGCD(p, Z2(r, 1));  // xi or xi^adj
                if (!Dividable(num, xi)) { xi.Adj2Inplace(); }
                assert(Dividable(num, xi));
                const auto u = SqrtMod(p - 1, p);
                const auto x = ToCD2(
                    EuclidGCD(ZOmega(xi.Int(), xi.Sqrt(), 0, -xi.Sqrt()), ZOmega(u, 0, 1, 0)));
                for (auto j = 0u; j < n; ++j) { t *= x; }
            } else if (r == 3) {
                assert(0 && "Unreachable");
            } else if (r == 5) {
                assert(0 && "Unreachable");
            } else if (r == 7) {
                assert(0 && "Unreachable");
            } else {
                // p = 2
                for (auto j = 0u; j < n; ++j) { t *= Delta; }
            }
        }

#ifdef QROT_VERBOSE
        std::cout << "g = " << g << ", norm = " << g.Norm() << std::endl;
        std::cout << "t = " << t << std::endl;
        std::cout << "t * t^adj = " << t * t.Adj() << ", norm = " << t.Norm().Norm() << std::endl;
        std::cout << "xi = " << num << std::endl;
#endif
        t.RealMut().IntMut() >>= (den_exp / 2);
        t.RealMut().SqrtMut() >>= (den_exp / 2);
        t.ImagMut().IntMut() >>= (den_exp / 2);
        t.ImagMut().SqrtMut() >>= (den_exp / 2);

        // Calculate unit `u` s.t. u * t is a solution of x^adj * x = g in Z[\omega]
        const auto u = SqrtOfUnit(CalcUnit(g, (t * t.Adj()).Real()));
        t *= CD2(D2(u.Int(), u.Sqrt()));
    } catch (std::exception& ex) {
        // Do nothing
    }

    return g == t.Norm();
}
void Diophantine::FactorizeIntoPrime(Integer n,
                                     std::unordered_map<Integer, std::uint32_t>& fac) const {
    // Factorize into small primes
    for (const auto& p : primes_) {
        auto exponent = std::uint32_t{0};
        auto q = Integer();
        auto r = Integer();
        mp::divide_qr(n, p, q, r);
        while (r == 0) {
            n = q;
            exponent++;
            mp::divide_qr(n, p, q, r);
        }
        if (exponent != 0) { fac[p] = exponent; }
    }
    if (n == 1) { return; }

    // Factorize into large primes
    auto queue = std::queue<Integer>();
    queue.push(n);
    while (!queue.empty()) {
        auto& n = queue.front();
        auto p = Integer{};
        const auto success = PollardRho(n, p);
        if (success) {
            n /= p;
            queue.push(p);
        } else {
            fac[n]++;
            queue.pop();
        }
    }
}
bool Diophantine::PollardRho(const Integer& n, Integer& p) const {
    // Assert: n is large (n >= 100)
    constexpr auto MaxLoops = std::uint32_t{10'000};
    // TODO: Use OpenMP
    for (auto offset = Integer{1}; offset < 100; ++offset) {
        const auto f = [&n, offset](const auto& x) { return ((x * x) + offset) % n; };

        auto x = Integer{2};
        auto y = Integer{2};
        for (auto i = std::uint32_t{0}; i < MaxLoops; ++i) {
            x = f(x);
            y = f(f(y));
            if (x == y) { break; }  // Find circuit
            p = mp::gcd(n, mp::abs(x - y));
            if (p > 1) { return true; }
        }
    }
    return false;
}
}  // namespace qrot
