#ifndef QROT_DIOPHANTINE_H
#define QROT_DIOPHANTINE_H

#include <vector>

#include "qrot/number.h"

namespace qrot {
class Diophantine {
public:
    Diophantine();

    /**
     * @brief Solve diophantine equation t^adj * t = g.
     *
     * @param g input
     * @paragraph t output
     * @return true Find solution
     * @return false Cannot find solution
     */
    bool Solve(const D2& g, CD2& t) const;
    /**
     * @brief Calculate prime factorization.
     *
     * @param n input
     * @param fac factorization map (key: prime, value: exponent)
     */
    void FactorizeIntoPrime(Integer n, std::unordered_map<Integer, std::uint32_t>& fac) const;

private:
    /**
     * @brief Implementation of Pollard-Rho algorithm.
     * @details See https://qiita.com/Kiri8128/items/eca965fe86ea5f4cbb98 for more information.
     *
     * @return true if n % p == 0
     * @return false if the algorithm cannot find divisor of n
     */
    bool PollardRho(const Integer& n, Integer& p) const;

    std::vector<Integer> primes_;
};
}  // namespace qrot

#endif  // QROT_DIOPHANTINE_H
