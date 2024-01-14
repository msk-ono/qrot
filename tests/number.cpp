#include "qrot/number.h"

#include <gtest/gtest.h>

using namespace qrot;

TEST(Number, RingConcept) {
    static_assert(RealRingConcept<Integer>);
    static_assert(RealRingConcept<DyadicFraction>);
    static_assert(RealRingConcept<Z2>);
    static_assert(RealRingConcept<D2>);
    static_assert(RingConcept<CZ>);
    static_assert(RingConcept<CD>);
    static_assert(RingConcept<CD2>);
    static_assert(RingConcept<CZ2>);
    static_assert(RingConcept<ZOmega>);
    static_assert(RingConcept<DOmega>);
}
TEST(Number, EuclidGCDInZ2) {
    // Case1: prime % 8 == 1
    {
        // p = 17, 6^2 = 2 (mod 17)
        const auto prime = 17;
        ASSERT_EQ(1, prime % 8);
        const auto u = SqrtMod(2, prime);
        const auto gcd = EuclidGCD(Z2(prime), Z2(u, 1));
        EXPECT_EQ(prime, mp::abs(gcd.Norm()));
    }
    {
        // p = 41, 15^2 = 2 (mod 41)
        const auto prime = 41;
        ASSERT_EQ(1, prime % 8);
        const auto u = SqrtMod(2, prime);
        const auto gcd = EuclidGCD(Z2(prime), Z2(u, 1));
        EXPECT_EQ(prime, mp::abs(gcd.Norm()));
    }
    {
        // p = 73, 32^2 = 2 (mod 73)
        const auto prime = 73;
        ASSERT_EQ(1, prime % 8);
        const auto u = SqrtMod(2, prime);
        const auto gcd = EuclidGCD(Z2(prime), Z2(u, 1));
        EXPECT_EQ(prime, mp::abs(gcd.Norm()));
    }
    // Case2: prime % 8 == 7
    {
        // p = 7, 3^2 = 2 (mod 7)
        const auto prime = 7;
        ASSERT_EQ(7, prime % 8);
        const auto u = SqrtMod(2, prime);
        const auto gcd = EuclidGCD(Z2(prime), Z2(u, 1));
        EXPECT_EQ(prime, mp::abs(gcd.Norm()));
    }
    {
        // p = 23, 5^2 = 2 (mod 23)
        const auto prime = 23;
        ASSERT_EQ(7, prime % 8);
        const auto u = SqrtMod(2, prime);
        const auto gcd = EuclidGCD(Z2(prime), Z2(u, 1));
        EXPECT_EQ(prime, mp::abs(gcd.Norm()));
    }
    {
        // p = 31, 8^2 = 2 (mod 31)
        const auto prime = 31;
        ASSERT_EQ(7, prime % 8);
        const auto u = SqrtMod(2, prime);
        const auto gcd = EuclidGCD(Z2(prime), Z2(u, 1));
        EXPECT_EQ(prime, mp::abs(gcd.Norm()));
    }
}
TEST(Number, EuclidGCDInZOmega) {
    // Case1: prime % 8 == 5
    {
        const auto prime = 13;
        ASSERT_EQ(5, prime % 8);
        const auto u = SqrtMod(prime - 1, prime);
        const auto gcd = EuclidGCD(ZOmega(prime), ZOmega(u, 0, 1, 0));
        const auto mul = gcd * gcd.Adj();
        EXPECT_EQ(prime, mul.Get(0));
        EXPECT_EQ(0, mul.Get(1));
        EXPECT_EQ(0, mul.Get(2));
        EXPECT_EQ(0, mul.Get(3));
    }
    {
        const auto prime = 29;
        ASSERT_EQ(5, prime % 8);
        const auto u = SqrtMod(prime - 1, prime);
        const auto gcd = EuclidGCD(ZOmega(prime), ZOmega(u, 0, 1, 0));
        const auto mul = gcd * gcd.Adj();
        EXPECT_EQ(prime, mul.Get(0));
        EXPECT_EQ(0, mul.Get(1));
        EXPECT_EQ(0, mul.Get(2));
        EXPECT_EQ(0, mul.Get(3));
    }
    {
        const auto prime = 37;
        ASSERT_EQ(5, prime % 8);
        const auto u = SqrtMod(prime - 1, prime);
        const auto gcd = EuclidGCD(ZOmega(prime), ZOmega(u, 0, 1, 0));
        const auto mul = gcd * gcd.Adj();
        EXPECT_EQ(prime, mul.Get(0));
        EXPECT_EQ(0, mul.Get(1));
        EXPECT_EQ(0, mul.Get(2));
        EXPECT_EQ(0, mul.Get(3));
    }
    // Case2: prime % 8 == 3
    {
        const auto prime = 11;
        ASSERT_EQ(3, prime % 8);
        const auto u = SqrtMod(prime - 2, prime);
        const auto gcd = EuclidGCD(ZOmega(prime), ZOmega(u, 1, 0, 1));
        const auto mul = gcd * gcd.Adj();
        EXPECT_EQ(prime, mul.Get(0));
        EXPECT_EQ(0, mul.Get(1));
        EXPECT_EQ(0, mul.Get(2));
        EXPECT_EQ(0, mul.Get(3));
    }
    {
        const auto prime = 19;
        ASSERT_EQ(3, prime % 8);
        const auto u = SqrtMod(prime - 2, prime);
        const auto gcd = EuclidGCD(ZOmega(prime), ZOmega(u, 1, 0, 1));
        const auto mul = gcd * gcd.Adj();
        EXPECT_EQ(prime, mul.Get(0));
        EXPECT_EQ(0, mul.Get(1));
        EXPECT_EQ(0, mul.Get(2));
        EXPECT_EQ(0, mul.Get(3));
    }
    {
        const auto prime = 43;
        ASSERT_EQ(3, prime % 8);
        const auto u = SqrtMod(prime - 2, prime);
        const auto gcd = EuclidGCD(ZOmega(prime), ZOmega(u, 1, 0, 1));
        const auto mul = gcd * gcd.Adj();
        EXPECT_EQ(prime, mul.Get(0));
        EXPECT_EQ(0, mul.Get(1));
        EXPECT_EQ(0, mul.Get(2));
        EXPECT_EQ(0, mul.Get(3));
    }
}
TEST(Number, SqrtMod) {
    // Case1: prime % 8 == 1
    {
        const auto prime = 17;
        ASSERT_EQ(1, prime % 8);
        const auto x = SqrtMod(2, prime);
        EXPECT_EQ(2, (x * x) % prime);
    }
    {
        const auto prime = 41;
        ASSERT_EQ(1, prime % 8);
        const auto x = SqrtMod(2, prime);
        EXPECT_EQ(2, (x * x) % prime);
    }
    {
        const auto prime = 73;
        ASSERT_EQ(1, prime % 8);
        const auto x = SqrtMod(2, prime);
        EXPECT_EQ(2, (x * x) % prime);
    }
}
TEST(Number, ModPow) {
    EXPECT_EQ(1, ModPow(10, 0, 7));
    EXPECT_EQ(3, ModPow(10, 1, 7));
    EXPECT_EQ(2, ModPow(10, 2, 7));
    EXPECT_EQ(6, ModPow(10, 3, 7));
    EXPECT_EQ(4, ModPow(10, 4, 7));
    EXPECT_EQ(5, ModPow(10, 5, 7));
    EXPECT_EQ(1, ModPow(10, 6, 7));

    EXPECT_EQ(1, ModPow(10, Integer("60000000000000"), 7));
    EXPECT_EQ(3, ModPow(10, Integer("60000000000001"), 7));
    EXPECT_EQ(2, ModPow(10, Integer("60000000000002"), 7));
    EXPECT_EQ(6, ModPow(10, Integer("60000000000003"), 7));
    EXPECT_EQ(4, ModPow(10, Integer("60000000000004"), 7));
    EXPECT_EQ(5, ModPow(10, Integer("60000000000005"), 7));
    EXPECT_EQ(1, ModPow(10, Integer("60000000000006"), 7));
}
TEST(Number, Pow) {
    EXPECT_EQ(1, Pow(10, 0));
    EXPECT_EQ(10, Pow(10, 1));
    EXPECT_EQ(100, Pow(10, 2));
    EXPECT_EQ(1000, Pow(10, 3));
    EXPECT_EQ(10000, Pow(10, 4));
    EXPECT_EQ(100000, Pow(10, 5));
    EXPECT_EQ(1000000, Pow(10, 6));
}
