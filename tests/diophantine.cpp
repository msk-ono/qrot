#include "qrot/diophantine.h"

#include <gtest/gtest.h>

using namespace qrot;

TEST(Diophantine, Factorize) {
    auto dio = Diophantine();
    const auto p1 = Integer(3257);
    const auto p2 = Integer(9277);
    const auto p3 = Integer("2586442777");
    const auto p4 = Integer("2586442787");
    const auto p5 = Integer("2586442813");
    {
        auto mp = std::unordered_map<Integer, std::uint32_t>();
        dio.FactorizeIntoPrime(p1 * p2 * p3, mp);
        EXPECT_TRUE(mp.contains(p1));
        EXPECT_EQ(1, mp.at(p1));
        EXPECT_TRUE(mp.contains(p2));
        EXPECT_EQ(1, mp.at(p2));
        EXPECT_TRUE(mp.contains(p3));
        EXPECT_EQ(1, mp.at(p3));
    }
    {
        auto mp = std::unordered_map<Integer, std::uint32_t>();
        dio.FactorizeIntoPrime(p1 * p2 * p3 * p3, mp);
        EXPECT_TRUE(mp.contains(p1));
        EXPECT_EQ(1, mp.at(p1));
        EXPECT_TRUE(mp.contains(p2));
        EXPECT_EQ(1, mp.at(p2));
        EXPECT_TRUE(mp.contains(p3));
        EXPECT_EQ(2, mp.at(p3));
    }
    {
        auto mp = std::unordered_map<Integer, std::uint32_t>();
        dio.FactorizeIntoPrime(p1 * p2 * p3 * p3 * p3, mp);
        EXPECT_TRUE(mp.contains(p1));
        EXPECT_EQ(1, mp.at(p1));
        EXPECT_TRUE(mp.contains(p2));
        EXPECT_EQ(1, mp.at(p2));
        EXPECT_TRUE(mp.contains(p3));
        EXPECT_EQ(3, mp.at(p3));
    }
    {
        auto mp = std::unordered_map<Integer, std::uint32_t>();
        dio.FactorizeIntoPrime(p1 * p2 * p3 * p3 * p3 * p4 * p4, mp);
        EXPECT_TRUE(mp.contains(p1));
        EXPECT_EQ(1, mp.at(p1));
        EXPECT_TRUE(mp.contains(p2));
        EXPECT_EQ(1, mp.at(p2));
        EXPECT_TRUE(mp.contains(p3));
        EXPECT_EQ(3, mp.at(p3));
        EXPECT_TRUE(mp.contains(p4));
        EXPECT_EQ(2, mp.at(p4));
    }
    {
        auto mp = std::unordered_map<Integer, std::uint32_t>();
        dio.FactorizeIntoPrime(p1 * p2 * p3 * p3 * p3 * p4 * p4 * p5 * p5 * p5, mp);
        EXPECT_TRUE(mp.contains(p1));
        EXPECT_EQ(1, mp.at(p1));
        EXPECT_TRUE(mp.contains(p2));
        EXPECT_EQ(1, mp.at(p2));
        EXPECT_TRUE(mp.contains(p3));
        EXPECT_EQ(3, mp.at(p3));
        EXPECT_TRUE(mp.contains(p4));
        EXPECT_EQ(2, mp.at(p4));
        EXPECT_TRUE(mp.contains(p5));
        EXPECT_EQ(3, mp.at(p5));
    }
}
TEST(Diophantine, SolveHard) {
    const auto u = DOmega(DyadicFraction(40727366, 26), DyadicFraction(10614512, 26),
                          DyadicFraction(10541729, 26), DyadicFraction(-26687414, 26));
    const auto t = DOmega(DyadicFraction(20133911, 26), DyadicFraction(2332111, 26),
                          DyadicFraction(-23432014, 26), DyadicFraction(30805761, 26));
    const auto x = u * u.Adj() + t * t.Adj();

    auto dio = Diophantine();
    const auto tmp = u * u.Adj();
    EXPECT_EQ(tmp.Get(1), -tmp.Get(3));
    EXPECT_EQ(0, tmp.Get(2));
    const auto g = D2(1) - D2(tmp.Get(0), tmp.Get(1));
    auto actual_t = CD2();
    dio.Solve(g, actual_t);
    EXPECT_EQ(g, (actual_t * actual_t.Adj()).Real());
}
TEST(Diophantine, SolveEasy) {
    const auto u = CD2(D2(DyadicFraction(-1, 2)), D2(DyadicFraction(-3, 2)));
    auto dio = Diophantine();
    const auto tmp = u * u.Adj();
    EXPECT_EQ(D2(0), tmp.Imag());
    const auto g = D2(1) - tmp.Real();
    auto actual_t = CD2();
    dio.Solve(g, actual_t);
    EXPECT_EQ(g, (actual_t * actual_t.Adj()).Real());
    EXPECT_EQ(CD2(1), u * u.Adj() + actual_t * actual_t.Adj());
}
