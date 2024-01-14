#include "qrot/matrix.h"

#include <gtest/gtest.h>

using namespace qrot;

MCD2 CalcMatrix(const std::string& str, bool inverse = false) {
    using namespace constant::mcd2;
    auto ret = MCD2::Identity();
    for (const auto& c : str) {
        if (c == 'T') {
            ret *= inverse ? TDag : T;
        } else if (c == 'S') {
            ret *= inverse ? S * S * S : S;
        } else if (c == 'H') {
            ret *= H;
        } else if (c == 'X') {
            ret *= X;
        } else if (c == 'W') {
            ret *= inverse ? W * W * W * W * W * W * W : W;
        }
    }
    return ret;
}
void Show(const CD2& complex) {
    std::cout << complex.Real().ToFloat() << "+" << complex.Imag().ToFloat() << " i" << std::endl;
}
void Show(const MCD2& mat) {
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "(0,0): ";
    Show(mat.Get(0, 0));
    std::cout << "(0,1): ";
    Show(mat.Get(0, 1));
    std::cout << "(1,0): ";
    Show(mat.Get(1, 0));
    std::cout << "(1,1): ";
    Show(mat.Get(1, 1));
}

TEST(Matrix, RingConcept) {
    static_assert(RingConcept<Matrix<Integer>>);
    static_assert(RingConcept<Matrix<Float>>);
    static_assert(RingConcept<Matrix<DyadicFraction>>);
    static_assert(RingConcept<Matrix<Z2>>);
    static_assert(RingConcept<Matrix<D2>>);
    static_assert(RingConcept<Matrix<CD>>);
    static_assert(RingConcept<Matrix<CZ2>>);
    static_assert(RingConcept<Matrix<CD2>>);
    static_assert(RingConcept<Matrix<ZOmega>>);
}
TEST(Matrix, Basic) {
    using Mat = Matrix<Integer>;
    const auto x = Mat(1, 2, 3, 4);
    const auto y = Mat(2, 0, 0, 2);
    const auto z = x * y;
    EXPECT_EQ(2, z.Get(0, 0));
    EXPECT_EQ(4, z.Get(0, 1));
    EXPECT_EQ(6, z.Get(1, 0));
    EXPECT_EQ(8, z.Get(1, 1));
}
TEST(Matrix, Multiplication) {
    using namespace constant::mcd2;
    const auto input = std::string("THTTTHTTTHTHTTHTHTT");
    const auto reverse = std::string(input.rbegin(), input.rend());

    EXPECT_EQ(I, T * TDag);
    EXPECT_EQ(I, TDag * T);
    EXPECT_EQ(I, H * H);
    EXPECT_EQ(I, X * X);
    EXPECT_EQ(I, Y * Y);
    EXPECT_EQ(I, Z * Z);
    EXPECT_EQ(Z, S * S);
    EXPECT_EQ(I, CalcMatrix(input, true) * CalcMatrix(reverse));
    EXPECT_EQ(I, CalcMatrix(reverse, true) * CalcMatrix(input));
    EXPECT_EQ(I, CalcMatrix(input) * CalcMatrix(reverse, true));
    EXPECT_EQ(I, CalcMatrix(reverse) * CalcMatrix(input, true));
}
TEST(Matrix, Pi128) {
    const auto input = std::string(
        "HTHTHTHTHTHSTHSTHTHSTHSTHTHSTHTHSTHSTHSTHSTHSTHTHTHSTHTHSTHSTHTHSTHTHTHTHTHTHTHSTHSTHS"
        "THSTHTHTHSTHSTHTHTHTHSTHSTHTHSTHSTHTHSTHTHSTHTHTHTHTHTHSTHSTHTHTHTHTHSTHSTHSTHTHSTHSTH"
        "STHTHTHSTHTHSTHSTHTHSTHSTHTHTHTHSTHSTHTHSTHTHSTHSTHSTHTHTHSTHTHSTHTHTHSTHSTHSTHSTHSTHT"
        "HSTHTHTHSTHTHSTHSTHSTHTHTHTHSTHTHSTHSTHTHTHTHTHSTHTHTHTHSTHTHTHTHSTHTHSTHSTHSTHSTHSTHT"
        "HTHSTHSTHTHTHSTHTHTHSTHSTHSTHTHSTHTHTHTHSTHSTHTHSTHTHSTHTHSTHSTHTHSTHTHSTHSTHSTHSTHTHS"
        "THTHSTHTHTHTHSTHTHTHTHTHTHSTHTHTHTHSTHSTHSTHSTHSTHTHSTHSTHSTHSTHTHSTHTHSTHTHSTHSTHTHTH"
        "STHSTHSTHTHTHTHSTHSTHTHTHTHTHTHSTHTHSTHTHSTHSTHTHSTHSTHTHSTHSTHTHSTHTHTHSTHSTHSTHSTHST"
        "HTHTHTHSTHSTHTHSTHSTHSTHTHTHTHSTHTHSTHSTHTHTHSTHSTHTHSTHSTHSTHTHSTHSTHSTHSTHTHSTHSTHST"
        "HSTHSTHTHTHTHSTHSTHSTHTHTHSTHSTHSTHTHSTHSTHTHTHSTHTHSTHSTHSTHSTHSTHSTHTHTHSTHSTHTHSTHS"
        "THTHSTHTHSTHTHSTHSTHSTHTHTHSTHTHTHTHTHTHSTHTHSTHSTHSTHSTHTHTHTHSTHSTHSTHTHSTHTHTHTHSTH"
        "THTHSTHTHSTHTHTHSTHTHSTHTHTHSTHSTHSTHSTHSTHTHTHSTHTHSTHSTHSTHTHTHSTHSTHTHTHTHTHTHSTHST"
        "HTHTHSTHTHTHSTHSTHSTHTHTHTHSTHTHTHSTHSTHSTHSTHSTHTHSTHSTHSTHTHTHTHTHTHSTHSTHTHSTHTHSTH"
        "THTHSTHTHSTHSTHTHTHSTHSTHTHSTHSTHSTHTHSTHTHTHSTHTHSTHSTHTHSTHSTHTHTHSTHSTHTHSTHTHSTHST"
        "HTHTHSTHSTHSTHSTHTHTHSTHTHSTHSTHSTHTHTHTHSTHTHTHTHSTHTHTHSTHTHSTHTHSTHSTHSTHTHSTHSTHST"
        "HTHSTHSTHTHSTHTHTHSTHSTHSTHSTHTHTHTHTHTHTHTHSTHSTHSTHTHSTHSTHTHSTHSTHTHTHSTHSTHTHTHSTH"
        "STHTHTHSTHSTHTHTHSTHSTHSTHSTHSTHTHTHTHTHTHSTHTHSTHTHTHSTHTHSTHTHTHTHSTHTHTHTHSTHSTHSTH"
        "STHSTHSTHSTHTHSTHTHTHSTHSTHTHSTHTHSTHSTHSTHTHTHSTHTHTHSTHTHSTHTHSTHSTHSTHTHSTHSTHTHSTH"
        "STHTHSTHTHTHTHSTHTHTHSTHSTHTHSTHTHTHTHTHSTHSTHTHTHSTHSTHSTHSTHSTHSTHTHTHTHSTHTHSTHTHTH"
        "THTHSTHTHSTHTHTHSTHSTHSTHSTHTHSTHSTHSTHTHSTHTHTHSTHTHSTHSTHSTHTHTHSTHSTHSTHTHSTHTHSTHS"
        "THSTHTHSTHSTHTHTHTHSTHTHTHTHSTHSTHTHSTHTHSTHSTHTHSTHSTHTHSTHSTHTHTHTHSTHSTHTHTHTHTHTHT"
        "HTHSTHTHSTHTHSTHTHSTHTHSTHSTHTHSTHSTHTHTHSTHS");
    const auto x = CalcMatrix(input);
    std::cout << ToMatC(x) << std::endl;
}
