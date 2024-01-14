#include "qrot/gate.h"

#include <gtest/gtest.h>

#include <queue>
#include <utility>
#include <vector>

#include "qrot/matrix.h"

using namespace qrot;

TEST(Gate, Product) {
    using namespace constant::gate;
    EXPECT_EQ("XY", (X * Y).ToString());
    EXPECT_EQ("YZ", (Y * Z).ToString());
}
TEST(Gate, HS) {
    using namespace constant::mcd2;
    EXPECT_EQ(X, Gate::FromString("HSSH").Mat());
    EXPECT_EQ(Z, Gate::FromString("SS").Mat());
    EXPECT_EQ(W, Gate::FromString("HSHSHS").Mat());
}
TEST(CliffordDatabase, Ctor) {
    using constant::mcd2::H, constant::mcd2::S;
    auto database = CliffordDatabase();
}
TEST(Gate, Normalize) {
    // Example from https://www.mathstat.dal.ca/~selinger/newsynth/
    const auto input = std::string(
        "SHTHTHTHTHTHTHTSHTHTHTHTSHTHTHTHTHTSHTSHTHTHTHTHTSHTHTHTSHTSHTHTSHTSHTSHTHTHTHTSHTHTHTHTHT"
        "SHTSHTHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTSHTHTHTSHTSHTSHTHTHTHTSHTHTSHTHTHTSHTHTHTHTSHTHTSHTHTS"
        "HTSHTSHTHTHTHTHTHTHTSHTHTSHTHTHTSHTSHTHTHTSHTSHTSHTHTSHTHTHTHTSHTSHTSHSSSWWWWWWW");
    auto gate = Gate::FromString(input);
    const auto matrix = gate.Mat();
    gate.Normalize();
    EXPECT_EQ(matrix, gate.Mat());
    EXPECT_EQ(Gate::FromString("SHSSSWWWWWWW").Mat(), Gate::FromString("SHSSXSSSXW").Mat());
}
