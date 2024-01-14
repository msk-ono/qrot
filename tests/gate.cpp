#include "qrot/gate.h"

#include <gtest/gtest.h>

using namespace qrot;

TEST(Gate, Product) {
    using namespace constant::gate;
    EXPECT_EQ("XY", (X * Y).ToString());
    EXPECT_EQ("YZ", (Y * Z).ToString());
}
