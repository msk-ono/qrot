#include "qrot/geometry.h"

#include <gtest/gtest.h>

#include <limits>

using namespace qrot;

static const auto EPS = Float{std::numeric_limits<double>::epsilon()};

#define MY_EXPECT_FLOAT_EQ(lhs, rhs) \
    EXPECT_TRUE(mp::abs((lhs) - (rhs)) < EPS) << "lhs = " << (lhs) << ", rhs = " << (rhs)

TEST(Ellipse, FromRectangleClockwise) {
    // 5 x^2 - 6xy + 5 y^2 = 8
    const auto ellipse = Ellipse::FromRectangle(Vec(1, 0), Vec(0, 1), Vec(2, 3), Vec(3, 2));
    MY_EXPECT_FLOAT_EQ(mp::sqrt(Float{2}), ellipse.Scale());
    MY_EXPECT_FLOAT_EQ(Float("1.5"), ellipse.Center().X());
    MY_EXPECT_FLOAT_EQ(Float("1.5"), ellipse.Center().Y());
    MY_EXPECT_FLOAT_EQ(Float("1.25"), ellipse.A());
    MY_EXPECT_FLOAT_EQ(Float("-0.75"), ellipse.B());
    MY_EXPECT_FLOAT_EQ(Float("1.25"), ellipse.D());
    MY_EXPECT_FLOAT_EQ(Float("1.25"), ellipse.ExponentFormat().first);
    MY_EXPECT_FLOAT_EQ(0, ellipse.ExponentFormat().second);
    MY_EXPECT_FLOAT_EQ(Float{1}, ellipse.A() * ellipse.D() - ellipse.B() * ellipse.B());
    const auto bbox = ellipse.CalcBBox();
    const auto c = Float("1.5");
    const auto x = mp::sqrt(Float{10}) / 2;
    const auto y = mp::sqrt(Float{10}) / 2;
    MY_EXPECT_FLOAT_EQ(c - x, bbox.x_min);
    MY_EXPECT_FLOAT_EQ(c + x, bbox.x_max);
    MY_EXPECT_FLOAT_EQ(c - y, bbox.y_min);
    MY_EXPECT_FLOAT_EQ(c + y, bbox.y_max);
}
TEST(Ellipse, FromRectangleAntiClockwise) {
    // 5 x^2 - 6xy + 5 y^2 = 8
    const auto ellipse = Ellipse::FromRectangle(Vec(0, 1), Vec(1, 0), Vec(3, 2), Vec(2, 3));
    MY_EXPECT_FLOAT_EQ(mp::sqrt(Float{2}), ellipse.Scale());
    MY_EXPECT_FLOAT_EQ(Float("1.5"), ellipse.Center().X());
    MY_EXPECT_FLOAT_EQ(Float("1.5"), ellipse.Center().Y());
    MY_EXPECT_FLOAT_EQ(Float("1.25"), ellipse.A());
    MY_EXPECT_FLOAT_EQ(Float("-0.75"), ellipse.B());
    MY_EXPECT_FLOAT_EQ(Float("1.25"), ellipse.D());
    MY_EXPECT_FLOAT_EQ(Float("1.25"), ellipse.ExponentFormat().first);
    MY_EXPECT_FLOAT_EQ(0, ellipse.ExponentFormat().second);
    MY_EXPECT_FLOAT_EQ(1, ellipse.A() * ellipse.D() - ellipse.B() * ellipse.B());
    const auto bbox = ellipse.CalcBBox();
    const auto c = Float("1.5");
    const auto x = mp::sqrt(Float{10}) / 2;
    const auto y = mp::sqrt(Float{10}) / 2;
    MY_EXPECT_FLOAT_EQ(c - x, bbox.x_min);
    MY_EXPECT_FLOAT_EQ(c + x, bbox.x_max);
    MY_EXPECT_FLOAT_EQ(c - y, bbox.y_min);
    MY_EXPECT_FLOAT_EQ(c + y, bbox.y_max);
}
TEST(Ellipse, FromCircle) {
    const auto ellipse = Ellipse::FromCircle(Vec(1, -5), 10);
    MY_EXPECT_FLOAT_EQ(10, ellipse.Scale());
    MY_EXPECT_FLOAT_EQ(1, ellipse.Center().X());
    MY_EXPECT_FLOAT_EQ(-5, ellipse.Center().Y());
    MY_EXPECT_FLOAT_EQ(1, ellipse.A());
    MY_EXPECT_FLOAT_EQ(0, ellipse.B());
    MY_EXPECT_FLOAT_EQ(1, ellipse.D());
    const auto bbox = ellipse.CalcBBox();
    MY_EXPECT_FLOAT_EQ(1 - 10, bbox.x_min);
    MY_EXPECT_FLOAT_EQ(1 + 10, bbox.x_max);
    MY_EXPECT_FLOAT_EQ(-5 - 10, bbox.y_min);
    MY_EXPECT_FLOAT_EQ(-5 + 10, bbox.y_max);
}
