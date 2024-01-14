#include "qrot/parser.h"

#include <gtest/gtest.h>

#include <numbers>

using namespace qrot;

TEST(Parser, Value) {
    const auto pi = std::numbers::pi;
    EXPECT_DOUBLE_EQ(pi / 128, static_cast<double>(AST::Parse("pi/128").Value()));
    EXPECT_DOUBLE_EQ(pi / 256, static_cast<double>(AST::Parse("pi/256").Value()));
    EXPECT_DOUBLE_EQ(-pi / 128, static_cast<double>(AST::Parse("-pi/128").Value()));
    EXPECT_DOUBLE_EQ(21, static_cast<double>(AST::Parse("5+20-4").Value()));
    EXPECT_DOUBLE_EQ(41, static_cast<double>(AST::Parse("12 + 34 - 5").Value()));
    EXPECT_DOUBLE_EQ(47, static_cast<double>(AST::Parse("5+6*7").Value()));
    EXPECT_DOUBLE_EQ(4, static_cast<double>(AST::Parse("(3+5)/2").Value()));
    EXPECT_DOUBLE_EQ(10, static_cast<double>(AST::Parse("-10+20").Value()));
    EXPECT_DOUBLE_EQ(-1.28, static_cast<double>(AST::Parse("-1.28").Value()));
}
TEST(Parser, Exceptions) {
    EXPECT_THROW(AST::Parse("pi 10"), std::runtime_error);
    EXPECT_THROW(AST::Parse("10 10"), std::runtime_error);
    EXPECT_THROW(AST::Parse("("), std::runtime_error);
    EXPECT_THROW(AST::Parse("()"), std::runtime_error);
    EXPECT_THROW(AST::Parse("((pi)"), std::runtime_error);
    EXPECT_THROW(AST::Parse("(pi))"), std::runtime_error);
    EXPECT_THROW(AST::Parse("- - 10"), std::runtime_error);
}
