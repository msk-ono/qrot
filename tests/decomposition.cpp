#include "qrot/decomposition.h"

#include <gtest/gtest.h>

using namespace qrot;

void TestUnitaryDecomposition(UnitaryDecomposer& decomposer, const MCD2& input,
                              const Gate& expected_output) {
    const auto actual_output = decomposer.Decompose(input);
    EXPECT_EQ(expected_output.Mat(), actual_output.Mat());
}

TEST(Decomposition, Unitary) {
    auto decomposer = UnitaryDecomposer();
    {
        auto input = Gate::FromString("T");
        TestUnitaryDecomposition(decomposer, input.Mat(), input);
    }
    {
        auto input = Gate::FromString("H");
        TestUnitaryDecomposition(decomposer, input.Mat(), input);
    }
    {
        auto input = Gate::FromString("TTT");
        TestUnitaryDecomposition(decomposer, input.Mat(), input);
    }
    {
        auto input = Gate::FromString("TH");
        TestUnitaryDecomposition(decomposer, input.Mat(), input);
    }
    {
        auto input = Gate::FromString("THT");
        TestUnitaryDecomposition(decomposer, input.Mat(), input);
    }
    {
        auto input = Gate::FromString("THTH");
        TestUnitaryDecomposition(decomposer, input.Mat(), input);
    }
    {
        auto input = Gate::FromString("THTTTHTHTTTHTHTHTHTTTHTHTTTHTTTHTHTHTHT");
        TestUnitaryDecomposition(decomposer, input.Mat(), input);
    }
    {
        auto input =
            Gate::FromString("TTTHTHTTTHTHTTTHTTTHTTTHTHTTTHTHTTTHTHTHTHTTTHTHTTTHTTTHTHTHTHT");
        TestUnitaryDecomposition(decomposer, input.Mat(), input);
    }
}
