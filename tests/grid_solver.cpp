#include "qrot/grid_solver.h"

#include <gtest/gtest.h>

#include <numbers>

#include "qrot/decomposition.h"
#include "qrot/diophantine.h"
#include "qrot/matrix.h"
#include "qrot/number.h"

using namespace qrot;

void TestOneDimGrid(double x0, double x1, double y0, double y1) {
    std::cout << "===---------------------===" << std::endl;
    std::cout << "Test 1-dim grid solver" << std::endl;
    auto solver = OneDimGridSolver(x0, x1, y0, y1);
    solver.EnumerateAllSolutions();
    const auto& solutions = solver.GetSolutions();
    EXPECT_TRUE(!solutions.empty());
    for (const auto& solution : solutions) {
        std::cout << solution << std::endl;
        EXPECT_LE(x0, solution.ToFloat());
        EXPECT_LE(solution.ToFloat(), x1);
        EXPECT_LE(y0, solution.Adj2().ToFloat());
        EXPECT_LE(solution.Adj2().ToFloat(), y1);
    }
}
void TestTwoDimGrid(double theta, double epsilon) {
    std::cout << "===---------------------===" << std::endl;
    std::cout << "Test 2-dim grid solver" << std::endl;
    auto solver = TwoDimGridSolver::New(theta, epsilon);
    solver.EnumerateAllSolutions();
    const auto& solutions = solver.GetSolutions();
    if (!solutions.empty()) {
        for (const auto& u : solutions) { std::cout << u << std::endl; }
    }
}

TEST(GridSolver, OneDim) {
    TestOneDimGrid(0.0, 1.1 + std::sqrt(2), 0, 1.1 + std::sqrt(2));
    TestOneDimGrid(1.0, 2.1 + std::sqrt(2), 0, 1.1 + std::sqrt(2));
    TestOneDimGrid(0.0, 1.1 + std::sqrt(2), 1, 28.1 + std::sqrt(2));
    TestOneDimGrid(1.0, 2.1 + std::sqrt(2), 1, 2.1 + std::sqrt(2));
}
TEST(GridSolver, TwoDim) {
    TestTwoDimGrid(std::numbers::pi / 128, 0.000001);
    EXPECT_EQ(0, 0);
}
