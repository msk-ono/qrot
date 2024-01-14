#ifndef QROT_GRID_SOLVER_H
#define QROT_GRID_SOLVER_H

#include <cstdint>
#include <vector>

#include "qrot/geometry.h"
#include "qrot/number.h"

namespace qrot {
#pragma region OneDimGridSolver
/**
 * @brief Solve one-dimensional grid-problems defined in 1403.2975.
 * @details Implementation of section 4 of 1403.2975.
 *
 * Find solutions a + \sqrt{2} b \in Z2 such that next two inequalities hold:
 * * x0 <= a + \sqrt{2} b <= x1
 * * y0 <= a - \sqrt{2} b <= y1
 */
class OneDimGridSolver {
public:
    OneDimGridSolver(Float x0, Float x1, Float y0, Float y1);

    void EnumerateAllSolutions();

    const std::vector<Z2>& GetSolutions() { return solutions_; }

private:
    struct Problem {
        Float x0, x1, y0, y1;
        enum class Conversion { DoLambda, DoInvLambda };
        std::vector<Conversion> history;

        bool IsValidSolution(const Float& a, const Float& b) const;
        void DoLambda();
        void DoInvLambda();
    };

    Problem problem_;
    std::vector<Z2> solutions_ = {};
};
#pragma endregion
#pragma region TwoDimGridSolver
/**
 * @brief Solve two-dimensional grid problems defined in 1403.2975.
 * @details Implementation of section 5 of 1403.2975.
 */
class TwoDimGridSolver {
public:
    static TwoDimGridSolver New(const Float& theta, const Float& epsilon);

    void EnumerateAllSolutions();
    void EnumerateNextLevelAllSolutions();

    const std::vector<CD2>& GetSolutions() { return solutions_; }

private:
    struct Problem {
        Float theta;
        Float epsilon;
        Float cos;         // cos(theta)
        Float sin;         // sin(theta)
        Ellipse orig_el1;  // original ellipse (covering circular segment)
        Ellipse orig_el2;  // original ellipse (unit circle)
        Ellipse el1;       // mapped ellipse (covering circular segment)
        Ellipse el2;       // mapped ellipse (unit circle)
        BBox bbox1;        // Bounding box of el1
        BBox bbox2;        // Bounding box of el2
        MD2 g1_;           // grid operator (orig -> mapped)
        MD2 g2_;           // grid operator (orig -> mapped)
        MD2 inv_g1_;       // inverse grid operator (mapped -> orig)
        MD2 inv_g2_;       // inverse grid operator (mapped -> orig)
    };

    TwoDimGridSolver(Problem&& problem) : problem_{std::move(problem)} {}
    void Solve();

    const Problem problem_;
    std::uint32_t level_ = 0;  //!< Search level
    std::vector<CD2> solutions_ = {};
};
#pragma endregion
}  // namespace qrot

#endif  // QROT_GRID_SOLVER_H
