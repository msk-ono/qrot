#include "qrot/grid_solver.h"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>

namespace qrot {
#pragma region OneDimGridSolver
bool OneDimGridSolver::Problem::IsValidSolution(const Float& a, const Float& b) const {
    using constant::f::Sqrt;
    static const Float Eps = std::numeric_limits<Float>::epsilon();
    // Check: x0 <= a + \sqrt{2} * b <= x1
    // Check: y0 <= a - \sqrt{2} * b <= y1
    const Float tmp = Sqrt * b;
    const Float tmp1 = a + tmp;
    const Float tmp2 = a - tmp;
    return x0 - Eps <= tmp1 && tmp1 <= x1 + Eps && y0 - Eps <= tmp2 && tmp2 <= y1 + Eps;
}
void OneDimGridSolver::Problem::DoLambda() {
    using constant::f::Lambda, constant::f::InvLambda;
    x0 *= Lambda;
    x1 *= Lambda;
    y0 *= -InvLambda;
    y1 *= -InvLambda;
    mp::swap(y0, y1);
    if (!history.empty() && history.back() == Conversion::DoInvLambda) {
        history.pop_back();
    } else {
        history.emplace_back(Conversion::DoLambda);
    }
}
void OneDimGridSolver::Problem::DoInvLambda() {
    using constant::f::Lambda, constant::f::InvLambda;
    x0 *= InvLambda;
    x1 *= InvLambda;
    y0 *= -Lambda;
    y1 *= -Lambda;
    mp::swap(y0, y1);
    if (!history.empty() && history.back() == Conversion::DoLambda) {
        history.pop_back();
    } else {
        history.emplace_back(Conversion::DoInvLambda);
    }
}
OneDimGridSolver::OneDimGridSolver(Float x0, Float x1, Float y0, Float y1)
    : problem_{std::move(x0), std::move(x1), std::move(y0), std::move(y1), {}} {
#ifdef QROT_VERBOSE
    if (problem_.x1 - problem_.x0 <= 0) { throw std::runtime_error("x1 must be larger than x0"); }
    if (problem_.y1 - problem_.y0 <= 0) { throw std::runtime_error("y1 must be larger than y0"); }
#endif
}
void OneDimGridSolver::EnumerateAllSolutions() {
    using namespace constant;
    using f::Sqrt, f::InvSqrt3, f::InvLambda;
    while (problem_.x1 - problem_.x0 >= 1) { problem_.DoInvLambda(); }
    while (problem_.x1 - problem_.x0 < InvLambda) { problem_.DoLambda(); }
    const Float min_b = mp::floor((problem_.x0 - problem_.y1) * InvSqrt3);
    const Float max_b = mp::ceil((problem_.x1 - problem_.y0) * InvSqrt3);
    for (Float tmp_b = min_b; tmp_b <= max_b; tmp_b += Float{1}) {
        const Float tmp_a = mp::floor(problem_.x1 - tmp_b * Sqrt);

        // Check & Push
        if (problem_.IsValidSolution(tmp_a, tmp_b)) {
            solutions_.emplace_back(Z2(static_cast<Integer>(tmp_a), static_cast<Integer>(tmp_b)));
        }
    }

    while (!problem_.history.empty()) {
        const auto& event = problem_.history.back();
        switch (event) {
            case Problem::Conversion::DoLambda: {
                problem_.DoInvLambda();
                for (auto&& solution : solutions_) { solution *= z2::InvLambda; }
                break;
            }
            case Problem::Conversion::DoInvLambda: {
                problem_.DoLambda();
                for (auto&& solution : solutions_) { solution *= z2::Lambda; }
                break;
            }
        }
    }
}
#pragma endregion
#pragma region FindGridOperation
namespace {
Float SinhL(const Float& x) {
    using constant::f::Lambda;
    return (mp::pow(Lambda, x) - mp::pow(Lambda, -x)) / 2;
}
Float CoshL(const Float& x) {
    using constant::f::Lambda;
    return (mp::pow(Lambda, x) + mp::pow(Lambda, -x)) / 2;
}
struct UnitGridOperation {
    enum class Type {
        Shift,
        R,
        K,
        A,
        B,
        Z,
        X,
    };

    static UnitGridOperation Shift(const Integer& n) { return {Type::Shift, n}; }
    static UnitGridOperation R() { return {Type::R}; }
    static UnitGridOperation K() { return {Type::K}; }
    static UnitGridOperation A(const Integer& n) { return {Type::A, n}; }
    static UnitGridOperation B(const Integer& n) { return {Type::B, n}; }
    static UnitGridOperation Z() { return {Type::Z}; }
    static UnitGridOperation X() { return {Type::X}; }

    MD2 GridOperator() const;

    Type type;
    Integer n = 0;
};
MD2 UnitGridOperation::GridOperator() const {
    static const auto HalfSqrt = D2{0, DyadicFraction{1, 1}};
    switch (type) {
        case UnitGridOperation::Type::Shift: throw std::runtime_error("Shift is not MD2 operator");
        case UnitGridOperation::Type::R: return MD2(HalfSqrt, -HalfSqrt, HalfSqrt, HalfSqrt);
        case UnitGridOperation::Type::K:
            return MD2(HalfSqrt - D2{1}, -HalfSqrt, HalfSqrt + D2{1}, HalfSqrt);
        case UnitGridOperation::Type::A: return MD2(1, D2(-2 * n), 0, 1);
        case UnitGridOperation::Type::B: return MD2(1, D2(0, n), 0, 1);
        case UnitGridOperation::Type::Z: return MD2(1, 0, 0, -1);
        case UnitGridOperation::Type::X: return MD2(0, 1, 1, 0);
    }
}
// std::ostream& operator<<(std::ostream& out, const UnitGridOperation& op) {
//     switch (op.type) {
//         case UnitGridOperation::Type::Shift: return out << "Shift(" << op.n << ')';
//         case UnitGridOperation::Type::R: return out << 'R';
//         case UnitGridOperation::Type::K: return out << 'K';
//         case UnitGridOperation::Type::A: return out << "A(" << op.n << ')';
//         case UnitGridOperation::Type::B: return out << "B(" << op.n << ')';
//         case UnitGridOperation::Type::Z: return out << 'Z';
//         case UnitGridOperation::Type::X: return out << 'X';
//     }
// }
struct EllipsePairState {
    Float e1, b1, z1;  // e1^2 - b1^2 = 1
    Float e2, b2, z2;  // e2^2 - b2^2 = 1

    Float Skew() const { return b1 * b1 + b2 * b2; }
    Float Bias() const { return z2 - z1; }
};
struct FindGridOperator {
    EllipsePairState state;
    std::vector<std::vector<UnitGridOperation>> history;

    void Find();
    MD2 GridOperator() const;

private:
    void Step();
    void R();
    void K();
    void A();
    void B();

    void Shift();
    void Z();
    void X();
};
void FindGridOperator::Find() {
    while (state.Skew() > Float{15}) { Step(); }
}
MD2 FindGridOperator::GridOperator() const {
    static const auto Lambda = Z2(1, 1);
    static const auto InvLambda = Z2(-1, 1);

    auto ret = MD2::Identity();
    for (auto itr1 = history.rbegin(); itr1 != history.rend(); ++itr1) {
        const auto& ops = *itr1;
        for (auto itr2 = ops.rbegin(); itr2 != ops.rend(); ++itr2) {
            const auto& op = *itr2;
            switch (op.type) {
                case UnitGridOperation::Type::Shift: {
                    assert(itr2 + 1 == ops.rend());
                    const auto x = op.n >= 0 ? Pow(Lambda, op.n) : Pow(InvLambda, -op.n);
                    const auto y = op.n >= 0 ? Pow(InvLambda, op.n) : Pow(Lambda, -op.n);
                    ret.GetMut(0, 0) *= D2(x.Int(), x.Sqrt());  // mul from left
                    ret.GetMut(1, 1) *= D2(y.Int(), y.Sqrt());  // mul from right
                    break;
                }
                default: ret.MulFromLeft(op.GridOperator());
            }
        }
    }
    return ret;
}
void FindGridOperator::Step() {
    static const Float m08 = Float{"-0.8"};
    static const Float m02 = Float{"-0.2"};
    static const Float p03 = Float{"0.3"};
    static const Float p08 = Float{"0.8"};

    history.emplace_back();

    Shift();
    Z();
    X();

    // Assert: -1 <= bias (state.z2 - state.z1) <= 1
    // Assert: 0 <= state.b2
    // Assert: 0 <= state.z1 + state.z2

    if (state.b1 >= 0) {
        // Case1
        if (m08 <= state.z1 && state.z1 <= p08 && m08 <= state.z2 && state.z2 <= p08) {
            R();
        } else if (state.z1 <= p03 && p08 <= state.z2) {
            K();
        } else if (p03 <= state.z1 && p03 <= state.z2) {
            A();
        } else if (p08 <= state.z1 && state.z2 <= p03) {
            K();
        } else {
            assert(0 && "Unreachable");
        }
    } else {
        // Case2
        if (m08 <= state.z1 && state.z1 <= p08 && m08 <= state.z2 && state.z2 <= p08) {
            R();
        } else if (m02 <= state.z1 && m02 <= state.z2) {
            B();
        } else {
            assert(0 && "Unreachable");
        }
    }
}
void FindGridOperator::R() {
    {
        const Float b = state.e1 * SinhL(state.z1);
        const Float x = state.e1 * CoshL(state.z1) + state.b1;
        const Float y = state.e1 * CoshL(state.z1) - state.b1;
        std::tie(state.e1, state.z1) = ToExponentFormat(x, y);
        state.b1 = b;
    }
    {
        const Float b = state.e2 * SinhL(state.z2);
        const Float x = state.e2 * CoshL(state.z2) + state.b2;
        const Float y = state.e2 * CoshL(state.z2) - state.b2;
        std::tie(state.e2, state.z2) = ToExponentFormat(x, y);
        state.b2 = b;
    }
    history.back().emplace_back(UnitGridOperation::R());
}
void FindGridOperator::K() {
    using constant::f::Sqrt;
    {
        const Float b = state.e1 * CoshL(state.z1 + 1) - Sqrt * state.b1;
        const Float x = state.e1 * CoshL(state.z1 + 2) - state.b1;
        const Float y = state.e1 * CoshL(state.z1) - state.b1;
        std::tie(state.e1, state.z1) = ToExponentFormat(x, y);
        state.b1 = b;
    }
    {
        const Float b = Sqrt * state.b2 - state.e2 * CoshL(state.z2 - 1);
        const Float x = state.e2 * CoshL(state.z2 - 2) - state.b2;
        const Float y = state.e2 * CoshL(state.z2) - state.b2;
        std::tie(state.e2, state.z2) = ToExponentFormat(x, y);
        state.b2 = b;
    }
    history.back().emplace_back(UnitGridOperation::K());
}
void FindGridOperator::A() {
    using constant::f::Lambda;
    Integer n = mp::max(
        Integer{1}, Integer{mp::floor(mp::pow(Lambda, std::min(state.z1, state.z2)) / Float{2})});
    const auto m = Float{n};
    {
        const Float x = state.e1 * mp::pow(Lambda, -state.z1);
        const Float b = state.b1 - 2 * m * x;
        const Float y = 4 * m * m * x - 4 * m * state.b1 + state.e1 * mp::pow(Lambda, state.z1);
        std::tie(state.e1, state.z1) = ToExponentFormat(x, y);
        state.b1 = b;
    }
    {
        const Float x = state.e2 * mp::pow(Lambda, -state.z2);
        const Float b = state.b2 - 2 * m * x;
        const Float y = 4 * m * m * x - 4 * m * state.b2 + state.e2 * mp::pow(Lambda, state.z2);
        std::tie(state.e2, state.z2) = ToExponentFormat(x, y);
        state.b2 = b;
    }
    history.back().emplace_back(UnitGridOperation::A(n));
}
void FindGridOperator::B() {
    using constant::f::Sqrt, constant::f::Lambda;
    Integer n = mp::max(Integer{1},
                        Integer{mp::floor(mp::pow(Lambda, std::min(state.z1, state.z2)) / Sqrt)});
    const auto m = Float{n};
    {
        const Float x = state.e1 * mp::pow(Lambda, -state.z1);
        const Float b = state.b1 + Sqrt * m * x;
        const Float y =
            2 * m * m * x + 2 * Sqrt * m * state.b1 + state.e1 * mp::pow(Lambda, state.z1);
        std::tie(state.e1, state.z1) = ToExponentFormat(x, y);
        state.b1 = b;
    }
    {
        const Float x = state.e2 * mp::pow(Lambda, -state.z2);
        const Float b = state.b2 - Sqrt * m * x;
        const Float y =
            2 * m * m * x - 2 * Sqrt * m * state.b2 + state.e2 * mp::pow(Lambda, state.z2);
        std::tie(state.e2, state.z2) = ToExponentFormat(x, y);
        state.b2 = b;
    }
    history.back().emplace_back(UnitGridOperation::B(n));
}
void FindGridOperator::Shift() {
    const auto bias = state.Bias();
    if (bias < -1 || 1 < bias) {
        Integer n = Integer{mp::floor((1 - bias) / 2)};
        state.z1 -= Float{n};
        state.z2 += Float{n};
        if (n % 2 != 0) { state.b2 *= -1; }
        history.back().emplace_back(UnitGridOperation::Shift(n));
    }
    assert(-1 <= state.Bias() && state.Bias() <= 1);
}
void FindGridOperator::Z() {
    if (state.b2 < 0) {
        state.b1 *= -1;
        state.b2 *= -1;
        history.back().emplace_back(UnitGridOperation::Z());
    }
}
void FindGridOperator::X() {
    if (state.z1 + state.z2 < 0) {
        state.z1 *= -1;
        state.z2 *= -1;
        history.back().emplace_back(UnitGridOperation::X());
    }
}
}  // namespace
#pragma endregion FindGridOperation
#pragma region TwoDimGridSolver
namespace {
void Mul(const MD2& m, CD2& p) {
    const D2 x = m.Get(0, 0) * p.Real() + m.Get(0, 1) * p.Imag();
    const D2 y = m.Get(1, 0) * p.Real() + m.Get(1, 1) * p.Imag();
    p = CD2(x, y);
}
void DivSqrt(CD2& p, std::uint32_t e) {
    // p /= \sqrt{2}^e
    p.RealMut().IntMut() >>= (e / 2);
    p.RealMut().SqrtMut() >>= (e / 2);
    p.ImagMut().IntMut() >>= (e / 2);
    p.ImagMut().SqrtMut() >>= (e / 2);
    if (e % 2 != 0) {
        p.RealMut().DivSqrt();
        p.ImagMut().DivSqrt();
    }
}
}  // namespace
TwoDimGridSolver TwoDimGridSolver::New(const Float& theta, const Float& epsilon) {
    // Calculate the edge coordinates of the rectangle
    const auto cos = mp::cos(theta);
    const auto sin = mp::sin(theta);
    const auto deps = epsilon * epsilon;
    const auto v1 = Vec(cos, sin);
    const auto v2 = Vec(sin, -cos);
    const auto t = epsilon * mp::sqrt(1 - epsilon / 4);
    const auto orig_el1 = Ellipse::FromRectangle((1 - deps / 2) * v1 - t * v2,  // a
                                                 (1 - deps / 2) * v1 + t * v2,  // b
                                                 v1 + t * v2,                   // c
                                                 v1 - t * v2                    // d
    );
    const auto orig_el2 = Ellipse::FromCircle(Vec(0, 0), 1);

    // Find grid operator
    auto finder = FindGridOperator();
    std::tie(finder.state.e1, finder.state.z1) = orig_el1.ExponentFormat();
    finder.state.b1 = orig_el1.B();
    std::tie(finder.state.e2, finder.state.z2) = orig_el2.ExponentFormat();
    finder.state.b2 = orig_el2.B();
    finder.Find();

    // Apply grid operator
    const auto inv_g1 = finder.GridOperator();
    const auto inv_g2 = Adj2(inv_g1);
    const auto g1 = inv_g1.Inv();
    const auto g2 = inv_g2.Inv();
    auto el1 = Ellipse::FromCircle(Vec(0, 0), 0);
    auto el2 = Ellipse::FromCircle(Vec(0, 0), 0);
    {
        // Calculate center of el1
        const auto g = ToMat(g1);
        const auto& p1 = orig_el1.Center();
        const auto mapped_p1 = g * p1;
        // Map ellipse pair
        const auto d1 = orig_el1.Mat();
        const auto d2 = orig_el2.Mat();
        const auto x1 = ToMat(inv_g1);
        const auto x2 = ToMat(inv_g2);
        const auto y1 = x1.Transpose() * d1 * x1;
        const auto y2 = x2.Transpose() * d2 * x2;
        el1 = Ellipse(mapped_p1, orig_el1.Scale(), y1.Get(0, 0), y1.Get(0, 1), y1.Get(1, 1));
        el2 = Ellipse(Vec(0, 0), orig_el2.Scale(), y2.Get(0, 0), y2.Get(0, 1), y2.Get(1, 1));
    }
    return Problem{theta,          epsilon,        cos, sin, orig_el1, orig_el2, el1, el2,
                   el1.CalcBBox(), el2.CalcBBox(), g1,  g2,  inv_g1,   inv_g2};
}
void TwoDimGridSolver::EnumerateAllSolutions() {
    using constant::f::Lambda, constant::f::InvLog2;
    static const Float Thresh = Lambda * Lambda;

    // Estimate good level
    const auto& p = problem_;
    const auto width =
        std::max(p.bbox1.XWidth() * p.bbox2.XWidth(), p.bbox1.YWidth() * p.bbox2.YWidth());
    level_ = static_cast<std::uint32_t>(mp::floor(mp::log(Thresh / width) * InvLog2));

    // Solve
    while (solutions_.empty()) {
        level_++;
        solutions_.clear();
        Solve();
    }
}
void TwoDimGridSolver::EnumerateNextLevelAllSolutions() {
    level_++;
    solutions_.clear();
    Solve();
}
void TwoDimGridSolver::Solve() {
    using constant::f::Sqrt, constant::f::InvSqrt, constant::cd2::Omega;

    // Solve upright 2-dim grid problem
    const auto& cos = problem_.cos;
    const auto& sin = problem_.sin;

    // Case1: a + b i
    {
        auto bbox1 = problem_.bbox1;
        auto bbox2 = problem_.bbox2;
        bbox1.Rescale(mp::pow(Sqrt, level_));
        bbox2.Rescale(mp::pow(-Sqrt, level_));
        auto x_solver = OneDimGridSolver(bbox1.x_min, bbox1.x_max, bbox2.x_min, bbox2.x_max);
        auto y_solver = OneDimGridSolver(bbox1.y_min, bbox1.y_max, bbox2.y_min, bbox2.y_max);
        x_solver.EnumerateAllSolutions();
        y_solver.EnumerateAllSolutions();
        for (const auto& x : x_solver.GetSolutions()) {
            for (const auto& y : y_solver.GetSolutions()) {
                auto p1 = CD2(ToD2(x), ToD2(y));
                auto p2 = CD2(ToD2(x.Adj2()), ToD2(y.Adj2()));
                // Rescale
                DivSqrt(p1, level_);
                DivSqrt(p2, level_);
                if (level_ % 2 != 0) { p2 = -p2; }
                // mapped -> original
                Mul(problem_.inv_g1_, p1);
                Mul(problem_.inv_g2_, p2);
                // Push valid solution
                auto is_valid = true;
                is_valid &= (p1.Norm() <= D2(1));
                is_valid &= (p1.Real().ToFloat() * cos + p1.Imag().ToFloat() * sin <= 1);
                is_valid &= (p2.Norm() <= D2(1));
                if (is_valid) { solutions_.emplace_back(p1); }
            }
        }
    }

    // a + b i + \omega
    {
        auto bbox1 = problem_.bbox1;
        auto bbox2 = problem_.bbox2;
        bbox1.Rescale(mp::pow(Sqrt, level_));
        bbox2.Rescale(mp::pow(-Sqrt, level_));
        bbox1.Translate(Vec(-InvSqrt, -InvSqrt));
        bbox2.Translate(Vec(InvSqrt, InvSqrt));
        auto x_solver = OneDimGridSolver(bbox1.x_min, bbox1.x_max, bbox2.x_min, bbox2.x_max);
        auto y_solver = OneDimGridSolver(bbox1.y_min, bbox1.y_max, bbox2.y_min, bbox2.y_max);
        x_solver.EnumerateAllSolutions();
        y_solver.EnumerateAllSolutions();
        for (const auto& x : x_solver.GetSolutions()) {
            for (const auto& y : y_solver.GetSolutions()) {
                auto p1 = CD2(ToD2(x), ToD2(y));
                auto p2 = CD2(ToD2(x.Adj2()), ToD2(y.Adj2()));
                // Translate
                p1 += Omega;
                p2 -= Omega;
                // Rescale
                DivSqrt(p1, level_);
                DivSqrt(p2, level_);
                if (level_ % 2 != 0) { p2 = -p2; }
                // mapped -> original
                Mul(problem_.inv_g1_, p1);
                Mul(problem_.inv_g2_, p2);
                // Push valid solutions
                auto is_valid = true;
                is_valid &= (p1.Norm() <= D2(1));
                is_valid &= (p1.Real().ToFloat() * cos + p1.Imag().ToFloat() * sin <= 1);
                is_valid &= (p2.Norm() <= D2(1));
                if (is_valid) { solutions_.emplace_back(p1); }
            }
        }
    }

#ifdef QROT_VERBOSE
    std::cout << "Solve level = " << level_ << " and found " << solutions_.size() << " solutions"
              << std::endl;
#endif
}
#pragma endregion
}  // namespace qrot
