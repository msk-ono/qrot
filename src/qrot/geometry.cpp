#include "qrot/geometry.h"

namespace qrot {
std::pair<Float, Float> ToExponentFormat(const Float& a, const Float& d) {
    using constant::f::InvLogLambda;
    const Float e = mp::sqrt(a * d);
    const Float z = mp::log(d / e) * InvLogLambda;
    return {e, z};
}
std::pair<Float, Float> ToNormalFormat(const Float& e, const Float& z) {
    using constant::f::Lambda;
    return {e * mp::pow(Lambda, -z), e * mp::pow(Lambda, z)};
}
Ellipse Ellipse::FromRectangle(Vec a, Vec b, Vec c, Vec d) {
    using constant::f::InvSqrt;
    const auto abs = [](const Vec& v) { return mp::sqrt(v.X() * v.X() + v.Y() * v.Y()); };

    const auto center = (a + b + c + d) * Float("0.25");

    // Move to center
    a -= center;
    b -= center;
    c -= center;
    d -= center;

    // If the rectangle is rotated by -theta, the rectangle will be upright
    const auto tangent = (b - a) * (Float{1} / abs(b - a));
    // const auto theta = mp::atan2(tangent.imag(), tangent.real());
    const Float& cos = tangent.X();   // cos(-theta)
    const Float& sin = -tangent.Y();  // sin(-theta)
    const Float width = abs(b - a);
    const Float height = abs(c - b);
    const Float x = width * InvSqrt;
    const Float y = height * InvSqrt;

    // Get scale
    const Float scale = mp::sqrt(x * y);

    // Rotate -theta to get D
    const auto tmp_d = ::qrot::Mat(y / x, 0, 0, x / y);
    const auto rotate = ::qrot::Mat(cos, -sin, sin, cos);
    const auto D = rotate.Transpose() * tmp_d * rotate;
    return Ellipse(center, scale, D.Get(0, 0), D.Get(0, 1), D.Get(1, 1));
}
BBox Ellipse::CalcBBox() const {
    using constant::f::Lambda;

    // Ellipse equation:
    // a x^2 + 2b xy + d y^2 = s^2
    // Assert: da - b^2 = 1, a > 0, d > 0
    const auto& a = a_;
    const auto& d = d_;
    const auto& s = s_;

    auto bbox = BBox();
    if (b_ == 0) {
        // not rotated
        const Float x = s * mp::sqrt(d);
        const Float y = s * mp::sqrt(a);
        bbox.x_min = c_.X() - x;
        bbox.x_max = c_.X() + x;
        bbox.y_min = c_.Y() - y;
        bbox.y_max = c_.Y() + y;
    } else {
        // Calc x-limit
        const Float x = s * mp::sqrt(d);
        bbox.x_min = c_.X() - x;
        bbox.x_max = c_.X() + x;
        // Calc y-limit
        const Float y = s * mp::sqrt(a);
        bbox.y_min = c_.Y() - y;
        bbox.y_max = c_.Y() + y;
    }
    return bbox;
}
}  // namespace qrot
