#ifndef QROT_GEOMETRY_H
#define QROT_GEOMETRY_H

#include <utility>

#include "qrot/matrix.h"
#include "qrot/number.h"

namespace qrot {
struct BBox {
    Float x_min;
    Float x_max;
    Float y_min;
    Float y_max;

    Float XWidth() const { return x_max - x_min; }
    Float YWidth() const { return y_max - y_min; }
    Float Area() const { return (x_max - x_min) * (y_max - y_min); }

    void Translate(const Vec& v) {
        x_min += v.X();
        x_max += v.X();
        y_min += v.Y();
        y_max += v.Y();
    }
    void Rescale(const Float& s) {
        x_min *= s;
        x_max *= s;
        y_min *= s;
        y_max *= s;
        if (s < 0) {
            mp::swap(x_min, x_max);
            mp::swap(y_min, y_max);
        }
    }
};
/**
 * @brief Calculate exponent format.
 * @details Ellipse matrix uses two formats.
 * 1. Normal format:   a, d
 * 2. Exponent format: e*lambda^{-z}, e*lambda^z
 *
 * e = \sqrt{a * d}
 * z = \log_{\lambda}{a * d / e}
 */
std::pair<Float, Float> ToExponentFormat(const Float& a, const Float& d);
/**
 * @brief Calculate normal format.
 * @details Ellipse matrix uses two formats.
 * 1. Normal format:   a, d
 * 1. Exponent format: e*lambda^{-z}, e*lambda^z
 *
 * a = e*lambda^{-z}
 * d = e*lambda^z
 */
std::pair<Float, Float> ToNormalFormat(const Float& e, const Float& z);
/**
 * @brief Ellipse
 * @details Define ellipse in the quadratic form:
 * (x-c)^dag D (x-c) <= s^2
 *
 * D is defined by the next matrix:
 * D = a b
 *     b d
 * Assert: det D = 1
 */
class Ellipse {
public:
    Ellipse(const Vec& c, const Float& s, const Float& a, const Float& b, const Float& d)
        : c_{c}, s_{s}, a_{a}, b_{b}, d_{d} {}

    static Ellipse FromCircle(const Vec& c, const Float& r) { return Ellipse(c, r, 1, 0, 1); }
    static Ellipse FromRectangle(Vec a, Vec b, Vec c, Vec d);

    const Vec& Center() const { return c_; }
    const Float& Scale() const { return s_; }
    Mat Mat() const { return ::qrot::Mat(a_, b_, b_, d_); }
    const Float& A() const { return a_; }
    const Float& B() const { return b_; }
    const Float& D() const { return d_; }
    std::pair<Float, Float> NormalFormat() const { return {a_, d_}; }
    std::pair<Float, Float> ExponentFormat() const { return ToExponentFormat(a_, d_); }

    BBox CalcBBox() const;

    void Translate(const Vec& v) { c_ += v; }
    void Rescale(const Float& s) { s_ *= s; }

private:
    Vec c_;            // center
    Float s_;          // scale
    Float a_, b_, d_;  // matrix
};
}  // namespace qrot

#endif  // QROT_GEOMETRY_H
