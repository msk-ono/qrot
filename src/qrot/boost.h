#ifndef QROT_BOOST_H
#define QROT_BOOST_H

#include "boost/multiprecision/cpp_complex.hpp"
#include "boost/multiprecision/cpp_dec_float.hpp"
#include "boost/multiprecision/cpp_int.hpp"

namespace qrot {
namespace mp = boost::multiprecision;
static constexpr auto FloatPrecision = std::size_t{1728};
using Integer = mp::cpp_int;
using FloatBackend =
    mp::cpp_bin_float<FloatPrecision, mp::backends::digit_base_2, void, std::int64_t>;
using Float = mp::number<FloatBackend, mp::et_off>;
using Complex = mp::number<mp::complex_adaptor<FloatBackend>, mp::et_off>;
namespace constant::f {
static inline const Float Pi = mp::default_ops::get_constant_pi<FloatBackend>();
static inline const Float Sqrt = mp::sqrt(Float{2});
static inline const Float InvSqrt = 1 / Sqrt;
static inline const Float Sqrt3 = Sqrt * Sqrt * Sqrt;
static inline const Float InvSqrt3 = 1 / Sqrt3;
static inline const Float Lambda = 1 + Sqrt;
static inline const Float InvLambda = -1 + Sqrt;
static inline const Float InvLog2 = 1 / mp::log(Float{2});
static inline const Float InvLogLambda = 1 / mp::log(Lambda);
}  // namespace constant::f
}  // namespace qrot

#endif  // QROT_BOOST_H
