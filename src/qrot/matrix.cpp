#include "qrot/matrix.h"

namespace qrot {
#pragma region Algorithm
MD2 Adj2(const MD2& m) {
    return {m.Get(0, 0).Adj2(), m.Get(0, 1).Adj2(), m.Get(1, 0).Adj2(), m.Get(1, 1).Adj2()};
}
Mat ToMat(const MD2& m) {
    return {m.Get(0, 0).ToFloat(), m.Get(0, 1).ToFloat(), m.Get(1, 0).ToFloat(),
            m.Get(1, 1).ToFloat()};
}
MatC ToMatC(const MCD2& m) {
    return {Complex(m.Get(0, 0).Real().ToFloat(), m.Get(0, 0).Imag().ToFloat()),
            Complex(m.Get(0, 1).Real().ToFloat(), m.Get(0, 1).Imag().ToFloat()),
            Complex(m.Get(1, 0).Real().ToFloat(), m.Get(1, 0).Imag().ToFloat()),
            Complex(m.Get(1, 1).Real().ToFloat(), m.Get(1, 1).Imag().ToFloat())};
}
#pragma endregion Algorithm
#pragma region Explicit Instantiation
template class Matrix<Integer>;
template class Matrix<Float>;  // Mat
template class Matrix<DyadicFraction>;
template class Matrix<Z2>;
template class Matrix<D2>;  // MD2
template class Matrix<CD>;
template class Matrix<CZ2>;
template class Matrix<CD2>;  // MCD2
template class Matrix<ZOmega>;
#pragma endregion Explicit Instantiation
}  // namespace qrot
