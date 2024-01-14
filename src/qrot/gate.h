#ifndef QROT_GATE_H
#define QROT_GATE_H

#include <compare>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "qrot/matrix.h"

namespace qrot {
#pragma region Atom
class Atom {
public:
    enum class Type { I, H, S, T, X, Y, Z, W };

    constexpr Atom() : t_{Type::I} {}
    explicit constexpr Atom(Type t) : t_{t} {}
    static constexpr Atom FromChar(char c);
    static constexpr Atom I() { return Atom{Type::I}; }
    static constexpr Atom H() { return Atom{Type::H}; }
    static constexpr Atom S() { return Atom{Type::S}; }
    static constexpr Atom T() { return Atom{Type::T}; }
    static constexpr Atom X() { return Atom{Type::X}; }
    static constexpr Atom Y() { return Atom{Type::Y}; }
    static constexpr Atom Z() { return Atom{Type::Z}; }
    static constexpr Atom W() { return Atom{Type::W}; }

    constexpr Type GetType() const { return t_; }
    constexpr char ToChar() const;
    constexpr bool IsClifford() const { return t_ != Type::T; }
    MCD2 Mat() const;

    constexpr auto operator<=>(const Atom&) const = default;

private:
    Type t_;
};
std::ostream& operator<<(std::ostream& out, Atom a);
namespace constant::gate {
static constexpr auto I = Atom::I();
static constexpr auto H = Atom::H();
static constexpr auto S = Atom::S();
static constexpr auto T = Atom::T();
static constexpr auto X = Atom::X();
static constexpr auto Y = Atom::Y();
static constexpr auto Z = Atom::Z();
static constexpr auto W = Atom::W();
}  // namespace constant::gate
#pragma endregion Atom
#pragma region Gate
class Gate {
public:
    using ConstIterator = std::vector<Atom>::const_iterator;

    Gate() = default;
    explicit Gate(Atom a) : atoms_{1, a} {}
    static Gate FromString(const std::string& s);

    bool Empty() const { return atoms_.empty(); }
    std::size_t Size() const { return atoms_.size(); }
    Atom Get(std::size_t idx) const { return atoms_[idx]; }
    std::size_t CountT() const;
    std::string ToString() const;
    bool IsClifford() const;
    MCD2 Mat() const;

    /**
     * @brief Normalize Cliffort+T gates.
     * @details See https://arxiv.org/abs/0806.3834 for more information.
     */
    void Normalize();

    ConstIterator begin() { return atoms_.begin(); }
    ConstIterator end() { return atoms_.end(); }
    ConstIterator begin() const { return atoms_.begin(); }
    ConstIterator end() const { return atoms_.end(); }

    Gate& operator*=(Atom a);
    Gate& operator*=(const Gate& g);

private:
    std::vector<Atom> atoms_;
};
std::ostream& operator<<(std::ostream& out, const Gate& g);
bool operator==(const Gate& lhs, const Gate& rhs);
bool operator==(const Gate& lhs, const Atom& rhs);
bool operator==(const Atom& lhs, const Gate& rhs);
bool operator!=(const Gate& lhs, const Gate& rhs);
bool operator!=(const Gate& lhs, const Atom& rhs);
bool operator!=(const Atom& lhs, const Gate& rhs);
Gate operator*(Atom l, Atom r);
Gate operator*(const Gate& l, Atom r);
Gate operator*(Atom l, const Gate& r);
Gate operator*(const Gate& l, const Gate& r);
#pragma endregion Gate
#pragma region CliffordDatabase
class CliffordDatabase {
public:
    enum class Type {
        CT,
        HCT,
        SHCT,
        NotClifford,
    };

    CliffordDatabase();

    static Type GetType(std::size_t idx);
    std::size_t SearchIndex(const MCD2& mat) const;
    const MCD2& GetMatrix(std::size_t idx) const { return c1_[idx].first; }
    const Gate& GetGate(std::size_t idx) const { return c1_[idx].second; }
    /**
     * @brief Get the index of T-moved gate.
     * @details
     * Let C be the gate of `idx`.
     * Calculate C' such that the equation C * T = P * T * C' holds true.
     * if      idx < 64,  P = I
     * else if idx < 128, P = H
     * else               P = SH
     *
     * Return the index of C'
     */
    std::size_t GetTMove(std::size_t idx) const;

private:
    static constexpr std::size_t NumElements = 192;
    static constexpr std::size_t NumCT = 64;

    // Aligned in C_T(C_1), H C_T(C_1), SH C_T(C_1) order
    std::vector<std::pair<MCD2, Gate>> c1_;
    std::vector<std::size_t> move_;  // T^dag C_T T
};
#pragma endregion CliffordDatabase
}  // namespace qrot

#endif  // QROT_GATE_H
