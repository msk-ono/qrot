#ifndef QROT_GATE_H
#define QROT_GATE_H

#include <compare>
#include <iostream>
#include <string>
#include <vector>

#include "qrot/matrix.h"

namespace qrot {
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
Gate operator*(Atom l, Atom r);
Gate operator*(const Gate& l, Atom r);
Gate operator*(Atom l, const Gate& r);
Gate operator*(const Gate& l, const Gate& r);
}  // namespace qrot

#endif  // QROT_GATE_H
