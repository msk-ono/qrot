#include "qrot/gate.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qrot {
constexpr Atom Atom::FromChar(char c) {
    switch (c) {
        case 'i':
        case 'I': return Atom{Type::I};
        case 'h':
        case 'H': return Atom{Type::H};
        case 's':
        case 'S': return Atom{Type::S};
        case 't':
        case 'T': return Atom{Type::T};
        case 'x':
        case 'X': return Atom{Type::X};
        case 'y':
        case 'Y': return Atom{Type::Y};
        case 'z':
        case 'Z': return Atom{Type::Z};
        case 'w':
        case 'W': return Atom{Type::W};
    }
    throw std::runtime_error("Unknown char: " + std::to_string(c));
}
constexpr char Atom::ToChar() const {
    switch (t_) {
        case Type::I: return 'I';
        case Type::H: return 'H';
        case Type::S: return 'S';
        case Type::T: return 'T';
        case Type::X: return 'X';
        case Type::Y: return 'Y';
        case Type::Z: return 'Z';
        case Type::W: return 'W';
    }
}
MCD2 Atom::Mat() const {
    switch (t_) {
        case Type::I: return constant::mcd2::I;
        case Type::H: return constant::mcd2::H;
        case Type::S: return constant::mcd2::S;
        case Type::T: return constant::mcd2::T;
        case Type::X: return constant::mcd2::X;
        case Type::Y: return constant::mcd2::Y;
        case Type::Z: return constant::mcd2::Z;
        case Type::W: return constant::mcd2::W;
    }
}
std::ostream& operator<<(std::ostream& out, Atom a) { return out << a.ToChar(); }
Gate Gate::FromString(const std::string& s) {
    auto gate = Gate();
    gate.atoms_.resize(s.size());
    std::transform(s.begin(), s.end(), gate.atoms_.begin(),
                   [](const char c) { return Atom::FromChar(c); });
    return gate;
}
std::size_t Gate::CountT() const {
    using constant::gate::T;
    auto ret = std::size_t{0};
    for (const auto atom : atoms_) { ret += (atom == T) ? 1 : 0; }
    return ret;
}
std::string Gate::ToString() const {
    if (Empty()) { return "I"; }
    auto ret = std::string();
    ret.resize(atoms_.size());
    std::transform(atoms_.begin(), atoms_.end(), ret.begin(),
                   [](const Atom a) { return a.ToChar(); });
    return ret;
}
bool Gate::IsClifford() const {
    using constant::gate::T;
    return std::all_of(atoms_.begin(), atoms_.end(), [](const Atom a) { return a.IsClifford(); });
}
MCD2 Gate::Mat() const {
    auto ret = MCD2::Identity();
    for (const auto a : atoms_) { ret *= a.Mat(); }
    return ret;
}
void Gate::Normalize() {
    // TODO: Implement https://arxiv.org/abs/0806.3834 !!!!!
    using namespace constant::gate;
    auto ret = std::vector<Atom>();
    ret.reserve(atoms_.size());
    for (const auto a : atoms_) {
        if (ret.empty()) {
            ret.emplace_back(a);
            continue;
        }

        if (ret.back() == H && a == H) {
            ret.pop_back();
        } else if (ret.back() == S && a == S) {
            ret.back() = Z;
        } else if (ret.back() == T && a == T) {
            ret.back() = S;
        } else if (ret.back() == X && a == X) {
            ret.pop_back();
        } else if (ret.back() == Y && a == Y) {
            ret.pop_back();
        } else if (ret.back() == Z && a == Z) {
            ret.pop_back();
        } else {
            ret.emplace_back(a);
        }
    }
    atoms_ = ret;
}
Gate& Gate::operator*=(Atom a) {
    atoms_.emplace_back(a);
    return *this;
}
Gate& Gate::operator*=(const Gate& g) {
    const auto tmp_size = atoms_.size();
    atoms_.resize(tmp_size + g.Size());
    std::copy(g.begin(), g.end(), atoms_.begin() + tmp_size);
    return *this;
}
std::ostream& operator<<(std::ostream& out, const Gate& g) { return out << g.ToString(); }
Gate operator*(Atom l, Atom r) {
    auto ret = Gate(l);
    return ret *= r;
}
Gate operator*(const Gate& l, Atom r) {
    // TODO: Implement more efficiently
    auto ret = Gate(l);
    return ret *= r;
}
Gate operator*(Atom l, const Gate& r) {
    // TODO: Implement more efficiently
    auto ret = Gate(l);
    return ret *= r;
}
Gate operator*(const Gate& l, const Gate& r) {
    // TODO: Implement more efficiently
    auto ret = Gate(l);
    return ret *= r;
}
}  // namespace qrot
