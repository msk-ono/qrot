#include "qrot/gate.h"

#include <algorithm>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace qrot {
#pragma region Atom
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
#pragma endregion Atom
#pragma region Gate
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
    using namespace constant;
    static CliffordDatabase database;

    auto normal = std::vector<Gate>();
    auto clifford = mcd2::I;
    for (const auto a : atoms_) {
        if (a == gate::T) {
            const auto index = database.SearchIndex(clifford);
            const auto move = database.GetTMove(index);
            switch (CliffordDatabase::GetType(index)) {
                case CliffordDatabase::Type::CT: {
                    if (!normal.empty() && normal.back() == gate::T) {
                        normal.pop_back();
                        if (!normal.empty()) {
                            clifford = normal.back().Mat() * mcd2::S * database.GetMatrix(move);
                            normal.pop_back();
                        } else {
                            clifford = mcd2::S * database.GetMatrix(move);
                        }
                    } else {
                        normal.emplace_back(gate::T);
                        clifford = database.GetMatrix(move);
                    }
                    break;
                }
                case CliffordDatabase::Type::HCT: {
                    normal.emplace_back(gate::H);
                    normal.emplace_back(gate::T);
                    clifford = database.GetMatrix(move);
                    break;
                }
                case CliffordDatabase::Type::SHCT: {
                    normal.emplace_back(gate::S * gate::H);
                    normal.emplace_back(gate::T);
                    clifford = database.GetMatrix(move);
                    break;
                }
                case CliffordDatabase::Type::NotClifford:
                    throw std::logic_error("clifford is not Clifford");
            }
        } else {
            clifford *= a.Mat();
        }
    }
    if (clifford != mcd2::I) {
        normal.emplace_back(database.GetGate(database.SearchIndex(clifford)));
    }

    atoms_.clear();
    for (const auto& gate : normal) {
        for (const auto atom : gate) { atoms_.emplace_back(atom); }
    }
}
bool operator==(const Gate& lhs, const Gate& rhs) {
    if (lhs.Size() != rhs.Size()) { return false; }
    for (auto i = std::size_t{0}; i < lhs.Size(); ++i) {
        if (lhs.Get(i) != rhs.Get(i)) { return false; }
    }
    return true;
}
bool operator==(const Gate& lhs, const Atom& rhs) { return lhs.Size() == 1 && lhs.Get(0) == rhs; }
bool operator==(const Atom& lhs, const Gate& rhs) { return rhs.Size() == 1 && rhs.Get(0) == lhs; }
bool operator!=(const Gate& lhs, const Gate& rhs) { return !(lhs == rhs); }
bool operator!=(const Gate& lhs, const Atom& rhs) { return !(lhs == rhs); }
bool operator!=(const Atom& lhs, const Gate& rhs) { return !(lhs == rhs); }
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
#pragma endregion Gate
#pragma region CliffordDatabase
CliffordDatabase::CliffordDatabase() {
    using namespace constant;

    // Calculate C_T of 0806.3834
    auto queue = std::queue<std::pair<MCD2, Gate>>();
    auto database = std::vector<std::pair<MCD2, Gate>>();
    const auto pushed = [&database](const MCD2& m) {
        return std::find_if(database.begin(), database.end(),
                            [&m](const auto& entry) { return entry.first == m; }) != database.end();
    };
    // Initialize
    queue.push({mcd2::I, Gate()});
    database.emplace_back(mcd2::I, Gate());
    // BFS
    while (!queue.empty()) {
        const auto [m, g] = queue.front();
        queue.pop();
        {
            const auto n = m * mcd2::S;
            if (!pushed(n)) {
                queue.push({n, g * gate::S});
                database.emplace_back(n, g * gate::S);
            }
        }
        {
            const auto n = m * mcd2::X;
            if (!pushed(n)) {
                queue.push({n, g * gate::X});
                database.emplace_back(n, g * gate::X);
            }
        }
        {
            const auto n = m * mcd2::W;
            if (!pushed(n)) {
                queue.push({n, g * gate::W});
                database.emplace_back(n, g * gate::W);
            }
        }
    }

    // Add prefix H
    for (auto i = std::size_t{0}; i < NumCT; ++i) {
        const auto& [m, g] = database[i];
        assert(!pushed(mcd2::H * m));
        database.emplace_back(mcd2::H * m, gate::H * g);
    }
    // Add prefix SH
    for (auto i = std::size_t{0}; i < NumCT; ++i) {
        const auto& [m, g] = database[i];
        assert(!pushed(mcd2::S * mcd2::H * m));
        database.emplace_back(mcd2::S * mcd2::H * m, gate::S * gate::H * g);
    }
    c1_.swap(database);

    // Calculate move of TDag C_T T
    move_.resize(NumCT);
    for (auto i = std::size_t{0}; i < NumCT; ++i) {
        const auto& m = GetMatrix(i);
        const auto n = mcd2::TDag * m * mcd2::T;
        move_[i] = SearchIndex(n);
    }
}
CliffordDatabase::Type CliffordDatabase::GetType(std::size_t idx) {
    if (idx < NumCT) {
        return Type::CT;
    } else if (idx < 2 * NumCT) {
        return Type::HCT;
    } else if (idx < 3 * NumCT) {
        return Type::SHCT;
    }
    return Type::NotClifford;
}
std::size_t CliffordDatabase::SearchIndex(const MCD2& mat) const {
    for (auto i = std::size_t{0}; i < c1_.size(); ++i) {
        if (mat == c1_[i].first) { return i; }
    }
    return std::numeric_limits<std::size_t>::max();
}
std::size_t CliffordDatabase::GetTMove(std::size_t idx) const { return move_[idx % NumCT]; }
#pragma endregion CliffordDatabase
}  // namespace qrot
