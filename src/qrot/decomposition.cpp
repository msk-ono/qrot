#include "qrot/decomposition.h"

#include <execution>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>

namespace qrot {
namespace {
static inline const auto StepGate = constant::mcd2::H * constant::mcd2::TDag * constant::mcd2::H;
/**
 * @brief Calculate the smallest denominator exponent.
 */
std::int32_t SDE(const D2& x) {
    auto ret = std::int32_t{0};
    ret = std::max(ret, 2 * x.Int().DenExp());
    ret = std::max(ret, 2 * x.Sqrt().DenExp() - 1);
    return ret;
}
}  // namespace
UnitaryDecomposer::UnitaryDecomposer() {
    auto ss = std::stringstream(
#include "qrot/s3.txt"
    );
    LoadS3Impl(ss);
}
void UnitaryDecomposer::InitializeStorage() {
    const auto path = std::string("s3.txt");
    if (LoadS3(path)) { return; }
    InitializeStorageImpl();
    if (!StoreS3(path)) {
        std::cerr << "Failed to store s3" << std::endl;
        return;
    }
}
void UnitaryDecomposer::InitializeStorageImpl() {
    using namespace constant;
    using constant::mcd2::I, constant::mcd2::H, constant::mcd2::T, constant::mcd2::TDag;
    constexpr auto MaxSDE = 4;
    constexpr auto MaxDepth = 30;

    // Queues
    auto i_queue = std::queue<std::pair<MCD2, Gate>>();
    auto o_queue = std::queue<std::pair<MCD2, Gate>>();

    // Initialize
    i_queue.push({I, Gate()});
    s3_.emplace_back(I, Gate());

    // Cache
    auto cache = std::vector<MCD2>();
    const auto searched = [&cache](const MCD2& mat) -> bool {
        // return std::find(std::execution::par, cache.begin(), cache.end(), mat) != cache.end();
        return std::find(cache.begin(), cache.end(), mat) != cache.end();
    };

    // Search depth
    auto depth = 0;
    while (true) {
        const auto size_before = s3_.size();
        while (!i_queue.empty()) {
            const auto [top, gate] = i_queue.front();
            i_queue.pop();

            // Test: H * U
            const auto hx = H * top;
            if (!searched(hx)) {
                const auto s = SDE(hx.Get(0, 0).Norm());
                if (s <= MaxSDE) {
                    cache.emplace_back(hx);
                    o_queue.push({hx, gate::H * gate});
                }
                if (s <= 3) { s3_.emplace_back(hx, gate::H * gate); }
            }
            // Test: T * U
            const auto tx = T * top;
            if (!searched(tx)) {
                const auto s = SDE(tx.Get(0, 0).Norm());
                if (s <= MaxSDE) {
                    cache.emplace_back(tx);
                    o_queue.push({tx, gate::T * gate});
                }
                if (s <= 3) { s3_.emplace_back(tx, gate::T * gate); }
            }
        }
        const auto size_after = s3_.size();
        std::cout << depth << " " << size_before << " " << size_after << " " << o_queue.size()
                  << " " << cache.size() << std::endl;

        // Next depth
        std::swap(i_queue, o_queue);
        if (depth++ >= MaxDepth) { break; }
        if (i_queue.empty()) { break; }
    }
}
bool UnitaryDecomposer::LoadS3(const std::string& path) {
    auto ifs = std::ifstream(path);
    if (!ifs.good()) { return false; }
    return LoadS3Impl(ifs);
}
bool UnitaryDecomposer::LoadS3Impl(std::istream& is) {
    // TODO: Implement in more efficient way
    std::int32_t num;
    is >> num;
    for (auto i = 0; i < num; ++i) {
        Integer rin, rsn, iin, isn;
        std::int32_t rid, rsd, iid, isd;
        std::string s;
        const auto to_cd2 = [&rin, &rsn, &iin, &isn, &rid, &rsd, &iid, &isd]() {
            return CD2(D2(DyadicFraction(rin, rid), DyadicFraction(rsn, rsd)),
                       D2(DyadicFraction(iin, iid), DyadicFraction(isn, isd)));
        };
        is >> rin >> rid >> rsn >> rsd >> iin >> iid >> isn >> isd;
        const auto x00 = to_cd2();
        is >> rin >> rid >> rsn >> rsd >> iin >> iid >> isn >> isd;
        const auto x01 = to_cd2();
        is >> rin >> rid >> rsn >> rsd >> iin >> iid >> isn >> isd;
        const auto x10 = to_cd2();
        is >> rin >> rid >> rsn >> rsd >> iin >> iid >> isn >> isd;
        const auto x11 = to_cd2();
        is >> s;
        s3_.emplace_back(MCD2(x00, x01, x10, x11), Gate::FromString(s));
    }
    return true;
}
bool UnitaryDecomposer::StoreS3(const std::string& path) {
    auto ofs = std::ofstream(path);
    if (!ofs.good()) { return false; }
    ofs << s3_.size() << '\n';
    for (const auto& [mat, gate] : s3_) {
        for (auto i = 0; i < 2; ++i) {
            for (auto j = 0; j < 2; ++j) {
                const auto& x = mat.Get(i, j);
                ofs << x.Real().Int().Num() << ' ';
                ofs << x.Real().Int().DenExp() << ' ';
                ofs << x.Real().Sqrt().Num() << ' ';
                ofs << x.Real().Sqrt().DenExp() << ' ';
                ofs << x.Imag().Int().Num() << ' ';
                ofs << x.Imag().Int().DenExp() << ' ';
                ofs << x.Imag().Sqrt().Num() << ' ';
                ofs << x.Imag().Sqrt().DenExp() << ' ';
            }
        }
        ofs << gate << '\n';
    }
    return true;
}
Gate UnitaryDecomposer::LookUpS3(const MCD2& x) const {
    for (const auto& [mat, s] : s3_) {
        if (x == mat) { return s; }
    }
    throw std::logic_error("Cannot find unitary for input matrix in s3 database");
    return Gate();
}
Gate UnitaryDecomposer::Decompose(const MCD2& input) {
    using namespace constant;
    using constant::mcd2::H;

    auto unitary = input;
    auto s = SDE(unitary.Get(0, 0).Norm());
    auto output = Gate();

    while (s > 3) {
        auto tmp = H * unitary;
        auto found_state = false;
        for (auto i = 0; i < 4; ++i) {
            const auto tmp_s = SDE(tmp.Get(0, 0).Norm());
            if (tmp_s == s - 1) {
                s = tmp_s;
                for (auto j = 0; j < i; ++j) { output *= gate::T; }
                output *= gate::H;
                unitary = tmp;
                found_state = true;
                break;
            }
            tmp.MulFromLeft(StepGate);
        }
        if (!found_state) {
            assert(0 && "Unreachable: not found state for unitary");
            break;
        }
    }

    // Look up
    output *= LookUpS3(unitary);

    output.Normalize();
    return output;
}
}  // namespace qrot
