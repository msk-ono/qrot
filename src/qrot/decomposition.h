#ifndef QROT_DECOMPOSITION_H
#define QROT_DECOMPOSITION_H

#include <string>
#include <utility>
#include <vector>

#include "qrot/gate.h"
#include "qrot/matrix.h"
#include "qrot/number.h"

namespace qrot {
/**
 * @brief Decompose unitary matrix to quantum gates.
 * @details Implementation of Algorithm 1 in 1206.5236.
 */
class UnitaryDecomposer {
public:
    UnitaryDecomposer();

    Gate Decompose(const MCD2& input);

private:
    void InitializeStorage();
    void InitializeStorageImpl();
    bool LoadS3(const std::string& path);
    bool LoadS3Impl(std::istream& is);
    bool StoreS3(const std::string& path);
    Gate LookUpS3(const MCD2& mat) const;

    std::vector<std::pair<MCD2, Gate>> s3_;
};
}  // namespace qrot

#endif  // QROT_DECOMPOSITION_H
