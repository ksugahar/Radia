#include "radia_mex_sparse_symmetry.h"

#include <complex>
#include <iostream>
#include <utility>
#include <vector>

template <class Scalar>
bool IsSymmetric(const std::vector<std::vector<std::pair<int, Scalar>>>& rows) {
    const auto result = radia::matlab::CheckSparseTransposeSymmetry<Scalar>(
        static_cast<int>(rows.size()),
        [&](int row) { return static_cast<int>(rows[row].size()); },
        [&](int row, int entry) { return rows[row][entry].first; },
        [&](int row, int entry) { return rows[row][entry].second; });
    return result.symmetric;
}

int main() {
    using Complex = std::complex<double>;
    const std::vector<std::vector<std::pair<int, double>>> real_symmetric{
        {{0, 4.0}, {1, -2.0}},
        {{0, -2.0}, {1, 3.0}},
    };
    const std::vector<std::vector<std::pair<int, Complex>>> complex_symmetric{
        {{0, {4.0, 1.0}}, {1, {-2.0, 0.5}}},
        {{0, {-2.0, 0.5}}, {1, {3.0, -1.0}}},
    };
    const std::vector<std::vector<std::pair<int, double>>> real_nonsymmetric{
        {{0, 4.0}, {1, -2.0}},
        {{0, -1.0}, {1, 3.0}},
    };
    const std::vector<std::vector<std::pair<int, Complex>>> complex_hermitian_not_symmetric{
        {{0, {4.0, 0.0}}, {1, {-2.0, 0.5}}},
        {{0, {-2.0, -0.5}}, {1, {3.0, 0.0}}},
    };

    if (!IsSymmetric(real_symmetric)) {
        std::cerr << "real symmetric matrix was rejected\n";
        return 1;
    }
    if (!IsSymmetric(complex_symmetric)) {
        std::cerr << "complex symmetric matrix was rejected\n";
        return 1;
    }
    if (IsSymmetric(real_nonsymmetric)) {
        std::cerr << "real nonsymmetric matrix was accepted\n";
        return 1;
    }
    if (IsSymmetric(complex_hermitian_not_symmetric)) {
        std::cerr << "Hermitian but non-complex-symmetric matrix was accepted\n";
        return 1;
    }
    return 0;
}
