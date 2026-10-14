/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file sparse_matrix_csr.hpp
/// @brief Owned CSR matrix storage (IC L factor, etc.)

#ifndef SPARSESOLV_CORE_SPARSE_MATRIX_CSR_HPP
#define SPARSESOLV_CORE_SPARSE_MATRIX_CSR_HPP

#include "types.hpp"
#include "parallel.hpp"
#include <vector>

namespace sparsesolv {

/// Sparse matrix in CSR format with owned storage
template<typename Scalar>
struct SparseMatrixCSR {
    std::vector<index_t> row_ptr;
    std::vector<index_t> col_idx;
    std::vector<Scalar> values;
    index_t rows = 0;
    index_t cols = 0;

    index_t nnz() const { return static_cast<index_t>(values.size()); }

    void clear() {
        row_ptr.clear();
        col_idx.clear();
        values.clear();
        rows = cols = 0;
    }

    /// y = A * x (parallel)
    void multiply(const Scalar* x, Scalar* y) const {
        parallel_for(rows, [&](index_t i) {
            Scalar sum = Scalar(0);
            for (index_t k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                sum += values[k] * x[col_idx[k]];
            }
            y[i] = sum;
        });
    }

    /// Compute transpose: B = A^T
    SparseMatrixCSR transpose() const {
        SparseMatrixCSR result;
        result.rows = cols;
        result.cols = rows;
        const index_t nnz_total = static_cast<index_t>(values.size());
        result.row_ptr.assign(cols + 1, 0);
        result.col_idx.resize(nnz_total);
        result.values.resize(nnz_total);

        for (index_t k = 0; k < nnz_total; ++k)
            result.row_ptr[col_idx[k] + 1]++;
        for (index_t i = 0; i < cols; ++i)
            result.row_ptr[i + 1] += result.row_ptr[i];

        std::vector<index_t> counter(cols, 0);
        for (index_t i = 0; i < rows; ++i) {
            for (index_t k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                index_t j = col_idx[k];
                index_t pos = result.row_ptr[j] + counter[j];
                result.col_idx[pos] = i;
                result.values[pos] = values[k];
                counter[j]++;
            }
        }
        return result;
    }
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_SPARSE_MATRIX_CSR_HPP
