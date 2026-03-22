/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file sparse_matrix_coo.hpp
/// @brief COO accumulator for element-by-element assembly, converts to CSR

#ifndef SPARSESOLV_CORE_SPARSE_MATRIX_COO_HPP
#define SPARSESOLV_CORE_SPARSE_MATRIX_COO_HPP

#include "types.hpp"
#include "sparse_matrix_csr.hpp"
#include "dense_matrix.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

namespace sparsesolv {

/// COO sparse matrix accumulator. Duplicate entries are summed on to_csr().
template<typename Scalar>
class SparseMatrixCOO {
public:
    SparseMatrixCOO(index_t rows, index_t cols)
        : rows_(rows), cols_(cols) {}

    /// Reserve capacity for entries
    void reserve(size_t n) { entries_.reserve(n); }

    /// Add a single entry
    void add(index_t row, index_t col, Scalar value) {
        if (row >= rows_ || col >= cols_)
            throw std::out_of_range("SparseMatrixCOO::add: index (" + std::to_string(row)
                + ", " + std::to_string(col) + ") out of range (" + std::to_string(rows_)
                + " x " + std::to_string(cols_) + ")");
        entries_.push_back({row, col, value});
    }

    /// Add dense submatrix at specified global DOF positions
    void add_submatrix(const index_t* row_dofs, index_t nrow,
                       const index_t* col_dofs, index_t ncol,
                       const DenseMatrix<Scalar>& mat) {
        for (index_t i = 0; i < nrow; ++i) {
            for (index_t j = 0; j < ncol; ++j) {
                Scalar val = mat(i, j);
                if (detail::dense_abs(val) > 0) {
                    entries_.push_back({row_dofs[i], col_dofs[j], val});
                }
            }
        }
    }

    /// Merge entries from another COO accumulator (move)
    void merge_from(SparseMatrixCOO&& other) {
        entries_.insert(entries_.end(),
                        std::make_move_iterator(other.entries_.begin()),
                        std::make_move_iterator(other.entries_.end()));
        other.entries_.clear();
    }

    /// Convert to CSR: sort by (row,col), merge duplicates, build row_ptr
    SparseMatrixCSR<Scalar> to_csr() const {
        SparseMatrixCSR<Scalar> csr;
        csr.rows = rows_;
        csr.cols = cols_;

        if (entries_.empty()) {
            csr.row_ptr.assign(rows_ + 1, 0);
            return csr;
        }

        // Sort entries by (row, col)
        auto sorted = entries_;
        std::sort(sorted.begin(), sorted.end(),
            [](const Entry& a, const Entry& b) {
                return (a.row < b.row) || (a.row == b.row && a.col < b.col);
            });

        // Merge duplicates and build CSR
        std::vector<index_t> merged_row, merged_col;
        std::vector<Scalar> merged_val;
        merged_row.reserve(sorted.size());
        merged_col.reserve(sorted.size());
        merged_val.reserve(sorted.size());

        merged_row.push_back(sorted[0].row);
        merged_col.push_back(sorted[0].col);
        merged_val.push_back(sorted[0].value);

        for (size_t k = 1; k < sorted.size(); ++k) {
            if (sorted[k].row == merged_row.back() &&
                sorted[k].col == merged_col.back()) {
                merged_val.back() += sorted[k].value;
            } else {
                merged_row.push_back(sorted[k].row);
                merged_col.push_back(sorted[k].col);
                merged_val.push_back(sorted[k].value);
            }
        }

        // Build row_ptr
        csr.row_ptr.assign(rows_ + 1, 0);
        for (size_t k = 0; k < merged_row.size(); ++k) {
            csr.row_ptr[merged_row[k] + 1]++;
        }
        for (index_t i = 1; i <= rows_; ++i) {
            csr.row_ptr[i] += csr.row_ptr[i - 1];
        }

        csr.col_idx = std::move(merged_col);
        csr.values = std::move(merged_val);

        return csr;
    }

private:
    struct Entry {
        index_t row, col;
        Scalar value;
    };

    std::vector<Entry> entries_;
    index_t rows_, cols_;
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_SPARSE_MATRIX_COO_HPP
