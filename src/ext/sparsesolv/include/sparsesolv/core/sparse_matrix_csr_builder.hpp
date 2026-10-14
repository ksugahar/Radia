/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file sparse_matrix_csr_builder.hpp
/// @brief Row-based CSR builder: O(sum d_i log d_i) per-row sort instead of O(N log N) global sort

#ifndef SPARSESOLV_CORE_SPARSE_MATRIX_CSR_BUILDER_HPP
#define SPARSESOLV_CORE_SPARSE_MATRIX_CSR_BUILDER_HPP

#include "types.hpp"
#include "sparse_matrix_csr.hpp"
#include "dense_matrix.hpp"
#include "parallel.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

namespace sparsesolv {

/// Row-based sparse matrix builder. Drop-in replacement for SparseMatrixCOO.
/// Entries are stored per-row, enabling per-row sort + dedup in build().
template<typename Scalar>
class SparseMatrixCSRBuilder {
public:
    SparseMatrixCSRBuilder(index_t rows, index_t cols)
        : rows_(rows), cols_(cols), row_entries_(rows) {}

    /// Add a single entry
    void add(index_t row, index_t col, Scalar value) {
        row_entries_[row].push_back({col, value});
    }

    /// Add dense submatrix at specified global DOF positions
    void add_submatrix(const index_t* row_dofs, index_t nrow,
                       const index_t* col_dofs, index_t ncol,
                       const DenseMatrix<Scalar>& mat) {
        for (index_t i = 0; i < nrow; ++i) {
            for (index_t j = 0; j < ncol; ++j) {
                Scalar val = mat(i, j);
                if (detail::dense_abs(val) > 0) {
                    row_entries_[row_dofs[i]].push_back({col_dofs[j], val});
                }
            }
        }
    }

    /// Merge entries from another builder (move, per-row append)
    void merge_from(SparseMatrixCSRBuilder&& other) {
        for (index_t i = 0; i < rows_; ++i) {
            auto& mine = row_entries_[i];
            auto& theirs = other.row_entries_[i];
            if (theirs.empty()) continue;
            mine.insert(mine.end(),
                        std::make_move_iterator(theirs.begin()),
                        std::make_move_iterator(theirs.end()));
            theirs.clear();
        }
    }

    /// Build CSR: per-row sort + dedup + flatten. O(sum d_i log d_i).
    SparseMatrixCSR<Scalar> build() {
        SparseMatrixCSR<Scalar> csr;
        csr.rows = rows_;
        csr.cols = cols_;
        csr.row_ptr.resize(rows_ + 1, 0);

        // Phase 1: Per-row sort and dedup (parallel)
        parallel_for(rows_, [&](index_t i) {
            auto& entries = row_entries_[i];
            if (entries.size() <= 1) return;

            std::sort(entries.begin(), entries.end(),
                [](const Entry& a, const Entry& b) { return a.col < b.col; });

            // In-place merge of duplicates
            size_t w = 0;
            for (size_t k = 1; k < entries.size(); ++k) {
                if (entries[k].col == entries[w].col)
                    entries[w].val += entries[k].val;
                else
                    entries[++w] = entries[k];
            }
            entries.resize(w + 1);
        });

        // Phase 2: Build row_ptr (prefix sum)
        for (index_t i = 0; i < rows_; ++i)
            csr.row_ptr[i + 1] = csr.row_ptr[i]
                + static_cast<index_t>(row_entries_[i].size());

        // Phase 3: Flatten to col_idx and values
        index_t nnz = csr.row_ptr[rows_];
        csr.col_idx.resize(nnz);
        csr.values.resize(nnz);

        parallel_for(rows_, [&](index_t i) {
            index_t base = csr.row_ptr[i];
            const auto& entries = row_entries_[i];
            for (size_t k = 0; k < entries.size(); ++k) {
                csr.col_idx[base + k] = entries[k].col;
                csr.values[base + k] = entries[k].val;
            }
        });

        return csr;
    }

private:
    struct Entry {
        index_t col;
        Scalar val;
    };

    index_t rows_, cols_;
    std::vector<std::vector<Entry>> row_entries_;
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_SPARSE_MATRIX_CSR_BUILDER_HPP
