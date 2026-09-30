#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>

namespace radia::matlab {

struct SparseSymmetryCheck {
    bool symmetric = true;
    int row = -1;
    int column = -1;
    double difference = 0.0;
    double tolerance = 0.0;
};

// Check A == transpose(A), including complex-symmetric matrices.  Sparse
// Cholesky uses transpose symmetry; conjugating the counterpart here would
// incorrectly reject the complex-symmetric eddy-current systems it supports.
template <class Scalar, class RowSize, class Column, class Value>
SparseSymmetryCheck CheckSparseTransposeSymmetry(
    int size, RowSize row_size, Column column, Value value,
    double relative_tolerance = 64.0 * std::numeric_limits<double>::epsilon()) {
    SparseSymmetryCheck result;
    double scale = 0.0;
    for (int row = 0; row < size; ++row)
        for (int entry = 0; entry < row_size(row); ++entry)
            scale = std::max(scale, static_cast<double>(std::abs(value(row, entry))));
    result.tolerance = relative_tolerance * scale;

    for (int row = 0; row < size; ++row) {
        for (int entry = 0; entry < row_size(row); ++entry) {
            const int col = column(row, entry);
            if (col < 0 || col >= size) {
                result.symmetric = false;
                result.row = row;
                result.column = col;
                result.difference = std::numeric_limits<double>::infinity();
                return result;
            }
            Scalar counterpart{};
            for (int transpose_entry = 0;
                 transpose_entry < row_size(col); ++transpose_entry) {
                if (column(col, transpose_entry) == row) {
                    counterpart = value(col, transpose_entry);
                    break;
                }
            }
            const double difference = static_cast<double>(
                std::abs(value(row, entry) - counterpart));
            if (!std::isfinite(difference) || difference > result.tolerance) {
                result.symmetric = false;
                result.row = row;
                result.column = col;
                result.difference = difference;
                return result;
            }
        }
    }
    return result;
}

}  // namespace radia::matlab
