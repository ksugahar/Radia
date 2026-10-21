#pragma once
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace radia_hmatrix_cuda {
// Schema-2 describes a partition, rather than a sum of overlapping operators.
// Diagonal squares own only their upper triangle; off-diagonal blocks are mirrored.
template<class Leaves, class Allocator = std::allocator<std::uint64_t>>
void requireUpperLeafPartition(std::size_t n, const Leaves& leaves,
                               const Allocator& allocator = Allocator{}) {
    if (!n || n > 4608) throw std::invalid_argument("Leaf coverage dimension budget exceeded");
    const auto words = (n + 63) / 64;
    std::vector<std::uint64_t, Allocator> bits(n * words, 0, allocator);
    std::size_t covered = 0;
    for (const auto& leaf : leaves) {
        if (!leaf.rows || !leaf.columns || leaf.row >= n || leaf.column >= n ||
            leaf.rows > n - leaf.row || leaf.columns > n - leaf.column ||
            leaf.row > leaf.column ||
            (leaf.row == leaf.column && leaf.rows != leaf.columns) ||
            (leaf.row < leaf.column && leaf.rows > leaf.column - leaf.row))
            throw std::invalid_argument("Invalid upper leaf coverage range");
        for (auto row = leaf.row; row < leaf.row + leaf.rows; ++row) {
            for (auto col = std::max(row, leaf.column); col < leaf.column + leaf.columns; ++col) {
                auto& word = bits[row * words + col / 64];
                const auto mask = std::uint64_t(1) << (col % 64);
                if (word & mask) throw std::invalid_argument("Duplicate or overlapping upper leaf coverage");
                word |= mask;
                ++covered;
            }
        }
    }
    if (covered != n * (n + 1) / 2)
        throw std::invalid_argument("Incomplete upper leaf coverage");
}
}
