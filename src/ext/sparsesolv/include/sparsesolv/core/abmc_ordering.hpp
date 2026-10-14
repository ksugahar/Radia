/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file abmc_ordering.hpp
 * @brief Algebraic Block Multi-Color (ABMC) ordering for parallel triangular solves
 *
 * Implements the ABMC ordering algorithm from:
 *   T. Iwashita, H. Nakashima, Y. Takahashi,
 *   "Algebraic Block Multi-Color Ordering Method for Parallel Multi-Threaded
 *    Sparse Triangular Solver in ICCG Method", IPDPS 2012
 *
 * Ported from JP-MARs/SparseSolv (MPL 2.0, original author: Shingo Hiruma).
 *
 * The algorithm works in two stages:
 * 1. Algebraic blocking: BFS aggregation groups nearby rows into blocks
 * 2. Multi-coloring: Greedy graph coloring of the block adjacency graph
 *
 * The resulting schedule enables parallel triangular solves:
 *   - Colors are processed sequentially (inter-color dependencies)
 *   - Blocks within the same color are processed in parallel (independent)
 *   - Rows within a block are processed sequentially
 */

#ifndef SPARSESOLV_CORE_ABMC_ORDERING_HPP
#define SPARSESOLV_CORE_ABMC_ORDERING_HPP

#include "types.hpp"
#include <vector>
#include <queue>
#include <algorithm>
#include <map>
#include <cassert>


namespace sparsesolv {

/**
 * @brief ABMC schedule for parallel triangular solves
 *
 * Stores the two-level hierarchy (color -> block -> rows) and the
 * permutation arrays for reordering matrices and vectors.
 */
struct ABMCSchedule {
    /// Flat CSR-like storage for color -> block mapping
    std::vector<index_t> color_offsets;   // size: num_colors + 1
    std::vector<index_t> color_blocks;    // flat list of block indices per color

    /// CSR-like storage for block -> row mapping
    /// Rows within block b are block_offsets[b] .. block_offsets[b+1]-1
    /// (contiguous new indices, no separate row array needed)
    std::vector<index_t> block_offsets;   // size: num_blocks + 1

    /// ordering[old_row] = new_row
    std::vector<index_t> ordering;

    /// reverse_ordering[new_row] = old_row
    std::vector<index_t> reverse_ordering;

    /// Check if the schedule has been built
    bool is_built() const { return !color_offsets.empty(); }

    /// Number of colors
    index_t num_colors() const {
        return color_offsets.empty() ? 0 : static_cast<index_t>(color_offsets.size()) - 1;
    }

    /// Number of blocks
    index_t num_blocks() const {
        return block_offsets.empty() ? 0 : static_cast<index_t>(block_offsets.size()) - 1;
    }

    /**
     * @brief Build ABMC schedule from sparse matrix pattern
     *
     * @param row_ptr CSR row pointer array (size n+1)
     * @param col_idx CSR column index array
     * @param n Number of rows
     * @param block_size Target number of rows per block
     * @param target_colors Target number of colors for graph coloring
     */
    void build(const index_t* row_ptr, const index_t* col_idx, index_t n,
               int block_size, int target_colors) {
        if (n <= 0) {
            clear();
            return;
        }

        // Clamp parameters
        if (block_size < 1) block_size = 1;
        if (target_colors < 2) target_colors = 2;

        // Stage 1: BFS aggregation into blocks
        const index_t max_blocks = (n + block_size - 1) / block_size;
        std::vector<std::vector<index_t>> raw_block_list;
        std::vector<index_t> block_assign(n, -1);

        make_algebraic_blocks(row_ptr, col_idx, n, block_size,
                              max_blocks, raw_block_list, block_assign);

        const index_t actual_num_blocks = static_cast<index_t>(raw_block_list.size());

        // Stage 2: Build block-level adjacency graph (CSR)
        std::vector<index_t> blk_row_ptr;
        std::vector<index_t> blk_col_idx;
        build_block_graph(row_ptr, col_idx, actual_num_blocks,
                          raw_block_list, block_assign,
                          blk_row_ptr, blk_col_idx);

        // Stage 3: Multi-color the block graph
        std::vector<std::vector<index_t>> blk_color_list;

        color_graph(blk_row_ptr.data(), blk_col_idx.data(), actual_num_blocks,
                    target_colors, blk_color_list);

        // Stage 4: Build final row ordering from block ordering
        build_row_ordering(n, actual_num_blocks, raw_block_list,
                           blk_color_list);

    }

    /// Clear all data
    void clear() {
        color_offsets.clear();
        color_blocks.clear();
        block_offsets.clear();
        ordering.clear();
        reverse_ordering.clear();
    }

private:
    /**
     * @brief BFS-based algebraic block aggregation
     *
     * Groups nearby rows into blocks of approximately block_size rows each.
     * Uses breadth-first search on the adjacency graph of the matrix.
     *
     * @param row_ptr CSR row pointer
     * @param col_idx CSR column indices
     * @param n Matrix dimension
     * @param block_size Target rows per block
     * @param max_blocks Maximum number of blocks
     * @param raw_block_list Output: raw_block_list[b] = {original row indices}
     * @param block_assign Output: block_assign[row] = block id
     */
    static void make_algebraic_blocks(
        const index_t* row_ptr, const index_t* col_idx, index_t n,
        int block_size, index_t max_blocks,
        std::vector<std::vector<index_t>>& raw_block_list,
        std::vector<index_t>& block_assign)
    {
        // block_assign values: -1 = unassigned, -2 = in BFS queue, >= 0 = assigned
        raw_block_list.clear();
        raw_block_list.reserve(max_blocks);

        std::queue<index_t> bfs_queue;
        index_t next_unassigned = 0;

        for (index_t blk = 0; blk < max_blocks; ++blk) {
            // Find the first unassigned row
            while (next_unassigned < n && block_assign[next_unassigned] >= 0) {
                ++next_unassigned;
            }
            if (next_unassigned >= n) break;

            // Start a new block with this seed row
            raw_block_list.emplace_back();
            auto& current_block = raw_block_list.back();
            const index_t blk_id = static_cast<index_t>(raw_block_list.size() - 1);

            block_assign[next_unassigned] = blk_id;
            current_block.push_back(next_unassigned);

            // Enqueue neighbors of the seed (mark as -2 to avoid duplicates)
            for (index_t k = row_ptr[next_unassigned]; k < row_ptr[next_unassigned + 1]; ++k) {
                index_t j = col_idx[k];
                if (block_assign[j] == -1) {
                    bfs_queue.push(j);
                    block_assign[j] = -2;
                }
            }

            // BFS expansion until block is full
            while (!bfs_queue.empty() &&
                   static_cast<int>(current_block.size()) < block_size) {
                index_t row = bfs_queue.front();
                bfs_queue.pop();

                if (block_assign[row] >= 0) continue;

                block_assign[row] = blk_id;
                current_block.push_back(row);

                if (static_cast<int>(current_block.size()) >= block_size) break;

                // Enqueue unassigned neighbors (only those not already queued)
                for (index_t k = row_ptr[row]; k < row_ptr[row + 1]; ++k) {
                    index_t j = col_idx[k];
                    if (block_assign[j] == -1) {
                        bfs_queue.push(j);
                        block_assign[j] = -2;
                    }
                }
            }

            // Reset queued-but-unassigned nodes back to -1 for next block
            while (!bfs_queue.empty()) {
                index_t row = bfs_queue.front();
                bfs_queue.pop();
                if (block_assign[row] == -2) {
                    block_assign[row] = -1;
                }
            }
        }

        // Handle any remaining unassigned rows (shouldn't happen normally)
        for (index_t i = 0; i < n; ++i) {
            if (block_assign[i] < 0) {
                if (raw_block_list.empty()) {
                    raw_block_list.emplace_back();
                }
                block_assign[i] = static_cast<index_t>(raw_block_list.size() - 1);
                raw_block_list.back().push_back(i);
            }
        }
    }

    /**
     * @brief Build block-level adjacency graph in CSR format
     *
     * Two blocks are adjacent if any row in one block has a nonzero entry
     * connecting to a row in the other block.
     */
    static void build_block_graph(
        const index_t* row_ptr, const index_t* col_idx,
        index_t num_blocks,
        const std::vector<std::vector<index_t>>& raw_block_list,
        const std::vector<index_t>& block_assign,
        std::vector<index_t>& blk_row_ptr,
        std::vector<index_t>& blk_col_idx)
    {
        // Iterate block-by-block so that last_seen[bj] = bi correctly deduplicates
        std::vector<std::vector<index_t>> neighbors(num_blocks);
        std::vector<index_t> last_seen(num_blocks, static_cast<index_t>(-1));

        for (index_t bi = 0; bi < num_blocks; ++bi) {
            for (index_t row : raw_block_list[bi]) {
                for (index_t k = row_ptr[row]; k < row_ptr[row + 1]; ++k) {
                    index_t bj = block_assign[col_idx[k]];
                    if (bi != bj && last_seen[bj] != bi) {
                        last_seen[bj] = bi;
                        neighbors[bi].push_back(bj);
                    }
                }
            }
        }

        // Convert to CSR
        blk_row_ptr.resize(num_blocks + 1);
        blk_row_ptr[0] = 0;
        for (index_t b = 0; b < num_blocks; ++b) {
            blk_row_ptr[b + 1] = blk_row_ptr[b] + static_cast<index_t>(neighbors[b].size());
        }

        blk_col_idx.resize(blk_row_ptr[num_blocks]);
        for (index_t b = 0; b < num_blocks; ++b) {
            index_t pos = blk_row_ptr[b];
            for (index_t j : neighbors[b]) {
                blk_col_idx[pos++] = j;
            }
        }
    }

    /**
     * @brief Greedy multi-color ordering of a graph using forbidden-color-set
     *
     * Assigns colors to nodes such that no two adjacent nodes in the
     * lower-triangular part share the same color. Uses a sentinel-based
     * forbidden array for O(degree) coloring per node.
     *
     * @param row_ptr CSR row pointer of the graph
     * @param col_idx CSR column indices
     * @param n Number of nodes
     * @param target_colors Suggested number of colors (lower bound)
     * @param out_color_list Output: color_list[c] = {node indices with color c}
     */
    static void color_graph(
        const index_t* row_ptr, const index_t* col_idx, index_t n,
        int target_colors,
        std::vector<std::vector<index_t>>& out_color_list)
    {
        int num_colors = target_colors;

        // Determine minimum number of colors needed
        // (max lower-triangular degree + 1)
        for (index_t i = 0; i < n; ++i) {
            int lower_count = 0;
            for (index_t k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                if (col_idx[k] < i) {
                    ++lower_count;
                }
            }
            if (lower_count >= num_colors) {
                num_colors = lower_count + 1;
            }
        }

        // Greedy coloring with forbidden-color-set approach: O(degree) per node
        std::vector<index_t> node_color(n, -1);
        out_color_list.resize(num_colors);
        for (auto& cl : out_color_list) cl.clear();

        // Sentinel-based forbidden array: forbidden[c] == i means color c
        // is used by a lower-triangular neighbor of node i
        std::vector<index_t> forbidden(num_colors, static_cast<index_t>(-1));

        for (index_t i = 0; i < n; ++i) {
            // Mark forbidden colors in one pass over neighbors
            for (index_t k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                index_t j = col_idx[k];
                if (j < i && node_color[j] >= 0) {
                    forbidden[node_color[j]] = i;
                }
            }

            // Pick first available color
            int chosen = -1;
            for (int c = 0; c < num_colors; ++c) {
                if (forbidden[c] != i) {
                    chosen = c;
                    break;
                }
            }

            // All existing colors conflict: create a new color
            if (chosen < 0) {
                chosen = num_colors;
                ++num_colors;
                forbidden.push_back(static_cast<index_t>(-1));
                out_color_list.emplace_back();
            }

            node_color[i] = chosen;
            out_color_list[chosen].push_back(i);
        }

        // Remove empty color groups
        out_color_list.erase(
            std::remove_if(out_color_list.begin(), out_color_list.end(),
                           [](const std::vector<index_t>& v) { return v.empty(); }),
            out_color_list.end());
    }

    /**
     * @brief Build final row-level ordering from block ordering
     *
     * Enumerates rows in the order: color 0 blocks -> color 1 blocks -> ...
     * Within each block, rows are sorted by original index for locality.
     *
     * The new row indices within each block are contiguous:
     * block b contains rows block_offsets[b] .. block_offsets[b+1]-1.
     * No separate block_rows array is needed (it would be the identity).
     */
    void build_row_ordering(
        index_t n, index_t actual_num_blocks,
        std::vector<std::vector<index_t>>& raw_block_list,
        const std::vector<std::vector<index_t>>& blk_color_list)
    {
        const index_t nc = static_cast<index_t>(blk_color_list.size());
        ordering.resize(n);
        reverse_ordering.resize(n);

        // Build flat color -> block mapping
        color_offsets.resize(nc + 1);
        color_offsets[0] = 0;
        for (index_t c = 0; c < nc; ++c) {
            color_offsets[c + 1] = color_offsets[c] +
                static_cast<index_t>(blk_color_list[c].size());
        }
        color_blocks.resize(color_offsets[nc]);

        // Build flat block -> row mapping
        block_offsets.resize(actual_num_blocks + 1);

        index_t row_idx = 0;
        index_t global_blk_id = 0;
        for (index_t c = 0; c < nc; ++c) {
            for (index_t bi = 0; bi < static_cast<index_t>(blk_color_list[c].size()); ++bi) {
                index_t orig_blk = blk_color_list[c][bi];
                color_blocks[color_offsets[c] + bi] = global_blk_id;

                // Sort rows within block in-place for better locality
                auto& rows = raw_block_list[orig_blk];
                std::sort(rows.begin(), rows.end());

                block_offsets[global_blk_id] = row_idx;
                for (index_t orig_row : rows) {
                    ordering[orig_row] = row_idx;
                    reverse_ordering[row_idx] = orig_row;
                    ++row_idx;
                }
                ++global_blk_id;
            }
        }
        block_offsets[actual_num_blocks] = row_idx;

        assert(row_idx == n);
        assert(global_blk_id == actual_num_blocks);
    }
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_ABMC_ORDERING_HPP
