/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file rcm_ordering.hpp
 * @brief Reverse Cuthill-McKee (RCM) bandwidth reduction ordering
 *
 * RCM reorders rows/columns to reduce matrix bandwidth, improving
 * cache locality for SpMV and triangular solves.
 *
 * Algorithm:
 * 1. Find a pseudo-peripheral node (BFS, pick farthest)
 * 2. BFS from that node, ordering neighbors by ascending degree
 * 3. Reverse the ordering
 */

#ifndef SPARSESOLV_CORE_RCM_ORDERING_HPP
#define SPARSESOLV_CORE_RCM_ORDERING_HPP

#include "types.hpp"
#include <vector>
#include <queue>
#include <algorithm>
#include <numeric>

namespace sparsesolv {

/**
 * @brief Compute RCM (Reverse Cuthill-McKee) ordering
 *
 * @param row_ptr CSR row pointers (size n+1)
 * @param col_idx CSR column indices
 * @param n Matrix dimension
 * @param ordering Output: ordering[old_row] = new_row
 * @param reverse_ordering Output: reverse_ordering[new_row] = old_row
 */
inline void compute_rcm_ordering(
    const index_t* row_ptr, const index_t* col_idx, index_t n,
    std::vector<index_t>& ordering,
    std::vector<index_t>& reverse_ordering)
{
    ordering.resize(n);
    reverse_ordering.resize(n);

    if (n <= 0) return;

    // Compute degree for each node
    std::vector<index_t> degree(n);
    for (index_t i = 0; i < n; ++i) {
        degree[i] = row_ptr[i + 1] - row_ptr[i];
    }

    // Find pseudo-peripheral node using BFS
    // Start from node with minimum degree
    index_t start = 0;
    index_t min_deg = degree[0];
    for (index_t i = 1; i < n; ++i) {
        if (degree[i] < min_deg) {
            min_deg = degree[i];
            start = i;
        }
    }

    // Two BFS passes to find a better pseudo-peripheral node
    std::vector<index_t> dist(n, -1);
    for (int pass = 0; pass < 2; ++pass) {
        std::fill(dist.begin(), dist.end(), -1);
        std::queue<index_t> q;
        q.push(start);
        dist[start] = 0;
        index_t farthest = start;
        index_t max_dist = 0;

        while (!q.empty()) {
            index_t u = q.front();
            q.pop();
            for (index_t k = row_ptr[u]; k < row_ptr[u + 1]; ++k) {
                index_t v = col_idx[k];
                if (dist[v] < 0) {
                    dist[v] = dist[u] + 1;
                    q.push(v);
                    if (dist[v] > max_dist) {
                        max_dist = dist[v];
                        farthest = v;
                    }
                }
            }
        }
        start = farthest;
    }

    // BFS from pseudo-peripheral node, ordering neighbors by degree
    std::vector<bool> visited(n, false);
    std::vector<index_t> cm_order; // Cuthill-McKee order (will be reversed)
    cm_order.reserve(n);

    // Handle potentially disconnected graphs
    for (index_t seed = 0; seed < n; ++seed) {
        index_t actual_seed = (seed == 0) ? start : seed;
        if (visited[actual_seed]) continue;

        std::queue<index_t> bfs;
        bfs.push(actual_seed);
        visited[actual_seed] = true;
        cm_order.push_back(actual_seed);

        while (!bfs.empty()) {
            index_t u = bfs.front();
            bfs.pop();

            // Collect unvisited neighbors
            std::vector<index_t> nbrs;
            for (index_t k = row_ptr[u]; k < row_ptr[u + 1]; ++k) {
                index_t v = col_idx[k];
                if (!visited[v]) {
                    nbrs.push_back(v);
                    visited[v] = true;
                }
            }

            // Sort neighbors by ascending degree (Cuthill-McKee heuristic)
            std::sort(nbrs.begin(), nbrs.end(), [&](index_t a, index_t b) {
                return degree[a] < degree[b];
            });

            for (index_t v : nbrs) {
                cm_order.push_back(v);
                bfs.push(v);
            }
        }
    }

    // Reverse to get RCM ordering
    for (index_t i = 0; i < n; ++i) {
        index_t old_row = cm_order[n - 1 - i]; // reverse
        reverse_ordering[i] = old_row;
        ordering[old_row] = i;
    }
}

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_RCM_ORDERING_HPP
