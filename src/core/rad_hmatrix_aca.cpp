/*-------------------------------------------------------------------------*\
*
*   Radia H-Matrix with ACA Compression - Implementation
*
*   Hierarchical matrix implementation for O(N log N) solver acceleration
*   Uses Adaptive Cross Approximation (ACA) for low-rank block compression
*   OpenMP parallelization for multi-threaded performance
*
*   Based on HACApK concepts for hierarchical matrix operations
*
*   Copyright (C): 2024 by European Synchrotron Radiation Facility, France
*
*-------------------------------------------------------------------------*/

#include "rad_hmatrix_aca.h"
#include "rad_interaction.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

//-------------------------------------------------------------------------
// Constructor
//-------------------------------------------------------------------------

radTHMatrixACA::radTHMatrixACA(radTInteraction* interaction, double eps, int max_rank)
    : m_interaction(interaction)
    , m_eps(eps)
    , m_max_rank(max_rank)
    , m_min_cluster_size(16)  // Minimum elements per leaf cluster
    , m_eta(2.0)              // Admissibility parameter (larger = more compression)
    , m_n_elem(0)
    , m_n_dense_blocks(0)
    , m_n_lowrank_blocks(0)
    , m_dense_memory(0)
    , m_lowrank_memory(0)
    , m_is_built(false)
{
}

//-------------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------------

radTHMatrixACA::~radTHMatrixACA()
{
    // unique_ptr handles cleanup automatically
}

//-------------------------------------------------------------------------
// Build H-matrix structure
//-------------------------------------------------------------------------

int radTHMatrixACA::Build()
{
    if(m_interaction == nullptr) return -1;

    // Get number of elements from interaction matrix dimensions
    m_n_elem = m_interaction->AmOfMainElem;
    if(m_n_elem <= 0) return -1;

    // User has full control over H-matrix usage
    // No automatic problem size threshold (per CLAUDE.md policy)

    // Adjust min_cluster_size based on problem size
    if(m_n_elem < 1000) {
        m_min_cluster_size = 16;
    } else if(m_n_elem < 5000) {
        m_min_cluster_size = 32;
    } else {
        m_min_cluster_size = 64;
    }

    // Initialize ordering (identity permutation)
    m_ordering.resize(m_n_elem);
    for(int i = 0; i < m_n_elem; ++i) {
        m_ordering[i] = i;
    }

    // Build cluster tree
    m_cluster_tree.reset(BuildClusterTree(0, m_n_elem, 0));
    if(!m_cluster_tree) return -2;

    // Clear existing blocks
    m_blocks.clear();
    m_n_dense_blocks = 0;
    m_n_lowrank_blocks = 0;
    m_dense_memory = 0;
    m_lowrank_memory = 0;

    // Phase 1: Enumerate block structure (sequential)
    // Collect leaf block pairs that need to be computed
    std::vector<std::pair<radTCluster*, radTCluster*>> leaf_pairs;
    EnumerateLeafPairs(m_cluster_tree.get(), m_cluster_tree.get(), leaf_pairs);

    // Phase 2: Compute blocks in parallel (HACApK style)
    int n_pairs = (int)leaf_pairs.size();
    m_blocks.resize(n_pairs);

    // Thread-local statistics
    int total_dense = 0;
    int total_lowrank = 0;
    int64_t total_dense_mem = 0;
    int64_t total_lowrank_mem = 0;

    #ifdef _OPENMP
    #pragma omp parallel reduction(+:total_dense, total_lowrank, total_dense_mem, total_lowrank_mem)
    {
        #pragma omp for schedule(dynamic)
    #endif
        for(int p = 0; p < n_pairs; ++p) {
            radTCluster* row_cluster = leaf_pairs[p].first;
            radTCluster* col_cluster = leaf_pairs[p].second;

            // Try low-rank approximation if admissible
            if(IsAdmissible(row_cluster, col_cluster)) {
                radTHBlock block;
                block.type = radTHBlockType::LOW_RANK;
                block.lr_block = std::make_unique<radTLowRankBlock>();

                block.lr_block->row_start = row_cluster->start;
                block.lr_block->row_size = row_cluster->size;
                block.lr_block->col_start = col_cluster->start;
                block.lr_block->col_size = col_cluster->size;

                bool success = CreateLowRankBlock(
                    row_cluster->start, row_cluster->size,
                    col_cluster->start, col_cluster->size,
                    *block.lr_block);

                if(success) {
                    total_lowrank++;
                    total_lowrank_mem += (int64_t)(block.lr_block->row_size * 3 + block.lr_block->col_size * 3) *
                                        block.lr_block->rank * sizeof(double);
                    m_blocks[p] = std::move(block);
                    continue;
                }
            }

            // Create dense block
            radTHBlock block;
            block.type = radTHBlockType::DENSE;
            block.dense_block = std::make_unique<radTDenseBlock>();

            CreateDenseBlock(
                row_cluster->start, row_cluster->size,
                col_cluster->start, col_cluster->size,
                *block.dense_block);

            total_dense++;
            total_dense_mem += (int64_t)row_cluster->size * col_cluster->size * 9 * sizeof(double);
            m_blocks[p] = std::move(block);
        }
    #ifdef _OPENMP
    }
    #endif

    m_n_dense_blocks = total_dense;
    m_n_lowrank_blocks = total_lowrank;
    m_dense_memory = total_dense_mem;
    m_lowrank_memory = total_lowrank_mem;

    m_is_built = true;
    return 0;
}

//-------------------------------------------------------------------------
// Build cluster tree recursively
//-------------------------------------------------------------------------

radTCluster* radTHMatrixACA::BuildClusterTree(int start, int size, int depth)
{
    radTCluster* cluster = new radTCluster();
    cluster->start = start;
    cluster->size = size;

    // Compute bounding box
    ComputeBoundingBox(start, size, cluster->bbox);

    // Check if we should create a leaf
    if(size <= m_min_cluster_size || depth >= 25) {
        cluster->left = nullptr;
        cluster->right = nullptr;
        return cluster;
    }

    // Partition elements by longest axis
    int split = PartitionByLongestAxis(start, size);

    // Check if partition succeeded
    if(split <= start || split >= start + size) {
        cluster->left = nullptr;
        cluster->right = nullptr;
        return cluster;
    }

    // Create children
    int left_size = split - start;
    int right_size = size - left_size;

    cluster->left = BuildClusterTree(start, left_size, depth + 1);
    cluster->right = BuildClusterTree(split, right_size, depth + 1);

    return cluster;
}

//-------------------------------------------------------------------------
// Compute bounding box for a range of elements
//-------------------------------------------------------------------------

void radTHMatrixACA::ComputeBoundingBox(int start, int size, radTBoundingBox& bbox) const
{
    if(size <= 0 || m_interaction == nullptr) {
        bbox = radTBoundingBox();
        return;
    }

    bool first = true;

    for(int i = start; i < start + size; ++i) {
        int elem_idx = m_ordering[i];

        if(elem_idx < (int)m_interaction->MainElemCenPointsVect.size()) {
            const TVector3d& center = m_interaction->MainElemCenPointsVect[elem_idx];

            if(first) {
                bbox.xmin = bbox.xmax = center.x;
                bbox.ymin = bbox.ymax = center.y;
                bbox.zmin = bbox.zmax = center.z;
                first = false;
            } else {
                bbox.xmin = std::min(bbox.xmin, center.x);
                bbox.xmax = std::max(bbox.xmax, center.x);
                bbox.ymin = std::min(bbox.ymin, center.y);
                bbox.ymax = std::max(bbox.ymax, center.y);
                bbox.zmin = std::min(bbox.zmin, center.z);
                bbox.zmax = std::max(bbox.zmax, center.z);
            }
        }
    }
}

//-------------------------------------------------------------------------
// Partition elements by longest bounding box axis
//-------------------------------------------------------------------------

int radTHMatrixACA::PartitionByLongestAxis(int start, int size)
{
    if(size <= 1) return start;

    radTBoundingBox bbox;
    ComputeBoundingBox(start, size, bbox);

    double dx = bbox.xmax - bbox.xmin;
    double dy = bbox.ymax - bbox.ymin;
    double dz = bbox.zmax - bbox.zmin;

    int axis = 0;
    if(dy > dx && dy > dz) axis = 1;
    else if(dz > dx && dz > dy) axis = 2;

    auto compare = [this, axis](int a, int b) {
        if(a >= (int)m_interaction->MainElemCenPointsVect.size() ||
           b >= (int)m_interaction->MainElemCenPointsVect.size()) {
            return a < b;
        }
        const TVector3d& ca = m_interaction->MainElemCenPointsVect[a];
        const TVector3d& cb = m_interaction->MainElemCenPointsVect[b];
        switch(axis) {
            case 0: return ca.x < cb.x;
            case 1: return ca.y < cb.y;
            case 2: return ca.z < cb.z;
            default: return ca.x < cb.x;
        }
    };

    std::sort(m_ordering.begin() + start, m_ordering.begin() + start + size, compare);

    return start + size / 2;
}

//-------------------------------------------------------------------------
// Enumerate leaf block pairs for parallel processing
// This collects all leaf cluster pairs that need block computation
//-------------------------------------------------------------------------

void radTHMatrixACA::EnumerateLeafPairs(radTCluster* row_cluster, radTCluster* col_cluster,
                                         std::vector<std::pair<radTCluster*, radTCluster*>>& pairs)
{
    if(row_cluster == nullptr || col_cluster == nullptr) return;

    // If admissible, this is a leaf pair (will try ACA)
    if(IsAdmissible(row_cluster, col_cluster)) {
        pairs.push_back({row_cluster, col_cluster});
        return;
    }

    // Check if both clusters are leaves
    if(row_cluster->isLeaf() && col_cluster->isLeaf()) {
        // This is a dense leaf pair
        pairs.push_back({row_cluster, col_cluster});
        return;
    }

    // Recursive subdivision
    if(!row_cluster->isLeaf() && !col_cluster->isLeaf()) {
        EnumerateLeafPairs(row_cluster->left, col_cluster->left, pairs);
        EnumerateLeafPairs(row_cluster->left, col_cluster->right, pairs);
        EnumerateLeafPairs(row_cluster->right, col_cluster->left, pairs);
        EnumerateLeafPairs(row_cluster->right, col_cluster->right, pairs);
    } else if(!row_cluster->isLeaf()) {
        EnumerateLeafPairs(row_cluster->left, col_cluster, pairs);
        EnumerateLeafPairs(row_cluster->right, col_cluster, pairs);
    } else {
        EnumerateLeafPairs(row_cluster, col_cluster->left, pairs);
        EnumerateLeafPairs(row_cluster, col_cluster->right, pairs);
    }
}

//-------------------------------------------------------------------------
// Check admissibility condition (eta criterion)
//-------------------------------------------------------------------------

bool radTHMatrixACA::IsAdmissible(const radTCluster* c1, const radTCluster* c2) const
{
    // Admissibility criterion based on HACApK (ELF_MAGIC implementation)
    // Standard criterion: dist(c1, c2) >= eta * min(diam(c1), diam(c2))
    //
    // Parameters from ELF_MAGIC m_HACApK_base.f90:
    //   param(51) = 2.0  (eta - admissibility parameter)
    //   param(21) = 15   (minimum cluster size)

    // Require minimum cluster sizes for ACA to be effective
    // Small clusters have too few elements for meaningful low-rank approximation
    const int min_cluster_size = 15;  // From HACApK param(21)
    if(c1->size < min_cluster_size || c2->size < min_cluster_size) {
        return false;
    }

    // Compute distance between cluster bounding boxes
    double dist = c1->bbox.distance(c2->bbox);

    // Compute minimum diameter of the two clusters
    double diam1 = c1->bbox.diameter();
    double diam2 = c2->bbox.diameter();
    double min_diam = std::min(diam1, diam2);

    // Admissibility: well-separated clusters can be approximated as low-rank
    // Using eta = 2.0 from HACApK param(51)
    return dist >= m_eta * min_diam;
}

//-------------------------------------------------------------------------
// Create low-rank block using ACA
//-------------------------------------------------------------------------

bool radTHMatrixACA::CreateLowRankBlock(int row_start, int row_size,
                                        int col_start, int col_size,
                                        radTLowRankBlock& block)
{
    double approx_norm = 0.0;

    int rank = ACAPlus(row_start, row_size, col_start, col_size,
                       block.U, block.V, approx_norm);

    if(rank <= 0 || rank >= m_max_rank) {
        return false;
    }

    // Check if low-rank is actually beneficial
    int dense_size = row_size * 3 * col_size * 3;
    int lowrank_size = (row_size * 3 + col_size * 3) * rank;

    if(lowrank_size >= dense_size * 0.8) {
        // Not enough compression, use dense
        return false;
    }

    block.rank = rank;
    return true;
}

//-------------------------------------------------------------------------
// Create dense block (OpenMP parallelized)
//-------------------------------------------------------------------------

void radTHMatrixACA::CreateDenseBlock(int row_start, int row_size,
                                      int col_start, int col_size,
                                      radTDenseBlock& block)
{
    block.row_start = row_start;
    block.row_size = row_size;
    block.col_start = col_start;
    block.col_size = col_size;

    int dense_rows = row_size * 3;
    int dense_cols = col_size * 3;
    block.data.resize(dense_rows * dense_cols, 0.0);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(row_size > 8)
    #endif
    for(int i = 0; i < row_size; ++i) {
        double entry[9];
        for(int j = 0; j < col_size; ++j) {
            int row_elem = m_ordering[row_start + i];
            int col_elem = m_ordering[col_start + j];

            ComputeEntry(row_elem, col_elem, entry);

            for(int ii = 0; ii < 3; ++ii) {
                for(int jj = 0; jj < 3; ++jj) {
                    int row_idx = i * 3 + ii;
                    int col_idx = j * 3 + jj;
                    block.data[row_idx * dense_cols + col_idx] = entry[ii * 3 + jj];
                }
            }
        }
    }
}

//-------------------------------------------------------------------------
// Compute single matrix entry N[i][j] (3x3 block)
//-------------------------------------------------------------------------

void radTHMatrixACA::ComputeEntry(int i, int j, double* entry_3x3) const
{
    std::memset(entry_3x3, 0, 9 * sizeof(double));

    if(m_interaction == nullptr || m_interaction->IntrcMat == nullptr) return;

    int n = m_interaction->AmOfMainElem;
    double* mat = m_interaction->IntrcMat;

    int base_idx = (i * n + j) * 9;
    for(int k = 0; k < 9; ++k) {
        entry_3x3[k] = mat[base_idx + k];
    }
}

//-------------------------------------------------------------------------
// ACA+ (Adaptive Cross Approximation Plus) - HACApK style
// Based on cHACApK_acaplus from HACApK library
// Computes low-rank approximation A = U * V^T
//-------------------------------------------------------------------------

// Helper: Compute 2-norm of a vector
static double compute_norm(const double* vec, int n) {
    double sum = 0.0;
    for(int i = 0; i < n; ++i) {
        sum += vec[i] * vec[i];
    }
    return std::sqrt(sum);
}

// Helper: Find index and value of maximum absolute element
static void find_max_abs(const double* vec, int n, double& max_val, int& max_idx) {
    max_val = 0.0;
    max_idx = 0;
    for(int i = 0; i < n; ++i) {
        double abs_val = std::abs(vec[i]);
        if(abs_val > max_val) {
            max_val = abs_val;
            max_idx = i;
        }
    }
}

int radTHMatrixACA::ACAPlus(int row_start, int row_size,
                            int col_start, int col_size,
                            std::vector<double>& U, std::vector<double>& V,
                            double& approx_norm)
{
    // HACApK ACA+ algorithm adapted for Radia's 3x3 block structure
    // Reference: cHACApK_base.c, cHACApK_acaplus function

    const int ndl = row_size * 3;  // Number of DOF rows
    const int ndt = col_size * 3;  // Number of DOF columns
    const int kmax = std::min(m_max_rank, std::min(ndl, ndt));

    // HACApK parameters (from ELF_MAGIC m_HACApK_base.f90)
    const double ACA_EPS = m_eps;          // ACA tolerance
    const double za_ACA_EPS = 1.0e-30;     // Absolute zero threshold

    U.clear();
    V.clear();

    // Row and column masks (track used pivots)
    std::vector<int> lrow_msk(ndl, 0);
    std::vector<int> lcol_msk(ndt, 0);

    // Reference vectors for pivot selection (HACApK style)
    std::vector<double> pa_ref(ndl, 0.0);  // Reference column
    std::vector<double> pb_ref(ndt, 0.0);  // Reference row

    // Working vectors
    std::vector<double> zaa(ndl * kmax, 0.0);  // U matrix storage
    std::vector<double> zab(ndt * kmax, 0.0);  // V matrix storage

    int k = 0;
    int ntries = ndl + ndt + 1;
    int ntries_row = 6;
    int ntries_col = 6;

    // Initialize with arbitrary column j_ref = 0
    int j_ref = 0;

    // Compute initial reference column pa_ref = A[:, j_ref]
    {
        int col_elem_idx = j_ref / 3;
        int col_dof = j_ref % 3;
        int col_elem = m_ordering[col_start + col_elem_idx];

        for(int i = 0; i < row_size; ++i) {
            int row_elem = m_ordering[row_start + i];
            double entry[9];
            ComputeEntry(row_elem, col_elem, entry);

            for(int ii = 0; ii < 3; ++ii) {
                pa_ref[i * 3 + ii] = entry[ii * 3 + col_dof];
            }
        }
    }

    double colnorm = compute_norm(pa_ref.data(), ndl);

    // Find initial row pivot i_ref = argmax |pa_ref|
    double max_val;
    int i_ref;
    find_max_abs(pa_ref.data(), ndl, max_val, i_ref);

    // Compute reference row pb_ref = A[i_ref, :]
    {
        int row_elem_idx = i_ref / 3;
        int row_dof = i_ref % 3;
        int row_elem = m_ordering[row_start + row_elem_idx];

        for(int j = 0; j < col_size; ++j) {
            int col_elem = m_ordering[col_start + j];
            double entry[9];
            ComputeEntry(row_elem, col_elem, entry);

            for(int jj = 0; jj < 3; ++jj) {
                pb_ref[j * 3 + jj] = entry[row_dof * 3 + jj];
            }
        }
    }

    double rownorm = compute_norm(pb_ref.data(), ndt);
    double apxnorm = 0.0;
    int lstop_aca = 0;

    while((k < kmax) && (ntries_row > 0 || ntries_col > 0) && (ntries > 0)) {
        ntries--;

        double* pcol = &zaa[ndl * k];  // k-th column of U
        double* prow = &zab[ndt * k];  // k-th column of V (row of V^T)

        double col_maxval, row_maxval;
        int i, j;
        find_max_abs(pa_ref.data(), ndl, col_maxval, i);
        find_max_abs(pb_ref.data(), ndt, row_maxval, j);

        if(row_maxval > col_maxval) {
            // Row pivot dominant: compute column first
            if(j != j_ref) {
                // Compute column A[:, j]
                int col_elem_idx = j / 3;
                int col_dof = j % 3;
                int col_elem = m_ordering[col_start + col_elem_idx];

                for(int ii = 0; ii < row_size; ++ii) {
                    int row_elem = m_ordering[row_start + ii];
                    double entry[9];
                    ComputeEntry(row_elem, col_elem, entry);

                    for(int iii = 0; iii < 3; ++iii) {
                        pcol[ii * 3 + iii] = entry[iii * 3 + col_dof];
                    }
                }

                // Subtract approximation: pcol -= sum_t(zaa[:, t] * zab[j, t])
                for(int t = 0; t < k; ++t) {
                    double vt = zab[j + ndt * t];
                    for(int ii = 0; ii < ndl; ++ii) {
                        pcol[ii] -= zaa[ii + ndl * t] * vt;
                    }
                }

                // Apply row mask
                for(int ii = 0; ii < ndl; ++ii) {
                    if(lrow_msk[ii] == 1) pcol[ii] = 0.0;
                }
            } else {
                // Use reference column
                for(int ii = 0; ii < ndl; ++ii) {
                    pcol[ii] = pa_ref[ii];
                }
            }

            find_max_abs(pcol, ndl, col_maxval, i);

            if(col_maxval < ACA_EPS && k >= 1) {
                lstop_aca = 1;
            } else {
                // Compute row A[i, :]
                int row_elem_idx = i / 3;
                int row_dof = i % 3;
                int row_elem = m_ordering[row_start + row_elem_idx];

                for(int jj = 0; jj < col_size; ++jj) {
                    int col_elem = m_ordering[col_start + jj];
                    double entry[9];
                    ComputeEntry(row_elem, col_elem, entry);

                    for(int jjj = 0; jjj < 3; ++jjj) {
                        prow[jj * 3 + jjj] = entry[row_dof * 3 + jjj];
                    }
                }

                // Subtract approximation
                for(int t = 0; t < k; ++t) {
                    double ut = zaa[i + ndl * t];
                    for(int jj = 0; jj < ndt; ++jj) {
                        prow[jj] -= ut * zab[jj + ndt * t];
                    }
                }

                // Apply column mask
                for(int jj = 0; jj < ndt; ++jj) {
                    if(lcol_msk[jj] == 1) prow[jj] = 0.0;
                }

                if(std::abs(pcol[i]) > 1.0e-20) {
                    double zinvmax = 1.0 / pcol[i];
                    for(int ii = 0; ii < ndl; ++ii) {
                        pcol[ii] *= zinvmax;
                    }
                } else {
                    k = std::max(k - 1, 0);
                    break;
                }
            }
        } else {
            // Column pivot dominant: compute row first
            if(i != i_ref) {
                // Compute row A[i, :]
                int row_elem_idx = i / 3;
                int row_dof = i % 3;
                int row_elem = m_ordering[row_start + row_elem_idx];

                for(int jj = 0; jj < col_size; ++jj) {
                    int col_elem = m_ordering[col_start + jj];
                    double entry[9];
                    ComputeEntry(row_elem, col_elem, entry);

                    for(int jjj = 0; jjj < 3; ++jjj) {
                        prow[jj * 3 + jjj] = entry[row_dof * 3 + jjj];
                    }
                }

                // Subtract approximation
                for(int t = 0; t < k; ++t) {
                    double ut = zaa[i + ndl * t];
                    for(int jj = 0; jj < ndt; ++jj) {
                        prow[jj] -= ut * zab[jj + ndt * t];
                    }
                }

                // Apply column mask
                for(int jj = 0; jj < ndt; ++jj) {
                    if(lcol_msk[jj] == 1) prow[jj] = 0.0;
                }
            } else {
                // Use reference row
                for(int jj = 0; jj < ndt; ++jj) {
                    prow[jj] = pb_ref[jj];
                }
            }

            find_max_abs(prow, ndt, row_maxval, j);

            if(row_maxval < ACA_EPS && k >= 1) {
                lstop_aca = 1;
            } else {
                // Compute column A[:, j]
                int col_elem_idx = j / 3;
                int col_dof = j % 3;
                int col_elem = m_ordering[col_start + col_elem_idx];

                for(int ii = 0; ii < row_size; ++ii) {
                    int row_elem = m_ordering[row_start + ii];
                    double entry[9];
                    ComputeEntry(row_elem, col_elem, entry);

                    for(int iii = 0; iii < 3; ++iii) {
                        pcol[ii * 3 + iii] = entry[iii * 3 + col_dof];
                    }
                }

                // Subtract approximation
                for(int t = 0; t < k; ++t) {
                    double vt = zab[j + ndt * t];
                    for(int ii = 0; ii < ndl; ++ii) {
                        pcol[ii] -= zaa[ii + ndl * t] * vt;
                    }
                }

                // Apply row mask
                for(int ii = 0; ii < ndl; ++ii) {
                    if(lrow_msk[ii] == 1) pcol[ii] = 0.0;
                }

                if(std::abs(prow[j]) > 1.0e-20) {
                    double zinvmax = 1.0 / prow[j];
                    for(int jj = 0; jj < ndt; ++jj) {
                        prow[jj] *= zinvmax;
                    }
                } else {
                    k = std::max(k - 1, 0);
                    break;
                }
            }
        }

        // Mark pivots as used
        lrow_msk[i] = 1;
        lcol_msk[j] = 1;

        // Update reference vectors
        if(i != i_ref) {
            double zinvmax = -pcol[i_ref];
            for(int jj = 0; jj < ndt; ++jj) {
                pb_ref[jj] += prow[jj] * zinvmax;
            }
            rownorm = compute_norm(pb_ref.data(), ndt);
        }

        // Try to find new reference row
        if(i == i_ref || rownorm < ACA_EPS) {
            if(i == i_ref) ntries_row++;
            if(ntries_row > 0) {
                rownorm = 0.0;
                int ii = i_ref;
                while(ii != (i_ref + ndl - 1) % ndl && rownorm < za_ACA_EPS && ntries_row > 0) {
                    if(lrow_msk[ii] == 0) {
                        // Compute new reference row
                        int row_elem_idx = ii / 3;
                        int row_dof = ii % 3;
                        int row_elem = m_ordering[row_start + row_elem_idx];

                        for(int jj = 0; jj < col_size; ++jj) {
                            int col_elem = m_ordering[col_start + jj];
                            double entry[9];
                            ComputeEntry(row_elem, col_elem, entry);

                            for(int jjj = 0; jjj < 3; ++jjj) {
                                pb_ref[jj * 3 + jjj] = entry[row_dof * 3 + jjj];
                            }
                        }

                        // Subtract approximation
                        for(int t = 0; t <= k; ++t) {
                            double ut = zaa[ii + ndl * t];
                            for(int jj = 0; jj < ndt; ++jj) {
                                pb_ref[jj] -= ut * zab[jj + ndt * t];
                            }
                        }

                        // Apply column mask
                        for(int jj = 0; jj < ndt; ++jj) {
                            if(lcol_msk[jj] == 1) pb_ref[jj] = 0.0;
                        }

                        rownorm = compute_norm(pb_ref.data(), ndt);
                        if(rownorm < ACA_EPS) lrow_msk[ii] = 1;
                        ntries_row--;
                    } else {
                        rownorm = 0.0;
                    }
                    ii = (ii + 1) % ndl;
                }
                i_ref = (ii + ndl - 1) % ndl;
            }
        }

        // Update reference column
        if(j != j_ref) {
            double zinvmax = -prow[j_ref];
            for(int ii = 0; ii < ndl; ++ii) {
                pa_ref[ii] += pcol[ii] * zinvmax;
            }
            colnorm = compute_norm(pa_ref.data(), ndl);
        }

        // Try to find new reference column
        if(j == j_ref || colnorm < ACA_EPS) {
            if(j == j_ref) ntries_col++;
            if(ntries_col > 0) {
                colnorm = 0.0;
                int jj = j_ref;
                while(jj != (j_ref + ndt - 1) % ndt && colnorm < za_ACA_EPS && ntries_col > 0) {
                    if(lcol_msk[jj] == 0) {
                        // Compute new reference column
                        int col_elem_idx = jj / 3;
                        int col_dof = jj % 3;
                        int col_elem = m_ordering[col_start + col_elem_idx];

                        for(int ii = 0; ii < row_size; ++ii) {
                            int row_elem = m_ordering[row_start + ii];
                            double entry[9];
                            ComputeEntry(row_elem, col_elem, entry);

                            for(int iii = 0; iii < 3; ++iii) {
                                pa_ref[ii * 3 + iii] = entry[iii * 3 + col_dof];
                            }
                        }

                        // Subtract approximation
                        for(int t = 0; t <= k; ++t) {
                            double vt = zab[jj + ndt * t];
                            for(int ii = 0; ii < ndl; ++ii) {
                                pa_ref[ii] -= zaa[ii + ndl * t] * vt;
                            }
                        }

                        // Apply row mask
                        for(int ii = 0; ii < ndl; ++ii) {
                            if(lrow_msk[ii] == 1) pa_ref[ii] = 0.0;
                        }

                        colnorm = compute_norm(pa_ref.data(), ndl);
                        if(colnorm < ACA_EPS) lcol_msk[jj] = 1;
                        ntries_col--;
                    } else {
                        colnorm = 0.0;
                    }
                    jj = (jj + 1) % ndt;
                }
                j_ref = (jj + ndt - 1) % ndt;
            }
        }

        // Check convergence
        if(colnorm < ACA_EPS && rownorm < ACA_EPS && k >= 1) {
            lstop_aca = 1;
            k++;
        }

        if(lstop_aca == 0) {
            double blknorm = compute_norm(pcol, ndl) * compute_norm(prow, ndt);
            if(k == 0) {
                apxnorm = blknorm;
            } else {
                if(blknorm < apxnorm * m_eps &&
                   rownorm < apxnorm * m_eps &&
                   colnorm < apxnorm * m_eps &&
                   k >= 1) {
                    lstop_aca = 1;
                }
            }
        }

        if(lstop_aca == 1 && k >= 1) break;
        k++;
    }

    approx_norm = apxnorm;

    if(k == 0) return 0;

    // Pack U and V matrices
    U.resize(ndl * k);
    V.resize(ndt * k);

    for(int t = 0; t < k; ++t) {
        for(int ii = 0; ii < ndl; ++ii) {
            U[ii * k + t] = zaa[ii + ndl * t];
        }
        for(int jj = 0; jj < ndt; ++jj) {
            V[jj * k + t] = zab[jj + ndt * t];
        }
    }

    return k;
}

//-------------------------------------------------------------------------
// Matrix-vector product: y = A * x (OpenMP parallelized)
//-------------------------------------------------------------------------

void radTHMatrixACA::MatVec(const double* x, double* y) const
{
    if(!m_is_built || m_interaction == nullptr) return;

    int n_dof = m_n_elem * 3;

    // Zero output
    std::memset(y, 0, n_dof * sizeof(double));

    int n_blocks = (int)m_blocks.size();

    #ifdef _OPENMP
    // Parallel version with thread-local storage
    int n_threads = omp_get_max_threads();
    std::vector<std::vector<double>> y_local(n_threads, std::vector<double>(n_dof, 0.0));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        double* y_thread = y_local[tid].data();

        #pragma omp for schedule(dynamic)
        for(int b = 0; b < n_blocks; ++b) {
            const auto& block = m_blocks[b];
            if(block.type == radTHBlockType::LOW_RANK && block.lr_block) {
                LowRankMatVec(*block.lr_block, x, y_thread);
            } else if(block.type == radTHBlockType::DENSE && block.dense_block) {
                DenseMatVec(*block.dense_block, x, y_thread);
            }
        }
    }

    // Reduce thread-local results (OpenMP parallelized)
    #pragma omp parallel for
    for(int i = 0; i < n_dof; ++i) {
        for(int t = 0; t < n_threads; ++t) {
            y[i] += y_local[t][i];
        }
    }
    #else
    // Serial version
    for(int b = 0; b < n_blocks; ++b) {
        const auto& block = m_blocks[b];
        if(block.type == radTHBlockType::LOW_RANK && block.lr_block) {
            LowRankMatVec(*block.lr_block, x, y);
        } else if(block.type == radTHBlockType::DENSE && block.dense_block) {
            DenseMatVec(*block.dense_block, x, y);
        }
    }
    #endif
}

//-------------------------------------------------------------------------
// Low-rank block matrix-vector product: y += U * V^T * x
//-------------------------------------------------------------------------

void radTHMatrixACA::LowRankMatVec(const radTLowRankBlock& block, const double* x, double* y) const
{
    int m = block.row_size * 3;
    int n = block.col_size * 3;
    int rank = block.rank;

    if(rank <= 0) return;

    // Temporary storage for V^T * x
    std::vector<double> temp(rank, 0.0);

    // Compute temp = V^T * x (using permuted indices)
    for(int j = 0; j < block.col_size; ++j) {
        int col_elem = m_ordering[block.col_start + j];
        for(int jj = 0; jj < 3; ++jj) {
            int col_idx = j * 3 + jj;
            int x_idx = col_elem * 3 + jj;

            for(int r = 0; r < rank; ++r) {
                temp[r] += block.V[col_idx * rank + r] * x[x_idx];
            }
        }
    }

    // Compute y += U * temp (using permuted indices)
    for(int i = 0; i < block.row_size; ++i) {
        int row_elem = m_ordering[block.row_start + i];
        for(int ii = 0; ii < 3; ++ii) {
            int row_idx = i * 3 + ii;
            int y_idx = row_elem * 3 + ii;

            for(int r = 0; r < rank; ++r) {
                y[y_idx] += block.U[row_idx * rank + r] * temp[r];
            }
        }
    }
}

//-------------------------------------------------------------------------
// Dense block matrix-vector product: y += A * x
//-------------------------------------------------------------------------

void radTHMatrixACA::DenseMatVec(const radTDenseBlock& block, const double* x, double* y) const
{
    int dense_cols = block.col_size * 3;

    for(int i = 0; i < block.row_size; ++i) {
        int row_elem = m_ordering[block.row_start + i];

        for(int j = 0; j < block.col_size; ++j) {
            int col_elem = m_ordering[block.col_start + j];

            for(int ii = 0; ii < 3; ++ii) {
                int row_idx = i * 3 + ii;
                int y_idx = row_elem * 3 + ii;

                for(int jj = 0; jj < 3; ++jj) {
                    int col_idx = j * 3 + jj;
                    int x_idx = col_elem * 3 + jj;

                    y[y_idx] += block.data[row_idx * dense_cols + col_idx] * x[x_idx];
                }
            }
        }
    }
}

//-------------------------------------------------------------------------
// Get statistics
//-------------------------------------------------------------------------

void radTHMatrixACA::GetStatistics(int& n_dense, int& n_lowrank, double& compression_ratio) const
{
    n_dense = m_n_dense_blocks;
    n_lowrank = m_n_lowrank_blocks;

    int64_t full_memory = (int64_t)m_n_elem * m_n_elem * 9 * sizeof(double);
    int64_t hmat_memory = m_dense_memory + m_lowrank_memory;

    if(full_memory > 0) {
        compression_ratio = (double)hmat_memory / (double)full_memory;
    } else {
        compression_ratio = 1.0;
    }
}

//-------------------------------------------------------------------------
// Print statistics
//-------------------------------------------------------------------------

void radTHMatrixACA::PrintStatistics() const
{
    int64_t full_memory = (int64_t)m_n_elem * m_n_elem * 9 * sizeof(double);
    int64_t hmat_memory = m_dense_memory + m_lowrank_memory;

    std::cout << "H-Matrix Statistics:" << std::endl;
    std::cout << "  Elements: " << m_n_elem << std::endl;
    std::cout << "  Dense blocks: " << m_n_dense_blocks << std::endl;
    std::cout << "  Low-rank blocks: " << m_n_lowrank_blocks << std::endl;
    std::cout << "  Full matrix memory: " << std::fixed << std::setprecision(2)
              << (double)full_memory / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "  H-matrix memory: " << std::fixed << std::setprecision(2)
              << (double)hmat_memory / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "  Compression ratio: " << std::fixed << std::setprecision(2)
              << 100.0 * (double)hmat_memory / (double)full_memory << "%" << std::endl;
}

//-------------------------------------------------------------------------
