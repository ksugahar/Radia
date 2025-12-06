/*-------------------------------------------------------------------------
*
* File name:      rad_hmatrix.cpp
*
* Project:        RADIA
*
* Description:    H-matrix (Hierarchical Matrix) implementation
*                 Based on HACApK library (ppOpen-HPC project, MIT License)
*
* Author(s):      Claude Code AI Assistant
*
* First release:  2025-12
*
-------------------------------------------------------------------------*/

#include "rad_hmatrix.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

//-------------------------------------------------------------------------
// Global H-matrix enable flag
//-------------------------------------------------------------------------

bool g_UseHMatrix = false;

void RadSolverHMatrixEnable()
{
    g_UseHMatrix = true;
}

void RadSolverHMatrixDisable()
{
    g_UseHMatrix = false;
}

bool RadSolverHMatrixIsEnabled()
{
    return g_UseHMatrix;
}

//-------------------------------------------------------------------------
// BoundingBox methods
//-------------------------------------------------------------------------

double radTBoundingBox::distance(const radTBoundingBox& other) const
{
    // Compute distance between two bounding boxes
    double dx = 0.0, dy = 0.0, dz = 0.0;

    if (xmax < other.xmin) dx = other.xmin - xmax;
    else if (other.xmax < xmin) dx = xmin - other.xmax;

    if (ymax < other.ymin) dy = other.ymin - ymax;
    else if (other.ymax < ymin) dy = ymin - other.ymax;

    if (zmax < other.zmin) dz = other.zmin - zmax;
    else if (other.zmax < zmin) dz = zmin - other.zmax;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

//-------------------------------------------------------------------------
// radTHMatrix implementation
//-------------------------------------------------------------------------

radTHMatrix::radTHMatrix()
    : m_initialized(false)
    , m_n_elem(0)
    , m_ndof(0)
    , m_IntrcMat(nullptr)
    , m_dense_elements(0)
    , m_lowrank_elements(0)
{
}

radTHMatrix::~radTHMatrix()
{
    // Smart pointers handle cleanup
}

void radTHMatrix::Initialize(TMatrix3df** IntrcMat, int n_elem,
                             const std::vector<std::vector<double>>& centers,
                             const radTHMatrixParams& params)
{
    if (IntrcMat == nullptr || n_elem <= 0) {
        m_initialized = false;
        return;
    }

    m_IntrcMat = IntrcMat;
    m_n_elem = n_elem;
    m_ndof = 3 * n_elem;
    m_params = params;

    // Initialize permutation array (identity initially)
    m_perm.resize(n_elem);
    std::iota(m_perm.begin(), m_perm.end(), 0);

    // Build cluster tree
    m_root_cluster = std::make_unique<radTCluster>();
    m_root_cluster->start = 0;
    m_root_cluster->size = n_elem;

    BuildClusterTree(m_root_cluster.get(), m_perm, centers);

    // Compute inverse permutation
    m_inv_perm.resize(n_elem);
    for (int i = 0; i < n_elem; i++) {
        m_inv_perm[m_perm[i]] = i;
    }

    // Clear existing leaves
    m_leaves.clear();
    m_dense_elements = 0;
    m_lowrank_elements = 0;

    // Build H-matrix structure
    BuildHMatrix(m_root_cluster.get(), m_root_cluster.get());

    // Fill leaf matrices with data
    for (auto& leaf : m_leaves) {
        if (leaf.type == radTLeafMatrix::DENSE) {
            FillDenseLeaf(leaf);
        } else {
            FillLowRankLeaf(leaf);
        }
    }

    m_initialized = true;
}

void radTHMatrix::BuildClusterTree(radTCluster* cluster,
                                   std::vector<int>& perm,
                                   const std::vector<std::vector<double>>& centers)
{
    if (cluster == nullptr) return;

    int start = cluster->start;
    int size = cluster->size;

    // Compute bounding box for this cluster
    cluster->bbox.xmin = cluster->bbox.xmax = centers[perm[start]][0];
    cluster->bbox.ymin = cluster->bbox.ymax = centers[perm[start]][1];
    cluster->bbox.zmin = cluster->bbox.zmax = centers[perm[start]][2];

    for (int i = start; i < start + size; i++) {
        int idx = perm[i];
        cluster->bbox.xmin = std::min(cluster->bbox.xmin, centers[idx][0]);
        cluster->bbox.xmax = std::max(cluster->bbox.xmax, centers[idx][0]);
        cluster->bbox.ymin = std::min(cluster->bbox.ymin, centers[idx][1]);
        cluster->bbox.ymax = std::max(cluster->bbox.ymax, centers[idx][1]);
        cluster->bbox.zmin = std::min(cluster->bbox.zmin, centers[idx][2]);
        cluster->bbox.zmax = std::max(cluster->bbox.zmax, centers[idx][2]);
    }

    // Stop if cluster is small enough
    if (size <= m_params.min_leaf_size) {
        return;
    }

    // Find longest axis for splitting
    double dx = cluster->bbox.xmax - cluster->bbox.xmin;
    double dy = cluster->bbox.ymax - cluster->bbox.ymin;
    double dz = cluster->bbox.zmax - cluster->bbox.zmin;

    int split_axis = 0;
    if (dy >= dx && dy >= dz) split_axis = 1;
    else if (dz >= dx && dz >= dy) split_axis = 2;

    // Sort elements by coordinate along split axis
    std::sort(perm.begin() + start, perm.begin() + start + size,
              [&centers, split_axis](int a, int b) {
                  return centers[a][split_axis] < centers[b][split_axis];
              });

    // Split in half
    int mid = size / 2;

    // Create child clusters
    cluster->left = std::make_unique<radTCluster>();
    cluster->left->start = start;
    cluster->left->size = mid;

    cluster->right = std::make_unique<radTCluster>();
    cluster->right->start = start + mid;
    cluster->right->size = size - mid;

    // Recursively build children
    BuildClusterTree(cluster->left.get(), perm, centers);
    BuildClusterTree(cluster->right.get(), perm, centers);
}

void radTHMatrix::BuildHMatrix(radTCluster* row_cluster, radTCluster* col_cluster)
{
    if (row_cluster == nullptr || col_cluster == nullptr) return;

    // Check if both clusters are leaves
    if (row_cluster->isLeaf() && col_cluster->isLeaf()) {
        // Create leaf matrix
        radTLeafMatrix leaf;
        leaf.row_start = row_cluster->start;
        leaf.row_size = row_cluster->size;
        leaf.col_start = col_cluster->start;
        leaf.col_size = col_cluster->size;

        // Decide if admissible (low-rank) or dense
        if (IsAdmissible(row_cluster, col_cluster)) {
            leaf.type = radTLeafMatrix::LOW_RANK;
        } else {
            leaf.type = radTLeafMatrix::DENSE;
        }

        m_leaves.push_back(std::move(leaf));
        return;
    }

    // Check admissibility for non-leaf blocks
    if (IsAdmissible(row_cluster, col_cluster)) {
        // Create low-rank leaf
        radTLeafMatrix leaf;
        leaf.type = radTLeafMatrix::LOW_RANK;
        leaf.row_start = row_cluster->start;
        leaf.row_size = row_cluster->size;
        leaf.col_start = col_cluster->start;
        leaf.col_size = col_cluster->size;
        m_leaves.push_back(std::move(leaf));
        return;
    }

    // Recursively subdivide
    if (!row_cluster->isLeaf() && !col_cluster->isLeaf()) {
        // Both have children - split both
        BuildHMatrix(row_cluster->left.get(), col_cluster->left.get());
        BuildHMatrix(row_cluster->left.get(), col_cluster->right.get());
        BuildHMatrix(row_cluster->right.get(), col_cluster->left.get());
        BuildHMatrix(row_cluster->right.get(), col_cluster->right.get());
    } else if (!row_cluster->isLeaf()) {
        // Only row cluster has children
        BuildHMatrix(row_cluster->left.get(), col_cluster);
        BuildHMatrix(row_cluster->right.get(), col_cluster);
    } else {
        // Only col cluster has children
        BuildHMatrix(row_cluster, col_cluster->left.get());
        BuildHMatrix(row_cluster, col_cluster->right.get());
    }
}

bool radTHMatrix::IsAdmissible(const radTCluster* c1, const radTCluster* c2) const
{
    // For small problems (< 500 elements), don't use low-rank approximation
    // The overhead is not worth it and can cause accuracy issues
    if(m_n_elem < 500) return false;

    // Standard admissibility: dist(c1, c2) >= eta * min(diam(c1), diam(c2))
    double dist = c1->bbox.distance(c2->bbox);
    double min_diam = std::min(c1->bbox.diameter(), c2->bbox.diameter());

    // Require significant cluster sizes for ACA to be effective
    if(c1->size < 16 || c2->size < 16) return false;

    return dist >= m_params.eta * min_diam;
}

double radTHMatrix::GetEntry(int i, int j, int comp_i, int comp_j) const
{
    // Get matrix entry from Radia's interaction matrix
    // i, j are element indices (after permutation)
    // comp_i, comp_j are component indices (0, 1, 2)
    int orig_i = m_perm[i];
    int orig_j = m_perm[j];

    TMatrix3df& Nij = m_IntrcMat[orig_i][orig_j];

    // Return -N[i][j][comp_i][comp_j] (the base matrix without chi)
    double val = 0.0;
    if (comp_i == 0) {
        if (comp_j == 0) val = Nij.Str0.x;
        else if (comp_j == 1) val = Nij.Str0.y;
        else val = Nij.Str0.z;
    } else if (comp_i == 1) {
        if (comp_j == 0) val = Nij.Str1.x;
        else if (comp_j == 1) val = Nij.Str1.y;
        else val = Nij.Str1.z;
    } else {
        if (comp_j == 0) val = Nij.Str2.x;
        else if (comp_j == 1) val = Nij.Str2.y;
        else val = Nij.Str2.z;
    }

    return -val;  // Return -N (base matrix)
}

void radTHMatrix::FillDenseLeaf(radTLeafMatrix& leaf)
{
    // Fill dense block with matrix entries
    // Each element is a 3x3 block, so total size is 3*row_size x 3*col_size
    int nrows = 3 * leaf.row_size;
    int ncols = 3 * leaf.col_size;

    leaf.dense_data.resize(nrows * ncols);
    m_dense_elements += nrows * ncols;

    for (int i = 0; i < leaf.row_size; i++) {
        for (int j = 0; j < leaf.col_size; j++) {
            int elem_i = leaf.row_start + i;
            int elem_j = leaf.col_start + j;

            for (int ci = 0; ci < 3; ci++) {
                for (int cj = 0; cj < 3; cj++) {
                    int row = 3 * i + ci;
                    int col = 3 * j + cj;
                    leaf.dense_data[row * ncols + col] = GetEntry(elem_i, elem_j, ci, cj);
                }
            }
        }
    }
}

void radTHMatrix::FillLowRankLeaf(radTLeafMatrix& leaf)
{
    // Use ACA+ to compute low-rank approximation
    int rank = ACAPlus(leaf, m_params.eps, m_params.max_rank);
    leaf.rank = rank;

    m_lowrank_elements += (int64_t)leaf.row_size * 3 * rank +
                          (int64_t)leaf.col_size * 3 * rank;
}

int radTHMatrix::ACAPlus(radTLeafMatrix& leaf, double eps, int max_rank)
{
    // ACA+ (Adaptive Cross Approximation with recompression)
    // Reference: Bebendorf, M. "Approximation of boundary element matrices"

    int nrows = 3 * leaf.row_size;
    int ncols = 3 * leaf.col_size;

    // Limit max_rank based on matrix dimensions
    max_rank = std::min(max_rank, std::min(nrows, ncols));

    std::vector<double> U_temp(nrows * max_rank);
    std::vector<double> V_temp(ncols * max_rank);

    std::vector<bool> row_used(nrows, false);
    std::vector<bool> col_used(ncols, false);

    // Compute Frobenius norm estimate for tolerance
    double fnorm2 = 0.0;
    int rank = 0;

    // Start with random row
    int pivot_row = 0;

    for (rank = 0; rank < max_rank; rank++) {
        // Find pivot column in current row
        double max_val = 0.0;
        int pivot_col = -1;

        for (int j = 0; j < ncols; j++) {
            if (col_used[j]) continue;

            int elem_j = leaf.col_start + j / 3;
            int comp_j = j % 3;
            int elem_i = leaf.row_start + pivot_row / 3;
            int comp_i = pivot_row % 3;

            double val = GetEntry(elem_i, elem_j, comp_i, comp_j);

            // Subtract current approximation
            for (int k = 0; k < rank; k++) {
                val -= U_temp[pivot_row * max_rank + k] * V_temp[j * max_rank + k];
            }

            if (std::abs(val) > std::abs(max_val)) {
                max_val = val;
                pivot_col = j;
            }
        }

        if (pivot_col < 0 || std::abs(max_val) < 1.0e-15) {
            break;  // No suitable pivot found
        }

        col_used[pivot_col] = true;

        // Compute V column (entire column at pivot_col)
        double inv_pivot = 1.0 / max_val;
        for (int j = 0; j < ncols; j++) {
            int elem_j = leaf.col_start + j / 3;
            int comp_j = j % 3;
            int elem_i = leaf.row_start + pivot_row / 3;
            int comp_i = pivot_row % 3;

            double val = GetEntry(elem_i, elem_j, comp_i, comp_j);

            for (int k = 0; k < rank; k++) {
                val -= U_temp[pivot_row * max_rank + k] * V_temp[j * max_rank + k];
            }

            V_temp[j * max_rank + rank] = val * inv_pivot;
        }

        // Compute U column (entire row at pivot_col, transposed)
        for (int i = 0; i < nrows; i++) {
            int elem_i = leaf.row_start + i / 3;
            int comp_i = i % 3;
            int elem_j = leaf.col_start + pivot_col / 3;
            int comp_j = pivot_col % 3;

            double val = GetEntry(elem_i, elem_j, comp_i, comp_j);

            for (int k = 0; k < rank; k++) {
                val -= U_temp[i * max_rank + k] * V_temp[pivot_col * max_rank + k];
            }

            U_temp[i * max_rank + rank] = val;
        }

        // Update Frobenius norm estimate
        double u_norm2 = 0.0, v_norm2 = 0.0;
        for (int i = 0; i < nrows; i++) {
            u_norm2 += U_temp[i * max_rank + rank] * U_temp[i * max_rank + rank];
        }
        for (int j = 0; j < ncols; j++) {
            v_norm2 += V_temp[j * max_rank + rank] * V_temp[j * max_rank + rank];
        }

        double delta_norm2 = u_norm2 * v_norm2;
        fnorm2 += delta_norm2;

        // Check convergence
        if (delta_norm2 < eps * eps * fnorm2) {
            rank++;
            break;
        }

        // Find next pivot row
        max_val = 0.0;
        int next_row = -1;
        for (int i = 0; i < nrows; i++) {
            if (row_used[i]) continue;
            if (std::abs(U_temp[i * max_rank + rank]) > std::abs(max_val)) {
                max_val = U_temp[i * max_rank + rank];
                next_row = i;
            }
        }

        if (next_row < 0) {
            rank++;
            break;
        }

        row_used[pivot_row] = true;
        pivot_row = next_row;
    }

    // Copy U and V to leaf with final rank
    leaf.U.resize(nrows * rank);
    leaf.V.resize(ncols * rank);

    for (int i = 0; i < nrows; i++) {
        for (int k = 0; k < rank; k++) {
            leaf.U[i * rank + k] = U_temp[i * max_rank + k];
        }
    }

    for (int j = 0; j < ncols; j++) {
        for (int k = 0; k < rank; k++) {
            leaf.V[j * rank + k] = V_temp[j * max_rank + k];
        }
    }

    return rank;
}

void radTHMatrix::LeafMatVec(const radTLeafMatrix& leaf,
                             const std::vector<double>& x,
                             std::vector<double>& y) const
{
    // Matrix-vector product for a single leaf block
    // Note: x and y are in permuted order

    int nrows = 3 * leaf.row_size;
    int ncols = 3 * leaf.col_size;

    if (leaf.type == radTLeafMatrix::DENSE) {
        // Dense matvec: y[row_start:row_end] += dense * x[col_start:col_end]
        for (int i = 0; i < nrows; i++) {
            int row = 3 * leaf.row_start + i;
            double sum = 0.0;
            for (int j = 0; j < ncols; j++) {
                int col = 3 * leaf.col_start + j;
                sum += leaf.dense_data[i * ncols + j] * x[col];
            }
            #pragma omp atomic
            y[row] += sum;
        }
    } else {
        // Low-rank matvec: y += U * (V^T * x)
        int rank = leaf.rank;
        if (rank <= 0) return;

        // First: temp = V^T * x[col_start:col_end]
        std::vector<double> temp(rank, 0.0);
        for (int j = 0; j < ncols; j++) {
            int col = 3 * leaf.col_start + j;
            for (int k = 0; k < rank; k++) {
                temp[k] += leaf.V[j * rank + k] * x[col];
            }
        }

        // Second: y[row_start:row_end] += U * temp
        for (int i = 0; i < nrows; i++) {
            int row = 3 * leaf.row_start + i;
            double sum = 0.0;
            for (int k = 0; k < rank; k++) {
                sum += leaf.U[i * rank + k] * temp[k];
            }
            #pragma omp atomic
            y[row] += sum;
        }
    }
}

void radTHMatrix::MatVec(const std::vector<double>& x, std::vector<double>& y,
                         const std::vector<double>& inv_chi)
{
    if (!m_initialized) {
        throw std::runtime_error("H-matrix not initialized");
    }

    // Permute input vector
    std::vector<double> x_perm(m_ndof);
    for (int i = 0; i < m_n_elem; i++) {
        int perm_i = m_inv_perm[i];  // Where does original element i go?
        x_perm[3 * perm_i + 0] = x[3 * i + 0];
        x_perm[3 * perm_i + 1] = x[3 * i + 1];
        x_perm[3 * perm_i + 2] = x[3 * i + 2];
    }

    // Initialize output (permuted)
    std::vector<double> y_perm(m_ndof, 0.0);

    // Process all leaf blocks
    #pragma omp parallel for if(m_leaves.size() > 10)
    for (size_t leaf_idx = 0; leaf_idx < m_leaves.size(); leaf_idx++) {
        LeafMatVec(m_leaves[leaf_idx], x_perm, y_perm);
    }

    // Add diagonal (1/chi) contribution and permute back
    for (int i = 0; i < m_n_elem; i++) {
        int perm_i = m_inv_perm[i];

        // y = (base_matrix + diag(1/chi)) * x
        // base_matrix part is already in y_perm
        // Add diagonal 1/chi * x
        y[3 * i + 0] = y_perm[3 * perm_i + 0] + inv_chi[3 * i + 0] * x[3 * i + 0];
        y[3 * i + 1] = y_perm[3 * perm_i + 1] + inv_chi[3 * i + 1] * x[3 * i + 1];
        y[3 * i + 2] = y_perm[3 * perm_i + 2] + inv_chi[3 * i + 2] * x[3 * i + 2];
    }
}

int radTHMatrix::GetNumDenseLeaves() const
{
    int count = 0;
    for (const auto& leaf : m_leaves) {
        if (leaf.type == radTLeafMatrix::DENSE) count++;
    }
    return count;
}

int radTHMatrix::GetNumLowRankLeaves() const
{
    int count = 0;
    for (const auto& leaf : m_leaves) {
        if (leaf.type == radTLeafMatrix::LOW_RANK) count++;
    }
    return count;
}

double radTHMatrix::GetCompressionRatio() const
{
    if (m_n_elem <= 0) return 1.0;

    int64_t full_size = (int64_t)m_ndof * m_ndof;
    int64_t compressed_size = m_dense_elements + m_lowrank_elements;

    if (compressed_size <= 0) return 1.0;

    return (double)full_size / (double)compressed_size;
}
