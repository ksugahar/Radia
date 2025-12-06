/*-------------------------------------------------------------------------
*
* File name:      rad_hmatrix.h
*
* Project:        RADIA
*
* Description:    H-matrix (Hierarchical Matrix) for fast matrix-vector products
*                 Used with BiCGSTAB solver for O(N log N) complexity
*
* Author(s):      Claude Code AI Assistant
*
* First release:  2025-12
*
* Reference:      HACApK library (ppOpen-HPC project, MIT License)
*                 Ida, A., et al. "HACApK: An H-matrix library for Krylov methods"
*
-------------------------------------------------------------------------*/

#ifndef __RAD_HMATRIX_H
#define __RAD_HMATRIX_H

#include "gmvectf.h"
#include <vector>
#include <memory>
#include <functional>

//-------------------------------------------------------------------------
// Global control: Enable/Disable H-matrix acceleration
//-------------------------------------------------------------------------

// Global H-matrix enable flag (set via Python API)
extern bool g_UseHMatrix;

// API functions to control H-matrix usage
void RadSolverHMatrixEnable();
void RadSolverHMatrixDisable();
bool RadSolverHMatrixIsEnabled();

//-------------------------------------------------------------------------
// H-matrix parameters
//-------------------------------------------------------------------------

struct radTHMatrixParams {
    double eta;           // Admissibility parameter (default: 2.0)
    double eps;           // ACA+ tolerance (default: 1e-4)
    int min_leaf_size;    // Minimum leaf size (default: 32)
    int max_rank;         // Maximum rank for low-rank blocks (default: 100)

    radTHMatrixParams()
        : eta(2.0), eps(1.0e-4), min_leaf_size(32), max_rank(100) {}
};

//-------------------------------------------------------------------------
// Axis-aligned bounding box for clustering
//-------------------------------------------------------------------------

struct radTBoundingBox {
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;

    radTBoundingBox() : xmin(0), xmax(0), ymin(0), ymax(0), zmin(0), zmax(0) {}

    // Compute diameter (maximum extent)
    double diameter() const {
        double dx = xmax - xmin;
        double dy = ymax - ymin;
        double dz = zmax - zmin;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    // Compute center
    void center(double& cx, double& cy, double& cz) const {
        cx = 0.5 * (xmin + xmax);
        cy = 0.5 * (ymin + ymax);
        cz = 0.5 * (zmin + zmax);
    }

    // Compute distance to another bounding box
    double distance(const radTBoundingBox& other) const;
};

//-------------------------------------------------------------------------
// Cluster structure for hierarchical partitioning
//-------------------------------------------------------------------------

struct radTCluster {
    int start;          // Start index in permutation array
    int size;           // Number of elements in this cluster
    radTBoundingBox bbox;
    std::unique_ptr<radTCluster> left;   // Left child
    std::unique_ptr<radTCluster> right;  // Right child

    radTCluster() : start(0), size(0), left(nullptr), right(nullptr) {}

    bool isLeaf() const { return left == nullptr && right == nullptr; }
};

//-------------------------------------------------------------------------
// Leaf matrix block (either dense or low-rank)
//-------------------------------------------------------------------------

struct radTLeafMatrix {
    enum Type { DENSE, LOW_RANK };

    Type type;
    int row_start, row_size;    // Row cluster range
    int col_start, col_size;    // Column cluster range

    // For DENSE blocks
    std::vector<double> dense_data;  // row_size * col_size elements

    // For LOW_RANK blocks: A approx= U * V^T
    int rank;
    std::vector<double> U;  // row_size * rank
    std::vector<double> V;  // col_size * rank

    radTLeafMatrix() : type(DENSE), row_start(0), row_size(0),
                       col_start(0), col_size(0), rank(0) {}
};

//-------------------------------------------------------------------------
// H-matrix class for Radia interaction matrix
//-------------------------------------------------------------------------

class radTHMatrix {
public:
    // Constructor
    radTHMatrix();
    ~radTHMatrix();

    // Initialize from Radia interaction matrix
    // IntrcMat: TMatrix3df** interaction matrix (3x3 blocks)
    // n_elem: number of elements
    // centers: element center coordinates [n_elem][3]
    // params: H-matrix parameters
    void Initialize(TMatrix3df** IntrcMat, int n_elem,
                    const std::vector<std::vector<double>>& centers,
                    const radTHMatrixParams& params);

    // Matrix-vector product: y = A * x
    // System matrix A = -N + diag(inv_chi)
    // x: input vector [ndof]
    // y: output vector [ndof]
    // inv_chi: diagonal 1/chi values [ndof]
    void MatVec(const std::vector<double>& x, std::vector<double>& y,
                const std::vector<double>& inv_chi);

    // Check if H-matrix is initialized
    bool IsInitialized() const { return m_initialized; }

    // Get statistics
    int GetNumElements() const { return m_n_elem; }
    int GetNumLeaves() const { return static_cast<int>(m_leaves.size()); }
    int GetNumDenseLeaves() const;
    int GetNumLowRankLeaves() const;
    double GetCompressionRatio() const;

private:
    // Build cluster tree using bisection
    void BuildClusterTree(radTCluster* cluster,
                          std::vector<int>& perm,
                          const std::vector<std::vector<double>>& centers);

    // Build H-matrix structure (determine which blocks are admissible)
    void BuildHMatrix(radTCluster* row_cluster, radTCluster* col_cluster);

    // Check admissibility criterion
    bool IsAdmissible(const radTCluster* c1, const radTCluster* c2) const;

    // Fill leaf matrix with dense data
    void FillDenseLeaf(radTLeafMatrix& leaf);

    // Fill leaf matrix with low-rank ACA+ approximation
    void FillLowRankLeaf(radTLeafMatrix& leaf);

    // ACA+ algorithm for low-rank approximation
    int ACAPlus(radTLeafMatrix& leaf, double eps, int max_rank);

    // Get matrix entry (i, j) from interaction matrix
    // i, j are 3x3 block indices, comp_i, comp_j are component indices (0-2)
    double GetEntry(int i, int j, int comp_i, int comp_j) const;

    // Matrix-vector product for a single leaf
    void LeafMatVec(const radTLeafMatrix& leaf,
                    const std::vector<double>& x,
                    std::vector<double>& y) const;

    // Member variables
    bool m_initialized;
    int m_n_elem;       // Number of elements (3x3 blocks)
    int m_ndof;         // Number of DOF (= 3 * n_elem)

    TMatrix3df** m_IntrcMat;  // Pointer to Radia's interaction matrix

    radTHMatrixParams m_params;
    std::vector<int> m_perm;           // Permutation array
    std::vector<int> m_inv_perm;       // Inverse permutation
    std::unique_ptr<radTCluster> m_root_cluster;
    std::vector<radTLeafMatrix> m_leaves;

    // Statistics
    int64_t m_dense_elements;
    int64_t m_lowrank_elements;
};

#endif // __RAD_HMATRIX_H
