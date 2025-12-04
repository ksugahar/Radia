/*-------------------------------------------------------------------------*
*
*   Radia H-Matrix with ACA Compression
*
*   Hierarchical matrix implementation for O(N log N) solver acceleration
*   Uses Adaptive Cross Approximation (ACA+) for low-rank block compression
*
*   Based on HACApK concepts but simplified for single-process Radia usage
*
*   Copyright (C): 2024 by European Synchrotron Radiation Facility, France
*
*-------------------------------------------------------------------------*/

#ifndef RAD_HMATRIX_ACA_H
#define RAD_HMATRIX_ACA_H

#include <vector>
#include <memory>
#include <functional>
#include <utility>
#include <cmath>

//-------------------------------------------------------------------------
// Forward declarations
//-------------------------------------------------------------------------

class radTInteraction;

//-------------------------------------------------------------------------
// Bounding Box for clustering
//-------------------------------------------------------------------------

struct radTBoundingBox {
    double xmin, xmax, ymin, ymax, zmin, zmax;

    radTBoundingBox() : xmin(0), xmax(0), ymin(0), ymax(0), zmin(0), zmax(0) {}

    double diameter() const {
        double dx = xmax - xmin;
        double dy = ymax - ymin;
        double dz = zmax - zmin;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    void center(double& cx, double& cy, double& cz) const {
        cx = 0.5 * (xmin + xmax);
        cy = 0.5 * (ymin + ymax);
        cz = 0.5 * (zmin + zmax);
    }

    double distance(const radTBoundingBox& other) const {
        double cx1, cy1, cz1, cx2, cy2, cz2;
        center(cx1, cy1, cz1);
        other.center(cx2, cy2, cz2);
        double dx = cx1 - cx2;
        double dy = cy1 - cy2;
        double dz = cz1 - cz2;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
};

//-------------------------------------------------------------------------
// Cluster Tree Node
//-------------------------------------------------------------------------

struct radTCluster {
    int start;          // Start index in element ordering
    int size;           // Number of elements in cluster
    radTBoundingBox bbox;
    radTCluster* left;   // Left child (null for leaf)
    radTCluster* right;  // Right child (null for leaf)

    radTCluster() : start(0), size(0), left(nullptr), right(nullptr) {}
    ~radTCluster() {
        delete left;
        delete right;
    }

    bool isLeaf() const { return left == nullptr && right == nullptr; }
};

//-------------------------------------------------------------------------
// Low-Rank Block (A = U * V^T)
//-------------------------------------------------------------------------

struct radTLowRankBlock {
    int row_start, row_size;    // Row cluster indices
    int col_start, col_size;    // Column cluster indices
    int rank;                   // Rank of approximation
    std::vector<double> U;      // U matrix (row_size * 3 x rank)
    std::vector<double> V;      // V matrix (col_size * 3 x rank)
};

//-------------------------------------------------------------------------
// Dense Block
//-------------------------------------------------------------------------

struct radTDenseBlock {
    int row_start, row_size;
    int col_start, col_size;
    std::vector<double> data;   // Dense matrix data (row_size*3 x col_size*3)
};

//-------------------------------------------------------------------------
// H-Matrix Block Types
//-------------------------------------------------------------------------

enum class radTHBlockType {
    DENSE,      // Full dense storage
    LOW_RANK    // ACA compressed
};

//-------------------------------------------------------------------------
// H-Matrix Block
//-------------------------------------------------------------------------

struct radTHBlock {
    radTHBlockType type;
    std::unique_ptr<radTLowRankBlock> lr_block;
    std::unique_ptr<radTDenseBlock> dense_block;
};

//-------------------------------------------------------------------------
// H-Matrix ACA Class
//-------------------------------------------------------------------------

class radTHMatrixACA {
public:
    // Constructor
    // Parameters tuned for MMM (Magnetic Moment Method):
    //   eps = 1e-5 (tighter than HACApK default for better accuracy)
    //   max_rank = 100 (reasonable limit for MMM)
    radTHMatrixACA(radTInteraction* interaction, double eps = 1e-5, int max_rank = 100);

    // Destructor
    ~radTHMatrixACA();

    // Build H-matrix structure
    int Build();

    // Matrix-vector product: y = A * x (where A is the interaction matrix N)
    // Input: x is magnetization vector (3*n_elem)
    // Output: y is field vector (3*n_elem)
    void MatVec(const double* x, double* y) const;

    // Get statistics
    void GetStatistics(int& n_dense, int& n_lowrank, double& compression_ratio) const;

    // Print statistics
    void PrintStatistics() const;

private:
    // Build cluster tree
    radTCluster* BuildClusterTree(int start, int size, int depth);

    // Enumerate leaf block pairs for parallel processing
    void EnumerateLeafPairs(radTCluster* row_cluster, radTCluster* col_cluster,
                            std::vector<std::pair<radTCluster*, radTCluster*>>& pairs);

    // Check admissibility condition (eta criterion)
    bool IsAdmissible(const radTCluster* c1, const radTCluster* c2) const;

    // Create low-rank block using ACA+
    bool CreateLowRankBlock(int row_start, int row_size,
                            int col_start, int col_size,
                            radTLowRankBlock& block);

    // Create dense block
    void CreateDenseBlock(int row_start, int row_size,
                          int col_start, int col_size,
                          radTDenseBlock& block);

    // Compute single matrix entry N[i][j] (3x3 block)
    void ComputeEntry(int i, int j, double* entry_3x3) const;

    // ACA+ algorithm
    int ACAPlus(int row_start, int row_size,
                int col_start, int col_size,
                std::vector<double>& U, std::vector<double>& V,
                double& approx_norm);

    // Compute bounding box for elements
    void ComputeBoundingBox(int start, int size, radTBoundingBox& bbox) const;

    // Partition elements by longest axis
    int PartitionByLongestAxis(int start, int size);

    // Matrix-vector product for low-rank block
    void LowRankMatVec(const radTLowRankBlock& block, const double* x, double* y) const;

    // Matrix-vector product for dense block
    void DenseMatVec(const radTDenseBlock& block, const double* x, double* y) const;

private:
    radTInteraction* m_interaction;     // Pointer to interaction object
    double m_eps;                       // ACA tolerance
    int m_max_rank;                     // Maximum rank for low-rank blocks
    int m_min_cluster_size;             // Minimum cluster size (leaf threshold)
    double m_eta;                       // Admissibility parameter

    int m_n_elem;                       // Number of elements
    std::vector<int> m_ordering;        // Element ordering for clustering

    std::unique_ptr<radTCluster> m_cluster_tree;  // Cluster tree
    std::vector<radTHBlock> m_blocks;             // H-matrix blocks

    // Statistics
    int m_n_dense_blocks;
    int m_n_lowrank_blocks;
    int64_t m_dense_memory;
    int64_t m_lowrank_memory;
    bool m_is_built;
};

//-------------------------------------------------------------------------

#endif // RAD_HMATRIX_ACA_H
