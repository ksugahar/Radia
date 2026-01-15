/*
 * rad_peec_matrices.h
 *
 * PEEC (Partial Element Equivalent Circuit) Matrix Construction
 *
 * Loop-Star decomposition for quasi-static electromagnetic analysis.
 * Uses Darwin approximation: Laplace kernel G(r) = 1/(4*pi*r)
 *
 * Matrix Types:
 *   L: Inductance matrix (Loop-Loop) - Neumann formula
 *   P: Potential coefficient matrix (Star-Star) - Capacitance = P^{-1}
 *   R: Resistance matrix (diagonal)
 *   M_LS: Loop-Star coupling matrix
 *
 * PEEC System Equation:
 *   [Z_LL   Z_LS] [I_L]   [V_L]
 *   [Z_SL   Z_SS] [I_S] = [V_S]
 *
 *   Z_LL = R + jw*L + Z_s (with optional skin effect)
 *   Z_SS = P / jw
 *   Z_LS = jw*M_LS
 *   Z_SL = Z_LS^T (reciprocity)
 *
 * Valid Frequency Range: DC to ~100 MHz (Darwin approximation)
 *
 * References:
 * [1] A. E. Ruehli, "Equivalent Circuit Models...", IEEE MTT, 1974
 * [2] G. Vecchi, "Loop-Star decomposition...", IEEE TAP, 1999
 *
 * Part of Radia project
 */

#ifndef RAD_PEEC_MATRICES_H
#define RAD_PEEC_MATRICES_H

#include "gmvect.h"
#include "rad_constants.h"
#include <vector>
#include <complex>
#include <memory>

namespace radia {

// Physical constants
constexpr double PEEC_MU_0 = RadConst::MU_0;         // 4*pi*1e-7 H/m
constexpr double PEEC_EPS_0 = 8.854187817e-12;       // F/m
constexpr double PEEC_INV_FOUR_PI = RadConst::INV_FOUR_PI;

/**
 * @brief PEEC conductor segment (Loop element)
 *
 * Represents a filament segment for PEEC analysis.
 */
struct PEECSegment {
    TVector3d center;       // Segment center [m]
    TVector3d direction;    // Unit direction vector
    double length;          // Segment length [m]
    double width;           // Cross-section width [m]
    double height;          // Cross-section height [m]
    double sigma;           // Conductivity [S/m]

    PEECSegment()
        : length(0), width(0), height(0), sigma(5.8e7) {
        center = TVector3d(0, 0, 0);
        direction = TVector3d(1, 0, 0);
    }

    PEECSegment(const TVector3d& c, const TVector3d& d, double l,
                double w, double h, double s = 5.8e7)
        : center(c), direction(d), length(l), width(w), height(h), sigma(s) {}

    // Cross-section area [m^2]
    double area() const { return width * height; }

    // DC resistance [Ohm]
    double resistance() const {
        double a = area();
        return (a > 0 && sigma > 0) ? length / (sigma * a) : 0.0;
    }
};

/**
 * @brief PEEC node (Star element)
 *
 * Represents a charge node for capacitive effects.
 */
struct PEECNode {
    TVector3d position;     // Node position [m]
    double area;            // Associated area [m^2] (for self-potential)

    PEECNode() : area(0) {
        position = TVector3d(0, 0, 0);
    }

    PEECNode(const TVector3d& p, double a) : position(p), area(a) {}
};

/**
 * @brief PEEC matrix set
 *
 * Contains all matrices for Loop-Star PEEC analysis.
 */
struct PEECMatrices {
    std::vector<double> L;      // Inductance matrix [H] (n_loop x n_loop, row-major)
    std::vector<double> P;      // Potential coefficient [1/F] (n_star x n_star, row-major)
    std::vector<double> R;      // Resistance [Ohm] (n_loop, diagonal only)
    std::vector<double> M_LS;   // Loop-Star coupling (n_loop x n_star, row-major)
    int n_loop;                 // Number of loop elements
    int n_star;                 // Number of star elements

    PEECMatrices() : n_loop(0), n_star(0) {}

    // Access L(i,j)
    double& L_at(int i, int j) { return L[i * n_loop + j]; }
    double L_at(int i, int j) const { return L[i * n_loop + j]; }

    // Access P(i,j)
    double& P_at(int i, int j) { return P[i * n_star + j]; }
    double P_at(int i, int j) const { return P[i * n_star + j]; }

    // Access M_LS(i,j)
    double& M_LS_at(int i, int j) { return M_LS[i * n_star + j]; }
    double M_LS_at(int i, int j) const { return M_LS[i * n_star + j]; }
};

/**
 * @brief PEEC Matrix Builder
 *
 * Constructs L, P, R, M_LS matrices from segments and nodes.
 */
class PEECMatrixBuilder {
public:
    PEECMatrixBuilder();
    ~PEECMatrixBuilder();

    // ========== Input geometry ==========

    /**
     * @brief Add conductor segment
     */
    void AddSegment(const PEECSegment& segment);

    /**
     * @brief Add multiple segments
     */
    void AddSegments(const std::vector<PEECSegment>& segments);

    /**
     * @brief Add charge node
     */
    void AddNode(const PEECNode& node);

    /**
     * @brief Add multiple nodes
     */
    void AddNodes(const std::vector<PEECNode>& nodes);

    /**
     * @brief Auto-generate nodes from segment endpoints
     *
     * Creates nodes at segment start/end points, merging coincident points.
     */
    void AutoGenerateNodes();

    /**
     * @brief Clear all geometry
     */
    void Clear();

    // ========== Matrix computation ==========

    /**
     * @brief Build all PEEC matrices
     * @param includeStar If true, compute P and M_LS matrices
     * @return PEECMatrices structure
     */
    PEECMatrices Build(bool includeStar = true);

    /**
     * @brief Compute L matrix only (inductance)
     */
    void ComputeL(PEECMatrices& matrices);

    /**
     * @brief Compute P matrix only (potential coefficient)
     */
    void ComputeP(PEECMatrices& matrices);

    /**
     * @brief Compute R matrix only (resistance)
     */
    void ComputeR(PEECMatrices& matrices);

    /**
     * @brief Compute M_LS matrix (Loop-Star coupling)
     */
    void ComputeM_LS(PEECMatrices& matrices);

    // ========== Query ==========

    int NumSegments() const { return static_cast<int>(segments_.size()); }
    int NumNodes() const { return static_cast<int>(nodes_.size()); }

    const std::vector<PEECSegment>& GetSegments() const { return segments_; }
    const std::vector<PEECNode>& GetNodes() const { return nodes_; }

private:
    std::vector<PEECSegment> segments_;
    std::vector<PEECNode> nodes_;

    // Compute self-inductance using GMD approximation
    double SelfInductance(const PEECSegment& seg) const;

    // Compute mutual inductance via Neumann formula
    double MutualInductance(const PEECSegment& seg_i, const PEECSegment& seg_j) const;

    // Compute self-potential coefficient (disk approximation)
    double SelfPotential(const PEECNode& node) const;

    // Compute mutual potential coefficient
    double MutualPotential(const PEECNode& node_i, const PEECNode& node_j) const;
};

/**
 * @brief PEEC Solver
 *
 * Solves the Loop-Star PEEC system at a given frequency.
 */
class PEECSolver {
public:
    PEECSolver();
    ~PEECSolver();

    /**
     * @brief Set PEEC matrices
     */
    void SetMatrices(const PEECMatrices& matrices);

    /**
     * @brief Set frequency
     */
    void SetFrequency(double freq_hz);

    /**
     * @brief Set surface impedance (skin effect)
     * @param Zs_diag Diagonal elements of surface impedance matrix
     */
    void SetSurfaceImpedance(const std::vector<std::complex<double>>& Zs_diag);

    /**
     * @brief Build impedance matrix at current frequency
     * @param Z Output impedance matrix (row-major, size = (n_loop+n_star)^2)
     */
    void BuildImpedanceMatrix(std::vector<std::complex<double>>& Z);

    /**
     * @brief Compute port impedance
     * @param portVector Port excitation vector (size = n_loop)
     * @return Port impedance [Ohm]
     */
    std::complex<double> ComputePortImpedance(const std::vector<double>& portVector);

    /**
     * @brief Solve for currents given voltage
     * @param V Voltage vector (size = n_loop + n_star)
     * @param I Output current vector (size = n_loop + n_star)
     */
    void Solve(const std::vector<std::complex<double>>& V,
               std::vector<std::complex<double>>& I);

private:
    PEECMatrices matrices_;
    double frequency_;
    double omega_;
    std::vector<std::complex<double>> Zs_;  // Surface impedance (diagonal)
    bool hasSurfaceImpedance_;
};

// ========== Utility functions ==========

/**
 * @brief Create wire segments along a straight path
 */
std::vector<PEECSegment> CreateWireSegments(
    const TVector3d& start,
    const TVector3d& end,
    double width,
    double height,
    int n_segments,
    double sigma = 5.8e7);

/**
 * @brief Create circular loop segments
 */
std::vector<PEECSegment> CreateLoopSegments(
    const TVector3d& center,
    double radius,
    const TVector3d& normal,
    double width,
    double height,
    int n_segments,
    double sigma = 5.8e7);

} // namespace radia

#endif // RAD_PEEC_MATRICES_H
