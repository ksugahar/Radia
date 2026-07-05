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
 *   Z_LL = R_dc + jw*L + Z_s  (Z_s computed in Python: Bessel/Dowell/ESIM)
 *   Z_SS = P / jw
 *   Z_LS = jw*M_LS
 *   Z_SL = Z_LS^T (reciprocity)
 *
 * Note: build() returns frequency-independent matrices (L, R_dc, P, M_LS).
 *       Frequency-dependent surface impedance Z_s is computed in Python
 *       using scipy.special.jv (Bessel) or other cross-section models.
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
 * @brief Cross-section type for conductor segments
 */
enum class CrossSectionType {
    RECTANGULAR = 0,  // Rectangular cross-section (default)
    CIRCULAR = 1      // Circular cross-section
};

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
    CrossSectionType cross_section_type;  // Cross-section shape

    // Connectivity (topology mode)
    int node_from;          // Source node ID (-1 = legacy mode, no topology)
    int node_to;            // Destination node ID (-1 = legacy mode)

    // Multi-filament subdivision (nwinc/nhinc)
    int nwinc;              // Width subdivisions (default: 1 = no subdivision)
    int nhinc;              // Height subdivisions (default: 1 = no subdivision)
    int parent_segment;     // Parent segment index (-1 if not a sub-filament)

    PEECSegment()
        : length(0), width(0), height(0), sigma(5.8e7),
          cross_section_type(CrossSectionType::RECTANGULAR),
          node_from(-1), node_to(-1),
          nwinc(1), nhinc(1), parent_segment(-1) {
        center = TVector3d(0, 0, 0);
        direction = TVector3d(1, 0, 0);
    }

    PEECSegment(const TVector3d& c, const TVector3d& d, double l,
                double w, double h, double s = 5.8e7,
                CrossSectionType type = CrossSectionType::RECTANGULAR)
        : center(c), direction(d), length(l), width(w), height(h), sigma(s),
          cross_section_type(type), node_from(-1), node_to(-1),
          nwinc(1), nhinc(1), parent_segment(-1) {}

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
 * @brief PEEC Panel (Star element - 2D surface)
 *
 * Represents a 2D surface panel for accurate capacitive effects.
 * Uses analytical edge integration (FastImp approach).
 */
struct PEECPanel {
    std::vector<TVector3d> vertices;  // 3 (triangle) or 4 (quad) vertices
    TVector3d center;                 // Panel centroid
    TVector3d normal;                 // Outward normal vector
    double area;                      // Panel area [m^2]

    enum Type { Triangle, Quadrilateral } type;

    PEECPanel() : area(0), type(Triangle) {
        center = TVector3d(0, 0, 0);
        normal = TVector3d(0, 0, 1);
    }

    PEECPanel(const std::vector<TVector3d>& verts)
        : vertices(verts) {
        type = (verts.size() == 3) ? Triangle : Quadrilateral;
        ComputeGeometry();
    }

    // Compute centroid, normal, and area from vertices
    void ComputeGeometry();
};

/**
 * @brief PEEC port definition
 *
 * Defines a port between two nodes for impedance extraction.
 */
struct PEECPort {
    int node_positive;      // Positive terminal node ID
    int node_negative;      // Negative terminal node ID
    int port_id;            // Port identifier

    PEECPort() : node_positive(-1), node_negative(-1), port_id(0) {}
    PEECPort(int pos, int neg, int id) : node_positive(pos), node_negative(neg), port_id(id) {}
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

    // Incidence matrix A (CSR format): n_junction x n_filament
    // A[m,n] = +1 if filament n leaves junction m
    // A[m,n] = -1 if filament n enters junction m
    std::vector<int> incidence_indptr;    // Row pointers (size: n_junction+1)
    std::vector<int> incidence_indices;   // Column indices (filament index)
    std::vector<double> incidence_data;   // Values (+1 or -1)
    int n_junction;                       // Number of internal junction nodes

    // Port definitions
    std::vector<PEECPort> ports;

    // Gathered capacitance for proper MNA (n_nodes x n_nodes, row-major, real)
    // C_gathered = G × (P/eps_eff)⁻¹ × G^T
    // Used in MNA: Y_node = A × Y_branch⁻¹ × A^T + jω × C_gathered
    std::vector<double> C_gathered;
    int n_nodes_gathered;               // Number of nodes in C_gathered

    PEECMatrices() : n_loop(0), n_star(0), n_junction(0), n_nodes_gathered(0) {}

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
     * @brief Add 2D surface panel
     */
    void AddPanel(const PEECPanel& panel);

    /**
     * @brief Add multiple panels
     */
    void AddPanels(const std::vector<PEECPanel>& panels);

    /**
     * @brief Auto-generate nodes from segment endpoints
     *
     * Creates nodes at segment start/end points, merging coincident points.
     */
    void AutoGenerateNodes();

    /**
     * @brief Auto-generate face panels from existing segments
     *
     * Creates quadrilateral panels on segment surfaces for accurate
     * capacitance extraction (panel-based P matrix).
     *
     * @param mode 0=all 6 faces, 1=top+bottom (2 per seg), 2=4 long faces (no end caps)
     * @param eps_r Relative permittivity for dielectric substrate (default: 1.0 = vacuum)
     */
    void GenerateFacePanels(int mode = 2, double eps_r = 1.0);

    /**
     * @brief Clear all geometry
     */
    void Clear();

    // ========== Topology-aware input (node-segment model) ==========

    /**
     * @brief Add a node at a given position
     * @param position Node position [m]
     * @param area Associated area [m^2] (for capacitive effects)
     * @return Node ID (auto-incrementing)
     */
    int AddNodeAt(const TVector3d& position, double area = 0.0);

    /**
     * @brief Add a segment connecting two nodes
     *
     * Computes center, direction, and length from node positions.
     * If nwinc*nhinc > 1, sub-filaments are created during Build().
     *
     * @param node_from Source node ID
     * @param node_to Destination node ID
     * @param width Cross-section width [m]
     * @param height Cross-section height [m]
     * @param sigma Conductivity [S/m]
     * @param type Cross-section type
     * @param nwinc Width subdivisions (default: 1)
     * @param nhinc Height subdivisions (default: 1)
     */
    void AddConnectedSegment(int node_from, int node_to,
                             double width, double height,
                             double sigma = 5.8e7,
                             CrossSectionType type = CrossSectionType::RECTANGULAR,
                             int nwinc = 1, int nhinc = 1);

    /**
     * @brief Add a port between two nodes
     * @param node_positive Positive terminal node ID
     * @param node_negative Negative terminal node ID
     * @return Port ID (auto-incrementing)
     */
    int AddPort(int node_positive, int node_negative);

    /**
     * @brief Expand multi-filament segments into sub-filaments
     *
     * For segments with nwinc*nhinc > 1, creates sub-filaments by
     * subdividing the cross-section. Sub-filaments share the same
     * node_from/node_to (parallel connection). The original parent
     * segment is replaced by its sub-filaments.
     */
    void ExpandFilaments();

    /**
     * @brief Build incidence matrix from segment connectivity
     *
     * Identifies junction nodes (internal nodes not used as port terminals
     * with 2+ connections) and builds the CSR incidence matrix.
     */
    void BuildIncidenceMatrix(PEECMatrices& matrices);

    // ========== Matrix computation ==========

    /**
     * @brief Build all PEEC matrices (frequency-independent: L, R_dc, P, M_LS)
     * @param includeStar If true, compute P and M_LS matrices
     * @param eps_eff Effective permittivity for gathered capacitance (default 1.0).
     *               If > 0 and panels present, also computes C_gathered.
     *               For PCB half-space model: eps_eff = (1 + eps_r) / 2.
     * @return PEECMatrices structure
     */
    PEECMatrices Build(bool includeStar = true, double eps_eff = 1.0);

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

    /**
     * @brief Compute gathered capacitance C_eff = G × (P/eps_eff)⁻¹ × G^T
     *
     * G is the gathering matrix (n_nodes × n_star) that maps panel charges to nodes.
     * Each panel's charge is split 50/50 between its parent segment's endpoint nodes.
     *
     * @param matrices Must have P already computed. C_gathered is written here.
     * @param eps_eff Effective permittivity for half-space model (default 1.0).
     *               For PCB: eps_eff = (1 + eps_r) / 2
     */
    void ComputeGatheredCapacitance(PEECMatrices& matrices, double eps_eff = 1.0);

    // ========== Query ==========

    int NumSegments() const { return static_cast<int>(segments_.size()); }
    int NumNodes() const { return static_cast<int>(nodes_.size()); }
    int NumPanels() const { return static_cast<int>(panels_.size()); }

    const std::vector<PEECSegment>& GetSegments() const { return segments_; }
    const std::vector<PEECNode>& GetNodes() const { return nodes_; }
    const std::vector<PEECPanel>& GetPanels() const { return panels_; }

    /**
     * @brief Return L(i, j) entry of the inductance matrix.
     *
     * Dispatches to SelfInductance(i) when i == j and MutualInductance(i, j)
     * otherwise. Exposed publicly for on-demand consumption by the HACApK
     * PEEC adapter (see RadHACApKPEECManager::GetInteractionMatrixElement)
     * without forcing O(N^2) assembly of the full L matrix.
     */
    double InductanceEntry(int i, int j) const;

private:
    std::vector<PEECSegment> segments_;
    std::vector<PEECNode> nodes_;
    std::vector<PEECPanel> panels_;
    std::vector<int> panel_segment_ids_;  // panel index -> parent segment index
    std::vector<PEECPort> ports_;
    int nextPortId_;
    double eps_r_;  // Relative permittivity for P matrix scaling (default: 1.0)

    // Compute self-inductance (dispatches to rectangular or circular)
    double SelfInductance(const PEECSegment& seg) const;

    // Compute self-inductance for rectangular cross-section (Grover formula)
    double SelfInductanceRectangular(const PEECSegment& seg) const;

    // Compute self-inductance for circular cross-section (Grover formula)
    double SelfInductanceCircular(const PEECSegment& seg) const;

    // Compute mutual inductance via Neumann formula
    double MutualInductance(const PEECSegment& seg_i, const PEECSegment& seg_j) const;

    /** Exact mutual inductance between two parallel rectangular bars.
     *
     * Uses Gauss quadrature over both cross-sections (n_gauss x n_gauss
     * points per bar, total n_gauss^4 filament-pair evaluations) to
     * compute the volume-volume integral that the filamentary Neumann
     * formula approximates as a line-line integral.
     *
     * Reference: A. E. Ruehli, IBM J. Res. Dev. 16(5), 1972, Section 6
     *            (filament decomposition of the mutual inductance integral).
     *
     * @param seg_i First segment (parallel to seg_j)
     * @param seg_j Second segment
     * @param n_gauss Gauss points per dimension per cross-section (default 3)
     */
    double MutualInductanceRectBar(const PEECSegment& seg_i,
                                    const PEECSegment& seg_j,
                                    int n_gauss = 3) const;

    // Compute self-potential coefficient (disk approximation)
    double SelfPotential(const PEECNode& node) const;

    // Compute mutual potential coefficient
    double MutualPotential(const PEECNode& node_i, const PEECNode& node_j) const;

    // Analytical panel integration (Hess-Smith / Wilton)
    double SelfPotentialPanelTriangle(const PEECPanel& panel) const;
    double SelfPotentialPanelQuad(const PEECPanel& panel) const;
    double MutualPotentialPanelTriangle(const PEECPanel& panel_i, const PEECPanel& panel_j) const;

    // Generalized mutual potential between any panel types (triangle or quad)
    double MutualPotentialPanel(const PEECPanel& panel_i, const PEECPanel& panel_j) const;

    // Hess-Smith analytical single-layer potential from a flat polygon
    // evaluated at a given observation point (Hess & Smith, 1967; Newman, 1986)
    double HessSmithPotential(const std::vector<TVector3d>& vertices,
                              const TVector3d& normal,
                              const TVector3d& obs_point) const;

    // Recursive fourfil mutual inductance (Romberg subdivision)
    double MutualInductanceFourfil(const PEECSegment& seg_i,
                                   const PEECSegment& seg_j,
                                   int depth) const;
};

/**
 * @brief PEEC Solver
 *
 * Solves the Loop-Star PEEC system at a given frequency.
 *
 * Supports two modes:
 * 1. Legacy mode: Direct impedance solve (Solve, ComputePortImpedance)
 * 2. MNA mode: Multi-port Z-parameter extraction via Modified Nodal Analysis
 *
 * MNA solver methods:
 *   0 = LU (LAPACK zgesv_/zgetrf_/zgetrs_) - default, O(n^3)
 *   1 = BiCGSTAB (rad_bicgstab.h templated solver) - iterative, O(k*n^2)
 *
 * Both share MKL BLAS/LAPACK infrastructure with Radia solver paths.
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
     * @brief Set segment connectivity for MNA topology
     * @param seg_nodes Vector of (node_from, node_to) pairs per filament
     * @param n_nodes Total number of nodes
     */
    void SetSegmentNodes(const std::vector<std::pair<int,int>>& seg_nodes, int n_nodes);

    /**
     * @brief Set port definitions for multi-port extraction
     * @param ports Vector of PEECPort
     */
    void SetPorts(const std::vector<PEECPort>& ports);

    // ========== MNA Multi-Port Solver (uses LAPACK zgetrf_/zgetrs_) ==========

    /**
     * @brief Compute multi-port Z-parameter matrix at a single frequency
     *
     * Z[i,j] = V_port_i when I_port_j = 1A (all other ports open-circuited).
     *
     * @param freq Frequency in Hz
     * @param Zs Per-filament surface impedance (n_loop complex values), or nullptr
     * @param n_Zs Number of Zs elements
     * @param Z_out Output Z matrix (n_ports x n_ports, row-major complex)
     * @param n_ports_out Number of ports (for validation)
     */
    void ComputeZMatrix(double freq,
                        const std::complex<double>* Zs, int n_Zs,
                        std::complex<double>* Z_out, int n_ports_out);

    /**
     * @brief Compute coupling coefficient from Z-parameters
     *
     * k_ij = L_ij / sqrt(L_ii * L_jj) where L_ij = Im(Z_ij) / omega
     *
     * @param freq Frequency in Hz
     * @param Zs Per-filament surface impedance, or nullptr
     * @param n_Zs Number of Zs elements
     * @param k_out Coupling coefficient matrix (n_ports x n_ports, row-major)
     * @param L_out Inductance matrix (n_ports x n_ports, row-major) [H]
     * @param n_ports_out Number of ports
     */
    void ComputeCouplingCoefficient(double freq,
                                    const std::complex<double>* Zs, int n_Zs,
                                    double* k_out, double* L_out,
                                    int n_ports_out);

    /**
     * @brief Batch frequency sweep for multi-port Z-parameters
     *
     * @param freqs Frequency array [Hz]
     * @param n_freq Number of frequencies
     * @param Zs_all Per-frequency surface impedance (n_freq x n_loop), or nullptr
     * @param Z_out Output Z matrices (n_freq x n_ports x n_ports, row-major)
     * @param n_ports_out Number of ports
     */
    void FrequencySweep(const double* freqs, int n_freq,
                        const std::complex<double>* Zs_all,
                        std::complex<double>* Z_out, int n_ports_out);

    /**
     * @brief Compute per-filament complex branch currents under user-specified
     *        port current injection.
     *
     * Solves MNA with current-source port excitation (not voltage-source)
     * and returns I_branch = Y_branch * A^T * V_node (n_loop complex entries,
     * signed so that positive = current flowing from node_from to node_to).
     *
     * Use case: "All filaments in parallel between port+ and port-" -> set
     * port_currents[0] = I_total, then the returned vector is the per-
     * filament AC current distribution including skin/proximity coupling.
     *
     * @param freq Frequency in Hz
     * @param port_currents Complex current injected at each port (size n_ports)
     * @param Zs Optional per-filament surface impedance (n_loop), or nullptr
     * @param n_Zs Number of Zs elements
     * @param I_branch_out Output branch currents (n_loop complex, row-major)
     * @param n_loop_out Size of I_branch_out (for validation)
     */
    void ComputeBranchCurrents(double freq,
                                const std::complex<double>* port_currents,
                                int n_ports_in,
                                const std::complex<double>* Zs, int n_Zs,
                                std::complex<double>* I_branch_out,
                                int n_loop_out);

    // ========== Legacy API (upgraded to LAPACK zgesv_) ==========

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

    // ========== Solver method selection ==========

    /**
     * @brief Set solver method for MNA Z_branch inversion
     * @param method 0 = LU (default), 1 = BiCGSTAB
     */
    void SetSolverMethod(int method) { solverMethod_ = method; }
    int GetSolverMethod() const { return solverMethod_; }

    /**
     * @brief Set BiCGSTAB parameters
     * @param tol Convergence tolerance (default 1e-10)
     * @param max_iter Maximum iterations (default 1000)
     */
    void SetBiCGSTABParams(double tol, int max_iter) {
        bicgstab_tol_ = tol;
        bicgstab_max_iter_ = max_iter;
    }

    // ========== Gathered Capacitance (proper MNA) ==========

    /**
     * @brief Set pre-computed gathered capacitance for proper MNA
     *
     * When set, MNA uses: Y_node = A * Y_branch^{-1} * A^T + jw * C_gathered
     * instead of the Schur complement (which cannot represent resonance).
     *
     * @param C_gath Gathered capacitance matrix (n_nodes x n_nodes, row-major, real)
     * @param n_nodes Number of nodes
     */
    void SetGatheredCapacitance(const std::vector<double>& C_gath, int n_nodes);
    bool HasGatheredCapacitance() const { return hasGatheredCapacitance_; }

    // ========== Query ==========

    int NumPorts() const { return static_cast<int>(ports_.size()); }
    bool HasTopology() const { return hasTopology_; }

private:
    PEECMatrices matrices_;
    double frequency_;
    double omega_;
    std::vector<std::complex<double>> Zs_;
    bool hasSurfaceImpedance_;

    // MNA topology data
    std::vector<std::pair<int,int>> segment_nodes_;
    int n_nodes_;
    std::vector<PEECPort> ports_;
    bool hasTopology_;

    // Gathered capacitance for proper MNA
    bool hasGatheredCapacitance_;
    std::vector<double> C_gathered_;    // n_nodes x n_nodes, real, row-major

    // Solver method: 0 = LU, 1 = BiCGSTAB
    int solverMethod_;
    double bicgstab_tol_;
    int bicgstab_max_iter_;

    // ========== MNA internal methods ==========

    /**
     * @brief Build full incidence matrix A_full (n_nodes x n_loop, dense, real)
     *
     * A_full[node, fil] = +1 if filament leaves node (node_from)
     * A_full[node, fil] = -1 if filament enters node (node_to)
     */
    void BuildFullIncidenceMatrix(std::vector<double>& A_full);

    /**
     * @brief Find connected components of the node graph via BFS
     * @param components Output: list of node sets (one per component)
     */
    void FindConnectedComponents(std::vector<std::vector<int>>& components);

    /**
     * @brief Select ground nodes (one per connected component)
     *
     * Prefers negative terminal of a port in each component.
     */
    void SelectGroundNodes(const std::vector<std::vector<int>>& components,
                           std::vector<int>& ground_nodes);

    /**
     * @brief Build Z_branch = diag(R + Zs) + jw*L
     */
    void BuildZBranch(double freq,
                      const std::complex<double>* Zs, int n_Zs,
                      std::vector<std::complex<double>>& Z_branch);

    /**
     * @brief Build Z_eff = Z_branch - Z_LS * Z_SS^{-1} * Z_SL (Schur complement)
     *
     * When panels are present, eliminates Star unknowns.
     * Falls back to Z_branch if no panels.
     */
    void BuildZEff(double freq,
                   const std::complex<double>* Zs, int n_Zs,
                   std::vector<std::complex<double>>& Z_eff);

    /**
     * @brief Full MNA multi-port solve
     *
     * 1. Invert Z_eff -> Y_branch (via zgesv_)
     * 2. Y_node = A * Y_branch * A^T (via cblas_zgemm)
     * 2b. If gathered capacitance: Y_node += jw * C_gathered
     * 3. Ground selection (one per connected component)
     * 4. LU factorize Y_reduced (zgetrf_)
     * 5. Solve for each port excitation (zgetrs_)
     * 6. Extract Z_mat from voltage solutions
     */
    void MNASolveMultiPort(const std::vector<std::complex<double>>& Z_eff,
                           double omega,
                           std::complex<double>* Z_out, int n_ports_out);
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
