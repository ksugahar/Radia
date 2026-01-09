/*
 * rad_rwg_basis.h
 *
 * RWG (Rao-Wilton-Glisson) Basis Functions for EFIE
 *
 * RWG basis functions are vector basis functions defined on edges of
 * triangular meshes. They are divergence-conforming and well-suited
 * for the Electric Field Integral Equation (EFIE).
 *
 * Key features:
 * - Edge-based vector basis functions
 * - Normal component continuity across edges
 * - Proper handling of surface charge via div(J)
 *
 * References:
 * [1] S.M. Rao, D.R. Wilton, A.W. Glisson, "Electromagnetic scattering
 *     by surfaces of arbitrary shape", IEEE TAP, 1982
 * [2] D.R. Wilton et al., "Potential integrals for uniform and linear
 *     source distributions on polygonal and polyhedral domains",
 *     IEEE TAP, 1984
 *
 * Part of Radia project
 */

#ifndef RAD_RWG_BASIS_H
#define RAD_RWG_BASIS_H

#include "gmvect.h"
#include "rad_constants.h"
#include <vector>
#include <array>
#include <complex>
#include <memory>
#include <unordered_map>

namespace radia {

// Physical constants
constexpr double RWG_MU_0 = RadConst::MU_0;
constexpr double RWG_EPS_0 = 8.854187817e-12;
constexpr double RWG_INV_FOUR_PI = RadConst::INV_FOUR_PI;

/**
 * @brief Triangle element for RWG mesh
 */
struct RWGTriangle {
    int vertices[3];          // Global vertex indices
    TVector3d center;         // Triangle centroid
    TVector3d normal;         // Unit outward normal
    double area;              // Triangle area

    // Edge information (local to triangle)
    // edge i is opposite to vertex i
    int edges[3];             // Global edge indices (-1 if boundary)
    int edgeSign[3];          // +1 if edge direction matches, -1 otherwise

    // Precomputed vectors for RWG evaluation
    TVector3d rhoVec[3];      // rho_i = r - r_i (vertex to point vectors)

    RWGTriangle() : area(0) {
        vertices[0] = vertices[1] = vertices[2] = -1;
        edges[0] = edges[1] = edges[2] = -1;
        edgeSign[0] = edgeSign[1] = edgeSign[2] = 1;
    }
};

/**
 * @brief Edge element for RWG basis
 *
 * Each internal edge supports one RWG basis function.
 * The basis function spans two adjacent triangles.
 */
struct RWGEdge {
    int vertex1, vertex2;     // Edge vertices (v1 < v2 for consistency)
    int triPlus, triMinus;    // Adjacent triangles (T+, T-)
    int localIdxPlus;         // Local edge index in T+ (0,1,2)
    int localIdxMinus;        // Local edge index in T-
    double length;            // Edge length
    TVector3d midpoint;       // Edge midpoint
    TVector3d direction;      // Unit vector from v1 to v2

    // Precomputed for integration
    TVector3d rhoPlusC;       // rho_c^+ = center of T+ to opposite vertex
    TVector3d rhoMinusC;      // rho_c^- = center of T- to opposite vertex

    bool isBoundary;          // True if on mesh boundary

    RWGEdge() : vertex1(-1), vertex2(-1), triPlus(-1), triMinus(-1),
                localIdxPlus(-1), localIdxMinus(-1), length(0), isBoundary(true) {}
};

/**
 * @brief RWG Mesh for EFIE discretization
 *
 * Manages triangular mesh with RWG basis function support.
 */
class RWGMesh {
public:
    RWGMesh();
    ~RWGMesh();

    /**
     * @brief Create mesh from vertex and triangle lists
     * @param vertices List of vertex coordinates
     * @param triangles List of triangle vertex indices (3 per triangle)
     */
    void CreateFromTriangles(const std::vector<TVector3d>& vertices,
                             const std::vector<std::array<int, 3>>& triangles);

    /**
     * @brief Create mesh by triangulating quadrilateral panels
     * Each quad is split into 2 triangles
     */
    void CreateFromQuadPanels(const std::vector<TVector3d>& vertices,
                              const std::vector<std::array<int, 4>>& quads);

    /**
     * @brief Create mesh for wire conductor (circular/rectangular cross-section)
     */
    void CreateWireMesh(const std::vector<TVector3d>& path,
                        double wireRadius,
                        int numAroundWire = 8,
                        int numAlongSegment = 1);

    /**
     * @brief Create rectangular plate mesh (open surface for workpiece)
     * @param center Plate center
     * @param size Plate size [Lx, Ly]
     * @param normal Plate normal direction
     * @param nx Number of divisions in x direction
     * @param ny Number of divisions in y direction
     */
    void CreateRectangularPlate(const TVector3d& center,
                                double Lx, double Ly,
                                const TVector3d& normal,
                                int nx = 10, int ny = 10);

    /**
     * @brief Create circular disk mesh (open surface)
     * @param center Disk center
     * @param radius Disk radius
     * @param normal Disk normal direction
     * @param nr Number of radial divisions
     * @param ntheta Number of angular divisions
     */
    void CreateCircularDisk(const TVector3d& center,
                            double radius,
                            const TVector3d& normal,
                            int nr = 5, int ntheta = 16);

    /**
     * @brief Create cylindrical shell mesh (open surface)
     * @param center Cylinder center
     * @param radius Cylinder radius
     * @param height Cylinder height
     * @param axis Cylinder axis direction
     * @param nz Number of divisions along axis
     * @param ntheta Number of angular divisions
     */
    void CreateCylindricalShell(const TVector3d& center,
                                double radius,
                                double height,
                                const TVector3d& axis,
                                int nz = 10, int ntheta = 16);

    /**
     * @brief Merge another mesh into this mesh
     * Used to combine coil and workpiece meshes
     * @param other Mesh to merge
     */
    void MergeMesh(const RWGMesh& other);

    /**
     * @brief Get boundary edges (edges with only one adjacent triangle)
     * @return Vector of boundary edge indices
     */
    std::vector<int> GetBoundaryEdges() const;

    /**
     * @brief Check if this mesh is a closed surface
     * @return True if no boundary edges exist
     */
    bool IsClosed() const { return numInteriorEdges_ == static_cast<int>(edges_.size()); }

    // Access
    int NumVertices() const { return static_cast<int>(vertices_.size()); }
    int NumTriangles() const { return static_cast<int>(triangles_.size()); }
    int NumEdges() const { return static_cast<int>(edges_.size()); }
    int NumInteriorEdges() const { return numInteriorEdges_; }

    const std::vector<TVector3d>& Vertices() const { return vertices_; }
    const std::vector<RWGTriangle>& Triangles() const { return triangles_; }
    const std::vector<RWGEdge>& Edges() const { return edges_; }

    const TVector3d& Vertex(int i) const { return vertices_[i]; }
    const RWGTriangle& Triangle(int i) const { return triangles_[i]; }
    const RWGEdge& Edge(int i) const { return edges_[i]; }

    /**
     * @brief Evaluate RWG basis function at a point
     *
     * The RWG basis function on edge n is:
     *   f_n(r) = (l_n / 2*A_+) * rho_n^+(r)  in T+
     *          = (l_n / 2*A_-) * rho_n^-(r)  in T-
     *          = 0  outside T+ and T-
     *
     * where:
     *   l_n = edge length
     *   A_+, A_- = areas of T+ and T-
     *   rho_n^+(r) = r - r_n^+  (vector from free vertex of T+ to point r)
     *   rho_n^-(r) = r_n^- - r  (vector from point r to free vertex of T-)
     *
     * @param edgeIdx Edge index (basis function index)
     * @param triIdx Triangle index where point lies
     * @param point Evaluation point
     * @return Vector value of basis function
     */
    TVector3d EvaluateRWG(int edgeIdx, int triIdx, const TVector3d& point) const;

    /**
     * @brief Evaluate divergence of RWG basis function
     *
     * div(f_n) = l_n / A_+  in T+
     *          = -l_n / A_-  in T-
     *
     * @param edgeIdx Edge index
     * @param triIdx Triangle index
     * @return Divergence value (scalar)
     */
    double EvaluateRWGDiv(int edgeIdx, int triIdx) const;

    /**
     * @brief Get free vertex of edge in given triangle
     * @param edgeIdx Edge index
     * @param triIdx Triangle index (must be triPlus or triMinus of edge)
     * @return Vertex position of free vertex (not on edge)
     */
    TVector3d GetFreeVertex(int edgeIdx, int triIdx) const;

private:
    std::vector<TVector3d> vertices_;
    std::vector<RWGTriangle> triangles_;
    std::vector<RWGEdge> edges_;
    int numInteriorEdges_;

    void BuildEdgeConnectivity();
    void ComputeTriangleProperties();
    void ComputeEdgeProperties();

    // Hash function for edge lookup
    static inline uint64_t EdgeHash(int v1, int v2) {
        if (v1 > v2) std::swap(v1, v2);
        return (static_cast<uint64_t>(v1) << 32) | static_cast<uint64_t>(v2);
    }
};

/**
 * @brief Analytical integration routines for RWG-EFIE
 *
 * Implements the formulas from Wilton et al. (1984) for
 * potential integrals over triangular domains.
 */
class RWGIntegration {
public:
    /**
     * @brief Compute potential integral over triangle
     *
     * I = integral_T { 1 / |r - r'| } dA'
     *
     * Uses analytical formula from Wilton et al.
     *
     * @param point Observation point r
     * @param v1, v2, v3 Triangle vertices
     * @return Potential integral value
     */
    static double PotentialIntegral(const TVector3d& point,
                                    const TVector3d& v1,
                                    const TVector3d& v2,
                                    const TVector3d& v3);

    /**
     * @brief Compute gradient of potential integral
     *
     * grad_I = integral_T { (r - r') / |r - r'|^3 } dA'
     *
     * @param point Observation point r
     * @param v1, v2, v3 Triangle vertices
     * @return Gradient vector
     */
    static TVector3d GradPotentialIntegral(const TVector3d& point,
                                           const TVector3d& v1,
                                           const TVector3d& v2,
                                           const TVector3d& v3);

    /**
     * @brief Compute vector potential integral for RWG
     *
     * A_mn = integral_Tm integral_Tn { f_m(r) . f_n(r') / |r-r'| } dA dA'
     *
     * This is the L matrix element for EFIE.
     *
     * @param mesh RWG mesh
     * @param edgeM Source edge index
     * @param edgeN Test edge index
     * @return Complex integral value
     */
    static std::complex<double> VectorPotentialIntegral(
        const RWGMesh& mesh, int edgeM, int edgeN);

    /**
     * @brief Compute scalar potential integral for RWG
     *
     * P_mn = integral_Tm integral_Tn { div(f_m) . div(f_n) / |r-r'| } dA dA'
     *
     * This is the P matrix element for EFIE.
     *
     * @param mesh RWG mesh
     * @param edgeM Source edge index
     * @param edgeN Test edge index
     * @return Complex integral value
     */
    static std::complex<double> ScalarPotentialIntegral(
        const RWGMesh& mesh, int edgeM, int edgeN);

    /**
     * @brief Self-term integral (m == n, same triangle)
     *
     * Uses analytical formula to handle 1/R singularity.
     *
     * @param mesh RWG mesh
     * @param edgeIdx Edge index
     * @param triIdx Triangle index
     * @param isVectorPotential True for A integral, false for Phi integral
     * @return Self-term integral value
     */
    static std::complex<double> SelfTermIntegral(
        const RWGMesh& mesh, int edgeIdx, int triIdx, bool isVectorPotential);

    // Gauss-Legendre quadrature for triangle (public for field computation)
    static constexpr int NQUAD_TRI = 7;  // 7-point rule
    static const double triQuadPts_[NQUAD_TRI][2];   // Barycentric coordinates
    static const double triQuadWts_[NQUAD_TRI];      // Weights

private:
    // Helper: signed solid angle
    static double SignedSolidAngle(const TVector3d& point,
                                   const TVector3d& v1,
                                   const TVector3d& v2,
                                   const TVector3d& v3);

    // Helper: line integral term for Wilton formula
    static double LineIntegralTerm(const TVector3d& point,
                                   const TVector3d& v1,
                                   const TVector3d& v2,
                                   double& R_plus, double& R_minus);
};

/**
 * @brief Loop-Star Decomposition for low-frequency stability
 *
 * The Loop-Star decomposition transforms the RWG basis into a combination
 * of Loop (solenoidal/curl) and Star (irrotational/divergence) functions.
 *
 * Benefits:
 * - Eliminates low-frequency breakdown in EFIE
 * - Separates inductive (Loop) and capacitive (Star) effects
 * - For MQS problems (div J = 0), only Loop equations are needed
 * - Reduces DOF and improves convergence
 *
 * Loop functions: Divergence-free, represent inductive effects
 * Star functions: Curl-free, represent capacitive effects
 *
 * Reference:
 *   G. Vecchi, "Loop-Star Decomposition of Basis Functions in the
 *   Discretization of the EFIE", IEEE TAP, 1999
 */
class LoopStarDecomposition {
public:
    LoopStarDecomposition();
    ~LoopStarDecomposition();

    /**
     * @brief Build Loop-Star decomposition from RWG mesh
     * @param mesh RWG mesh
     */
    void Build(const RWGMesh& mesh);

    /**
     * @brief Get number of Loop (solenoidal) basis functions
     * For closed surfaces: N_loop = N_edges - N_vertices + 1 (Euler)
     * For open surfaces: N_loop = N_edges - N_vertices + 1 - N_boundary_loops
     */
    int NumLoops() const { return numLoops_; }

    /**
     * @brief Get number of Star (irrotational) basis functions
     * N_star = N_vertices - 1 (one vertex grounded)
     */
    int NumStars() const { return numStars_; }

    /**
     * @brief Get Loop transformation matrix T_loop
     * I_rwg = T_loop * I_loop + T_star * I_star
     * T_loop: N_edges x N_loop
     */
    const std::vector<double>& LoopMatrix() const { return T_loop_; }

    /**
     * @brief Get Star transformation matrix T_star
     * T_star: N_edges x N_star
     */
    const std::vector<double>& StarMatrix() const { return T_star_; }

    /**
     * @brief Transform matrix from RWG to Loop-Star basis
     * Z_LS = T^T * Z_RWG * T
     *
     * @param Z_RWG Input matrix in RWG basis (N_edges x N_edges)
     * @param Z_LL Output Loop-Loop block (N_loop x N_loop)
     * @param Z_LS Output Loop-Star block (N_loop x N_star)
     * @param Z_SL Output Star-Loop block (N_star x N_loop)
     * @param Z_SS Output Star-Star block (N_star x N_star)
     */
    void TransformMatrix(const std::vector<std::complex<double>>& Z_RWG,
                         std::vector<std::complex<double>>& Z_LL,
                         std::vector<std::complex<double>>& Z_LS,
                         std::vector<std::complex<double>>& Z_SL,
                         std::vector<std::complex<double>>& Z_SS) const;

    /**
     * @brief Transform RHS vector from RWG to Loop-Star basis
     * V_L = T_loop^T * V_RWG
     * V_S = T_star^T * V_RWG
     */
    void TransformRHS(const std::vector<std::complex<double>>& V_RWG,
                      std::vector<std::complex<double>>& V_L,
                      std::vector<std::complex<double>>& V_S) const;

    /**
     * @brief Transform solution from Loop-Star to RWG basis
     * I_RWG = T_loop * I_L + T_star * I_S
     */
    void TransformSolution(const std::vector<std::complex<double>>& I_L,
                           const std::vector<std::complex<double>>& I_S,
                           std::vector<std::complex<double>>& I_RWG) const;

    /**
     * @brief Check if decomposition is valid
     */
    bool IsValid() const { return isValid_; }

private:
    int numEdges_;
    int numVertices_;
    int numLoops_;
    int numStars_;
    bool isValid_;

    // Transformation matrices (row-major storage)
    // T_loop: numEdges_ x numLoops_
    // T_star: numEdges_ x numStars_
    std::vector<double> T_loop_;
    std::vector<double> T_star_;

    // Build spanning tree for Loop-Star construction
    void BuildSpanningTree(const RWGMesh& mesh,
                           std::vector<int>& parent,
                           std::vector<int>& treeEdges,
                           std::vector<int>& coTreeEdges);
};

/**
 * @brief EFIE Solver using RWG basis functions
 */
class RWGEFIESolver {
public:
    RWGEFIESolver();
    ~RWGEFIESolver();

    /**
     * @brief Set mesh
     */
    void SetMesh(std::shared_ptr<RWGMesh> mesh) { mesh_ = mesh; }

    /**
     * @brief Set material properties
     */
    void SetConductivity(double sigma) { conductivity_ = sigma; }
    void SetFrequency(double freq);

    /**
     * @brief Enable/disable MQS (Magneto-Quasi-Static) mode
     * In MQS mode, only Loop equations are solved (div J = 0).
     * MQS is ENABLED by default for power electronics applications.
     * @param enable True for MQS (default), false for full EFIE
     *
     * NOTE: Loop-Star decomposition is ALWAYS used (mandatory).
     * The legacy direct RWG solver has been removed.
     */
    void SetUseMQS(bool enable) { useMQS_ = enable; }
    bool GetUseMQS() const { return useMQS_; }

    /**
     * @brief Define port for impedance calculation
     * @param edgeIndices Edges that form the port gap
     */
    void DefinePort(const std::vector<int>& edgeIndices);

    /**
     * @brief Set voltage excitation at port
     */
    void SetVoltageExcitation(std::complex<double> V);

    /**
     * @brief Assemble and solve EFIE system
     *
     * EFIE: (jωL + R + P/(jω)) J = V
     *
     * where:
     *   L_mn = μ₀/(4π) ∫∫ f_m · f_n / R dA dA'  (vector potential)
     *   P_mn = 1/(4πε₀) ∫∫ div(f_m) div(f_n) / R dA dA'  (scalar potential)
     *   R = surface resistance (from skin effect)
     */
    void Solve();

    /**
     * @brief Get port impedance after solve
     */
    std::complex<double> GetImpedance() const { return impedance_; }

    /**
     * @brief Get total current through port
     */
    std::complex<double> GetCurrent() const { return totalCurrent_; }

    /**
     * @brief Get solution coefficients (current on each edge)
     */
    const std::vector<std::complex<double>>& GetSolution() const { return solution_; }

    /**
     * @brief Compute B field at a point from solved currents
     */
    void ComputeB(const TVector3d& point,
                  std::complex<double>& Bx,
                  std::complex<double>& By,
                  std::complex<double>& Bz) const;

private:
    std::shared_ptr<RWGMesh> mesh_;

    double conductivity_;
    double frequency_;
    double omega_;
    std::complex<double> jOmega_;

    // Port definition
    std::vector<int> portEdges_;
    std::complex<double> voltageExcitation_;

    // System matrix and solution
    std::vector<std::complex<double>> L_matrix_;  // Vector potential matrix
    std::vector<std::complex<double>> P_matrix_;  // Scalar potential matrix
    std::vector<std::complex<double>> R_matrix_;  // Resistance matrix
    std::vector<std::complex<double>> Z_matrix_;  // Combined: jωL + R + P/(jω)
    std::vector<std::complex<double>> rhs_;
    std::vector<std::complex<double>> solution_;

    // Results
    std::complex<double> impedance_;
    std::complex<double> totalCurrent_;

    // Loop-Star decomposition (always enabled for low-frequency stability)
    // Policy: Loop-Star is mandatory for power electronics applications (DC - 1 MHz)
    bool useMQS_;        // True for MQS mode (Loop equations only, default)
    std::unique_ptr<LoopStarDecomposition> loopStar_;

    void AssembleMatrices();
    void AssembleRHS();
    void SolveLoopStar();  // Loop-Star based solve (mandatory)
    void ComputeImpedance();
};

} // namespace radia

#endif // RAD_RWG_BASIS_H
