/* rad_hdiv_vim.h -- Symmetric HDiv-type VIM (Volume Integral Method) demag operator.
 *
 * The HDiv-type VIM is the Galerkin complement to the multipole-moment MMM MSC kernel: a SYMMETRIC
 * demag operator N = B^T G B with
 *   B = charge map   M |-> (rho = -div M per cell [P0],  sigma = M.n per boundary face [P0])
 *   G = Coulomb Gram (charge-charge interaction; symmetric since 1/r is symmetric)
 * The magnetisation M lives in lowest-order HDiv (RT0): one signed normal-flux DOF per face.
 * Loops = ker(B) (charge-free fields) are FIELD-NULL BY CONSTRUCTION: B.loop = 0 => N.loop =
 * B^T G (B loop) = 0 for ANY G.  The material system is A = (1/chi) M_mass - N (symmetric
 * INDEFINITE), solved iteratively by MINRES (symmetric Krylov); mu_r-independent by construction, so
 * no direct factorization is needed (a symmetric H-LDL^T was prototyped and removed 2026-06-08).
 *
 * This first increment builds the structured-hex RT0 topology + charge map B + a DENSE Coulomb
 * Gram G + assembles N, by hand (no NGSolve), validated against the retired NGSolve prototype
 * golden recorded in docs/hdiv_vim/vim_examples_retirement_results.json
 * (regular 3x3x3 -> ndof=108, n_loop=28).  Later
 * phases swap the dense G for a HACApK symmetric H-matrix (rad_hacapk_hdiv, reusing the Wilton
 * 1/r face integral in rad_poly_analytical.cpp) and add the sparse HDiv mass M_mass.
 *
 * Conventions (CLAUDE.md): row-major [target][source]; +N physical sign; 1/(4pi) in G's kernel.
 */
#ifndef RAD_HDIV_VIM_H
#define RAD_HDIV_VIM_H

#include <vector>
#include <array>

namespace rad_hdiv {

typedef std::array<double, 3> Vec3;

/* One RT0 face of a structured hex grid.  Global face normals are +x / +y / +z (ax = 0/1/2).
 * The face separates cell `lo` (on the -normal side; this face is that cell's HI face) from cell
 * `hi` (on the +normal side; this face is that cell's LO face); -1 == domain boundary on that side. */
struct Face {
    int    ax;        // normal axis: 0=x, 1=y, 2=z (global orientation; topological role of the face)
    double area;      // |F| (true area, distortion-aware)
    Vec3   c;         // centroid (mean of the 4 corner vertices)
    Vec3   v[4];      // 4 corner vertices (for bilinear sub-points on distorted faces)
    int    lo, hi;    // adjacent cell ids (-1 if boundary)
    bool   bnd;       // on the domain boundary
};

struct Mesh {
    int nx, ny, nz;
    std::vector<Face>   faces;     // size = ndof (one normal-flux DOF per face)
    int                 n_cell;
    std::vector<Vec3>   cell_c;    // cell centroids
    std::vector<double> cell_V;    // cell volumes (true volume, distortion-aware)
    std::vector<std::array<Vec3, 8>> cell_verts;  // 8 corner vertices per cell (NGSolve hex order:
                                   // v0(000) v1(001) v2(011) v3(010) v4(100) v5(101) v6(111) v7(110))
    int n_face() const { return (int)faces.size(); }
};

/* Build the structured nx*ny*nz hex grid (cell size h), enumerating x-, y-, z-faces in that order
 * (matches the Python design spec hdiv_vim_structured.py).  For 3x3x3: ndof=108, n_cell=27.
 * distort != 0 applies a smooth node displacement (x + d sin(pi y/L) z, y + 0.83 d x z,
 * z + 0.67 d x sin(pi y/L); L = domain size) so the cells become distorted hexes -- cell volumes,
 * face areas/centroids, and vertices are then derived from the displaced nodes (vertex-based). */
Mesh BuildStructuredRT0(int nx, int ny, int nz, double h = 1.0, double distort = 0.0);

/* Dense charge map B, row-major [charge][face], shape (n_charge x n_face).  Charge rows are the
 * n_cell volume charges (rho = -div M, normalized per unit volume) followed by the n_bnd boundary
 * charges (sigma = M.n, per unit area).  Matches the prototype Bv/V, Bb/area normalization. */
void AssembleChargeMap(const Mesh& m, std::vector<double>& B, int& n_charge, int& n_bnd);

/* Dense symmetric Coulomb Gram G, row-major (n_charge x n_charge).
 *   nsub == 0 : centroid-monopole off-diagonal meas_a meas_b/(4pi r) + cube/square self-energy
 *               (fast, but ~7% off the true demag factor and slightly NON-PD; matches the NGSolve
 *               prototype).  Use for the structural tests (symmetry / loop-nullity are G-agnostic).
 *   nsub >= 1 : ACCURATE -- each charge cell -> nsub^3 (volume) / nsub^2 (boundary face) sub-points,
 *               G_ab = sum_{p,q} w_p w_q/(4pi r_pq), coincident sub-point -> analytic sub-cell self.
 *               Converges to the true Coulomb integral; makes N positive semi-definite and drives
 *               the uniform demag factor to the true 1/3.  Cost is O(n_charge^2 * nsub^6) (golden
 *               sizes only; the production path is the Wilton analytic face integral + H-matrix). */
void AssembleCoulombGram(const Mesh& m, std::vector<double>& G, int& n_charge, int nsub = 0);

/* Assemble the symmetric demag operator N = B^T G B, row-major (n_face x n_face).  nsub passed to
 * AssembleCoulombGram (0 = centroid-monopole, >=1 = accurate sub-point quadrature). */
void AssembleN(const Mesh& m, std::vector<double>& N, int nsub = 0);

/* Lowest-order RT0 HDiv mass matrix M_mass = int phi_i . phi_j, row-major (n_face x n_face).
 * Unit-flux basis -> per-cell per-axis 2x2 block (1/h)[[1/3,1/6],[1/6,1/3]] on (lo_face, hi_face);
 * shared faces accumulate from both adjacent cells.  The material system is A = (1/chi) M_mass - N
 * (symmetric indefinite); the generalized eigenvalues of (N, M_mass) are the demag factors. */
void AssembleMass(const Mesh& m, std::vector<double>& M_mass);

/* Sparse (COO triplet) form of the same RT0 mass M_mass, for the scalable system apply
 * y += inv_chi * M_mass * x (the material system A = (1/chi) M_mass - N is applied as the
 * H-matvec N*x plus this sparse local mass).  Appends (I[k], J[k], V[k]); diag is the per-face
 * accumulated diagonal (length n_face).  Same per-cell per-axis 2x2 blocks as AssembleMass. */
void BuildMassCOO(const Mesh& m, std::vector<int>& I, std::vector<int>& J,
                  std::vector<double>& V, std::vector<double>& diag);

//-------------------------------------------------------------------------
// Reusable charge quadrature + per-pair Gram + per-face charge support.
// These are the building blocks SHARED by the dense AssembleCoulombGram / AssembleN
// and the ON-DEMAND HACApK H-matrix entry function (rad_hacapk_hdiv): the H-matrix
// never forms the dense G/N, it evaluates N[i][j] = sum_a sum_b B[a][i] G[a][b] B[b][j]
// per (face i, face j) using exactly these helpers, so the H-matrix and the dense
// reference agree entry-by-entry (verified in tests/feec/test_hdiv_vim_hmatrix.py).
//-------------------------------------------------------------------------

/* Charge quadrature: centroid + measure per charge (n_cell volume charges then n_bnd
 * boundary charges, in the SAME row order as AssembleChargeMap), plus -- when nsub>=1 --
 * the per-charge sub-point cloud (trilinear hex / bilinear quad) with |detJ| / |x_u x x_v|
 * sub-weights.  nsub<=0 leaves sp/sw empty (centroid-monopole mode). */
struct ChargeQuad {
    int n_charge = 0;
    int n_cell   = 0;
    std::vector<Vec3>                cent;  // [n_charge] centroid
    std::vector<double>              meas;  // [n_charge] measure (cell volume or face area)
    std::vector<std::vector<Vec3>>   sp;    // [n_charge] sub-points (empty if nsub<=0)
    std::vector<std::vector<double>> sw;    // [n_charge] sub-weights (empty if nsub<=0)
};
void BuildChargeQuad(const Mesh& m, int nsub, ChargeQuad& q);

/* G[a][b]: the exact Coulomb Gram entry between charges a, b under quadrature q
 * (== AssembleCoulombGram(nsub)[a][b]).  a==b -> self-energy; centroid-monopole when
 * q has no sub-points, accurate sub-point quadrature otherwise. */
double CoulombGramEntry(const ChargeQuad& q, int a, int b);

/* Per-face charge support (CSC by face): each RT0 face feeds exactly two charges --
 * interior face -> (lo cell, hi cell); boundary face -> (its one cell, its sigma charge).
 * rows[f] = the (<=2) charge-row indices (-1 if unused), coef[f] = the matching B entries
 * (-1/V_lo, +1/V_hi, out_sign/area).  This is AssembleChargeMap stored column-wise. */
struct ChargeMapCSC {
    int n_charge = 0;
    std::vector<std::array<int, 2> >    rows;  // [n_face]
    std::vector<std::array<double, 2> > coef;  // [n_face]
};
void BuildChargeMapCSC(const Mesh& m, ChargeMapCSC& csc);

//-------------------------------------------------------------------------
// M2: analytic charge-Gram potentials (the accurate Wilton surface + phi_tet volume integrals,
// for the C++ scalable ChargeGram H-matrix entry; reused from the dense Python reference).
//-------------------------------------------------------------------------

/* Exact INT_T 1/|r-r'| dA' over a flat triangle (V = 3 verts, [3][3]) at obs r -- the Wilton/Graglia
 * analytic triangle potential (pure 1/r integral, NO 1/4pi).  Matches the dense Python
 * radia.vim._core.tri_potential to ~1e-12 and rad_poly_analytical's I_scalar Wilton edge formula. */
double TriPotential(const double V[3][3], const double r[3]);

/* Newtonian potential INT_tet 1/|P-r'| dV' of a uniform tetrahedron (V = 4 verts, [4][3]) at P, via the
 * divergence theorem reusing TriPotential.  Matches radia.vim._core.phi_tet. */
double PhiTet(const double V[4][3], const double P[3]);

/* Field INT_T (r-r')/|r-r'|^3 dA' of a uniform (sigma=1) flat triangle = -grad TriPotential (the Wilton
 * triangle FIELD; out[3], NO 1/4pi).  Matches radia.vim.flat_triangle_charge_field. */
void TriField(const double V[3][3], const double r[3], double out[3]);

/* Field INT_tet (P-r')/|P-r'|^3 dV' of a uniform tetrahedron = -grad PhiTet, via the divergence theorem
 * (n_face*TriPotential + d_face*TriField).  out[3], NO 1/4pi.  Matches radia.vim.tet_self_volume_field*4pi. */
void TetField(const double V[4][3], const double P[3], double out[3]);

/* ---- degree-1/2 polynomial-charge field kernels (the order<=2 fast path; NO 1/4pi) ----
 * Ports of radia.vim._field, validated entry-by-entry vs the Python (machine precision). */
void TriMoment1(const double V[3][3], const double r[3], double out[3]);          /* INT_T r'/R dS' */
void TriMoment2(const double V[3][3], const double r[3], double out[3][3]);       /* INT_T r'(x)r'/R dS' */
void TetMoment1(const double V[4][3], const double r[3], double out[3]);          /* INT_V r'/R dV' */
void TetVolFieldLinear(const double V[4][3], const double r[3], double rho0,
                       const double g[3], double out[3]);                          /* linear volume charge */
void TetVolFieldQuadratic(const double V[4][3], const double r[3], double rho0,
                          const double g[3], const double Q[3][3], double out[3]); /* quadratic volume charge */
void LinTriField(const double V[3][3], const double r[3], double sigma0,
                 const double s[3], double out[3]);                                /* linear surface charge */
void QuadTriField(const double V[3][3], const double r[3], double sigma0,
                  const double s[3], const double S[3][3], double out[3]);         /* quadratic surface charge */

/* Closest point on a flat triangle / in a tetrahedron to p -- the Duffy singularity origin x0 for the
 * order>=3 / curved singular-quadrature charge-potential path (PhiAtHO_Duffy).  out[3]. */
void ClosestPointTriangle(const double p[3], const double a[3], const double b[3], const double c[3],
                          double out[3]);
void ClosestPointTet(const double V[4][3], const double p[3], double out[3]);

/* CURVED (isoparametric P2) panel support.  ClosestRefTri: xi0 = argmin|X(xi)-p|^2 in reference coords
 * (out[2]).  CurvedTriPotential: INT_curvedtri xi^e0 eta^e1 / |p-X(xi)| dA_curved via the reference Duffy from
 * xi0, evaluating the curved P2 map X(xi) + curved area element at each point.  nodes = 6x3 P2 nodes (corners
 * 0,1,2 ; mid-edges 3=(0-1),4=(1-2),5=(2-0)); gl/gw = an nq-point Gauss-Legendre rule on [0,1]. */
void ClosestRefTri(const double nodes[6][3], const double p[3], double xi0[2]);
double CurvedTriPotential(const double nodes[6][3], int e0, int e1, const double p[3],
                          const double* gl, const double* gw, int nq);

/* CURVED (isoparametric P2) tetrahedron VOLUME-charge support.  nodes = 10x3 P2 nodes (corners 0,1,2,3 ;
 * mid-edges 4=(0-1),5=(1-2),6=(2-0),7=(0-3),8=(1-3),9=(2-3)).  CurvedTetPotential: INT_curvedtet
 * xi^e0 eta^e1 zeta^e2 / |p-X(xi)| dV_curved via the reference Duffy from xi0=ClosestRefTet. */
void ClosestRefTet(const double nodes[10][3], const double p[3], double xi0[3]);
double CurvedTetPotential(const double nodes[10][3], int e0, int e1, int e2, const double p[3],
                          const double* gl, const double* gw, int nq);

/* OUTER-quadrature helpers for the curved charge Gram: curved physical point X(xi) + curved MEASURE (area
 * element dA = |Xu x Xv| for a P2 face, volume element dV = |det dX/dxi| for a P2 tet) at a reference point. */
void CurvedTriMapMeasure(const double nodes[6][3], double xi, double eta, double X[3], double& dA);
void CurvedTetMapMeasure(const double nodes[10][3], double xi, double eta, double zeta, double X[3], double& dV);

} // namespace rad_hdiv
#endif
