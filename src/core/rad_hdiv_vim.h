/* rad_hdiv_vim.h -- Symmetric HDiv-type VIM (Volume Integral Method) demag operator.
 *
 * The HDiv-type VIM is the Galerkin alternative to the collocation MSC kernel: a SYMMETRIC
 * demag operator N = B^T G B with
 *   B = charge map   M |-> (rho = -div M per cell [P0],  sigma = M.n per boundary face [P0])
 *   G = Coulomb Gram (charge-charge interaction; symmetric since 1/r is symmetric)
 * The magnetisation M lives in lowest-order HDiv (RT0): one signed normal-flux DOF per face.
 * Loops = ker(B) (charge-free fields) are FIELD-NULL BY CONSTRUCTION: B.loop = 0 => N.loop =
 * B^T G (B loop) = 0 for ANY G.  The material system is A = (1/chi) M_mass - N (symmetric
 * INDEFINITE), solved by MINRES (symmetric Krylov) and the future symmetric H-LDL^T.
 *
 * This first increment builds the structured-hex RT0 topology + charge map B + a DENSE Coulomb
 * Gram G + assembles N, by hand (no NGSolve), validated against the NGSolve prototype golden
 * (examples/feec_vim/hdiv_demag_quad_self.json: regular 3x3x3 -> ndof=108, n_loop=28).  Later
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
    int    ax;        // normal axis: 0=x, 1=y, 2=z
    double area;      // |F|
    Vec3   c;         // centroid
    int    lo, hi;    // adjacent cell ids (-1 if boundary)
    bool   bnd;       // on the domain boundary
};

struct Mesh {
    int nx, ny, nz;
    std::vector<Face>   faces;     // size = ndof (one normal-flux DOF per face)
    int                 n_cell;
    std::vector<Vec3>   cell_c;    // cell centroids
    std::vector<double> cell_V;    // cell volumes
    int n_face() const { return (int)faces.size(); }
};

/* Build the structured nx*ny*nz hex grid (cell size h), enumerating x-, y-, z-faces in that order
 * (matches the Python design spec hdiv_vim_structured.py).  For 3x3x3: ndof=108, n_cell=27. */
Mesh BuildStructuredRT0(int nx, int ny, int nz, double h = 1.0);

/* Dense charge map B, row-major [charge][face], shape (n_charge x n_face).  Charge rows are the
 * n_cell volume charges (rho = -div M, normalized per unit volume) followed by the n_bnd boundary
 * charges (sigma = M.n, per unit area).  Matches the prototype Bv/V, Bb/area normalization. */
void AssembleChargeMap(const Mesh& m, std::vector<double>& B, int& n_charge, int& n_bnd);

/* Dense symmetric Coulomb Gram G, row-major (n_charge x n_charge).  Off-diagonal = centroid
 * monopole meas_a meas_b/(4pi r); diagonal = cube/square self-energy (placeholder; the accurate
 * distorted-hex self-energy / Wilton self-integral is a later phase). */
void AssembleCoulombGram(const Mesh& m, std::vector<double>& G, int& n_charge);

/* Assemble the symmetric demag operator N = B^T G B, row-major (n_face x n_face). */
void AssembleN(const Mesh& m, std::vector<double>& N);

/* Lowest-order RT0 HDiv mass matrix M_mass = int phi_i . phi_j, row-major (n_face x n_face).
 * Unit-flux basis -> per-cell per-axis 2x2 block (1/h)[[1/3,1/6],[1/6,1/3]] on (lo_face, hi_face);
 * shared faces accumulate from both adjacent cells.  The material system is A = (1/chi) M_mass - N
 * (symmetric indefinite); the generalized eigenvalues of (N, M_mass) are the demag factors. */
void AssembleMass(const Mesh& m, std::vector<double>& M_mass);

} // namespace rad_hdiv
#endif
