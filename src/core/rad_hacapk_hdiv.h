/* rad_hacapk_hdiv.h -- HACApK H-matrix for the symmetric HDiv-type VIM demag operator.
 *
 * The HDiv-type VIM demag operator N = B^T G B (rad_hdiv_vim.h) is SYMMETRIC and has the
 * loops field-null by construction.  This manager builds N as a HACApK H-matrix so the
 * operator SCALES: the DOFs are the RT0 faces (one normal-flux DOF per face), the cluster
 * tree is over the face centroids, and the ACA entry function evaluates N[i][j] ON DEMAND
 * as the sparse charge-cluster Coulomb interaction
 *
 *     N[i][j] = sum_{a in supp(face i)} sum_{b in supp(face j)} B[a][i] G[a][b] B[b][j]
 *
 * where each face's charge support has <= 2 entries (interior: lo/hi cell; boundary: cell +
 * sigma charge) and G[a][b] is the per-pair Coulomb Gram (rad_hdiv::CoulombGramEntry, nsub
 * controlled).  No dense G or N is ever formed -- the H-matrix is the storage.  The matvec
 * (O(N log N)) is all a symmetric Krylov solver (MINRES) needs; we measured the plain-Jacobi
 * MINRES iteration count to be mu_r-INDEPENDENT (~8), so this is a scalable symmetric solver.
 *
 * Stores +N (ComputeSystemEntry left at the base default).  The material system
 * A = (1/chi) M_mass - N is applied by the caller (N via this H-matvec + the sparse local
 * M_mass).  It is solved iteratively (SolveLinearMaterial Jacobi-PCG / SolveNonlinearPicard); the
 * operator is mu_r-independent by construction, so no direct factorization is used.
 *
 * Conventions (CLAUDE.md): row-major [target][source]; +N physical sign; HACApK ACA+ only.
 */
#ifndef __RAD_HACAPK_HDIV_H
#define __RAD_HACAPK_HDIV_H

#include "rad_hacapk.h"     // RadHACApKBase
#include "rad_hdiv_vim.h"   // rad_hdiv::Mesh / ChargeMapCSC / ChargeQuad

//-------------------------------------------------------------------------
// RadHACApKHDivManager: builds N = B^T G B (HDiv-type VIM) as a HACApK H-matrix.
//-------------------------------------------------------------------------

class RadHACApKHDivManager : public RadHACApKBase {
public:
    // Structured nx*ny*nz RT0 hex grid (cell size h, optional distortion).  nsub controls the
    // per-pair Gram quadrature in the entry function (0 = centroid-monopole, >=1 = accurate
    // sub-point, matching rad_hdiv::AssembleCoulombGram(nsub) entry-by-entry).
    RadHACApKHDivManager(int nx, int ny, int nz, double h, double distort, int nsub);
    ~RadHACApKHDivManager() override {}   // base dtor frees the H-matrix resources

    // +N(i,j) for 0-based ORIGINAL face indices (the charge-cluster Coulomb sum).
    double GetInteractionMatrixElement(int dof_i, int dof_j) const override;

    // Material system apply: y = A x = inv_chi * (M_mass x) - (N x), where N x is the O(N log N)
    // H-matvec and M_mass is the sparse local RT0 mass.  This is the operator a symmetric Krylov
    // solver (MINRES) iterates on; A is symmetric indefinite (its generalized eigenvalues vs
    // M_mass are inv_chi - demag_factor, demag_factor in [0,1]).
    void ApplySystem(const std::vector<double>& x, double inv_chi, std::vector<double>& y);

    // Jacobi preconditioner diagonal: diag(A)_f = inv_chi * M_mass_ff - N_ff (N_ff via the
    // O(1) entry function).  For preconditioned MINRES (M = |diag(A)|).
    std::vector<double> DiagSystem(double inv_chi) const;

    const rad_hdiv::Mesh& GetMesh() const { return m_mesh; }

protected:
    void ExtractCoordinates() override;   // mesh + face-centroid coordinates
    void OnBeforeBuild() override;        // sparse charge map + charge quadrature
    void InitializeInvChi() override;     // zeros (the +N H-matrix does not use 1/chi)
    bool IsVariableDOF() const override { return false; }
    int  GetUniformNFFC() const override { return 1; }   // one normal-flux DOF per face

private:
    int    m_nx, m_ny, m_nz, m_nsub;
    double m_h, m_distort;
    rad_hdiv::Mesh        m_mesh;
    rad_hdiv::ChargeMapCSC m_csc;    // per-face charge support (B columns)
    rad_hdiv::ChargeQuad   m_quad;   // centroids/measures (+ sub-points) for on-demand G
    // sparse RT0 mass M_mass (COO) + its per-face diagonal, for the scalable system apply
    std::vector<int>     m_mI, m_mJ;
    std::vector<double>  m_mV, m_mass_diag;
};

//-------------------------------------------------------------------------
// RadHACApKChargeGram: the charge-charge Coulomb Gram G as a HACApK H-matrix.
//-------------------------------------------------------------------------

/* The UNSTRUCTURED / general-mesh production path.  Charges = volume cells (rho = -div M) +
 * boundary faces (sigma = M.n) extracted from ANY RT0 mesh (e.g. NGSolve HDiv(order=0) on a tet
 * mesh -- see examples/feec_vim/hdiv_demag_tet.py).  This manager builds the n_charge x n_charge
 * Coulomb Gram G as a HACApK H-matrix (a clean 1/r kernel over the charge centroids):
 *   G[a][b] = meas_a meas_b / (4pi |c_a - c_b|)   (a != b, centroid monopole)
 *   G[a][a] = self_energy[a]                       (the accurate sub-divided self, computed by the
 *                                                   caller per element shape -- tet/tri/hex/quad).
 * The demag operator N = B^T G B is then applied matrix-free as B^T (G-Hmatvec (B m)) with B the
 * sparse charge map; the H-matrix gives the O(N log N) Gram matvec that makes the solve scalable
 * on real geometry.  Stores +G (ComputeSystemEntry = default); MatVec/stats from RadHACApKBase. */
class RadHACApKChargeGram : public RadHACApKBase {
public:
    // MONOPOLE mode: centroids [n*3] charge centroids; measures [n] cell volumes / face areas;
    // self_energy [n] the diagonal G[a][a] (caller-computed accurate self-energy, element-shape-aware).
    RadHACApKChargeGram(std::vector<double> centroids,
                        std::vector<double> measures,
                        std::vector<double> self_energy);

    // ANALYTIC mode (M2b): the EXACT charge Gram from per-charge GEOMETRY -- matches the dense Python
    // build_demag(analytic_gram=True) entry-by-entry.  cell_verts [n_el*12] (tets, 4 verts) then
    // face_verts [n_bf*9] (triangles, 3 verts); the n_charge charges are the n_el volume cells
    // (rho = -div M) followed by the n_bf boundary faces (sigma = M.n).  Entry
    //   G[a][b] = (1/4pi) INT_a Phi_b   (Phi_b = PhiTet/TriPotential of source b, exact analytic),
    // the outer INT_a by tet barycentric sub-points (cells) / Dunavant-5 (faces), symmetrized; the
    // diagonal is the analytic self (the Wilton/phi_tet potential is exact through the 1/r singularity).
    RadHACApKChargeGram(std::vector<double> cell_verts,
                        std::vector<double> face_verts,
                        int n_el);
    ~RadHACApKChargeGram() override {}

    double GetInteractionMatrixElement(int a, int b) const override;

    // M3 (the iterative-solve hot kernel in C++): solve the SPD HDiv-VIM linear material system
    //   ((1/chi) M_mass + B^T G B) m = rhs
    // by Jacobi-preconditioned conjugate gradients, with G applied as THIS charge-Gram H-matvec
    // (O(N log N)) -- no dense N, no Python per-iteration glue.  This is the linear soft-iron demag
    // solve AND the symmetric Picard warmstart of the nonlinear Newton.  Sparse inputs are caller-
    // provided: B as CSR over charges (B_indptr [n_charge+1], B_indices/B_data = face columns, so
    // (B x)[charge] = sum data*x[face]); M_mass as COO (mI,mJ,mV) on the n_face DOFs; prec = the
    // Jacobi diagonal of the system (length n_face).  Returns m (length n_face); iters_out = CG iters.
    std::vector<double> SolveLinearMaterial(
        const std::vector<int>& B_indptr, const std::vector<int>& B_indices,
        const std::vector<double>& B_data, int n_face,
        const std::vector<int>& mI, const std::vector<int>& mJ, const std::vector<double>& mV,
        double inv_chi, const std::vector<double>& prec, const std::vector<double>& rhs,
        double tol, int maxit, int& iters_out);

    // M3 (the NONLINEAR solve in C++): scalar-chi Picard for the isotropic nonlinear demag.
    // Each Picard step is a SolveLinearMaterial solve of ((1/chi) M_mass + B^T G B) m = H0*(M_mass mu),
    // then chi <- 0.5 chi + 0.5*chi_sec(|H|) with the closed-form saturating curve
    //   M(H) = chi0 H / (1 + chi0 |H|/Msat)   ->   chi_sec(|H|) = chi0/(1 + chi0|H|/Msat),
    // and the scalar self-consistent field H = H0 - Dscal*M_avg, Dscal = mu.(B^T G B mu)/denom,
    // M_avg = mu.(M_mass m)/denom.  Converges to the scalar fixed point M_avg = M(H0 - Dscal*M_avg)
    // -- the full nonlinear physics for an isotropic body, with NO NGSolve per iteration (the
    // per-element tensor-tangent refinement for non-uniform M stays NGSolve).  All sparse inputs as in
    // SolveLinearMaterial; Mmass_diag + N_diag build the per-chi Jacobi preconditioner.
    struct PicardResult { std::vector<double> m; double Mavg; double chi; double Dscal; int iters; };
    PicardResult SolveNonlinearPicard(
        const std::vector<int>& B_indptr, const std::vector<int>& B_indices,
        const std::vector<double>& B_data, int n_face,
        const std::vector<int>& mI, const std::vector<int>& mJ, const std::vector<double>& mV,
        const std::vector<double>& Mmass_diag, const std::vector<double>& N_diag,
        const std::vector<double>& mu, double denom,
        double chi0, double Msat, double H0,
        int picard_iters, double cg_tol, int cg_maxit);

protected:
    void ExtractCoordinates() override;
    void OnBeforeBuild() override {}
    void InitializeInvChi() override { m_inv_chi.assign(m_ndof, 0.0); }
    bool IsVariableDOF() const override { return false; }
    int  GetUniformNFFC() const override { return 1; }

private:
    double PhiAt(int src, const double p[3]) const;   // exact analytic potential of source charge src at p
    double QuadDot(int tgt, int src) const;            // (1/4pi) sum_p w_p PhiAt(src, p) over tgt's outer quad

    std::vector<double> m_cent, m_meas, m_self;        // monopole mode (m_cent also = the cluster-tree points)
    int  m_n = 0;
    // analytic mode (M2b)
    bool m_analytic = false;
    int  m_n_el = 0;
    std::vector<double> m_cellV, m_faceV;              // [n_el*12], [n_bf*9]
    std::vector<std::vector<rad_hdiv::Vec3>> m_qp;     // [n] outer-quad points per charge
    std::vector<std::vector<double>>          m_qw;    // [n] outer-quad weights per charge
};

#endif // __RAD_HACAPK_HDIV_H
