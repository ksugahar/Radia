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
#include <unordered_map>

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

    // SYSTEM-A mode for H-LU: when enabled (chi > 0), the H-matrix stores the form-1 SOFT-IRON
    // MATERIAL SYSTEM  A = M_mass + chi*N  (the uniform-chi material system ((1/chi)M_mass + N)
    // scaled by chi -- well-conditioned: loop modes see only M_mass), so the HACApK H-LU
    // (cHACApK_hlu_*) factors A directly = a scalable DIRECT solve / strong preconditioner for the
    // HDiv-VIM demag.  Default OFF (ComputeSystemEntry = +N for the matvec path).  Call BEFORE
    // BuildHMatrix.  M_mass[i][j] is read from the sparse RT0 mass via m_mass_map (built in
    // OnBeforeBuild).
    void SetSystemMode(double chi) { m_system_chi = chi; m_system_mode = (chi > 0.0); }
    double ComputeSystemEntry(int dof_i, int dof_j) const override;

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
    // SYSTEM-A H-LU mode: store A = M_mass + chi*N via ComputeSystemEntry (vs the default +N).
    bool   m_system_mode = false;
    double m_system_chi  = 0.0;
    std::unordered_map<long long, double> m_mass_map;  // (i*ndof+j) -> M_mass[i][j], O(1) lookup
};

//-------------------------------------------------------------------------
// RadHACApKChargeGram: the charge-charge Coulomb Gram G as a HACApK H-matrix.
//-------------------------------------------------------------------------

/* The UNSTRUCTURED / general-mesh production path.  Charges = volume cells (rho = -div M) +
 * boundary faces (sigma = M.n) extracted from ANY RT0 mesh (e.g. NGSolve HDiv(order=0) on a tet
 * mesh -- see examples/vim/hdiv_demag_tet.py).  This manager builds the n_charge x n_charge
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
    //
    // near_factor: the NEAR/FAR split that makes the BUILD fast.  A pair (a,b) is NEAR when
    // |c_a - c_b| <= near_factor*(size_a + size_b) and uses the expensive analytic entry; FAR pairs use
    // the cheap centroid-monopole meas_a*meas_b/(4pi r).  Physically correct: the analytic entry only
    // matters for the NEAR (non-uniform-M, div M != 0) interaction; the far field is monopole to
    // O((size/r)^2) (the validated monopole+near-correction split).  near_factor defaults to a HUGE
    // value => ALL pairs near => all-analytic (matches the dense build_demag(analytic_gram=True) golden);
    // pass near_factor ~= 2 for the fast split.
    // The CELL outer-quad uses a built-in 4-pt Gauss-Duffy tet rule (64 nodes) -- the order-0 charge is
    // constant so the inner is the EXACT analytic PhiTet and INT_T PhiTet dx is smooth, integrated to
    // ~machine precision.  (The old hardcoded equal-weight _bary_tet(3) under-integrated the volume self-
    // energy by ~6.5% -- golden-invisible because every uniform-M demag golden has div M = 0.)
    RadHACApKChargeGram(std::vector<double> cell_verts,
                        std::vector<double> face_verts,
                        int n_el, double near_factor = 1e30);

    // HIGH-ORDER (order-p) mode: POLYNOMIAL charges (a monomial basis on each host element), the order-p
    // extension validated against the dense Python build_demag_highorder.  charge_host[c] = host element
    // index (a cell when charge_kind[c]==0, a boundary face when ==1); charge_expo[3*c+{0,1,2}] = the
    // monomial exponents in the host's REFERENCE barycentric coords (tet cell: lam1^i lam2^j lam3^k;
    // face: lam1^i lam2^j, the 3rd exponent ignored).  The SAME reference convention is used by the
    // Python charge-density map B (so B and G share the basis; N = B^T G B is then basis-invariant and
    // matches the NGSolve-L2-basis dense reference's demag).  Entry
    //   G[a][b] = (1/4pi) INT_ha INT_hb m_a(x) m_b(y)/|x-y|
    // = the monomial-WEIGHTED outer quad (m_a folded into the outer weights) x the polynomial-charge
    // inner potential PhiAt(b,.) by singularity SUBTRACTION reusing the exact PhiTet/TriPotential through
    // the 1/r singularity.  ref_tet_pts[nqt*3]/ref_tet_w[nqt] (weights sum to 1/6) + ref_tri_pts[nqr*2]/
    // ref_tri_w[nqr] (sum to 1/2) are the REFERENCE-element Gauss-Duffy rules (Python-supplied), mapped per
    // host and used for BOTH the monomial-weighted outer quad and the FIXED inner-potential subtraction
    // table.
    //
    // NEAR/FAR adaptive quadrature (build speedup, ACCURACY-PRESERVING): the per-entry cost is the nested
    // outer x inner quadrature (~O(quad^6) for vol-vol).  The expensive HIGH-quad subtraction is only needed
    // for NEAR/self pairs (through the 1/r singularity); for a well-separated (FAR) pair the kernel 1/|x-y|
    // is SMOOTH, so a CHEAP LOW-quad PLAIN double-Gauss (no PhiTet, no subtraction) is already accurate.  If
    // the optional ref_*_lo LOW-quad rules are supplied (and ho_far_factor < inf), a pair (a,b) with
    // |c_a-c_b| > ho_far_factor*(size_a+size_b) uses QuadDotFar (the low-quad plain double sum); NEAR/self
    // pairs keep the full high-quad subtraction.  This is NOT FMM/multipole (the zero-mean high-order modes
    // have zero monopole, so a monopole far is WRONG); it is just adaptive quadrature order, and the HACApK
    // ACA still compresses the well-separated low-rank blocks (now from cheap entries).  ho_far_factor
    // defaults to 1e30 (=> all pairs NEAR => the original all-high-quad behavior, golden-equivalent).
    RadHACApKChargeGram(std::vector<double> cell_verts, std::vector<double> face_verts, int n_el,
                        std::vector<int> charge_host, std::vector<int> charge_kind,
                        std::vector<int> charge_expo,
                        std::vector<double> ref_tet_pts, std::vector<double> ref_tet_w,
                        std::vector<double> ref_tri_pts, std::vector<double> ref_tri_w,
                        std::vector<double> ref_tet_pts_lo = {}, std::vector<double> ref_tet_w_lo = {},
                        std::vector<double> ref_tri_pts_lo = {}, std::vector<double> ref_tri_w_lo = {},
                        double ho_far_factor = 1e30,
                        std::vector<double> ref_tet_pts_in = {}, std::vector<double> ref_tet_w_in = {},
                        std::vector<double> ref_tri_pts_in = {}, std::vector<double> ref_tri_w_in = {});
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

    // The mu_r-INDEPENDENT production MATERIAL solve in C++: Jacobi-preconditioned MINRES for the
    // SYMMETRIC INDEFINITE system A m = rhs, A = inv_chi*M_mass - B^T G B (eigenvalues vs M_mass =
    // inv_chi - demag_factor, demag_factor in [0,1] -> indefinite -> MINRES, not CG).  G applied as
    // the analytic charge-Gram H-matvec.  This is the loop-field-null payoff: the iteration count is
    // mu_r-INDEPENDENT (golden test_hdiv_vim_solve.py: ~80 flat, even decreasing on distorted meshes),
    // unlike the +N Picard CG (SolveLinearMaterial) whose near-null loop modes ill-condition at high
    // mu_r.  TaskManager: wrapped in ngcore::RegionTaskManager so the HACApK H-matvec runs parallel
    // (multi-thread) under the caller's `with TaskManager():` (or stands up its own pool); the O(N)
    // vector ops are serial (negligible vs the O(N log N) matvec).  prec = the SPD Jacobi
    // preconditioner (length n_face, e.g. |inv_chi*M_mass_diag - N_diag|).
    std::vector<double> SolveMaterialMINRES(
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
    std::vector<double> m_size;                        // [n] characteristic size: vol^(1/3) / area^(1/2)
    double m_near_factor = 1e30;                       // near/far split: NEAR if |c_a-c_b| <= nf*(size_a+size_b)

    // HIGH-ORDER (polynomial-charge) mode
    bool m_highorder = false;
    std::vector<int> m_host, m_kind, m_expo;           // [n] host elem, [n] 0=cell/1=face, [n*3] monomial exponents
    std::vector<int> m_nmono;                          // [n] # co-located charges per (kind,host) group -- QuadDot memo gating (skip groups of 1)
    long long m_build_id = 0;                          // monotonic per-build id -> the QuadDot thread_local memo owner key (pointer-reuse-safe)
    std::vector<double> m_cellInv;                     // [n_el*9] physical->ref affine inverse per cell (row-major)
    std::vector<double> m_faceGinv;                    // [n_bf*4] 2x2 (a.a) Gram inverse per face (for 2D ref coords)
    std::vector<std::vector<rad_hdiv::Vec3>> m_inP;    // [n] FIXED inner-potential Gauss points per HOST (cell/face)
    std::vector<std::vector<double>>          m_inW;   // [n] inner-potential Gauss weights (sum = host measure)
    std::vector<std::vector<double>>          m_srcval; // [n] PRECOMPUTED m_src(y_q) at the FIXED m_inP points -- bit-exact hoist of EvalMono out of the PhiAtHO inner loop (value depends only on (src,q))
    // near/far adaptive quadrature: LOW-quad tables for the cheap FAR plain double-Gauss (empty => disabled)
    double m_ho_far_factor = 1e30;                     // FAR if |c_a-c_b| > m_ho_far_factor*(size_a+size_b)
    std::vector<std::vector<rad_hdiv::Vec3>> m_qp_lo;  // [n] LOW-quad outer points (m_a folded into m_qw_lo)
    std::vector<std::vector<double>>          m_qw_lo; // [n] LOW-quad outer weights (monomial-folded)
    std::vector<std::vector<rad_hdiv::Vec3>> m_inP_lo; // [n] LOW-quad inner points (plain)
    std::vector<std::vector<double>>          m_inW_lo;// [n] LOW-quad inner weights (plain, NOT monomial-folded)
    std::vector<std::vector<double>>          m_srcval_lo; // [n] PRECOMPUTED m_src(y_q) at the FIXED m_inP_lo points (for QuadDotFar)
    double EvalMono(int charge, const double p[3]) const;   // charge's monomial at physical p (host ref-coord map)
    double PhiAtHO(int src, const double p[3]) const;       // polynomial-charge inner potential (subtraction, NEAR)
    double QuadDotFar(int tgt, int src) const;              // cheap LOW-quad plain double-Gauss (FAR, no subtraction)
};

#endif // __RAD_HACAPK_HDIV_H
