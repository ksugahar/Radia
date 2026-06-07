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
 * M_mass).  The symmetric H-LDL^T factorization of the COMPRESSED H-matrix needs rk-aware
 * symmetric H-arith (cHACApK_harith.c calls that "future work") and is a separate layer.
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
};

#endif // __RAD_HACAPK_HDIV_H
