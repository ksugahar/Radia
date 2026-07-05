/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk_bem.h
*
* Project:        RADIA
*
* Description:    HACApK adapter for scalar Laplace BEM Galerkin matrices
*                 (single-layer SL, double-layer DL).  Implements the
*                 RadHACApKBase virtual contract with a pre-computed dense
*                 entry table as the kernel.  This lets Python build the
*                 dense Galerkin matrix once via our in-tree Sauter-Schwab
*                 quadrature (src/radia/bem/sibc_hacapk.py) and then hand
*                 it to HACApK for ACA compression and fast MatVec.
*
* First release:  2026-04-30
*
-------------------------------------------------------------------------*/

#ifndef __RAD_HACAPK_BEM_H
#define __RAD_HACAPK_BEM_H

#include "rad_hacapk.h"
#include <vector>

/**
 * Default HACApK parameters tuned for the scalar Laplace SL kernel
 * (Galerkin P1 on flat triangles).
 *
 * SL has weak (logarithmic) singular-value decay, so eps must be small
 * to retain accuracy; leaf_size is moderate (1 DOF per vertex, smaller
 * blocks fine).  eta=2.0 is the standard admissibility for Laplace.
 */
inline RadHACApKParams RadHACApKBEMDefaultParams() {
    RadHACApKParams p;
    p.aca_eps = 1.0e-6;   // Laplace SL: ACA needs tighter eps than compact magnetic kernels
    p.leaf_size = 64;     // 1 DOF/vertex; 64 vertices per leaf
    p.eta = 2.0;
    p.max_rank = 400;
    p.print_level = 0;
    return p;
}

/**
 * RadHACApKBEMManager adapts the HACApK H-matrix pipeline to a scalar
 * Galerkin BEM matrix (e.g. Laplace SL) on flat P1 triangles.
 *
 * Design:
 *   - 1 DOF per vertex (uniform nffc = 1)
 *   - m_coordinates holds vertex coordinates (used by HACApK clustering)
 *   - the entry table is a pre-computed (n_v * n_v) dense matrix passed
 *     in at construction; GetInteractionMatrixElement(i, j) returns
 *     dense[i, j] in O(1)
 *   - the matrix stored by HACApK is +M(i, j) directly (default
 *     ComputeSystemEntry on the base returns GetInteractionMatrixElement
 *     unchanged)
 *
 * Memory tradeoff: the dense table costs N^2 * 8 bytes during the build
 * (e.g. 470 MB for N=2477), but after BuildHMatrix the H-matrix uses
 * only O(N log N) memory; the dense table can then be discarded by the
 * Python caller.
 */
class RadHACApKBEMManager : public RadHACApKBase {
public:
    /**
     * @param coordinates n_v x 3 vertex coordinate array, row-major,
     *                    not owned (caller keeps the buffer alive until
     *                    BuildHMatrix returns).
     * @param dense_entries n_v * n_v row-major dense matrix, M[i, j] at
     *                      offset i * n_v + j, not owned (kept alive
     *                      until BuildHMatrix returns).
     * @param n_v number of vertices (= number of DOFs).
     */
    RadHACApKBEMManager(const double* coordinates,
                        const double* dense_entries,
                        int n_v);
    ~RadHACApKBEMManager() override = default;

    /**
     * Return M(i, j) from the dense entry table in O(1).
     */
    double GetInteractionMatrixElement(int dof_i, int dof_j) const override;

    // Base::ComputeSystemEntry default (return GetInteractionMatrixElement)
    // is correct for BEM -- no override needed.

protected:
    void ExtractCoordinates() override;
    void OnBeforeBuild() override;
    void InitializeInvChi() override;
    bool IsVariableDOF() const override { return false; }
    int GetUniformNFFC() const override { return 1; }

private:
    const double* m_coords_ext;   // not owned, n_v * 3
    const double* m_entries_ext;  // not owned, n_v * n_v row-major
    int m_n_v;
};

#endif // __RAD_HACAPK_BEM_H
