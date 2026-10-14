/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk_peec.h
*
* Project:        RADIA
*
* Description:    HACApK adapter for PEEC (Partial Element Equivalent
*                 Circuit) filament inductance matrices. Implements the
*                 RadHACApKBase virtual contract with the Ruehli 1972
*                 closed-form mutual inductance as the kernel.
*
* First release:  2026-04-16
*
* Reference:      A. E. Ruehli, "Inductance Calculations in a Complex
*                 Integrated Circuit Environment", IBM J. Res. Dev.
*                 16(5), pp. 470-481, 1972.
*
-------------------------------------------------------------------------*/

#ifndef __RAD_HACAPK_PEEC_H
#define __RAD_HACAPK_PEEC_H

#include "rad_hacapk.h"
#include "rad_peec_matrices.h"

/**
 * Recommended default HACApK parameters for the PEEC filament kernel.
 *
 * PEEC has 1 DOF per filament, so leaf_size is larger than the compact
 * magnetic defaults to preserve a useful "DOFs per leaf" scale. max_rank
 * is higher because the 1/r Neumann kernel between slender filaments has
 * slower singular-value decay. eta is slightly larger so adjacent-turn blocks
 * (near-field)
 * stay dense rather than being forced through ACA on a kernel where
 * ACA struggles.
 *
 * Values validated 2026-04-16 on a 50 mm / 6 mm pitch helical solenoid
 * sanity check (see RadHACApKPEECSanityCheck).
 */
inline RadHACApKParams RadHACApKPEECDefaultParams() {
    RadHACApKParams p;
    p.aca_eps = 1.0e-8;  // PEEC sweet spot: ACA rank jumps 14->96, error 4.8%->0.38%
                         // between 1e-7 and 1e-8 for Ruehli 1/r kernel (bench 2026-04-16)
    p.leaf_size = 128;   // 1 DOF/fil -> useful DOFs/leaf for PEEC
    p.eta = 3.0;         // keep near-field blocks dense (helix adjacent turns)
    p.max_rank = 400;    // 1/r kernel needs higher rank than compact magnetic kernels
    p.print_level = 0;
    return p;
}

/**
 * RadHACApKPEECManager adapts the HACApK H-matrix pipeline to the PEEC
 * filament inductance kernel.
 *
 * Layout:
 *   - one DOF per filament (uniform nffc = 1)
 *   - m_coordinates holds filament center points (used by HACApK clustering)
 *   - system matrix stored by HACApK is +L(i, j) directly (default
 *     ComputeSystemEntry on the base returns GetInteractionMatrixElement
 *     unchanged)
 *
 * Frequency-dependent terms (jωL + R + Zs) are combined outside HACApK
 * by the enclosing BiCGSTAB loop / PRIMA reduction. HACApK stores only
 * the real, positive, symmetric L; the matvec gives y = L * x.
 *
 * Usage:
 *   radia::PEECMatrixBuilder builder;
 *   // ... populate segments / nodes ...
 *   RadHACApKPEECManager mgr(builder);
 *   mgr.BuildHMatrix();
 *   std::vector<double> x(N), y(N);
 *   mgr.MatVec(x, y);   // y = L * x
 *
 * The builder must remain alive for the lifetime of the manager; the
 * manager reads filament geometry and evaluates Ruehli entries on
 * demand via builder.InductanceEntry(i, j).
 */
class RadHACApKPEECManager : public RadHACApKBase {
public:
    explicit RadHACApKPEECManager(const radia::PEECMatrixBuilder& builder);
    ~RadHACApKPEECManager() override = default;

    /**
     * Ruehli mutual inductance L(i, j) on demand. Diagonal is self-
     * inductance (Grover closed form); off-diagonal is the parallel
     * analytical formula, fourfil near-field, or 8x8 Gauss-Legendre
     * quadrature depending on geometry — dispatched by
     * PEECMatrixBuilder::InductanceEntry.
     */
    double GetInteractionMatrixElement(int dof_i, int dof_j) const override;

    // Base::ComputeSystemEntry default (return GetInteractionMatrixElement)
    // is correct for PEEC — no override needed.

protected:
    void ExtractCoordinates() override;
    void OnBeforeBuild() override;
    void InitializeInvChi() override;
    bool IsVariableDOF() const override { return false; }
    int GetUniformNFFC() const override { return 1; }

private:
    const radia::PEECMatrixBuilder& m_builder;
};

/**
 * Self-contained sanity check: build N parallel filaments, form the
 * dense L matrix via PEECMatrixBuilder::ComputeL, form the HACApK L
 * via RadHACApKPEECManager, run MatVec on a deterministic random
 * vector in both, and return the max relative error between them.
 *
 * Used from Python via radia._TestPEECHACApKSanity(n) to verify the
 * HACApK PEEC adapter pipeline end-to-end (ExtractCoordinates,
 * BuildHMatrix, MatVec permutation, GetInteractionMatrixElement).
 *
 * Returns max|y_hmat - y_dense| / max|y_dense|. Should be below
 * ~aca_eps (typ. 1e-8) for well-separated filaments.
 */
double RadHACApKPEECSanityCheck(int n_filaments);

#endif // __RAD_HACAPK_PEEC_H
