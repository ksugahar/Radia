/*-------------------------------------------------------------------------
*
* File name:      rad_relaxation_methods.h
*
* Project:        RADIA
*
* Description:    Relaxation methods
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RAD_RELAXATION_METHODS_H
#define __RAD_RELAXATION_METHODS_H

#include "rad_interaction.h"
#include "rad_math_methods.h"
#include <memory>

//-------------------------------------------------------------------------
// Solver method constants
//-------------------------------------------------------------------------

namespace RadSolverMethod {
	constexpr int LU = 0;  // LU direct solver. Non-symmetric Krylov relaxation methods were retired.
}

//-------------------------------------------------------------------------
// Nonlinear Iteration Context (shared across all solver methods)
//
// Phase 1 of unification: Encapsulates common state for Newton-Raphson iteration
// Reference: internal/design/NONLINEAR_ITERATION_UNIFICATION_PLAN.md
//-------------------------------------------------------------------------

// Forward declaration
class radTg3dRelax;
class radTPolyhedron;
class radTHysteresisMaterial;

/**
 * Shared context for nonlinear magnetostatic iteration.
 *
 * This struct encapsulates the nonlinear state used by the retained LU
 * relaxation route. The former non-symmetric Krylov/HACApK subclasses were
 * retired; current H-matrix solvers use their own HDiv, PEEC, and BEM APIs.
 *
 * The Newton-Raphson iteration solves:
 *   (I - chi*N) * M = chi * H_ext
 * where chi is updated each iteration based on the B-H curve.
 */
struct NonlinearContext {
	// Problem dimensions
	int totalDOF;           // Total degrees of freedom (sum of element DOFs)
	int AmOfMainElem;       // Number of magnetic elements

	// Pointers to interaction matrix data (owned by IntrctPtr, DO NOT delete)
	double* FlatMagn;       // Magnetization array [totalDOF]
	double* FlatField;      // Field array [totalDOF]
	double* FlatExtern;     // External field array [totalDOF]

	// State vectors for Newton-Raphson iteration
	std::vector<double> OldMagn;           // Previous iteration magnetization [totalDOF]
	std::vector<double> OldChi;            // Previous iteration chi [AmOfMainElem]
	std::vector<double> OldBnorm;          // Previous iteration |B| [AmOfMainElem]
	std::vector<double> CurrentChiArray;   // Current chi values [AmOfMainElem]
	std::vector<double> NewFieldArray;     // Computed field after linear solve [totalDOF]

	// Element cache (avoid repeated virtual calls)
	std::vector<radTPolyhedron*> polyCache;  // Polyhedron pointers [AmOfMainElem]

	// Flags
	bool all_materials_linear;  // True if all materials are linear (converge in 1 iteration)

	// Convergence tracking
	double B_sat;               // Saturation B for relative convergence check
	double max_B_rel_change;    // Maximum relative B change this iteration
	double nonlinear_tol;       // Outer nonlinear convergence tolerance

	// Under-relaxation parameter (0 = full step, 0.5 = 50% damping)
	double relax_param;

	// Base matrix for linear system (geometric part without chi)
	std::vector<double> BaseMatrix;

	// Newton-Raphson fields
	bool use_newton;                           // True to use differential chi
	std::vector<double> DifferentialChiArray;  // chi_d per element [AmOfMainElem]
	std::vector<double> OldSigma;              // sigma_old for Newton RHS correction [totalDOF]

	// Newton line search damping
	bool newton_damping_enabled;               // Enable adaptive line search damping
	int newton_ls_max_iter;                    // Max line search backtracks (default: 5)
	double newton_ls_min_omega;                // Minimum omega threshold (default: 0.01)
	int total_ls_backtracks;                   // Statistics: cumulative backtracks
	std::vector<double> accepted_omegas;       // ω values per iteration (for debugging)

	// B-input Newton-Raphson fields (energy-based hysteresis)
	bool use_b_input;                          // True for B-input Newton solver
	std::vector<std::vector<TVector3d>> saved_hys_states;  // Per-element saved Jk states

	// Constructor
	NonlinearContext()
		: totalDOF(0)
		, AmOfMainElem(0)
		, FlatMagn(nullptr)
		, FlatField(nullptr)
		, FlatExtern(nullptr)
		, all_materials_linear(true)
		, B_sat(1.0)
		, max_B_rel_change(0.0)
		, nonlinear_tol(0.0)
		, relax_param(0.0)
		, use_newton(false)
		, newton_damping_enabled(false)
		, newton_ls_max_iter(5)
		, newton_ls_min_omega(0.01)
		, total_ls_backtracks(0)
		, use_b_input(false)
	{}
};

//-------------------------------------------------------------------------
// Helper functions for the retained nonlinear LU iteration
//-------------------------------------------------------------------------

/**
 * Initialize NonlinearContext with problem dimensions, arrays, and initial chi.
 * Called once at the start of AutoRelax.
 *
 * @param ctx Output context to initialize
 * @param IntrctPtr Interaction data (provides geometry, materials, arrays)
 * @param MagnResetIsNotNeeded If false, reset magnetization before solving
 * @return true if initialization successful, false if problem is empty/invalid
 */
bool InitializeNonlinearContext(NonlinearContext& ctx, radTInteraction* IntrctPtr, bool MagnResetIsNotNeeded);

/**
 * Build base matrix with correct sign convention for 3-DOF relaxable elements.
 *
 * @param ctx Context with BaseMatrix to fill
 * @param IntrctPtr Interaction data
 * @return true on success, false on memory allocation failure
 */
bool BuildBaseMatrix(NonlinearContext& ctx, radTInteraction* IntrctPtr);

/**
 * Store old values (M, chi, B-norm) before linear solve for convergence check.
 * Called at the start of each nonlinear iteration.
 *
 * @param ctx Context with OldMagn, OldChi, OldBnorm to update
 * @param IntrctPtr Interaction data
 */
void StoreOldValuesAndComputeBnorm(NonlinearContext& ctx, radTInteraction* IntrctPtr);

/**
 * Update element magnetization from flat array and compute H = M/chi.
 * Called after linear solve to sync element state.
 *
 * @param ctx Context with FlatMagn (solution)
 * @param IntrctPtr Interaction data
 */
void UpdateMagnAndComputeH(NonlinearContext& ctx, radTInteraction* IntrctPtr);

/**
 * Update chi using ELF-style dual method and check convergence.
 * Returns max relative B-field change (ELF mucal2 criterion).
 *
 * @param ctx Context with CurrentChiArray to update
 * @param IntrctPtr Interaction data
 * @return Maximum relative B-field change across all elements
 */
double UpdateChiAndCheckConvergence(NonlinearContext& ctx, radTInteraction* IntrctPtr);

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

class radTIterativeRelaxMeth {
protected:
	radTInteraction* IntrctPtr;

public:
	radTIterativeRelaxMeth(radTInteraction* InIntrctPtr) { IntrctPtr = InIntrctPtr;}
	radTIterativeRelaxMeth() { IntrctPtr = 0;}

	virtual void DefineNewMagnetizations() {}

	void MakeN_iter(int);
	void ComputeRelaxStatusParam(const TVector3d*, const TVector3d*, const TVector3d*);

	/**
	 * Unified nonlinear iteration using helper functions.
	 * Calls virtual SolveLinearStep which is overridden by each solver.
	 *
	 * @return Number of nonlinear iterations performed
	 */
	int AutoRelax_Unified(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded = 0);

	/**
	 * B-input Newton-Raphson solver for energy-based hysteresis materials.
	 *
	 * Solves F(M) = M - Inverse(mu_0*(H_ext + N*M + M))/mu_0 = 0
	 * using Newton-Raphson with analytical Jacobian dJ/dB from ComputeJacobian().
	 *
	 * Requires ALL elements to have radTEnergyHysteresisMaterial.
	 * Converges in 2-4 iterations (vs hundreds for Hantila/Picard).
	 *
	 * @param PrecOnMagnetiz Convergence tolerance (||F||/||M||)
	 * @param MaxIterNumber Maximum Newton iterations
	 * @param MagnResetIsNotNeeded If 0, reset magnetization before solving
	 * @return Number of Newton iterations performed
	 */
	int AutoRelax_BInput_Newton(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded = 0);

	/**
	 * B-input Hantila solver for energy-based hysteresis materials.
	 *
	 * Uses constant LHS (I - alpha*N), LU-factored ONCE.
	 * Each iteration: O(N^2) back-substitution + Inverse(B) per element.
	 * More iterations than Newton, but cheaper per iteration for large N.
	 *
	 * @param PrecOnMagnetiz Convergence tolerance (max|dB|/B_sat)
	 * @param MaxIterNumber Maximum Hantila iterations
	 * @param alpha Polarization parameter (>= max dM/dH). 0 = auto-compute.
	 * @param relax Under-relaxation (0 = full step)
	 * @param MagnResetIsNotNeeded If 0, reset magnetization before solving
	 * @return Number of iterations performed
	 */
	int AutoRelax_BInput_Hantila(double PrecOnMagnetiz, int MaxIterNumber,
	                             double alpha = 0.0, double relax = 0.0,
	                             char MagnResetIsNotNeeded = 0);

protected:
	/**
	 * Build system matrix and RHS, then solve the linear system.
	 *
	 * @param ctx Nonlinear context with BaseMatrix, CurrentChiArray, FlatExtern
	 * @param iterCount Current nonlinear iteration number (0 for first)
	 * @return Number of linear iterations (0 for direct solvers like LU)
	 */
	virtual int SolveLinearStep(NonlinearContext& ctx, int iterCount) { return 0; }

	/**
	 * Solve the B-input Newton linear step: J_F * dM = -F
	 * where J_F = I - block_diag(dJ/dB) * (N + I)
	 *
	 * @param ctx Nonlinear context
	 * @param NpI Interaction matrix + identity (column-major, dof x dof)
	 * @param dJdB_blocks Block-diagonal dJ/dB (n_elem * 9 doubles, column-major 3x3)
	 * @param F Residual vector (dof)
	 * @param dM Output: Newton step (dof)
	 * @return 0 on success
	 */
	virtual int SolveBInputLinearStep(NonlinearContext& ctx,
	                                  const std::vector<double>& NpI,
	                                  const std::vector<double>& dJdB_blocks,
	                                  const std::vector<double>& F,
	                                  std::vector<double>& dM) { return 0; }

	/**
	 * Whether this solver requires the dense base matrix (ctx.BaseMatrix).
	 * LU and BiCGSTAB need it.
	 * Override to return false for matrix-free solvers.
	 */
	virtual bool NeedsDenseMatrix() const { return true; }
};

//-------------------------------------------------------------------------
// Note: Legacy relaxation methods (Methods 1-8) have been removed
// Newton-style M(H) update is now integrated into LU and BiCGSTAB solvers
//-------------------------------------------------------------------------


/**
 * Direct solver using LU decomposition (Gaussian elimination with partial pivoting)
 * Method number 0
 *
 * This solver is required for tetrahedral elements with high permeability materials
 * where iterative relaxation methods diverge due to large interaction coefficients.
 *
 * The system to solve is:
 *   (I - chi*N) * M = chi * H_ext
 * where N is the interaction matrix (demagnetization coefficients)
 */
class radTRelaxationMethNo_0 : public radTIterativeRelaxMeth {

public:
	radTRelaxationMethNo_0(radTInteraction* InInteractionPtr)
	: radTIterativeRelaxMeth(InInteractionPtr)
	{
		IntrctPtr = InInteractionPtr;
	}

	~radTRelaxationMethNo_0() {}

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

protected:
	// Override: LU direct solver for linear step
	int SolveLinearStep(NonlinearContext& ctx, int iterCount) override;

	// Override: LU reads the dense interaction/base matrix.
	bool NeedsDenseMatrix() const override;

	// Override: Dense Jacobian assembly + LAPACK dgesv_ for B-input Newton
	int SolveBInputLinearStep(NonlinearContext& ctx,
	                          const std::vector<double>& NpI,
	                          const std::vector<double>& dJdB_blocks,
	                          const std::vector<double>& F,
	                          std::vector<double>& dM) override;

private:
	// Solve linear system Ax=b using LU decomposition with partial pivoting
	// A is modified in place (contains LU factors after call)
	// b is overwritten with solution x
	// Returns 0 on success, non-zero on failure (singular matrix)
	int SolveLU(std::vector<std::vector<double>>& A, std::vector<double>& b, int n);

	// Flat matrix version for variable DOF
	int SolveLU_Flat(std::vector<double>& A, std::vector<double>& b, int n);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

/* Non-symmetric BiCGSTAB / HACApK relaxation methods retired. */

//-------------------------------------------------------------------------

#endif
