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

// HACApK support (Method 2) - currently disabled, requires C/C++ interface fix
// #define RADIA_USE_HACAPK 1
#ifdef RADIA_USE_HACAPK
#include "rad_hacapk.h"
#endif

//-------------------------------------------------------------------------
// Solver method constants
//-------------------------------------------------------------------------

namespace RadSolverMethod {
	constexpr int LU         = 0;  // LU direct solver
	constexpr int BICGSTAB   = 1;  // BiCGSTAB iterative solver (default)
	constexpr int BICGSTAB_HMATRIX = 2;  // BiCGSTAB with H-matrix (HACApK ACA+)
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

/**
 * Shared context for nonlinear magnetostatic iteration.
 *
 * This struct encapsulates all state vectors and parameters that are
 * IDENTICAL across LU, BiCGSTAB, and HACApK solvers. By extracting
 * this common state, we eliminate ~1200 lines of duplicated code.
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
	bool last_solve_was_moment_hacapk;  // True when the last linear step used the moment method-2 path
	double last_moment_linear_tol;      // Effective tolerance used by the last moment Krylov solve
	int last_moment_krylov_solver;      // 0=BiCGSTAB, 1=GMRES

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

	// Multipole-moment outer acceleration history (safeguarded Anderson depth 1).
	bool moment_anderson_have_prev;
	std::vector<double> MomentAndersonPrevResidual;
	std::vector<double> MomentAndersonPrevImage;
	int moment_anderson_accepted;
	int moment_anderson_rejected;

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
		, last_solve_was_moment_hacapk(false)
		, last_moment_linear_tol(0.0)
		, last_moment_krylov_solver(0)
		, relax_param(0.0)
		, use_newton(false)
		, newton_damping_enabled(false)
		, newton_ls_max_iter(5)
		, newton_ls_min_omega(0.01)
		, total_ls_backtracks(0)
		, use_b_input(false)
		, moment_anderson_have_prev(false)
		, moment_anderson_accepted(0)
		, moment_anderson_rejected(0)
	{}
};

//-------------------------------------------------------------------------
// Helper functions for unified nonlinear iteration
// These functions encapsulate common logic shared across LU/BiCGSTAB/HACApK
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
 * Build base matrix with correct sign convention for MMM/MSC elements.
 * MMM (3DOF): negate interaction matrix (stores N, need -N)
 * MSC (6DOF): use as-is (stores -K/(4pi))
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
	 * This is the ONLY part that differs between LU, BiCGSTAB, and HACApK.
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

	// Override: the all-moment MMM/MSC Picard path assembles its own per-element moment system and never
	// reads the dense interaction/base matrix, so no dense BaseMatrix is needed (storage decoupling).
	// Returns true for the genuine dense LU / B-input Newton-Hantila paths (which DO read the dense BaseMatrix).
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

/**
 * BiCGSTAB iterative solver
 * Method number 1
 *
 * This solver uses BiCGSTAB (Biconjugate Gradient Stabilized) with
 * Jacobi (diagonal) preconditioning for faster convergence.
 * Stable for high permeability materials.
 *
 * Recommended for N > 100 elements where direct solver (Method 0)
 * becomes too slow due to O(N^3) complexity.
 *
 * Reference: van der Vorst, SIAM J. Sci. Stat. Comput. 13 (1992)
 */
class radTRelaxationMethNo_1 : public radTIterativeRelaxMeth {

public:
	radTRelaxationMethNo_1(radTInteraction* InInteractionPtr)
	: radTIterativeRelaxMeth(InInteractionPtr)
	{
		IntrctPtr = InInteractionPtr;
	}

	~radTRelaxationMethNo_1() {}

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

protected:
	// Override: BiCGSTAB iterative solver for linear step
	int SolveLinearStep(NonlinearContext& ctx, int iterCount) override;

	// Moment method-1 builds the parameter-free moment dense matrix inside
	// SolveLinearStep, so it does not need the legacy FlatInteract matrix.
	bool NeedsDenseMatrix() const override;

private:
	// Variable DOF version of BiCGSTAB
	// elemChiArray: chi for system matrix diagonal (chi_abs for Picard, chi_d for Newton)
	int SolveBiCGSTAB_VariableDOF(NonlinearContext& ctx,
	                              int totalDOF, double tol, int max_iter, double& residual,
	                              const std::vector<double>& elemChiArray,
	                              bool use_newton = false,
	                              const std::vector<double>* absChiArray = nullptr,
	                              const double* oldSigma = nullptr);

	// Matrix-vector product
	// Computes: y = A * x where A = -N - diag(1/chi) (ELF-compatible)
	void MatVec(const std::vector<double>& x, std::vector<double>& y,
	            const std::vector<double>& inv_chi, int ndof);

	// Dense matrix-vector product
	// IMPORTANT: Uses pre-computed inv_chi values that are FIXED during BiCGSTAB iterations
	void DenseMatVec(const std::vector<double>& x, std::vector<double>& y,
	                 const std::vector<double>& inv_chi, int ndof);

	// Build flat matrix for BLAS dgemv (row-major order)
	// Matrix A = -N - diag(1/chi) stored in column-major for BLAS (ELF-compatible)
	void BuildFlatMatrix(std::vector<double>& A_flat, const std::vector<double>& inv_chi, int ndof);

	// Dense matrix-vector product using BLAS dgemv
	void DenseMatVec_BLAS(const std::vector<double>& A_flat, const std::vector<double>& x,
	                      std::vector<double>& y, int ndof);

	// Variable DOF matrix-vector product using flat storage
	void MatVec_VariableDOF(const std::vector<double>& x, std::vector<double>& y,
	                        const std::vector<double>& inv_chi, int totalDOF);

	// Get diagonal elements for Jacobi preconditioner
	// IMPORTANT: Uses pre-computed inv_chi values that are FIXED during BiCGSTAB iterations
	void GetDiagonalElements(std::vector<double>& diag, const std::vector<double>& inv_chi, int n_elem);

	// Variable DOF version
	void GetDiagonalElements_VariableDOF(std::vector<double>& diag, const std::vector<double>& inv_chi, int totalDOF);

	// Block Jacobi preconditioner (for ill-conditioned matrices)
	bool BuildBlockJacobiPreconditioner_VariableDOF(std::vector<double>& blockInverse, std::vector<int>& blockOffsets,
	                                                 const std::vector<double>& inv_chi, int totalDOF);
	void ApplyBlockJacobiPreconditioner_VariableDOF(const std::vector<double>& x, std::vector<double>& y,
	                                                 const std::vector<double>& blockInverse, const std::vector<int>& blockOffsets);

	// BLAS-like operations
	double Dot(const std::vector<double>& a, const std::vector<double>& b, int n);
	double Norm2(const std::vector<double>& a, int n);
	void Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n);
	void Copy(const std::vector<double>& src, std::vector<double>& dst, int n);
	void Scale(double alpha, std::vector<double>& x, int n);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

#ifdef RADIA_USE_HACAPK
/**
 * BiCGSTAB iterative solver with H-matrix acceleration (HACApK ACA+)
 * Method number 2
 *
 * Uses HACApK library for O(N log N) matrix-vector products via ACA+ compression.
 * Effective when elements are spatially well-separated (multiple objects).
 * For single compact objects, Method 1 (dense BiCGSTAB) may be faster.
 *
 * Reference:
 *   - HACApK: ppOpen-HPC project (MIT License)
 *   - ACA+: Bebendorf & Rjasanow, Computing 70 (2003)
 */
class radTRelaxationMethNo_2 : public radTIterativeRelaxMeth {

public:
	radTRelaxationMethNo_2(radTInteraction* InInteractionPtr)
	: radTIterativeRelaxMeth(InInteractionPtr)
	, m_hacapk(nullptr)
	{
		IntrctPtr = InInteractionPtr;
	}

	~radTRelaxationMethNo_2() {
		if (m_hacapk) {
			delete m_hacapk;
			m_hacapk = nullptr;
		}
	}

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

	// Variable DOF version for 6DOF MSC hexahedra
	int AutoRelax_VariableDOF(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

	// Get H-matrix statistics (for debugging/analysis)
	const RadHACApKStats& GetHMatrixStats() const {
		static RadHACApKStats empty_stats;
		return m_hacapk ? m_hacapk->GetStats() : empty_stats;
	}

	// Set H-matrix parameters (call before AutoRelax)
	void SetHACApKParams(const RadHACApKParams& params) { m_hacapk_params = params; }

private:
	// HACApK manager (owns the H-matrix)
	RadHACApKMMMManager* m_hacapk;
	RadHACApKParams m_hacapk_params;

	// NOTE: Jacobi preconditioner is now recomputed every iteration (FIX 2025-12-27)
	// No longer cached. See SolveBiCGSTAB_HMatrix_VariableDOF for details.

	// BiCGSTAB with H-matrix for 6DOF MSC hexahedra
	// Returns number of iterations (0 on failure)
	// elemChiArray: chi for system matrix diagonal (chi_abs for Picard, chi_d for Newton)
	// Newton params (optional): diffChiArray=differential chi, absChiArray=absolute chi, oldSigma=prev solution
	int SolveBiCGSTAB_HMatrix_VariableDOF(NonlinearContext& ctx,
	                                       int totalDOF, double tol, int max_iter, double& residual,
	                                       const std::vector<double>& elemChiArray,
	                                       bool use_newton = false,
	                                       const std::vector<double>* absChiArray = nullptr,
	                                       const double* oldSigma = nullptr);

	// Matrix-vector product using H-matrix for 6DOF MSC hexahedra
	void MatVec_HMatrix_VariableDOF(const std::vector<double>& x, std::vector<double>& y, int totalDOF);

	// Get diagonal elements using H-matrix for 6DOF MSC hexahedra
	void GetDiagonalElements_HMatrix_VariableDOF(std::vector<double>& diag,
	                                              const std::vector<double>& inv_chi, int totalDOF);

	// Block Jacobi preconditioner for H-matrix BiCGSTAB
	bool BuildBlockJacobiPreconditioner_HMatrix(std::vector<double>& blockInverse,
	                                             std::vector<int>& blockOffsets,
	                                             const std::vector<double>& inv_chi, int totalDOF);
	void ApplyBlockJacobiPreconditioner_HMatrix(const std::vector<double>& x, std::vector<double>& y,
	                                             const std::vector<double>& blockInverse,
	                                             const std::vector<int>& blockOffsets);

	// BLAS-like operations (same as radTRelaxationMethNo_1)
	double Dot(const std::vector<double>& a, const std::vector<double>& b, int n);
	double Norm2(const std::vector<double>& a, int n);
	void Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n);
	void Copy(const std::vector<double>& src, std::vector<double>& dst, int n);
	void Scale(double alpha, std::vector<double>& x, int n);
};
#endif // RADIA_USE_HACAPK

//-------------------------------------------------------------------------

#endif
