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

	// Variable DOF version for hybrid MSC + standard element analysis
	int AutoRelax_VariableDOF(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

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
	, m_cachedDOF(0)
	, m_cachedDOFType(0)
	{
		IntrctPtr = InInteractionPtr;
	}

	~radTRelaxationMethNo_1() {}

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

	// Variable DOF version for hybrid MSC + standard element analysis
	int AutoRelax_VariableDOF(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

private:
	// Variable DOF version of BiCGSTAB
	// elemChiArray: isotropic chi for each element (3DOF elements use this, 6DOF uses poly->CurrentChi)
	int SolveBiCGSTAB_VariableDOF(int totalDOF, double tol, int max_iter, double& residual,
	                              const std::vector<double>& elemChiArray);

	// Matrix-vector product
	// Computes: y = A * x where A = -N + diag(1/chi)
	void MatVec(const std::vector<double>& x, std::vector<double>& y,
	            const std::vector<double>& inv_chi, int ndof);

	// Dense matrix-vector product
	// IMPORTANT: Uses pre-computed inv_chi values that are FIXED during BiCGSTAB iterations
	void DenseMatVec(const std::vector<double>& x, std::vector<double>& y,
	                 const std::vector<double>& inv_chi, int ndof);

	// Build flat matrix for BLAS dgemv (row-major order)
	// Matrix A = -N + diag(1/chi) stored in column-major for BLAS
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

	// BLAS-like operations (Intel MKL cblas)
	double Dot(const std::vector<double>& a, const std::vector<double>& b, int n);
	double Norm2(const std::vector<double>& a, int n);
	void Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n);
	void Copy(const std::vector<double>& src, std::vector<double>& dst, int n);
	void Scale(double alpha, std::vector<double>& x, int n);

	// Cached work vectors to avoid repeated allocation (performance optimization)
	// These are resized only when totalDOF changes
	int m_cachedDOF;
	std::vector<double> m_r, m_r0, m_p, m_v, m_s, m_t;
	std::vector<double> m_p_hat, m_s_hat, m_diag_inv;
	std::vector<double> m_inv_chi, m_rhs, m_sol;

	// Cached DOF type for fast MatVec path selection
	// 0 = not initialized, 1 = pure 3DOF (tetra), 2 = pure 6DOF (hex MSC), 3 = mixed
	int m_cachedDOFType;

	// Ensure work vectors are sized for given DOF
	void EnsureWorkVectors(int totalDOF);
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
	RadHACApKManager* m_hacapk;
	RadHACApKParams m_hacapk_params;

	// NOTE: Jacobi preconditioner is now recomputed every iteration (FIX 2025-12-27)
	// No longer cached. See SolveBiCGSTAB_HMatrix_VariableDOF for details.

	// BiCGSTAB with H-matrix for 6DOF MSC hexahedra
	// Returns number of iterations (0 on failure)
	// elemChiArray: isotropic chi for each element (3DOF elements use this, 6DOF uses poly->CurrentChi)
	int SolveBiCGSTAB_HMatrix_VariableDOF(int totalDOF, double tol, int max_iter, double& residual,
	                                       const std::vector<double>& elemChiArray);

	// Matrix-vector product using H-matrix for 6DOF MSC hexahedra
	void MatVec_HMatrix_VariableDOF(const std::vector<double>& x, std::vector<double>& y, int totalDOF);

	// Get diagonal elements using H-matrix for 6DOF MSC hexahedra
	void GetDiagonalElements_HMatrix_VariableDOF(std::vector<double>& diag,
	                                              const std::vector<double>& inv_chi, int totalDOF);

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
