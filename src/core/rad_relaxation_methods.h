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
#include "rad_hmatrix_aca.h"
#include <memory>

//-------------------------------------------------------------------------
// Solver method constants
//-------------------------------------------------------------------------

namespace RadSolverMethod {
	constexpr int LU        = 0;   // LU direct solver
	constexpr int BICGSTAB  = 1;  // BiCGSTAB iterative solver (default)
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

private:
	// Solve linear system Ax=b using LU decomposition with partial pivoting
	// A is modified in place (contains LU factors after call)
	// b is overwritten with solution x
	// Returns 0 on success, non-zero on failure (singular matrix)
	int SolveLU(std::vector<std::vector<double>>& A, std::vector<double>& b, int n);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

/**
 * BiCGSTAB iterative solver with H-matrix acceleration
 * Method number 1
 *
 * This solver uses BiCGSTAB (Biconjugate Gradient Stabilized) with
 * optional H-matrix (HACApK) for fast matrix-vector products, providing:
 * - O(N log N) per iteration instead of O(N^2) when H-matrix is enabled
 * - Jacobi (diagonal) preconditioning for faster convergence
 * - Stable for high permeability materials
 *
 * H-matrix is controlled via rad.SolverHMatrixEnable() / rad.SolverHMatrixDisable()
 *
 * Recommended for N > 100 elements where direct solver (Method 0)
 * becomes too slow due to O(N^3) complexity.
 *
 * Reference: van der Vorst, SIAM J. Sci. Stat. Comput. 13 (1992)
 */
class radTRelaxationMethNo_1 : public radTIterativeRelaxMeth {

public:
	radTRelaxationMethNo_1(radTInteraction* InInteractionPtr)
	: radTIterativeRelaxMeth(InInteractionPtr),
	  m_hmatrix(nullptr),
	  m_hmatrix_initialized(false)
	{
		IntrctPtr = InInteractionPtr;
	}

	~radTRelaxationMethNo_1() {}

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

private:
	// Initialize H-matrix if enabled
	bool InitializeHMatrix();

	// BiCGSTAB iterative solver
	// Solves: A*x = b using BiCGSTAB with Jacobi preconditioner
	// Returns number of iterations (0 on failure)
	int SolveBiCGSTAB(int ndof, double tol, int max_iter, double& residual);

	// Matrix-vector product (uses H-matrix if enabled, else dense)
	// Computes: y = A * x where A = -N + diag(1/chi)
	void MatVec(const std::vector<double>& x, std::vector<double>& y,
	            const std::vector<double>& inv_chi, int ndof);

	// Dense matrix-vector product
	// IMPORTANT: Uses pre-computed inv_chi values that are FIXED during BiCGSTAB iterations
	void DenseMatVec(const std::vector<double>& x, std::vector<double>& y,
	                 const std::vector<double>& inv_chi, int ndof);

	// H-matrix accelerated matrix-vector product
	// Uses H-matrix for -N*x, then adds (1/chi)*x
	void HMatrixMatVec(const std::vector<double>& x, std::vector<double>& y,
	                   const std::vector<double>& inv_chi, int ndof);

	// Get diagonal elements for Jacobi preconditioner
	// IMPORTANT: Uses pre-computed inv_chi values that are FIXED during BiCGSTAB iterations
	void GetDiagonalElements(std::vector<double>& diag, const std::vector<double>& inv_chi, int n_elem);

	// BLAS-like operations
	double Dot(const std::vector<double>& a, const std::vector<double>& b, int n);
	double Norm2(const std::vector<double>& a, int n);
	void Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n);
	void Copy(const std::vector<double>& src, std::vector<double>& dst, int n);
	void Scale(double alpha, std::vector<double>& x, int n);

	// H-matrix member
	std::unique_ptr<radTHMatrixACA> m_hmatrix;
	bool m_hmatrix_initialized;
};

//-------------------------------------------------------------------------

#endif
