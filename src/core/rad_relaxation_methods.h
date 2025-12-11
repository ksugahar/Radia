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
	constexpr int LU         = 0;  // LU direct solver
	constexpr int BICGSTAB   = 1;  // BiCGSTAB iterative solver (default)
	constexpr int IMPLICIT_SS = 2;  // Implicit Successive Substitution (Gauss-Seidel)
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
	{
		IntrctPtr = InInteractionPtr;
	}

	~radTRelaxationMethNo_1() {}

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

	// Variable DOF version for hybrid MSC + standard element analysis
	int AutoRelax_VariableDOF(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

private:
	// BiCGSTAB iterative solver
	// Solves: A*x = b using BiCGSTAB with Jacobi preconditioner
	// Returns number of iterations (0 on failure)
	int SolveBiCGSTAB(int ndof, double tol, int max_iter, double& residual);

	// Variable DOF version of BiCGSTAB
	int SolveBiCGSTAB_VariableDOF(int totalDOF, double tol, int max_iter, double& residual);

	// Matrix-vector product
	// Computes: y = A * x where A = -N + diag(1/chi)
	void MatVec(const std::vector<double>& x, std::vector<double>& y,
	            const std::vector<double>& inv_chi, int ndof);

	// Dense matrix-vector product
	// IMPORTANT: Uses pre-computed inv_chi values that are FIXED during BiCGSTAB iterations
	void DenseMatVec(const std::vector<double>& x, std::vector<double>& y,
	                 const std::vector<double>& inv_chi, int ndof);

	// Variable DOF matrix-vector product using flat storage
	void MatVec_VariableDOF(const std::vector<double>& x, std::vector<double>& y,
	                        const std::vector<double>& inv_chi, int totalDOF);

	// Get diagonal elements for Jacobi preconditioner
	// IMPORTANT: Uses pre-computed inv_chi values that are FIXED during BiCGSTAB iterations
	void GetDiagonalElements(std::vector<double>& diag, const std::vector<double>& inv_chi, int n_elem);

	// Variable DOF version
	void GetDiagonalElements_VariableDOF(std::vector<double>& diag, const std::vector<double>& inv_chi, int totalDOF);

	// BLAS-like operations
	double Dot(const std::vector<double>& a, const std::vector<double>& b, int n);
	double Norm2(const std::vector<double>& a, int n);
	void Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n);
	void Copy(const std::vector<double>& src, std::vector<double>& dst, int n);
	void Scale(double alpha, std::vector<double>& x, int n);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

/**
 * Implicit Successive Substitution (Gauss-Seidel) solver
 * Method number 2
 *
 * This solver implements the original Radia Method No. 3 algorithm:
 * - Sequential element updates (Gauss-Seidel style)
 * - For each element i: solve (I - chi*N[i][i]) * H[i] = chi*Mr + QuasiExtField[i]
 * - Immediately update M[i] = M(H[i]) before processing next element
 *
 * Benefits over BiCGSTAB:
 * - No need to build or store full system matrix
 * - O(N^2) per iteration instead of O(N^2) matvec per BiCGSTAB iteration
 * - Better convergence for nonlinear materials
 * - Under-relaxation support for stability
 *
 * Reference: Original Radia radTRelaxationMethNo_3 (ochubar/Radia)
 */
class radTRelaxationMethNo_2 : public radTIterativeRelaxMeth {

public:
	radTRelaxationMethNo_2(radTInteraction* InInteractionPtr)
	: radTIterativeRelaxMeth(InInteractionPtr), InstMisfitM(1.0e30), Omega(0.3)
	{
		IntrctPtr = InInteractionPtr;
	}

	~radTRelaxationMethNo_2() {}

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

	// Variable DOF version for hybrid MSC + standard element analysis
	int AutoRelax_VariableDOF(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

	// Set under-relaxation parameter (default 0.3, range [0.1, 1.0])
	void SetOmega(double omega) { Omega = (omega > 0.0 && omega <= 1.0) ? omega : 0.3; }

private:
	double InstMisfitM;  // Instantaneous misfit for convergence tracking
	double Omega;        // Under-relaxation parameter (default 0.3)

	// Perform one Gauss-Seidel sweep over all elements
	void DefineNewMagnetizations();

	// 3x3 matrix inversion
	int Matrix3d_inv(const double A[3][3], double Ainv[3][3]);
};

//-------------------------------------------------------------------------

#endif
