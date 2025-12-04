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

//-------------------------------------------------------------------------
// Solver method constants
//-------------------------------------------------------------------------

namespace RadSolverMethod {
	constexpr int NEWTON    = 8;   // Newton-Raphson for nonlinear materials
	constexpr int LU        = 9;   // LU direct solver
	constexpr int BICGSTAB  = 10;  // BiCGSTAB iterative solver (default)
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
// Note: radTSimpleRelaxation (Method 1) has been removed (deprecated)
// Note: radTRelaxationMethNo_2 (Method 2) has been removed (deprecated)
// Note: radTRelaxationMethNo_3 (Method 3) has been removed (deprecated)
// Note: radTRelaxationMethNo_4 (Method 4) has been removed (deprecated)
// Note: radTRelaxationMethNo_a5 (Method 5) has been removed (deprecated)
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Note: radTRelaxationMethNo_6 (Method 6) has been removed (deprecated)
// Note: radTRelaxationMethNo_7 (Method 7) has been removed (deprecated)
//-------------------------------------------------------------------------

/**
 * Newton-Raphson iterative solver for nonlinear materials
 * Method number 8
 *
 * This method uses Newton-Raphson iteration to solve nonlinear magnetic
 * problems with saturable materials (MatSatIso, MatSatIsoTab, etc.)
 *
 * Required for:
 * - Nonlinear (saturable) materials where M = M(H) is nonlinear
 * - Uses local Jacobian matrix for each element
 */
class radTRelaxationMethNo_8 : public radTIterativeRelaxMeth {

	double mInstMisfitMe2, mDesiredPrecOnMagnetizE2;

public:

	radTRelaxationMethNo_8(radTInteraction* InInteractionPtr) : radTIterativeRelaxMeth(InInteractionPtr)
	{
		double DesiredPrecOnMagnetiz = 1.E-03;
		mDesiredPrecOnMagnetizE2 = DesiredPrecOnMagnetiz*DesiredPrecOnMagnetiz;
		IntrctPtr = InInteractionPtr;
		mInstMisfitMe2 = 1.E+23;
	}

	~radTRelaxationMethNo_8() {}

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);
	void DefineNewMagnetizations();
};

//-------------------------------------------------------------------------

/**
 * Direct solver using LU decomposition (Gaussian elimination with partial pivoting)
 * Method number 9
 *
 * This solver is required for tetrahedral elements with high permeability materials
 * where iterative relaxation methods diverge due to large interaction coefficients.
 *
 * The system to solve is:
 *   (I - chi*N) * M = chi * H_ext
 * where N is the interaction matrix (demagnetization coefficients)
 */
class radTRelaxationMethNo_9 : public radTIterativeRelaxMeth {

public:
	radTRelaxationMethNo_9(radTInteraction* InInteractionPtr)
	: radTIterativeRelaxMeth(InInteractionPtr)
	{
		IntrctPtr = InInteractionPtr;
	}

	~radTRelaxationMethNo_9() {}

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
 * Method number 10
 *
 * This solver uses BiCGSTAB (Biconjugate Gradient Stabilized) with
 * H-matrix (HACApK) for fast matrix-vector products, providing:
 * - O(N log N) per iteration instead of O(N^2)
 * - Jacobi (diagonal) preconditioning for faster convergence
 * - Stable for high permeability materials
 *
 * Recommended for N > 100 elements where direct solver (Method 9)
 * becomes too slow due to O(N^3) complexity.
 *
 * Reference: van der Vorst, SIAM J. Sci. Stat. Comput. 13 (1992)
 */
class radTRelaxationMethNo_10 : public radTIterativeRelaxMeth {

public:
	radTRelaxationMethNo_10(radTInteraction* InInteractionPtr)
	: radTIterativeRelaxMeth(InInteractionPtr)
	{
		IntrctPtr = InInteractionPtr;
	}

	~radTRelaxationMethNo_10() {}

	int AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded=0);

private:
	// BiCGSTAB iterative solver
	// Solves: A*x = b using BiCGSTAB with Jacobi preconditioner
	// Returns number of iterations (0 on failure)
	int SolveBiCGSTAB(int ndof, double tol, int max_iter, double& residual);

	// Dense matrix-vector product
	void DenseMatVec(const std::vector<double>& x, std::vector<double>& y, int ndof);

	// Get diagonal elements for Jacobi preconditioner
	void GetDiagonalElements(std::vector<double>& diag, int n_elem);

	// BLAS-like operations
	double Dot(const std::vector<double>& a, const std::vector<double>& b, int n);
	double Norm2(const std::vector<double>& a, int n);
	void Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n);
	void Copy(const std::vector<double>& src, std::vector<double>& dst, int n);
	void Scale(double alpha, std::vector<double>& x, int n);
};

//-------------------------------------------------------------------------

#endif
