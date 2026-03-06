/*-------------------------------------------------------------------------
*
* File name:      radrlmet.cpp
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

#include "rad_relaxation_methods.h"
#include "rad_yield.h"
#include "rad_polyhedron.h"  // For radTPolyhedron in variable DOF solver
#include "rad_material_def.h"  // For radTNonlinearIsotropMaterial::GetHfromM
#include "rad_application.h"  // For radTApplication::NonlinearMethod

#include <time.h>
#include <chrono>   // For timing instrumentation
#include <cstring>  // For std::memcpy
#include <cstdio>   // For fprintf in debug logging
#include <cstdlib>  // For getenv

// External access to radTApplication for NonlinearMethod setting
extern radTApplication rad;

//-------------------------------------------------------------------------

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_LAPACK
#include "cblas.h"
extern "C" {
    // LAPACK dgesv: Solve A*x = b using LU factorization with partial pivoting
    // Parameters:
    //   n: order of matrix A
    //   nrhs: number of right-hand sides
    //   A: matrix (n x n), overwritten with LU factors
    //   lda: leading dimension of A
    //   ipiv: pivot indices
    //   b: right-hand side, overwritten with solution
    //   ldb: leading dimension of b
    //   info: 0 = success, < 0 = illegal arg, > 0 = singular
    void dgesv_(int* n, int* nrhs, double* A, int* lda, int* ipiv, double* b, int* ldb, int* info);
}
#endif

extern radTYield radYield;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTIterativeRelaxMeth::ComputeRelaxStatusParam(const TVector3d* NewMagnArray, const TVector3d* OldMagnArray, const TVector3d* NewFieldArray)
{
	double BufMisfitM, BufMaxModM, BufMaxModH, TestBufMaxModM, TestBufMaxModH;
	BufMisfitM=0.;
	BufMaxModM=BufMaxModH=TestBufMaxModM=TestBufMaxModH=1.E-17;
	TVector3d Mnew_mi_MoldVect;

	radTRelaxStatusParam& RelStatParR = IntrctPtr->RelaxStatusParam;

	#pragma omp parallel for reduction(+:BufMisfitM) if(IntrctPtr->AmOfMainElem > 100)
	for(int i=0; i<IntrctPtr->AmOfMainElem; i++)
	{
		double LocalTestBufMaxModM = 0., LocalTestBufMaxModH = 0.;
		if(RelStatParR.MisfitM >= 0. && OldMagnArray != nullptr)
		{
			Mnew_mi_MoldVect = NewMagnArray[i] - OldMagnArray[i];
			BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y
						+ Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;
		}
		if(RelStatParR.MaxModM >= 0.)
		{
			LocalTestBufMaxModM = sqrt(NewMagnArray[i].x*NewMagnArray[i].x
								+ NewMagnArray[i].y*NewMagnArray[i].y
								+ NewMagnArray[i].z*NewMagnArray[i].z);
			#pragma omp critical
			{
				if(LocalTestBufMaxModM > BufMaxModM) BufMaxModM = LocalTestBufMaxModM;
			}
		}
		if(RelStatParR.MaxModH >= 0.)
		{
			LocalTestBufMaxModH = sqrt(NewFieldArray[i].x*NewFieldArray[i].x
								+ NewFieldArray[i].y*NewFieldArray[i].y
								+ NewFieldArray[i].z*NewFieldArray[i].z);
			#pragma omp critical
			{
				if(LocalTestBufMaxModH > BufMaxModH) BufMaxModH = LocalTestBufMaxModH;
			}
		}
	}
	if(RelStatParR.MisfitM >= 0.) RelStatParR.MisfitM = sqrt(BufMisfitM/IntrctPtr->AmOfMainElem);
	if(RelStatParR.MaxModM >= 0.) RelStatParR.MaxModM = BufMaxModM;
	if(RelStatParR.MaxModH >= 0.) RelStatParR.MaxModH = BufMaxModH;
}

//-------------------------------------------------------------------------

void radTIterativeRelaxMeth::MakeN_iter(int IterNum)
{
	for(int i=0; i<(IterNum-1); i++)
	{
		DefineNewMagnetizations(); 

		if(radYield.Check()==0) return; // To allow multitasking on Mac: consider better places for this
	}

	//radTSend Send;
	std::vector<TVector3d> vOldMagnArray(IntrctPtr->AmOfMainElem);
	TVector3d* OldMagnArray = vOldMagnArray.data();

	for(int k=0; k<IntrctPtr->AmOfMainElem; k++) OldMagnArray[k] = (IntrctPtr->g3dRelaxPtrVect[k])->Magn;
	DefineNewMagnetizations();
	for(int q=0; q<IntrctPtr->AmOfMainElem; q++)
		IntrctPtr->NewMagnArray[q] = (IntrctPtr->g3dRelaxPtrVect[q])->Magn;

	ComputeRelaxStatusParam(IntrctPtr->NewMagnArray, OldMagnArray, IntrctPtr->NewFieldArray);
}

//-------------------------------------------------------------------------
// Note: Legacy relaxation methods (Methods 1-8) have been removed
// Newton-style M(H) update is now integrated into LU and BiCGSTAB solvers
//-------------------------------------------------------------------------

//=========================================================================
// Method 0: LU Direct Solver
//=========================================================================

int radTRelaxationMethNo_0::SolveLU(std::vector<std::vector<double>>& A, std::vector<double>& b, int n)
{
#ifdef HAVE_LAPACK
	// Use LAPACK dgesv for optimized LU decomposition (multi-threaded, SIMD optimized)
	// LAPACK uses column-major (Fortran) ordering, so we need to transpose
	// Or we can use the fact that solving A*x=b is equivalent to solving A^T*x=b
	// when we store A in row-major order and pass it to LAPACK as column-major

	// Create contiguous column-major array for LAPACK
	std::vector<double> A_col(n * n);
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			A_col[j * n + i] = A[i][j];  // Column-major: A_col[j][i] = A[i][j]
		}
	}

	std::vector<int> ipiv(n);
	int nrhs = 1;
	int info = 0;

	dgesv_(&n, &nrhs, A_col.data(), &n, ipiv.data(), b.data(), &n, &info);

	return (info == 0) ? 0 : -1;

#else
	// Fallback: Gaussian elimination with partial pivoting
	// Forward elimination
	for(int k = 0; k < n - 1; k++)
	{
		// Find pivot
		int maxRow = k;
		double maxVal = std::abs(A[k][k]);
		for(int i = k + 1; i < n; i++)
		{
			if(std::abs(A[i][k]) > maxVal)
			{
				maxVal = std::abs(A[i][k]);
				maxRow = i;
			}
		}

		// Check for singular matrix
		if(maxVal < 1.0e-15)
		{
			return -1;  // Singular matrix
		}

		// Swap rows if needed
		if(maxRow != k)
		{
			std::swap(A[k], A[maxRow]);
			std::swap(b[k], b[maxRow]);
		}

		// Eliminate below pivot
		for(int i = k + 1; i < n; i++)
		{
			double factor = A[i][k] / A[k][k];
			A[i][k] = 0.0;
			for(int j = k + 1; j < n; j++)
			{
				A[i][j] -= factor * A[k][j];
			}
			b[i] -= factor * b[k];
		}
	}

	// Check last diagonal element
	if(std::abs(A[n-1][n-1]) < 1.0e-15)
	{
		return -1;  // Singular matrix
	}

	// Back substitution
	for(int i = n - 1; i >= 0; i--)
	{
		double sum = b[i];
		for(int j = i + 1; j < n; j++)
		{
			sum -= A[i][j] * b[j];
		}
		b[i] = sum / A[i][i];
	}

	return 0;  // Success
#endif
}

int radTRelaxationMethNo_0::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	// Debug: write to log file to trace execution path
	FILE* path_log = std::fopen("S:/Radia/01_GitHub/radia_path.log", "a");
	if(path_log)
	{
		std::fprintf(path_log, "AutoRelax called: HasVariableDOF=%d, AmOfMainElem=%d\n",
		            IntrctPtr->HasVariableDOF() ? 1 : 0, IntrctPtr->AmOfMainElem);
		std::fflush(path_log);
		std::fclose(path_log);
	}

	// Check if variable DOF is active (e.g., 6 DOF MSC hexahedra)
	// If so, use the variable DOF solver for better convergence
	if(IntrctPtr->HasVariableDOF())
	{
		return AutoRelax_VariableDOF(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
	}

	// Reset magnetization if needed
	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	// Check pointers and sizes
	if(AmOfMainElem <= 0) return 0;
	if(IntrctPtr->InteractMatrix == nullptr) return 0;
	if(IntrctPtr->NewMagnArray == nullptr) return 0;
	if(IntrctPtr->ExternFieldArray == nullptr) return 0;
	if(IntrctPtr->NewFieldArray == nullptr) return 0;
	if(IntrctPtr->g3dRelaxPtrVect.empty()) return 0;
	int ndof = 3 * AmOfMainElem;  // 3 DOF per element (Mx, My, Mz)

	// Get access to interaction matrix and field arrays
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;

	// Store old magnetization and H field for convergence check
	std::vector<TVector3d> OldMagnArray(AmOfMainElem);
	std::vector<TVector3d> OldFieldArray(AmOfMainElem);  // H field from previous iteration

	double PrecOnMagnetizE2 = PrecOnMagnetiz * PrecOnMagnetiz;
	double MisfitE2 = 1.0e30;
	int iterCount = 0;

	// Initialize H field for first iteration
	// For nonlinear materials, we need a reasonable initial estimate of H
	// Using H = 0 gives chi=0 (at table origin), using H_ext can exceed saturation
	// Start with a moderate H value (1000 A/m) which is typically in the linear region
	const double H_init_mag = 1000.0;  // Initial H magnitude in A/m
	for(int i = 0; i < AmOfMainElem; i++)
	{
		// Use direction from external field but moderate magnitude
		double H_ext_mag = std::sqrt(ExternFieldAr[i].x*ExternFieldAr[i].x +
		                             ExternFieldAr[i].y*ExternFieldAr[i].y +
		                             ExternFieldAr[i].z*ExternFieldAr[i].z);
		if(H_ext_mag > 1.0e-10)
		{
			double scale = std::min(1.0, H_init_mag / H_ext_mag);
			NewFieldAr[i].x = ExternFieldAr[i].x * scale;
			NewFieldAr[i].y = ExternFieldAr[i].y * scale;
			NewFieldAr[i].z = ExternFieldAr[i].z * scale;
		}
		else
		{
			NewFieldAr[i] = ExternFieldAr[i];
		}
	}

	// Build base matrix (geometric part without chi) once - this is the -N matrix
	// OPTIMIZATION (2025-12-11): Use flat column-major array for LAPACK and fast copy
	// Only diagonal elements are updated each iteration
	//
	// MATRIX LAYOUT FIX (2025-12-24):
	// InteractMatrix[i][j] stores TMatrix3df where:
	//   Str0 = dH/dMx = (dHx/dMx, dHy/dMx, dHz/dMx) <- COLUMN vector (response to Mx)
	//   Str1 = dH/dMy = (dHx/dMy, dHy/dMy, dHz/dMy) <- COLUMN vector (response to My)
	//   Str2 = dH/dMz = (dHx/dMz, dHy/dMz, dHz/dMz) <- COLUMN vector (response to Mz)
	//
	// For row k of the 3x3 block, we need N[i][j] element (k, l) = dH_k/dM_l
	// This requires: A[3i+k, 3j+l] = Str_l.{xyz}[k] (transposed access)
	std::vector<double> BaseMatrix_flat(ndof * ndof, 0.0);
	#pragma omp parallel for if(AmOfMainElem > 50)
	for(int i = 0; i < AmOfMainElem; i++)
	{
		for(int j = 0; j < AmOfMainElem; j++)
		{
			TMatrix3df& Nij = IntrcMat[i][j];
			int row_base = 3 * i;
			int col_base = 3 * j;

			// Column-major storage for LAPACK: A[row + col*lda]
			// Row 0 (Hx response): A[3i+0, 3j+l] = -dHx/dM_l = -Str_l.x
			BaseMatrix_flat[(row_base + 0) + (col_base + 0)*ndof] = -Nij.Str0.x;  // -dHx/dMx
			BaseMatrix_flat[(row_base + 0) + (col_base + 1)*ndof] = -Nij.Str1.x;  // -dHx/dMy
			BaseMatrix_flat[(row_base + 0) + (col_base + 2)*ndof] = -Nij.Str2.x;  // -dHx/dMz

			// Row 1 (Hy response): A[3i+1, 3j+l] = -dHy/dM_l = -Str_l.y
			BaseMatrix_flat[(row_base + 1) + (col_base + 0)*ndof] = -Nij.Str0.y;  // -dHy/dMx
			BaseMatrix_flat[(row_base + 1) + (col_base + 1)*ndof] = -Nij.Str1.y;  // -dHy/dMy
			BaseMatrix_flat[(row_base + 1) + (col_base + 2)*ndof] = -Nij.Str2.y;  // -dHy/dMz

			// Row 2 (Hz response): A[3i+2, 3j+l] = -dHz/dM_l = -Str_l.z
			BaseMatrix_flat[(row_base + 2) + (col_base + 0)*ndof] = -Nij.Str0.z;  // -dHz/dMx
			BaseMatrix_flat[(row_base + 2) + (col_base + 1)*ndof] = -Nij.Str1.z;  // -dHz/dMy
			BaseMatrix_flat[(row_base + 2) + (col_base + 2)*ndof] = -Nij.Str2.z;  // -dHz/dMz
		}
	}

	// Store diagonal of base matrix (needed for efficient diagonal update)
	std::vector<double> BaseMatrix_diag(ndof);
	for(int i = 0; i < ndof; i++)
	{
		BaseMatrix_diag[i] = BaseMatrix_flat[i + i*ndof];
	}

	// Outer nonlinear iteration loop
	// For linear materials, this converges in 1 iteration
	// For nonlinear materials, chi(H) is updated each iteration
	std::vector<double> SystemMatrix_flat(ndof * ndof, 0.0);
	std::vector<double> RHS(ndof, 0.0);
	std::vector<int> ipiv(ndof);

	for(iterCount = 0; iterCount < MaxIterNumber; iterCount++)
	{
		// Store old magnetization and H field (for dB calculation)
		for(int i = 0; i < AmOfMainElem; i++)
		{
			OldMagnArray[i] = MagnAr[i];
			OldFieldArray[i] = NewFieldAr[i];  // Store H_old
		}

		// Update H field from current M using constitutive relation: H = M / chi
		// This approach is O(N) instead of O(N^2) for H = H_ext - N*M
		// Skip on first iteration (H already initialized above)
		if(iterCount > 0)
		{
			for(int i = 0; i < AmOfMainElem; i++)
			{
				radTg3dRelax* g3dRelaxPtr_i = IntrctPtr->g3dRelaxPtrVect[i];
				radTMaterial* MaterPtr_i = (radTMaterial*)(g3dRelaxPtr_i->MaterHandle.rep);

				// Get current chi value from previous iteration's H
				TMatrix3d KsiTensor;
				TVector3d MrVect;
				MaterPtr_i->DefineInstantKsiTensor(NewFieldAr[i], KsiTensor, MrVect);

				// Use average chi (isotropic approximation for H update)
				double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
				if(chi < 1.0e-6) chi = 1.0e-6;

				// H_int = M / chi (from M = chi * H)
				NewFieldAr[i].x = MagnAr[i].x / chi;
				NewFieldAr[i].y = MagnAr[i].y / chi;
				NewFieldAr[i].z = MagnAr[i].z / chi;
			}
		}

		// OPTIMIZATION (2025-12-11): Use BLAS dcopy for fast matrix copy
		// Copy base matrix to system matrix (flat arrays)
#ifdef HAVE_LAPACK
		cblas_dcopy(ndof * ndof, BaseMatrix_flat.data(), 1, SystemMatrix_flat.data(), 1);
#else
		std::copy(BaseMatrix_flat.begin(), BaseMatrix_flat.end(), SystemMatrix_flat.begin());
#endif

		// Update diagonal and RHS using current chi(H)
		// OPTIMIZATION: Store 1/chi values for diagonal update
		std::vector<double> inv_chi_diag(ndof);

		for(int i = 0; i < AmOfMainElem; i++)
		{
			radTg3dRelax* g3dRelaxPtr_i = IntrctPtr->g3dRelaxPtrVect[i];
			if(g3dRelaxPtr_i == nullptr) return 0;
			radTMaterial* MaterPtr_i = (radTMaterial*)(g3dRelaxPtr_i->MaterHandle.rep);
			if(MaterPtr_i == nullptr) return 0;

			// Get susceptibility tensor at current H (nonlinear: chi depends on H)
			TVector3d& InstH = NewFieldAr[i];  // Use current H field
			TMatrix3d KsiTensor;
			TVector3d MrVect;
			MaterPtr_i->DefineInstantKsiTensor(InstH, KsiTensor, MrVect);

			// Compute 1/chi values for this element
			double inv_chi_x = (KsiTensor.Str0.x > 1.0e-10) ? (1.0 / KsiTensor.Str0.x) : 1.0e10;
			double inv_chi_y = (KsiTensor.Str1.y > 1.0e-10) ? (1.0 / KsiTensor.Str1.y) : 1.0e10;
			double inv_chi_z = (KsiTensor.Str2.z > 1.0e-10) ? (1.0 / KsiTensor.Str2.z) : 1.0e10;

			inv_chi_diag[3*i + 0] = inv_chi_x;
			inv_chi_diag[3*i + 1] = inv_chi_y;
			inv_chi_diag[3*i + 2] = inv_chi_z;

			// Update diagonal: A[k,k] = BaseMatrix[k,k] + 1/chi[k]
			int k0 = 3*i + 0;
			int k1 = 3*i + 1;
			int k2 = 3*i + 2;
			SystemMatrix_flat[k0 + k0*ndof] = BaseMatrix_diag[k0] + inv_chi_x;
			SystemMatrix_flat[k1 + k1*ndof] = BaseMatrix_diag[k1] + inv_chi_y;
			SystemMatrix_flat[k2 + k2*ndof] = BaseMatrix_diag[k2] + inv_chi_z;

			// RHS = H_ext + Mr/chi
			double Mr_over_chi_x = (KsiTensor.Str0.x > 1.0e-10) ? (MrVect.x / KsiTensor.Str0.x) : 0.0;
			double Mr_over_chi_y = (KsiTensor.Str1.y > 1.0e-10) ? (MrVect.y / KsiTensor.Str1.y) : 0.0;
			double Mr_over_chi_z = (KsiTensor.Str2.z > 1.0e-10) ? (MrVect.z / KsiTensor.Str2.z) : 0.0;

			RHS[3*i + 0] = ExternFieldAr[i].x + Mr_over_chi_x;
			RHS[3*i + 1] = ExternFieldAr[i].y + Mr_over_chi_y;
			RHS[3*i + 2] = ExternFieldAr[i].z + Mr_over_chi_z;
		}

		// Solve the linear system using LAPACK dgesv directly
#ifdef HAVE_LAPACK
		int n = ndof;
		int nrhs = 1;
		int info = 0;
		dgesv_(&n, &nrhs, SystemMatrix_flat.data(), &n, ipiv.data(), RHS.data(), &n, &info);
		int ierr = info;
#else
		// Fallback to old SolveLU (requires conversion to 2D array)
		std::vector<std::vector<double>> SystemMatrix_2d(ndof, std::vector<double>(ndof));
		for(int row = 0; row < ndof; row++)
		{
			for(int col = 0; col < ndof; col++)
			{
				SystemMatrix_2d[row][col] = SystemMatrix_flat[row + col*ndof];
			}
		}
		int ierr = SolveLU(SystemMatrix_2d, RHS, ndof);
#endif

		if(ierr != 0)
		{
			// Solver failed - singular matrix
			return iterCount;
		}

		// Extract LU solution (M values) - this is the new magnetization from linearized system
		for(int i = 0; i < AmOfMainElem; i++)
		{
			MagnAr[i].x = RHS[3 * i + 0];
			MagnAr[i].y = RHS[3 * i + 1];
			MagnAr[i].z = RHS[3 * i + 2];
		}

		// Pure Newton-Raphson iteration:
		// Use the LU solution directly - do NOT apply Gauss-Seidel M(H) correction!
		// The Gauss-Seidel M(H) update was causing oscillations and non-convergence
		// for high-permeability nonlinear materials because it mixed two iteration schemes.

		// Compute convergence using relative change ||dM||/||M||
		// This is the standard Radia convergence criterion
		double M_diff_sq = 0.0;
		double M_norm_sq = 0.0;
		for(int i = 0; i < AmOfMainElem; i++)
		{
			// dM = M_new - M_old
			TVector3d dM;
			dM.x = MagnAr[i].x - OldMagnArray[i].x;
			dM.y = MagnAr[i].y - OldMagnArray[i].y;
			dM.z = MagnAr[i].z - OldMagnArray[i].z;

			M_diff_sq += dM.AmpE2();
			M_norm_sq += MagnAr[i].AmpE2();

			// Update the object's magnetization
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
			g3dRelaxPtr->Magn = MagnAr[i];
		}

		// Relative change: ||dM|| / ||M||
		double rel_change = 0.0;
		if(M_norm_sq > 1.0e-30)
		{
			rel_change = std::sqrt(M_diff_sq / M_norm_sq);
		}
		else
		{
			rel_change = std::sqrt(M_diff_sq);
		}
		MisfitE2 = rel_change * rel_change;  // For compatibility with status reporting

		// Check convergence using relative tolerance (PrecOnMagnetiz is the relative tolerance)
		if(rel_change <= PrecOnMagnetiz)
		{
			iterCount++;
			break;
		}

		// Allow multitasking
		if(radYield.Check() == 0) return iterCount;
	}

	// Update relaxation status
	IntrctPtr->RelaxStatusParam.MisfitM = std::sqrt(MisfitE2);
	ComputeRelaxStatusParam(MagnAr, OldMagnArray.data(), NewFieldAr);

	return iterCount;
}

//=========================================================================
// Method 1: BiCGSTAB Iterative Solver (default)
//=========================================================================

double radTRelaxationMethNo_1::Dot(const std::vector<double>& a, const std::vector<double>& b, int n)
{
#ifdef HAVE_LAPACK
	// Use OpenBLAS cblas_ddot for optimized dot product
	return cblas_ddot(n, a.data(), 1, b.data(), 1);
#else
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum) if(n > 100)
	for(int i = 0; i < n; i++)
	{
		sum += a[i] * b[i];
	}
	return sum;
#endif
}

double radTRelaxationMethNo_1::Norm2(const std::vector<double>& a, int n)
{
#ifdef HAVE_LAPACK
	// Use OpenBLAS cblas_dnrm2 for optimized norm
	return cblas_dnrm2(n, a.data(), 1);
#else
	return std::sqrt(Dot(a, a, n));
#endif
}

void radTRelaxationMethNo_1::Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n)
{
#ifdef HAVE_LAPACK
	// Use OpenBLAS cblas_daxpy: y = alpha*x + y
	cblas_daxpy(n, alpha, x.data(), 1, y.data(), 1);
#else
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		y[i] += alpha * x[i];
	}
#endif
}

void radTRelaxationMethNo_1::Copy(const std::vector<double>& src, std::vector<double>& dst, int n)
{
#ifdef HAVE_LAPACK
	// Use OpenBLAS cblas_dcopy for optimized copy
	cblas_dcopy(n, src.data(), 1, dst.data(), 1);
#else
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		dst[i] = src[i];
	}
#endif
}

void radTRelaxationMethNo_1::Scale(double alpha, std::vector<double>& x, int n)
{
#ifdef HAVE_LAPACK
	// Use OpenBLAS cblas_dscal for optimized scale
	cblas_dscal(n, alpha, x.data(), 1);
#else
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		x[i] *= alpha;
	}
#endif
}

void radTRelaxationMethNo_1::GetDiagonalElements(std::vector<double>& diag, const std::vector<double>& inv_chi, int n_elem)
{
	// Extract diagonal elements from interaction matrix for Jacobi preconditioner
	// Diagonal block [i][i] is a 3x3 matrix, we extract the diagonal of that
	// CRITICAL FIX: Use pre-computed 1/chi values that are FIXED for this BiCGSTAB solve
	// (chi is only updated in the outer nonlinear iteration loop)
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;

	for(int i = 0; i < n_elem; i++)
	{
		// Use pre-computed 1/chi values
		double inv_chi_x = inv_chi[3*i + 0];
		double inv_chi_y = inv_chi[3*i + 1];
		double inv_chi_z = inv_chi[3*i + 2];

		// Diagonal of system matrix: A = -N + 1/chi
		if(IntrcMat != nullptr)
		{
			// Nii is from InteractMatrix[i][i]
			TMatrix3df& Nii = IntrcMat[i][i];
			diag[3*i + 0] = -Nii.Str0.x + inv_chi_x;
			diag[3*i + 1] = -Nii.Str1.y + inv_chi_y;
			diag[3*i + 2] = -Nii.Str2.z + inv_chi_z;
		}
		else
		{
			// Fallback: just use 1/chi as diagonal (no N contribution)
			diag[3*i + 0] = inv_chi_x;
			diag[3*i + 1] = inv_chi_y;
			diag[3*i + 2] = inv_chi_z;
		}
	}
}

void radTRelaxationMethNo_1::DenseMatVec(const std::vector<double>& x, std::vector<double>& y,
                                         const std::vector<double>& inv_chi, int ndof)
{
	// Computes y = A * x where A = -N + 1/chi
	// Uses dense matrix-vector product
	// CRITICAL FIX: Use pre-computed 1/chi values that are FIXED for this BiCGSTAB solve
	// (chi is only updated in the outer nonlinear iteration loop)
	//
	// MATRIX LAYOUT FIX (2025-12-24):
	// InteractMatrix[i][j] stores TMatrix3df where:
	//   Str0 = dH/dMx = (dHx/dMx, dHy/dMx, dHz/dMx) <- COLUMN vector (response to Mx)
	//   Str1 = dH/dMy = (dHx/dMy, dHy/dMy, dHz/dMy) <- COLUMN vector (response to My)
	//   Str2 = dH/dMz = (dHx/dMz, dHy/dMz, dHz/dMz) <- COLUMN vector (response to Mz)
	//
	// For matrix-vector product H = N * M, we need:
	//   Hx = dHx/dMx*Mx + dHx/dMy*My + dHx/dMz*Mz = Str0.x*Mx + Str1.x*My + Str2.x*Mz
	//   Hy = dHy/dMx*Mx + dHy/dMy*My + dHy/dMz*Mz = Str0.y*Mx + Str1.y*My + Str2.y*Mz
	//   Hz = dHz/dMx*Mx + dHz/dMy*My + dHz/dMz*Mz = Str0.z*Mx + Str1.z*My + Str2.z*Mz

	int n_elem = ndof / 3;
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;

	// Initialize y to zero
	std::fill(y.begin(), y.end(), 0.0);

	if(IntrcMat != nullptr)
	{
		// Dense matrix-vector product: O(N^2)
		#pragma omp parallel for if(n_elem > 50)
		for(int i = 0; i < n_elem; i++)
		{
			// Use pre-computed 1/chi values
			double inv_chi_x = inv_chi[3*i + 0];
			double inv_chi_y = inv_chi[3*i + 1];
			double inv_chi_z = inv_chi[3*i + 2];

			// y[i] = -sum(N[i][j] * x[j]) + (1/chi) * x[i]
			double y0 = inv_chi_x * x[3*i + 0];
			double y1 = inv_chi_y * x[3*i + 1];
			double y2 = inv_chi_z * x[3*i + 2];

			for(int j = 0; j < n_elem; j++)
			{
				TMatrix3df& Nij = IntrcMat[i][j];
				double xj0 = x[3*j + 0];  // Mx
				double xj1 = x[3*j + 1];  // My
				double xj2 = x[3*j + 2];  // Mz

				// Correct transpose: use column k of N for output component k
				// y0 (Hx) -= dHx/dMx*Mx + dHx/dMy*My + dHx/dMz*Mz
				// y1 (Hy) -= dHy/dMx*Mx + dHy/dMy*My + dHy/dMz*Mz
				// y2 (Hz) -= dHz/dMx*Mx + dHz/dMy*My + dHz/dMz*Mz
				y0 -= Nij.Str0.x*xj0 + Nij.Str1.x*xj1 + Nij.Str2.x*xj2;
				y1 -= Nij.Str0.y*xj0 + Nij.Str1.y*xj1 + Nij.Str2.y*xj2;
				y2 -= Nij.Str0.z*xj0 + Nij.Str1.z*xj1 + Nij.Str2.z*xj2;
			}

			y[3*i + 0] = y0;
			y[3*i + 1] = y1;
			y[3*i + 2] = y2;
		}
	}
}

void radTRelaxationMethNo_1::BuildFlatMatrix(std::vector<double>& A_flat, const std::vector<double>& inv_chi, int ndof)
{
	// Build flat matrix A = -N + diag(1/chi) for BLAS dgemv
	// Stored in column-major order for BLAS compatibility
	//
	// MATRIX LAYOUT FIX (2025-12-24):
	// InteractMatrix[i][j] stores TMatrix3df where:
	//   Str0 = dH/dMx = (dHx/dMx, dHy/dMx, dHz/dMx) <- COLUMN vector (response to Mx)
	//   Str1 = dH/dMy = (dHx/dMy, dHy/dMy, dHz/dMy) <- COLUMN vector (response to My)
	//   Str2 = dH/dMz = (dHx/dMz, dHy/dMz, dHz/dMz) <- COLUMN vector (response to Mz)
	//
	// For row k of the 3x3 block, we need N[i][j] element (k, l) = dH_k/dM_l
	// This requires: A[3i+k, 3j+l] = Str_l.{xyz}[k] (transposed access)

	int n_elem = ndof / 3;
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;

	A_flat.resize(ndof * ndof, 0.0);

	if(IntrcMat != nullptr)
	{
		#pragma omp parallel for if(n_elem > 50)
		for(int i = 0; i < n_elem; i++)
		{
			for(int j = 0; j < n_elem; j++)
			{
				TMatrix3df& Nij = IntrcMat[i][j];

				// Column-major: A[row + col*lda]
				int row_base = 3*i;
				int col_base = 3*j;

				// Row 0 (Hx response): A[3i+0, 3j+l] = -dHx/dM_l = -Str_l.x
				A_flat[(row_base + 0) + (col_base + 0)*ndof] = -Nij.Str0.x;  // -dHx/dMx
				A_flat[(row_base + 0) + (col_base + 1)*ndof] = -Nij.Str1.x;  // -dHx/dMy
				A_flat[(row_base + 0) + (col_base + 2)*ndof] = -Nij.Str2.x;  // -dHx/dMz

				// Row 1 (Hy response): A[3i+1, 3j+l] = -dHy/dM_l = -Str_l.y
				A_flat[(row_base + 1) + (col_base + 0)*ndof] = -Nij.Str0.y;  // -dHy/dMx
				A_flat[(row_base + 1) + (col_base + 1)*ndof] = -Nij.Str1.y;  // -dHy/dMy
				A_flat[(row_base + 1) + (col_base + 2)*ndof] = -Nij.Str2.y;  // -dHy/dMz

				// Row 2 (Hz response): A[3i+2, 3j+l] = -dHz/dM_l = -Str_l.z
				A_flat[(row_base + 2) + (col_base + 0)*ndof] = -Nij.Str0.z;  // -dHz/dMx
				A_flat[(row_base + 2) + (col_base + 1)*ndof] = -Nij.Str1.z;  // -dHz/dMy
				A_flat[(row_base + 2) + (col_base + 2)*ndof] = -Nij.Str2.z;  // -dHz/dMz
			}

			// Add diagonal 1/chi terms
			A_flat[(3*i + 0) + (3*i + 0)*ndof] += inv_chi[3*i + 0];
			A_flat[(3*i + 1) + (3*i + 1)*ndof] += inv_chi[3*i + 1];
			A_flat[(3*i + 2) + (3*i + 2)*ndof] += inv_chi[3*i + 2];
		}
	}
}

void radTRelaxationMethNo_1::DenseMatVec_BLAS(const std::vector<double>& A_flat, const std::vector<double>& x,
                                              std::vector<double>& y, int ndof)
{
#ifdef HAVE_LAPACK
	// Use OpenBLAS cblas_dgemv: y = alpha*A*x + beta*y
	// A is stored in column-major order
	cblas_dgemv(CblasColMajor, CblasNoTrans, ndof, ndof,
	            1.0, A_flat.data(), ndof, x.data(), 1,
	            0.0, y.data(), 1);
#else
	// Fallback: manual matrix-vector multiply
	std::fill(y.begin(), y.end(), 0.0);
	for(int i = 0; i < ndof; i++)
	{
		for(int j = 0; j < ndof; j++)
		{
			y[i] += A_flat[i + j*ndof] * x[j];
		}
	}
#endif
}

void radTRelaxationMethNo_1::MatVec(const std::vector<double>& x, std::vector<double>& y,
                                    const std::vector<double>& inv_chi, int ndof)
{
	// Use dense matvec (H-matrix support removed 2025-12-06)
	DenseMatVec(x, y, inv_chi, ndof);
}

int radTRelaxationMethNo_1::SolveBiCGSTAB(int ndof, double tol, int max_iter, double& residual)
{
	// BiCGSTAB with Jacobi preconditioner
	// Reference: van der Vorst, SIAM J. Sci. Stat. Comput. 13 (1992)
	//
	// CRITICAL FIX:
	// - Pre-compute 1/chi values ONCE at the start of this solve
	// - Use these FIXED values throughout all BiCGSTAB iterations
	// - chi is only updated in the OUTER nonlinear iteration loop
	// - This ensures the linear system being solved is consistent
	//
	// OPTIMIZATION (2025-12-11):
	// - Build flat matrix ONCE for BLAS dgemv acceleration
	// - Use OpenBLAS cblas_dgemv for matrix-vector products

	int n_elem = ndof / 3;

	// Allocate work vectors
	std::vector<double> r(ndof), r0(ndof), p(ndof), v(ndof), s(ndof), t(ndof);
	std::vector<double> p_hat(ndof), s_hat(ndof), diag_inv(ndof);

	// Pre-compute 1/chi values for ALL elements ONCE
	// These values are FIXED for the entire BiCGSTAB solve
	std::vector<double> inv_chi(ndof);
	std::vector<double> rhs(ndof);

	if(IntrctPtr->ExternFieldArray == nullptr) return 0;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;

	for(int i = 0; i < n_elem; i++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
		radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		// Use current H for chi(H) computation (from outer loop)
		TVector3d InstH = (NewFieldAr != nullptr) ? NewFieldAr[i] : TVector3d(0., 0., 0.);
		TMatrix3d KsiTensor;
		TVector3d MrVect;
		MaterPtr->DefineInstantKsiTensor(InstH, KsiTensor, MrVect);

		// Pre-compute 1/chi values (FIXED for this solve)
		inv_chi[3*i + 0] = (KsiTensor.Str0.x > 1.0e-10) ? 1.0/KsiTensor.Str0.x : 1.0e10;
		inv_chi[3*i + 1] = (KsiTensor.Str1.y > 1.0e-10) ? 1.0/KsiTensor.Str1.y : 1.0e10;
		inv_chi[3*i + 2] = (KsiTensor.Str2.z > 1.0e-10) ? 1.0/KsiTensor.Str2.z : 1.0e10;

		// Compute RHS: b = H_ext + Mr/chi (using pre-computed 1/chi)
		double inv_chi_for_rhs_x = (KsiTensor.Str0.x > 1.0e-10) ? 1.0/KsiTensor.Str0.x : 0.0;
		double inv_chi_for_rhs_y = (KsiTensor.Str1.y > 1.0e-10) ? 1.0/KsiTensor.Str1.y : 0.0;
		double inv_chi_for_rhs_z = (KsiTensor.Str2.z > 1.0e-10) ? 1.0/KsiTensor.Str2.z : 0.0;

		rhs[3*i + 0] = IntrctPtr->ExternFieldArray[i].x + MrVect.x * inv_chi_for_rhs_x;
		rhs[3*i + 1] = IntrctPtr->ExternFieldArray[i].y + MrVect.y * inv_chi_for_rhs_y;
		rhs[3*i + 2] = IntrctPtr->ExternFieldArray[i].z + MrVect.z * inv_chi_for_rhs_z;
	}

	// Get initial guess (current magnetization)
	std::vector<double> sol(ndof);
	if(IntrctPtr->NewMagnArray == nullptr) return 0;
	for(int i = 0; i < n_elem; i++)
	{
		sol[3*i + 0] = IntrctPtr->NewMagnArray[i].x;
		sol[3*i + 1] = IntrctPtr->NewMagnArray[i].y;
		sol[3*i + 2] = IntrctPtr->NewMagnArray[i].z;
	}

	// Build Jacobi preconditioner: M^{-1} = diag(A)^{-1}
	// Uses pre-computed inv_chi values
	GetDiagonalElements(diag_inv, inv_chi, n_elem);
	for(int i = 0; i < ndof; i++)
	{
		if(std::abs(diag_inv[i]) > 1.0e-15)
		{
			diag_inv[i] = 1.0 / diag_inv[i];
		}
		else
		{
			diag_inv[i] = 1.0;  // Fallback for near-zero diagonal
		}
	}

#ifdef HAVE_LAPACK
	// Build flat matrix ONCE for BLAS dgemv (significant speedup for large problems)
	std::vector<double> A_flat;
	BuildFlatMatrix(A_flat, inv_chi, ndof);

	// Initialize: r0 = b - A*x0
	DenseMatVec_BLAS(A_flat, sol, v, ndof);  // v = A*x0 (uses BLAS dgemv)
	Copy(rhs, r, ndof);                       // r = rhs
	Axpy(-1.0, v, r, ndof);                   // r = r - v
#else
	// Initialize: r0 = b - A*x0
	MatVec(sol, v, inv_chi, ndof);  // v = A*x0 (uses pre-computed inv_chi)
	Copy(rhs, r, ndof);                   // r = rhs
	Axpy(-1.0, v, r, ndof);               // r = r - v
#endif

	// Choose r0* = r0
	Copy(r, r0, ndof);

	// Initialize BiCGSTAB parameters
	double rho = 1.0, alpha_bicg = 1.0, omega = 1.0;
	std::fill(p.begin(), p.end(), 0.0);
	std::fill(v.begin(), v.end(), 0.0);

	// Compute ||b|| for relative residual
	double rhs_norm = Norm2(rhs, ndof);
	if(rhs_norm < 1.0e-30) rhs_norm = 1.0;

	int iter;
	for(iter = 1; iter <= max_iter; iter++)
	{
		double rho_old = rho;
		rho = Dot(r0, r, ndof);

		// Check for breakdown
		if(std::abs(rho) < 1.0e-30)
		{
			residual = Norm2(r, ndof) / rhs_norm;
			break;
		}

		if(iter == 1)
		{
			Copy(r, p, ndof);
		}
		else
		{
			if(std::abs(rho_old * omega) < 1.0e-30)
			{
				residual = Norm2(r, ndof) / rhs_norm;
				break;
			}
			double beta = (rho / rho_old) * (alpha_bicg / omega);
			Axpy(-omega, v, p, ndof);
			Scale(beta, p, ndof);
			Axpy(1.0, r, p, ndof);
		}

		// Apply preconditioner: p_hat = M^{-1} * p
		#pragma omp parallel for if(ndof > 100)
		for(int i = 0; i < ndof; i++)
		{
			p_hat[i] = diag_inv[i] * p[i];
		}

		// v = A * p_hat
#ifdef HAVE_LAPACK
		DenseMatVec_BLAS(A_flat, p_hat, v, ndof);
#else
		MatVec(p_hat, v, inv_chi, ndof);
#endif

		// alpha_bicg = rho / (r0, v)
		double r0_dot_v = Dot(r0, v, ndof);
		if(std::abs(r0_dot_v) < 1.0e-30)
		{
			residual = Norm2(r, ndof) / rhs_norm;
			break;
		}
		alpha_bicg = rho / r0_dot_v;

		// s = r - alpha_bicg * v
		Copy(r, s, ndof);
		Axpy(-alpha_bicg, v, s, ndof);

		// Check if s is small enough
		double s_norm = Norm2(s, ndof);
		if(s_norm / rhs_norm < tol)
		{
			Axpy(alpha_bicg, p_hat, sol, ndof);
			residual = s_norm / rhs_norm;
			break;
		}

		// Apply preconditioner: s_hat = M^{-1} * s
		#pragma omp parallel for if(ndof > 100)
		for(int i = 0; i < ndof; i++)
		{
			s_hat[i] = diag_inv[i] * s[i];
		}

		// t = A * s_hat
#ifdef HAVE_LAPACK
		DenseMatVec_BLAS(A_flat, s_hat, t, ndof);
#else
		MatVec(s_hat, t, inv_chi, ndof);
#endif

		// omega = (t, s) / (t, t)
		double t_dot_s = Dot(t, s, ndof);
		double t_dot_t = Dot(t, t, ndof);
		if(std::abs(t_dot_t) < 1.0e-30)
		{
			Axpy(alpha_bicg, p_hat, sol, ndof);
			residual = s_norm / rhs_norm;
			break;
		}
		omega = t_dot_s / t_dot_t;

		// x = x + alpha_bicg * p_hat + omega * s_hat
		Axpy(alpha_bicg, p_hat, sol, ndof);
		Axpy(omega, s_hat, sol, ndof);

		// r = s - omega * t
		Copy(s, r, ndof);
		Axpy(-omega, t, r, ndof);

		// Check convergence
		double r_norm = Norm2(r, ndof);
		residual = r_norm / rhs_norm;
		if(residual < tol)
		{
			break;
		}

		// Check for stagnation
		if(std::abs(omega) < 1.0e-30)
		{
			break;
		}
	}

	// Copy solution back to NewMagnArray
	for(int i = 0; i < n_elem; i++)
	{
		IntrctPtr->NewMagnArray[i].x = sol[3*i + 0];
		IntrctPtr->NewMagnArray[i].y = sol[3*i + 1];
		IntrctPtr->NewMagnArray[i].z = sol[3*i + 2];
	}

	return iter;
}

int radTRelaxationMethNo_1::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	// Check if variable DOF is active (e.g., 6 DOF MSC hexahedra)
	// If so, use the variable DOF solver for better convergence
	if(IntrctPtr->HasVariableDOF())
	{
		return AutoRelax_VariableDOF(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
	}

	// Reset magnetization if needed
	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	int AmOfMainElem = IntrctPtr->AmOfMainElem;
	if(AmOfMainElem <= 0) return 0;

	// Check required arrays
	if(IntrctPtr->InteractMatrix == nullptr) return 0;
	if(IntrctPtr->NewMagnArray == nullptr) return 0;
	if(IntrctPtr->ExternFieldArray == nullptr) return 0;
	if(IntrctPtr->NewFieldArray == nullptr) return 0;
	if(IntrctPtr->g3dRelaxPtrVect.empty()) return 0;

	int ndof = 3 * AmOfMainElem;
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;

	// Store old magnetization and H field for convergence checking
	std::vector<TVector3d> OldMagnArray(AmOfMainElem);
	std::vector<TVector3d> OldFieldArray(AmOfMainElem);  // H field from previous iteration

	double MisfitE2 = 1.0e30;
	int totalIterCount = 0;
	int outerIter = 0;

	// Initialize H field for first iteration
	// For nonlinear materials, we need a reasonable initial estimate of H
	// Using H = 0 gives chi=0 (at table origin), using H_ext can exceed saturation
	// Start with a moderate H value (1000 A/m) which is typically in the linear region
	const double H_init_mag = 1000.0;  // Initial H magnitude in A/m
	for(int i = 0; i < AmOfMainElem; i++)
	{
		// Use direction from external field but moderate magnitude
		double H_ext_mag = std::sqrt(ExternFieldAr[i].x*ExternFieldAr[i].x +
		                             ExternFieldAr[i].y*ExternFieldAr[i].y +
		                             ExternFieldAr[i].z*ExternFieldAr[i].z);
		if(H_ext_mag > 1.0e-10)
		{
			double scale = std::min(1.0, H_init_mag / H_ext_mag);
			NewFieldAr[i].x = ExternFieldAr[i].x * scale;
			NewFieldAr[i].y = ExternFieldAr[i].y * scale;
			NewFieldAr[i].z = ExternFieldAr[i].z * scale;
		}
		else
		{
			NewFieldAr[i] = ExternFieldAr[i];
		}
	}

	// Outer nonlinear iteration loop
	// For linear materials, this converges in 1 iteration
	// For nonlinear materials, chi(H) is updated each iteration
	// Pure Newton iteration without Gauss-Seidel M(H) correction

	for(outerIter = 0; outerIter < MaxIterNumber; outerIter++)
	{
		// Store old magnetization and H field (for dB calculation)
		for(int i = 0; i < AmOfMainElem; i++)
		{
			OldMagnArray[i] = MagnAr[i];
			OldFieldArray[i] = NewFieldAr[i];  // Store H_old
		}

		// Update H field from current M using constitutive relation: H = M / chi
		// This approach is O(N) instead of O(N^2) for H = H_ext - N*M
		// Skip on first iteration (H already initialized above)
		if(outerIter > 0)
		{
			for(int i = 0; i < AmOfMainElem; i++)
			{
				radTg3dRelax* g3dRelaxPtr_i = IntrctPtr->g3dRelaxPtrVect[i];
				radTMaterial* MaterPtr_i = (radTMaterial*)(g3dRelaxPtr_i->MaterHandle.rep);

				// Get current chi value from previous iteration's H
				TMatrix3d KsiTensor;
				TVector3d MrVect;
				MaterPtr_i->DefineInstantKsiTensor(NewFieldAr[i], KsiTensor, MrVect);

				// Use average chi (isotropic approximation for H update)
				double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
				if(chi < 1.0e-6) chi = 1.0e-6;

				// H_int = M / chi (from M = chi * H)
				NewFieldAr[i].x = MagnAr[i].x / chi;
				NewFieldAr[i].y = MagnAr[i].y / chi;
				NewFieldAr[i].z = MagnAr[i].z / chi;
			}
		}

		// Solve linear system using BiCGSTAB with current chi(H)
		// System: A*M = b where A = -N + 1/chi(H), b = H_ext + Mr/chi
		//
		// NOTE: max_iter for BiCGSTAB inner loop should be FIXED (not user's MaxIterNumber)
		// User's MaxIterNumber controls OUTER nonlinear iterations, not inner BiCGSTAB
		// bicg_tol is set via rad.SetBiCGSTABTol() Python API (default: 1e-4, ELF-compatible)
		const double bicg_tol = rad.m_bicg_tol;
		const int bicg_max_iter = 10000;  // Inner BiCGSTAB max iterations (fixed)
		double residual = 0.0;
		int n_iter = SolveBiCGSTAB(ndof, bicg_tol, bicg_max_iter, residual);
		totalIterCount += n_iter;

		// Pure Newton-Raphson iteration:
		// Use the BiCGSTAB solution directly - do NOT apply Gauss-Seidel M(H) correction!
		// The Gauss-Seidel M(H) update was causing oscillations and non-convergence
		// for high-permeability nonlinear materials because it mixed two iteration schemes.

		// Compute convergence using relative change ||dM||/||M||
		// This is the standard Radia convergence criterion
		double M_diff_sq = 0.0;
		double M_norm_sq = 0.0;
		for(int i = 0; i < AmOfMainElem; i++)
		{
			// dM = M_new - M_old
			TVector3d dM;
			dM.x = MagnAr[i].x - OldMagnArray[i].x;
			dM.y = MagnAr[i].y - OldMagnArray[i].y;
			dM.z = MagnAr[i].z - OldMagnArray[i].z;

			M_diff_sq += dM.AmpE2();
			M_norm_sq += MagnAr[i].AmpE2();

			// Update the object's magnetization
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
			g3dRelaxPtr->Magn = MagnAr[i];
		}

		// Relative change: ||dM|| / ||M||
		double rel_change = 0.0;
		if(M_norm_sq > 1.0e-30)
		{
			rel_change = std::sqrt(M_diff_sq / M_norm_sq);
		}
		else
		{
			rel_change = std::sqrt(M_diff_sq);
		}
		MisfitE2 = rel_change * rel_change;  // For compatibility with status reporting

		// Check convergence using relative tolerance (PrecOnMagnetiz is the relative tolerance)
		if(rel_change <= PrecOnMagnetiz)
		{
			outerIter++;
			break;
		}

		// Allow multitasking
		if(radYield.Check() == 0) return outerIter;
	}

	// Update relaxation status
	IntrctPtr->RelaxStatusParam.MisfitM = std::sqrt(MisfitE2);
	ComputeRelaxStatusParam(MagnAr, OldMagnArray.data(), NewFieldAr);

	// Return nonlinear iteration count (outerIter), not BiCGSTAB total (totalIterCount)
	// Return nonlinear iteration count for consistent benchmarking
	return outerIter;
}

//=========================================================================
// Variable DOF Solver Methods for Hybrid MSC + Standard Element Analysis
// Reference: Yano & Sugahara, "MMM with MSC", J. Magn. Soc. Jpn., 2023
//=========================================================================

//-------------------------------------------------------------------------
// LU solver with flat matrix for variable DOF
//-------------------------------------------------------------------------

int radTRelaxationMethNo_0::SolveLU_Flat(std::vector<double>& A, std::vector<double>& b, int n)
{
#ifdef HAVE_LAPACK
	// Use LAPACK dgesv for optimized LU decomposition
	// A is stored row-major, need to transpose for LAPACK column-major
	std::vector<double> A_col(n * n);
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			A_col[j * n + i] = A[i * n + j];  // Transpose: A_col[j][i] = A[i][j]
		}
	}

	std::vector<int> ipiv(n);
	int nrhs = 1;
	int info = 0;

	dgesv_(&n, &nrhs, A_col.data(), &n, ipiv.data(), b.data(), &n, &info);

	return (info == 0) ? 0 : -1;
#else
	// Fallback: Gaussian elimination with partial pivoting
	for(int k = 0; k < n - 1; k++)
	{
		// Find pivot
		int maxRow = k;
		double maxVal = std::abs(A[k * n + k]);
		for(int i = k + 1; i < n; i++)
		{
			if(std::abs(A[i * n + k]) > maxVal)
			{
				maxVal = std::abs(A[i * n + k]);
				maxRow = i;
			}
		}

		if(maxVal < 1.0e-15) return -1;  // Singular

		// Swap rows
		if(maxRow != k)
		{
			for(int j = 0; j < n; j++)
			{
				std::swap(A[k * n + j], A[maxRow * n + j]);
			}
			std::swap(b[k], b[maxRow]);
		}

		// Eliminate
		for(int i = k + 1; i < n; i++)
		{
			double factor = A[i * n + k] / A[k * n + k];
			A[i * n + k] = 0.0;
			for(int j = k + 1; j < n; j++)
			{
				A[i * n + j] -= factor * A[k * n + j];
			}
			b[i] -= factor * b[k];
		}
	}

	if(std::abs(A[(n-1) * n + (n-1)]) < 1.0e-15) return -1;

	// Back substitution
	for(int i = n - 1; i >= 0; i--)
	{
		double sum = b[i];
		for(int j = i + 1; j < n; j++)
		{
			sum -= A[i * n + j] * b[j];
		}
		b[i] = sum / A[i * n + i];
	}

	return 0;
#endif
}

int radTRelaxationMethNo_0::AutoRelax_VariableDOF(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	// Check if variable DOF is active
	if(!IntrctPtr->HasVariableDOF())
	{
		// Fall back to standard solver
		return AutoRelax(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
	}

	// Reset if needed
	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	int AmOfMainElem = IntrctPtr->AmOfMainElem;
	if(AmOfMainElem <= 0) return 0;

	int totalDOF = IntrctPtr->GetTotalDOF();
	if(totalDOF <= 0) return 0;

	// Get flat arrays
	double* FlatInteract = IntrctPtr->GetFlatInteractMatrix();
	double* FlatMagn = IntrctPtr->GetFlatMagnArray();
	double* FlatField = IntrctPtr->GetFlatFieldArray();
	double* FlatExtern = IntrctPtr->GetFlatExternFieldArray();

	if(FlatInteract == nullptr || FlatMagn == nullptr || FlatField == nullptr || FlatExtern == nullptr)
		return 0;

	// Build base matrix (geometric part without chi)
	// For 3 DOF elements (MMM): store -N (need to negate from interaction matrix)
	// For 6 DOF elements (MSC): interaction matrix already stores -K/(4pi), copy directly
	std::vector<double> BaseMatrix(totalDOF * totalDOF);

	// First copy interaction matrix
	for(int i = 0; i < totalDOF * totalDOF; i++)
	{
		BaseMatrix[i] = FlatInteract[i];
	}

	// Sign convention for matrix equation:
	// MMM (3 DOF): FlatInteract stores N, equation is (-N + I/chi)*M = H_ext ↁEnegate
	// MSC (6 DOF): FlatInteract stores -K/(4pi), equation is (-K/(4pi) + I/chi)*sigma = H_ext.n ↁEuse as-is
	//
	// Only negate 3x3 blocks (MMM). 6x6 blocks (MSC) already have correct sign.
	for(int row_elem = 0; row_elem < AmOfMainElem; row_elem++)
	{
		int dof_row = IntrctPtr->GetElementDOF(row_elem);
		int offset_row = IntrctPtr->GetElementDOFOffset(row_elem);

		for(int col_elem = 0; col_elem < AmOfMainElem; col_elem++)
		{
			int dof_col = IntrctPtr->GetElementDOF(col_elem);
			int offset_col = IntrctPtr->GetElementDOFOffset(col_elem);

			// Only negate blocks involving 3 DOF elements (MMM)
			// 6x6 blocks (both MSC) keep their sign (-K/(4pi))
			// IMPORTANT: BaseMatrix is stored in COLUMN-MAJOR format (copied from FlatInteract)
			// Element at (row, col) is at index [col * totalDOF + row]
			if(dof_row == 3 || dof_col == 3)
			{
				for(int i = 0; i < dof_row; i++)
				{
					for(int j = 0; j < dof_col; j++)
					{
						// Column-major: (row, col) = [col * totalDOF + row]
						int idx = (offset_col + j) * totalDOF + (offset_row + i);
						BaseMatrix[idx] = -BaseMatrix[idx];
					}
				}
			}
		}
	}

	// Store old values for convergence check
	std::vector<double> OldMagn(totalDOF);
	// Store old chi values for MSC convergence check (ELF mucal1 style)
	std::vector<double> OldChi(AmOfMainElem, 0.0);
	// Store old B-field norms for ELF mucal2 (B-field) convergence check
	std::vector<double> OldBnorm(AmOfMainElem, 0.0);

	double PrecOnMagnetizE2 = PrecOnMagnetiz * PrecOnMagnetiz;
	double MisfitE2 = 1.0e30;
	int iterCount = 0;

	// Initialize H field (for chi(H) computation in nonlinear iteration)
	// ELF uses H = 100 A/m for first iteration (mucal0 style)
	const double H_init_mag = 100.0;
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		// For MSC elements (6 DOF), initialize differently than standard (3 DOF)
		if(dof == 3)
		{
			double H_ext_mag = 0.0;
			for(int k = 0; k < 3; k++)
			{
				H_ext_mag += FlatExtern[offset + k] * FlatExtern[offset + k];
			}
			H_ext_mag = std::sqrt(H_ext_mag);
			double scale = (H_ext_mag > 1.0e-10) ? std::min(1.0, H_init_mag / H_ext_mag) : 1.0;
			for(int k = 0; k < 3; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k] * scale;
			}
		}
		else if(dof == 6)
		{
			// MSC elements: initialize FlatField and NewFieldArray
			for(int k = 0; k < dof; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k];
			}
			// Initialize NewFieldArray with small H in z-direction (ELF mucal0 style)
			if(IntrctPtr->NewFieldArray != nullptr)
			{
				IntrctPtr->NewFieldArray[elem].x = 0.0;
				IntrctPtr->NewFieldArray[elem].y = 0.0;
				IntrctPtr->NewFieldArray[elem].z = H_init_mag;
			}
		}
		else
		{
			// Other elements: copy external field
			for(int k = 0; k < dof; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k];
			}
		}
	}

	// Work arrays - SystemMatrix needed because LU destroys the matrix
	// Use memcpy for fast O(n^2) copy
	std::vector<double> SystemMatrix(totalDOF * totalDOF);
	std::vector<double> RHS(totalDOF);
#ifdef HAVE_LAPACK
	// Pre-allocate LAPACK work arrays to avoid per-iteration allocation
	std::vector<double> A_col(totalDOF * totalDOF);  // Column-major for LAPACK
	std::vector<int> ipiv(totalDOF);
#endif

	// Cache polyhedron pointers to avoid repeated dynamic_cast
	std::vector<radTPolyhedron*> polyCache(AmOfMainElem, nullptr);
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		if(IntrctPtr->GetElementDOF(elem) == 6)
		{
			polyCache[elem] = dynamic_cast<radTPolyhedron*>(IntrctPtr->g3dRelaxPtrVect[elem]);
		}
	}

	// Debug logging: always write to hardcoded path
	FILE* debug_log = std::fopen("S:/Radia/01_GitHub/radia_iter.log", "w");
	if(debug_log)
	{
		std::fprintf(debug_log, "# Radia LU Solver AutoRelax_VariableDOF\n");
		std::fprintf(debug_log, "# AmOfMainElem=%d, totalDOF=%d\n", AmOfMainElem, totalDOF);
		std::fprintf(debug_log, "# iter, rel_change, M_avg_z, max_chi, min_chi\n");
		std::fflush(debug_log);
	}

	// Track previous iteration's average Mz for ELF-style convergence
	double prev_M_avg_z = 0.0;

	for(iterCount = 0; iterCount < MaxIterNumber; iterCount++)
	{
		// Store old values
		for(int i = 0; i < totalDOF; i++)
		{
			OldMagn[i] = FlatMagn[i];
		}
		// Store old chi and B-field for MSC convergence check
		const double MU_0_iter = 4.0 * 3.14159265358979323846 * 1.0e-7;
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			radTPolyhedron* poly = polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				OldChi[elem] = poly->CurrentChi;
				// Store old B-field norm: B = mu_0 * (H + M) = mu_0 * mu_r * H
				// For ELF mucal2 convergence check
				TVector3d M = poly->Magn;
				double chi = poly->CurrentChi;
				if(chi < 1.0e-6) chi = 1.0e-6;
				TVector3d H(M.x / chi, M.y / chi, M.z / chi);  // H = M / chi
				TVector3d B(MU_0_iter * (H.x + M.x), MU_0_iter * (H.y + M.y), MU_0_iter * (H.z + M.z));
				OldBnorm[elem] = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
			}
		}

		// NOTE: H field is updated AFTER solve (using H = M / chi from the solve's chi)
		// NewFieldArray is used as H_est for chi(H) computation in the next iteration

		// Copy base matrix (O(n^2)) and prepare for LU solve
		// Base matrix already contains -K/(4pi) for 6 DOF MSC elements
		// Equation: (-K/(4pi) + 1/chi * I) * sigma = H_ext_n
		std::memcpy(SystemMatrix.data(), BaseMatrix.data(), totalDOF * totalDOF * sizeof(double));

		// Update diagonal and RHS
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			int offset = IntrctPtr->GetElementDOFOffset(elem);

			if(dof == 3)
			{
				// Standard element (tetrahedron, 3 DOF MMM)
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

				TVector3d H_est(FlatField[offset], FlatField[offset+1], FlatField[offset+2]);
				TMatrix3d KsiTensor;
				TVector3d MrVect;
				MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);

				double chi_vals[3] = {KsiTensor.Str0.x, KsiTensor.Str1.y, KsiTensor.Str2.z};
				double Mr_vals[3] = {MrVect.x, MrVect.y, MrVect.z};

				for(int k = 0; k < 3; k++)
				{
					int row = offset + k;
					double chi = chi_vals[k];
					double inv_chi = (chi > 1.0e-10) ? (1.0 / chi) : 1.0e10;

					// Add 1/chi to diagonal
					SystemMatrix[row * totalDOF + row] += inv_chi;

					// RHS = H_ext + Mr/chi
					double Mr_over_chi = (chi > 1.0e-10) ? (Mr_vals[k] / chi) : 0.0;
					RHS[row] = FlatExtern[row] + Mr_over_chi;
				}
			}
			else if(dof == 6)
			{
				// MSC hexahedron (6 DOF): equation (-K/(4pi) + 1/chi * I) * sigma = H_ext_n
				// Base matrix already contains -K/(4pi), just add 1/chi to diagonal
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
				radTPolyhedron* poly = polyCache[elem];  // Use cached pointer

				double chi;
				if(iterCount == 0)
				{
					// First iteration: use initial chi from BH curve 2nd point (ELF mucal0 style)
					// This gives chi = B2/(mu0*H2) - 1 where (H2, B2) is the 2nd data point
					radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
					if(NonlinMater != nullptr)
					{
						double chi_init = NonlinMater->GetInitialChi_ELF_Style();
						if(chi_init > 0)
						{
							chi = chi_init;
						}
						else
						{
							// Fallback: use DefineInstantKsiTensor at H=100
							TVector3d H_est(0., 0., 100.0);
							TMatrix3d KsiTensor;
							TVector3d MrVect;
							MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);
							chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						}
					}
					else
					{
						// For non-nonlinear materials, use DefineInstantKsiTensor
						TVector3d H_est(0., 0., 100.0);
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);
						chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
					}
				}
				else
				{
					// Subsequent iterations: use chi from CurrentChi (updated after previous solve)
					// This ensures consistency between matrix build and chi update
					chi = (poly && poly->Use6DOF_MSC) ? poly->CurrentChi : 1.0;
				}

				if(chi < 1.0e-6) chi = 1.0e-6;
				double inv_chi = 1.0 / chi;

				// Store chi for post-solve H update
				if(poly && poly->Use6DOF_MSC)
				{
					poly->CurrentChi = chi;
				}

				for(int k = 0; k < 6; k++)
				{
					int row = offset + k;
					// Add 1/chi to diagonal for MSC constitutive relation
					SystemMatrix[row * totalDOF + row] += inv_chi;
					// RHS = H_ext dot n (already computed in m_flatExternFieldArray)
					RHS[row] = FlatExtern[row];
				}
			}
		}

		// Solve using pre-allocated LAPACK arrays
#ifdef HAVE_LAPACK
		// SystemMatrix is already in COLUMN-MAJOR format (copied from FlatInteract)
		// LAPACK dgesv expects column-major, so just copy directly
		std::memcpy(A_col.data(), SystemMatrix.data(), totalDOF * totalDOF * sizeof(double));
		int nrhs = 1;
		int info = 0;
		dgesv_(&totalDOF, &nrhs, A_col.data(), &totalDOF, ipiv.data(), RHS.data(), &totalDOF, &info);
		if(info != 0) return iterCount;
#else
		// SolveLU_Flat expects row-major, but SystemMatrix is column-major
		// Transpose in-place for the fallback solver
		for(int i = 0; i < totalDOF; i++)
		{
			for(int j = i + 1; j < totalDOF; j++)
			{
				std::swap(SystemMatrix[i * totalDOF + j], SystemMatrix[j * totalDOF + i]);
			}
		}
		int ierr = SolveLU_Flat(SystemMatrix, RHS, totalDOF);
		if(ierr != 0) return iterCount;
#endif

		// Extract solution
		for(int i = 0; i < totalDOF; i++)
		{
			FlatMagn[i] = RHS[i];
		}

		// Compute convergence (sigma-based for comparison, M-based computed later)
		double M_diff_sq = 0.0;
		double M_norm_sq = 0.0;
		for(int i = 0; i < totalDOF; i++)
		{
			double dM = FlatMagn[i] - OldMagn[i];
			M_diff_sq += dM * dM;
			M_norm_sq += FlatMagn[i] * FlatMagn[i];
		}

		// For 6DOF MSC: track actual M change (not sigma change)
		double M_sum_old_norm = 0.0;
		double M_sum_diff_sq = 0.0;
		double M_sum_new_z = 0.0;    // For average Mz
		double M_sum_diff_z = 0.0;   // For average dMz

		// Update element magnetization from flat array
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			int offset = IntrctPtr->GetElementDOFOffset(elem);

			if(dof == 3)
			{
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				g3dRelaxPtr->Magn.x = FlatMagn[offset + 0];
				g3dRelaxPtr->Magn.y = FlatMagn[offset + 1];
				g3dRelaxPtr->Magn.z = FlatMagn[offset + 2];
			}
			else if(dof == 6)
			{
				// MSC hexahedron: update sigma values and compute effective magnetization
				radTPolyhedron* poly = polyCache[elem];  // Use cached pointer
				if(poly && poly->Use6DOF_MSC)
				{
					// Store old M for convergence check
					TVector3d M_old = poly->Magn;

					// Store sigma values
					for(int k = 0; k < 6; k++)
					{
						poly->Sigma[k] = FlatMagn[offset + k];
					}

					// Compute effective magnetization from sigma using least-squares:
					// For each face: sigma_i = M dot n_i
					// Solve for M = [Mx, My, Mz] that best fits these constraints
					// Using weighted average based on face normals
					double Mx = 0.0, My = 0.0, Mz = 0.0;
					double wx = 0.0, wy = 0.0, wz = 0.0;
					for(int face = 0; face < 6; face++)
					{
						double sigma = poly->Sigma[face];
						TVector3d& n = poly->FaceNormal[face];
						double nx2 = n.x * n.x;
						double ny2 = n.y * n.y;
						double nz2 = n.z * n.z;
						// Weighted contribution: sigma_i * n_i decomposes M
						Mx += sigma * n.x;
						My += sigma * n.y;
						Mz += sigma * n.z;
						wx += nx2;
						wy += ny2;
						wz += nz2;
					}
					// Normalize (each direction has 2 opposing faces with |n|=1)
					if(wx > 1.0e-10) poly->Magn.x = Mx / wx;
					if(wy > 1.0e-10) poly->Magn.y = My / wy;
					if(wz > 1.0e-10) poly->Magn.z = Mz / wz;

					// Compute M change for convergence (ELF uses |dM|/|M_new|)
					double dMx = poly->Magn.x - M_old.x;
					double dMy = poly->Magn.y - M_old.y;
					double dMz = poly->Magn.z - M_old.z;
					double M_new_norm = poly->Magn.x*poly->Magn.x + poly->Magn.y*poly->Magn.y + poly->Magn.z*poly->Magn.z;
					M_sum_old_norm += M_new_norm;
					M_sum_diff_sq += dMx*dMx + dMy*dMy + dMz*dMz;
					// For average-based convergence (ELF style)
					M_sum_new_z += std::fabs(poly->Magn.z);
					M_sum_diff_z += std::fabs(dMz);

					// ELF-compatible H update: H = M / chi
					// Use the SAME chi that was used to build the matrix for this solve.
					// This is the ELF mucal1 algorithm.
					if(IntrctPtr->NewFieldArray != nullptr)
					{
						double chi_used = poly->CurrentChi;
						if(chi_used < 1.0e-6) chi_used = 1.0e-6;

						// H = M / chi (constitutive relation, same as ELF)
						IntrctPtr->NewFieldArray[elem].x = poly->Magn.x / chi_used;
						IntrctPtr->NewFieldArray[elem].y = poly->Magn.y / chi_used;
						IntrctPtr->NewFieldArray[elem].z = poly->Magn.z / chi_used;
					}
				}
			}
		}

		// Compute convergence for 6DOF MSC elements
		// NonlinearMethod: 0=mucal1 (chi-change), 1=mucal2 (B-change/Newton)
		// For 3DOF MMM elements: use M change
		double max_B_rel_change = 0.0;
		double max_chi_rel_change = 0.0;
		double sum_chi_rel_change = 0.0;  // For ELF-style average convergence
		double sum_dmu = 0.0;             // Sum of |mu_new - mu_old|
		double sum_mu = 0.0;              // Sum of mu_new
		int n_6dof_elements = 0;          // Count for average
		bool has_6dof_elements = false;
		const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;

		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			if(dof == 6)
			{
				has_6dof_elements = true;
				radTPolyhedron* poly = polyCache[elem];  // Use cached pointer
				if(poly && poly->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
				{
					TVector3d H_new = IntrctPtr->NewFieldArray[elem];
					radTMaterial* MaterPtr = (radTMaterial*)(IntrctPtr->g3dRelaxPtrVect[elem]->MaterHandle.rep);

					// Get chi used for this iteration's matrix build
					double chi_matrix = poly->CurrentChi;
					double mu_old = chi_matrix + 1.0;

					// Use ELF-style dual-method chi update:
					// Method 1: Standard mu = B(H)/(mu_0*H) with Newton correction
					// Method 2: H+B sum interpolation (physics-based)
					// Selects method with smaller |mu_new - mu_old|
					double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);

					// Try to cast to radTNonlinearIsotropMaterial to use ComputeChiDualMethod
					radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
					double chi_new;
					if(NonlinMater != nullptr)
					{
						// Use ELF-style dual-method chi update
						chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old);
						// Debug: write first element's chi update
						if(debug_log && elem == 0 && iterCount < 10)
						{
							std::fprintf(debug_log, "  iter%d elem0: H_mag=%.0f, chi_matrix=%.1f, mu_old=%.1f, chi_new=%.1f\n",
							            iterCount+1, H_mag, chi_matrix, mu_old, chi_new);
						}
					}
					else
					{
						// Fallback: use DefineInstantKsiTensor for other material types
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_new, KsiTensor, MrVect);
						chi_new = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						if(chi_new < 1.0e-6) chi_new = 1.0e-6;
						// Debug: write first element for fallback case
						if(debug_log && elem == 0 && iterCount < 10)
						{
							std::fprintf(debug_log, "  iter%d elem0 (FALLBACK): H_mag=%.0f, chi_matrix=%.1f, chi_new=%.1f\n",
							            iterCount+1, H_mag, chi_matrix, chi_new);
						}
					}

					// Update chi for next iteration with optional under-relaxation
					// relax=0: full step, relax>0: chi = chi_new*(1-relax) + chi_old*relax
					double relax = rad.m_relax;
					if(relax > 0.0 && relax <= 1.0)
					{
						chi_new = chi_new * (1.0 - relax) + chi_matrix * relax;
					}
					poly->CurrentChi = chi_new;

					// ELF mucal2 (Newton-Raphson) convergence: use B-field change
					// B = mu_0 * (H + M) = mu_0 * mu_r * H
					// ELF uses: rel_change = |B_new - B_old| / B_sat
					TVector3d M_new = poly->Magn;
					double chi_for_B = chi_new;
					if(chi_for_B < 1.0e-6) chi_for_B = 1.0e-6;
					TVector3d H_for_B(M_new.x / chi_for_B, M_new.y / chi_for_B, M_new.z / chi_for_B);
					TVector3d B_new_vec(MU_0 * (H_for_B.x + M_new.x),
					                    MU_0 * (H_for_B.y + M_new.y),
					                    MU_0 * (H_for_B.z + M_new.z));
					double B_new_norm = std::sqrt(B_new_vec.x*B_new_vec.x + B_new_vec.y*B_new_vec.y + B_new_vec.z*B_new_vec.z);

					// Get B_sat from BH curve (ELF uses B_max from last point)
					double B_sat = 1.0;  // fallback
					if(NonlinMater != nullptr)
					{
						B_sat = NonlinMater->GetBsaturation();
						if(B_sat < 1.0e-10) B_sat = 1.0;  // fallback
					}

					// B-field convergence: |B_new - B_old| / B_sat
					double B_old_norm = OldBnorm[elem];
					double B_rel_change = std::fabs(B_new_norm - B_old_norm) / B_sat;
					if(B_rel_change > max_B_rel_change)
						max_B_rel_change = B_rel_change;

					// Also track chi change for debugging
					double mu_matrix = chi_matrix + 1.0;
					double mu_new_val = chi_new + 1.0;
					double mu_rel_change = std::fabs(mu_new_val - mu_matrix) / mu_matrix;
					if(mu_rel_change > max_chi_rel_change)
						max_chi_rel_change = mu_rel_change;

					n_6dof_elements++;
				}
			}
		}

		// Convergence criterion depends on element type:
		// - 3DOF MMM elements: use M change (FlatMagn is M)
		// - 6DOF MSC elements: use B-field change (ELF mucal2/Newton-Raphson style)
		double rel_change;
		if(has_6dof_elements)
		{
			// For 6DOF MSC: use ELF-style B-field change (mucal2)
			// rel_change = MAX over all elements of |B_new - B_old| / B_sat
			rel_change = max_B_rel_change;
		}
		else
		{
			// For 3DOF MMM: use M change (original Radia style)
			rel_change = (M_norm_sq > 1.0e-30) ? std::sqrt(M_diff_sq / M_norm_sq) : std::sqrt(M_diff_sq);
		}
		MisfitE2 = rel_change * rel_change;

		// Debug logging - write to file
		if(debug_log)
		{
			double M_sum_z = 0.0;
			double max_chi = 0.0;
			double min_chi = 1.0e20;
			int n_6dof = 0;
			for(int elem = 0; elem < AmOfMainElem; elem++)
			{
				int dof = IntrctPtr->GetElementDOF(elem);
				if(dof == 6)
				{
					radTPolyhedron* poly = polyCache[elem];
					if(poly && poly->Use6DOF_MSC)
					{
						M_sum_z += poly->Magn.z;
						if(poly->CurrentChi > max_chi) max_chi = poly->CurrentChi;
						if(poly->CurrentChi < min_chi) min_chi = poly->CurrentChi;
						n_6dof++;
					}
				}
			}
			double M_avg_z = (n_6dof > 0) ? M_sum_z / n_6dof : 0.0;
			// Also compute M-based convergence (sqrt(M_diff_sq / M_norm_sq))
			double M_rel_change = (M_norm_sq > 1.0e-30) ? std::sqrt(M_diff_sq / M_norm_sq) : 0.0;
			// Log: rel_change is now B-field based for 6DOF, show both B and chi
			std::fprintf(debug_log, "%d, B_rel=%.6e, chi_rel=%.6e, M_avg_z=%.0f, max_chi=%.1f, min_chi=%.1f\n",
			            iterCount + 1, max_B_rel_change, max_chi_rel_change, M_avg_z, max_chi, min_chi);
			std::fflush(debug_log);
		}

		if(rel_change <= PrecOnMagnetiz)
		{
			iterCount++;
			break;
		}

		if(radYield.Check() == 0)
		{
			if(debug_log) std::fclose(debug_log);
			return iterCount;
		}
	}

	if(debug_log) std::fclose(debug_log);
	IntrctPtr->RelaxStatusParam.MisfitM = std::sqrt(MisfitE2);

	return iterCount;
}

//-------------------------------------------------------------------------
// BiCGSTAB solver variable DOF methods
//-------------------------------------------------------------------------

void radTRelaxationMethNo_1::MatVec_VariableDOF(const std::vector<double>& x, std::vector<double>& y,
                                                 const std::vector<double>& inv_chi, int totalDOF)
{
	// Computes y = A * x where A = (base matrix) + diag(1/chi)
	// Sign convention:
	// - MMM (3 DOF): FlatInteract stores N, equation is (-N + I/chi) -> negate
	// - MSC (6 DOF): FlatInteract stores -K/(4pi), equation is (-K/(4pi) + I/chi) -> use as-is
	//
	// IMPORTANT: FlatInteract is stored in COLUMN-MAJOR format (Fortran/LAPACK style)
	// Element at (row, col) is at index [col * totalDOF + row]

	const double* FlatInteract = IntrctPtr->GetFlatInteractMatrix();
	if(FlatInteract == nullptr) return;

	std::fill(y.begin(), y.end(), 0.0);

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	#pragma omp parallel for if(AmOfMainElem > 50)
	for(int row_elem = 0; row_elem < AmOfMainElem; row_elem++)
	{
		int dof_row = IntrctPtr->GetElementDOF(row_elem);
		int offset_row = IntrctPtr->GetElementDOFOffset(row_elem);

		// Diagonal contribution: (1/chi) * x
		for(int k = 0; k < dof_row; k++)
		{
			y[offset_row + k] = inv_chi[offset_row + k] * x[offset_row + k];
		}

		// Matrix-vector product
		for(int col_elem = 0; col_elem < AmOfMainElem; col_elem++)
		{
			int dof_col = IntrctPtr->GetElementDOF(col_elem);
			int offset_col = IntrctPtr->GetElementDOFOffset(col_elem);

			// Get block from flat matrix - COLUMN-MAJOR: block starts at [col * totalDOF + row]
			const double* block = &FlatInteract[offset_col * totalDOF + offset_row];

			// Sign: negate for 3 DOF blocks, use as-is for 6x6 MSC blocks
			double sign = (dof_row == 6 && dof_col == 6) ? 1.0 : -1.0;

			// Column-major block access: element (i, j) within block is at [j * totalDOF + i]
			for(int i = 0; i < dof_row; i++)
			{
				double sum = 0.0;
				for(int j = 0; j < dof_col; j++)
				{
					sum += block[j * totalDOF + i] * x[offset_col + j];
				}
				y[offset_row + i] += sign * sum;
			}
		}
	}
}

void radTRelaxationMethNo_1::GetDiagonalElements_VariableDOF(std::vector<double>& diag,
                                                              const std::vector<double>& inv_chi, int totalDOF)
{
	// Extract diagonal elements for Jacobi preconditioner
	// Sign convention matches MatVec_VariableDOF:
	// - MMM (3 DOF): negate stored values
	// - MSC (6 DOF): use as-is
	const double* FlatInteract = IntrctPtr->GetFlatInteractMatrix();
	if(FlatInteract == nullptr) return;

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		// Get diagonal block
		const double* diag_block = &FlatInteract[offset * totalDOF + offset];

		// Sign: negate for 3 DOF, use as-is for 6 DOF MSC
		double sign = (dof == 6) ? 1.0 : -1.0;

		for(int k = 0; k < dof; k++)
		{
			// Diagonal element: sign*matrix_ii + 1/chi
			diag[offset + k] = sign * diag_block[k * totalDOF + k] + inv_chi[offset + k];
		}
	}
}

int radTRelaxationMethNo_1::SolveBiCGSTAB_VariableDOF(int totalDOF, double tol, int max_iter, double& residual)
{
	// BiCGSTAB with Jacobi preconditioner for variable DOF systems
	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	// Allocate work vectors
	std::vector<double> r(totalDOF), r0(totalDOF), p(totalDOF), v(totalDOF), s(totalDOF), t(totalDOF);
	std::vector<double> p_hat(totalDOF), s_hat(totalDOF), diag_inv(totalDOF);
	std::vector<double> inv_chi(totalDOF);
	std::vector<double> rhs(totalDOF);
	std::vector<double> sol(totalDOF);

	double* FlatMagn = IntrctPtr->GetFlatMagnArray();
	double* FlatField = IntrctPtr->GetFlatFieldArray();
	double* FlatExtern = IntrctPtr->GetFlatExternFieldArray();

	if(FlatMagn == nullptr || FlatField == nullptr || FlatExtern == nullptr) return 0;

	// Pre-compute 1/chi and RHS for all elements
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		if(dof == 3)
		{
			// Standard element
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

			TVector3d InstH(FlatField[offset], FlatField[offset+1], FlatField[offset+2]);
			TMatrix3d KsiTensor;
			TVector3d MrVect;
			MaterPtr->DefineInstantKsiTensor(InstH, KsiTensor, MrVect);

			double chi_vals[3] = {KsiTensor.Str0.x, KsiTensor.Str1.y, KsiTensor.Str2.z};
			double Mr_vals[3] = {MrVect.x, MrVect.y, MrVect.z};

			for(int k = 0; k < 3; k++)
			{
				double chi = chi_vals[k];
				inv_chi[offset + k] = (chi > 1.0e-10) ? (1.0 / chi) : 1.0e10;
				double inv_chi_rhs = (chi > 1.0e-10) ? (1.0 / chi) : 0.0;
				rhs[offset + k] = FlatExtern[offset + k] + Mr_vals[k] * inv_chi_rhs;
			}
		}
		else if(dof == 6)
		{
			// MSC hexahedron (6 DOF): equation (-K/(4pi) + 1/chi * I) * sigma = H_ext_n
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);

			// Use CurrentChi which was initialized in AutoRelax_VariableDOF
			double chi = (poly && poly->Use6DOF_MSC && poly->CurrentChi > 1.0e-6)
			           ? poly->CurrentChi : 1.0;
			if(chi < 1.0e-6) chi = 1.0e-6;

			double inv_chi_val = 1.0 / chi;

			for(int k = 0; k < 6; k++)
			{
				inv_chi[offset + k] = inv_chi_val;
				// RHS = H_ext dot n (already computed in m_flatExternFieldArray)
				rhs[offset + k] = FlatExtern[offset + k];
			}
		}
		else
		{
			// Fallback for unknown DOF types
			for(int k = 0; k < dof; k++)
			{
				inv_chi[offset + k] = 1.0;
				rhs[offset + k] = FlatExtern[offset + k];
			}
		}
	}

	// Initial guess: use zero for better stability
	// (FlatMagn may contain garbage from previous failed iterations)
	for(int i = 0; i < totalDOF; i++)
	{
		sol[i] = 0.0;
	}

	// Build Jacobi preconditioner
	GetDiagonalElements_VariableDOF(diag_inv, inv_chi, totalDOF);
	for(int i = 0; i < totalDOF; i++)
	{
		diag_inv[i] = (std::abs(diag_inv[i]) > 1.0e-15) ? (1.0 / diag_inv[i]) : 1.0;
	}

	// Initialize: r0 = b - A*x0
	MatVec_VariableDOF(sol, v, inv_chi, totalDOF);
	Copy(rhs, r, totalDOF);
	Axpy(-1.0, v, r, totalDOF);
	Copy(r, r0, totalDOF);

	double rho = 1.0, alpha_bicg = 1.0, omega = 1.0;
	std::fill(p.begin(), p.end(), 0.0);
	std::fill(v.begin(), v.end(), 0.0);

	double rhs_norm = Norm2(rhs, totalDOF);
	if(rhs_norm < 1.0e-30) rhs_norm = 1.0;

	int iter;
	for(iter = 1; iter <= max_iter; iter++)
	{
		double rho_old = rho;
		rho = Dot(r0, r, totalDOF);

		if(std::abs(rho) < 1.0e-30)
		{
			residual = Norm2(r, totalDOF) / rhs_norm;
			break;
		}

		if(iter == 1)
		{
			Copy(r, p, totalDOF);
		}
		else
		{
			if(std::abs(rho_old * omega) < 1.0e-30)
			{
				residual = Norm2(r, totalDOF) / rhs_norm;
				break;
			}
			double beta = (rho / rho_old) * (alpha_bicg / omega);
			Axpy(-omega, v, p, totalDOF);
			Scale(beta, p, totalDOF);
			Axpy(1.0, r, p, totalDOF);
		}

		// Apply preconditioner
		#pragma omp parallel for if(totalDOF > 100)
		for(int i = 0; i < totalDOF; i++)
		{
			p_hat[i] = diag_inv[i] * p[i];
		}

		MatVec_VariableDOF(p_hat, v, inv_chi, totalDOF);

		double r0_dot_v = Dot(r0, v, totalDOF);
		if(std::abs(r0_dot_v) < 1.0e-30)
		{
			residual = Norm2(r, totalDOF) / rhs_norm;
			break;
		}
		alpha_bicg = rho / r0_dot_v;

		Copy(r, s, totalDOF);
		Axpy(-alpha_bicg, v, s, totalDOF);

		double s_norm = Norm2(s, totalDOF);
		if(s_norm / rhs_norm < tol)
		{
			Axpy(alpha_bicg, p_hat, sol, totalDOF);
			residual = s_norm / rhs_norm;
			break;
		}

		#pragma omp parallel for if(totalDOF > 100)
		for(int i = 0; i < totalDOF; i++)
		{
			s_hat[i] = diag_inv[i] * s[i];
		}

		MatVec_VariableDOF(s_hat, t, inv_chi, totalDOF);

		double t_dot_s = Dot(t, s, totalDOF);
		double t_dot_t = Dot(t, t, totalDOF);
		if(std::abs(t_dot_t) < 1.0e-30)
		{
			Axpy(alpha_bicg, p_hat, sol, totalDOF);
			residual = s_norm / rhs_norm;
			break;
		}
		omega = t_dot_s / t_dot_t;

		Axpy(alpha_bicg, p_hat, sol, totalDOF);
		Axpy(omega, s_hat, sol, totalDOF);

		Copy(s, r, totalDOF);
		Axpy(-omega, t, r, totalDOF);

		double r_norm = Norm2(r, totalDOF);
		double prev_residual = residual;
		residual = r_norm / rhs_norm;
		if(residual < tol) break;

		// Detect divergence: if residual increases by more than 10x, stop
		if(residual > 10.0 * prev_residual && prev_residual > 0.0 && iter > 5)
		{
			// BiCGSTAB is diverging, stop early
			break;
		}

		if(std::abs(omega) < 1.0e-30) break;
	}

	// Copy solution back to flat array
	for(int i = 0; i < totalDOF; i++)
	{
		FlatMagn[i] = sol[i];
	}

	return iter;
}

int radTRelaxationMethNo_1::AutoRelax_VariableDOF(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	// Check if variable DOF is active
	if(!IntrctPtr->HasVariableDOF())
	{
		return AutoRelax(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
	}

	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	int AmOfMainElem = IntrctPtr->AmOfMainElem;
	if(AmOfMainElem <= 0) return 0;

	int totalDOF = IntrctPtr->GetTotalDOF();
	if(totalDOF <= 0) return 0;

	double* FlatMagn = IntrctPtr->GetFlatMagnArray();
	double* FlatField = IntrctPtr->GetFlatFieldArray();
	double* FlatExtern = IntrctPtr->GetFlatExternFieldArray();

	if(FlatMagn == nullptr || FlatField == nullptr || FlatExtern == nullptr) return 0;

	std::vector<double> OldMagn(totalDOF);
	// Store old chi values for MSC convergence check (ELF mucal1 style)
	std::vector<double> OldChi_bicg(AmOfMainElem, 0.0);
	std::vector<double> OldBnorm_bicg(AmOfMainElem, 0.0);  // For ELF mucal2 B-field convergence
	double MisfitE2 = 1.0e30;
	int totalIterCount = 0;
	int outerIter = 0;

	// Cache polyhedron pointers to avoid repeated dynamic_cast
	// Cache ALL elements (not just DOF==6) so we can access them during chi initialization
	std::vector<radTPolyhedron*> polyCache(AmOfMainElem, nullptr);
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
		polyCache[elem] = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
	}

	// Initialize CurrentChi with ELF-style initial value (same as LU solver and HACApK)
	// This uses BH curve's 2nd point: chi = B2/(mu0*H2) - 1
	// Without this, CurrentChi starts at 1.0 causing slow convergence (18 iterations instead of 3-4)
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		if(dof == 6)
		{
			radTPolyhedron* poly = polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
				radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
				if(NonlinMater != nullptr)
				{
					double chi_init = NonlinMater->GetInitialChi_ELF_Style();
					if(chi_init > 0)
					{
						poly->CurrentChi = chi_init;
					}
				}
			}
		}
	}

	// Initialize H field and NewFieldArray for chi(H) computation
	// FIX (2025-12-24): Initialize NewFieldArray for 6DOF elements too
	const double H_init_mag = 1000.0;
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		if(dof == 3)
		{
			double H_ext_mag = 0.0;
			for(int k = 0; k < 3; k++)
			{
				H_ext_mag += FlatExtern[offset + k] * FlatExtern[offset + k];
			}
			H_ext_mag = std::sqrt(H_ext_mag);
			double scale = (H_ext_mag > 1.0e-10) ? std::min(1.0, H_init_mag / H_ext_mag) : 1.0;
			for(int k = 0; k < 3; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k] * scale;
			}
		}
		else if(dof == 6)
		{
			// For 6DOF MSC elements: initialize FlatField and NewFieldArray
			for(int k = 0; k < dof; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k];
			}
			// Initialize NewFieldArray with H_init_mag in z-direction for chi(H) computation
			// This ensures SolveBiCGSTAB_VariableDOF can detect first iteration properly
			if(IntrctPtr->NewFieldArray != nullptr)
			{
				IntrctPtr->NewFieldArray[elem].x = 0.0;
				IntrctPtr->NewFieldArray[elem].y = 0.0;
				IntrctPtr->NewFieldArray[elem].z = H_init_mag;
			}
		}
		else
		{
			for(int k = 0; k < dof; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k];
			}
		}
	}

	// Track previous iteration's average Mz for ELF-style convergence
	double prev_M_avg_z_bicg = 0.0;

	// Outer nonlinear iteration
	for(outerIter = 0; outerIter < MaxIterNumber; outerIter++)
	{
		// Store old values
		for(int i = 0; i < totalDOF; i++)
		{
			OldMagn[i] = FlatMagn[i];
		}
		// Store old chi for MSC convergence check
		const double MU_0_bicg = 4.0 * 3.14159265358979323846 * 1.0e-7;
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			radTPolyhedron* poly = polyCache[elem];  // Use cached pointer
			if(poly && poly->Use6DOF_MSC)
			{
				OldChi_bicg[elem] = poly->CurrentChi;
				// Store old B-field norm for ELF mucal2 convergence
				TVector3d M = poly->Magn;
				double chi = poly->CurrentChi;
				if(chi < 1.0e-6) chi = 1.0e-6;
				TVector3d H(M.x / chi, M.y / chi, M.z / chi);
				TVector3d B(MU_0_bicg * (H.x + M.x), MU_0_bicg * (H.y + M.y), MU_0_bicg * (H.z + M.z));
				OldBnorm_bicg[elem] = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
			}
		}

		// Update H from M
		if(outerIter > 0)
		{
			for(int elem = 0; elem < AmOfMainElem; elem++)
			{
				int dof = IntrctPtr->GetElementDOF(elem);
				int offset = IntrctPtr->GetElementDOFOffset(elem);

				if(dof == 3)
				{
					radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
					radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

					TVector3d H_est(FlatField[offset], FlatField[offset+1], FlatField[offset+2]);
					TMatrix3d KsiTensor;
					TVector3d MrVect;
					MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);

					double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
					if(chi < 1.0e-6) chi = 1.0e-6;

					for(int k = 0; k < 3; k++)
					{
						FlatField[offset + k] = FlatMagn[offset + k] / chi;
					}
				}
			}
		}

		// Solve with BiCGSTAB
		// NOTE: max_iter for BiCGSTAB inner loop should be FIXED (not user's MaxIterNumber)
		// User's MaxIterNumber controls OUTER nonlinear iterations, not inner BiCGSTAB
		// bicg_tol is set via rad.SetBiCGSTABTol() Python API (default: 1e-4, ELF-compatible)
		double residual = 0.0;
		const double bicg_tol = rad.m_bicg_tol;
		const int bicg_max_iter = 10000;  // Inner BiCGSTAB max iterations (fixed)
		int n_iter = SolveBiCGSTAB_VariableDOF(totalDOF, bicg_tol, bicg_max_iter, residual);
		totalIterCount += n_iter;

		// Update element magnetization
		double M_diff_sq = 0.0;
		double M_norm_sq = 0.0;
		for(int i = 0; i < totalDOF; i++)
		{
			double dM = FlatMagn[i] - OldMagn[i];
			M_diff_sq += dM * dM;
			M_norm_sq += FlatMagn[i] * FlatMagn[i];
		}

		// For 6DOF MSC: track actual Mz for ELF-style convergence
		double M_sum_new_z_bicg = 0.0;
		int n_6dof_elems = 0;

		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			int offset = IntrctPtr->GetElementDOFOffset(elem);

			if(dof == 3)
			{
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				g3dRelaxPtr->Magn.x = FlatMagn[offset + 0];
				g3dRelaxPtr->Magn.y = FlatMagn[offset + 1];
				g3dRelaxPtr->Magn.z = FlatMagn[offset + 2];
			}
			else if(dof == 6)
			{
				// MSC hexahedron: update sigma values and compute effective magnetization
				radTPolyhedron* poly = polyCache[elem];  // Use cached pointer
				if(poly && poly->Use6DOF_MSC)
				{
					// Store sigma values
					for(int k = 0; k < 6; k++)
					{
						poly->Sigma[k] = FlatMagn[offset + k];
					}

					// Compute effective magnetization from sigma using weighted average
					double Mx = 0.0, My = 0.0, Mz = 0.0;
					double wx = 0.0, wy = 0.0, wz = 0.0;
					for(int face = 0; face < 6; face++)
					{
						double sigma = poly->Sigma[face];
						TVector3d& n = poly->FaceNormal[face];
						double nx2 = n.x * n.x;
						double ny2 = n.y * n.y;
						double nz2 = n.z * n.z;
						Mx += sigma * n.x;
						My += sigma * n.y;
						Mz += sigma * n.z;
						wx += nx2;
						wy += ny2;
						wz += nz2;
					}
					if(wx > 1.0e-10) poly->Magn.x = Mx / wx;
					if(wy > 1.0e-10) poly->Magn.y = My / wy;
					if(wz > 1.0e-10) poly->Magn.z = Mz / wz;

					// Track Mz for ELF-style convergence
					M_sum_new_z_bicg += std::fabs(poly->Magn.z);
					n_6dof_elems++;

					// ELF-compatible H update: H = M / chi
					// Use the SAME chi that was used to build the matrix for this solve.
					// This is the ELF mucal1 algorithm.
					if(IntrctPtr->NewFieldArray != nullptr)
					{
						double chi_used = poly->CurrentChi;
						if(chi_used < 1.0e-6) chi_used = 1.0e-6;

						// H = M / chi (constitutive relation, same as ELF)
						IntrctPtr->NewFieldArray[elem].x = poly->Magn.x / chi_used;
						IntrctPtr->NewFieldArray[elem].y = poly->Magn.y / chi_used;
						IntrctPtr->NewFieldArray[elem].z = poly->Magn.z / chi_used;
					}
				}
			}
		}

		// Compute convergence for 6DOF MSC elements
		// NonlinearMethod: 0=mucal1 (chi-change), 1=mucal2 (B-change/Newton)
		// For 3DOF MMM elements: use M change
		double max_chi_rel_change = 0.0;
		double sum_chi_rel_change = 0.0;  // For ELF-style average convergence
		double sum_dmu = 0.0;             // Sum of |mu_new - mu_old|
		double sum_mu = 0.0;              // Sum of mu_new
		int n_6dof_elements = 0;          // Count for average
		double max_B_rel_change = 0.0;
		bool has_6dof_elements = false;
		const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;

		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			if(dof == 6)
			{
				has_6dof_elements = true;
				radTPolyhedron* poly = polyCache[elem];  // Use cached pointer
				if(poly && poly->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
				{
					TVector3d H_new = IntrctPtr->NewFieldArray[elem];
					radTMaterial* MaterPtr = (radTMaterial*)(IntrctPtr->g3dRelaxPtrVect[elem]->MaterHandle.rep);

					// Get chi used for this iteration's matrix build
					double chi_matrix = poly->CurrentChi;
					double mu_old = chi_matrix + 1.0;

					// Use ELF-style dual-method chi update
					double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);
					radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
					double chi_new;
					if(NonlinMater != nullptr)
					{
						// Use ELF-style dual-method chi update
						chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old);
					}
					else
					{
						// Fallback: use DefineInstantKsiTensor
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_new, KsiTensor, MrVect);
						chi_new = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						if(chi_new < 1.0e-6) chi_new = 1.0e-6;
					}

					// Update chi for next iteration with optional under-relaxation
					// relax=0: full step, relax>0: chi = chi_new*(1-relax) + chi_old*relax
					double relax_bicg = rad.m_relax;
					if(relax_bicg > 0.0 && relax_bicg <= 1.0)
					{
						chi_new = chi_new * (1.0 - relax_bicg) + chi_matrix * relax_bicg;
					}
					poly->CurrentChi = chi_new;

					// ELF mucal2 (Newton-Raphson) convergence: use B-field change
					TVector3d M_new = poly->Magn;
					double chi_for_B = chi_new;
					if(chi_for_B < 1.0e-6) chi_for_B = 1.0e-6;
					TVector3d H_for_B(M_new.x / chi_for_B, M_new.y / chi_for_B, M_new.z / chi_for_B);
					TVector3d B_new_vec(MU_0_bicg * (H_for_B.x + M_new.x),
					                    MU_0_bicg * (H_for_B.y + M_new.y),
					                    MU_0_bicg * (H_for_B.z + M_new.z));
					double B_new_norm = std::sqrt(B_new_vec.x*B_new_vec.x + B_new_vec.y*B_new_vec.y + B_new_vec.z*B_new_vec.z);

					// Get B_sat from BH curve
					double B_sat = 1.0;
					if(NonlinMater != nullptr)
					{
						B_sat = NonlinMater->GetBsaturation();
						if(B_sat < 1.0e-10) B_sat = 1.0;
					}

					// B-field convergence: |B_new - B_old| / B_sat
					double B_old_norm = OldBnorm_bicg[elem];
					double B_rel_change = std::fabs(B_new_norm - B_old_norm) / B_sat;
					if(B_rel_change > max_B_rel_change)
						max_B_rel_change = B_rel_change;

					// Also track chi change for debugging
					double mu_matrix = chi_matrix + 1.0;
					double mu_new_val = chi_new + 1.0;
					double mu_rel_change = std::fabs(mu_new_val - mu_matrix) / mu_matrix;
					if(mu_rel_change > max_chi_rel_change)
						max_chi_rel_change = mu_rel_change;

					n_6dof_elements++;
				}
			}
		}

		// Convergence criterion depends on element type:
		// - 3DOF MMM elements: use M change (FlatMagn is M)
		// - 6DOF MSC elements: use B-field change (ELF mucal2/Newton-Raphson style)
		double rel_change;
		if(has_6dof_elements)
		{
			// For 6DOF MSC: use ELF-style B-field change (mucal2)
			rel_change = max_B_rel_change;
		}
		else
		{
			// For 3DOF MMM: use M change (original Radia style)
			rel_change = (M_norm_sq > 1.0e-30) ? std::sqrt(M_diff_sq / M_norm_sq) : std::sqrt(M_diff_sq);
		}
		MisfitE2 = rel_change * rel_change;

		if(rel_change <= PrecOnMagnetiz)
		{
			outerIter++;
			break;
		}

		if(radYield.Check() == 0) return outerIter;
	}

	IntrctPtr->RelaxStatusParam.MisfitM = std::sqrt(MisfitE2);

	return outerIter;
}

//=========================================================================
// Method 2: BiCGSTAB with H-matrix (HACApK ACA+)
// Conditionally compiled when RADIA_USE_HACAPK is defined
//=========================================================================

#ifdef RADIA_USE_HACAPK

double radTRelaxationMethNo_2::Dot(const std::vector<double>& a, const std::vector<double>& b, int n)
{
#ifdef HAVE_LAPACK
	return cblas_ddot(n, a.data(), 1, b.data(), 1);
#else
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum) if(n > 100)
	for(int i = 0; i < n; i++)
	{
		sum += a[i] * b[i];
	}
	return sum;
#endif
}

double radTRelaxationMethNo_2::Norm2(const std::vector<double>& a, int n)
{
#ifdef HAVE_LAPACK
	return cblas_dnrm2(n, a.data(), 1);
#else
	return std::sqrt(Dot(a, a, n));
#endif
}

void radTRelaxationMethNo_2::Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n)
{
#ifdef HAVE_LAPACK
	cblas_daxpy(n, alpha, x.data(), 1, y.data(), 1);
#else
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		y[i] += alpha * x[i];
	}
#endif
}

void radTRelaxationMethNo_2::Copy(const std::vector<double>& src, std::vector<double>& dst, int n)
{
#ifdef HAVE_LAPACK
	cblas_dcopy(n, src.data(), 1, dst.data(), 1);
#else
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		dst[i] = src[i];
	}
#endif
}

void radTRelaxationMethNo_2::Scale(double alpha, std::vector<double>& x, int n)
{
#ifdef HAVE_LAPACK
	cblas_dscal(n, alpha, x.data(), 1);
#else
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		x[i] *= alpha;
	}
#endif
}

int radTRelaxationMethNo_2::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	// HACApK supports both 3DOF tetrahedra and 6DOF hexahedra
	return AutoRelax_VariableDOF(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
}

//-------------------------------------------------------------------------
// SolveBiCGSTAB_HMatrix_VariableDOF
// BiCGSTAB with H-matrix for both 3DOF tetrahedra and 6DOF hexahedra
//-------------------------------------------------------------------------

int radTRelaxationMethNo_2::SolveBiCGSTAB_HMatrix_VariableDOF(int totalDOF, double tol, int max_iter, double& residual)
{
	if (!m_hacapk || !m_hacapk->IsValid()) return 0;

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	// Allocate work vectors
	std::vector<double> r(totalDOF), r0(totalDOF), p(totalDOF), v(totalDOF), s(totalDOF), t(totalDOF);
	std::vector<double> p_hat(totalDOF), s_hat(totalDOF), diag_inv(totalDOF);
	std::vector<double> inv_chi(totalDOF);
	std::vector<double> rhs(totalDOF);
	std::vector<double> sol(totalDOF);

	double* FlatMagn = IntrctPtr->GetFlatMagnArray();
	double* FlatField = IntrctPtr->GetFlatFieldArray();
	double* FlatExtern = IntrctPtr->GetFlatExternFieldArray();

	if(FlatMagn == nullptr || FlatField == nullptr || FlatExtern == nullptr) return 0;

	// Pre-compute 1/chi and RHS for all elements
	// Supports both 3DOF tetrahedra and 6DOF hexahedra
	// CRITICAL FIX (2025-12-20): Use current H field estimate from FlatField, not FlatExtern
	// This matches the LU/BiCGSTAB solvers that use NewFieldAr[i] for chi computation
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		if(dof != 3 && dof != 6)
		{
			// HACApK supports 3DOF tetrahedra and 6DOF hexahedra
			std::cerr << "[HACApK] Error: Element " << elem << " has " << dof
			          << " DOF, expected 3 (tetrahedra) or 6 (hexahedra)" << std::endl;
			return 0;
		}

		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
		radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		// FIX (2025-12-21): For first iteration, use small H (like ELF's mucal0)
		// to get initial chi from the linear region of B-H curve.
		// Detect first iteration by checking if NewFieldArray equals ExternFieldArray
		TVector3d H_est;
		bool isFirstIter = false;
		if(IntrctPtr->NewFieldArray != nullptr && IntrctPtr->ExternFieldArray != nullptr)
		{
			TVector3d& newField = IntrctPtr->NewFieldArray[elem];
			TVector3d& extField = IntrctPtr->ExternFieldArray[elem];
			if(std::abs(newField.x - extField.x) < 1.0e-10 &&
			   std::abs(newField.y - extField.y) < 1.0e-10 &&
			   std::abs(newField.z - extField.z) < 1.0e-10)
			{
				isFirstIter = true;
			}
		}

		if(isFirstIter)
		{
			// First iteration: use H = 1000 A/m (same as LU solver line 284)
			H_est = TVector3d(0., 0., 1000.0);
		}
		else if(IntrctPtr->NewFieldArray != nullptr)
		{
			H_est = IntrctPtr->NewFieldArray[elem];
		}
		else
		{
			H_est = TVector3d(0., 0., FlatExtern[offset]);
		}
		TMatrix3d KsiTensor;
		TVector3d MrVect;
		MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);

		if(dof == 3)
		{
			// For 3DOF: use anisotropic chi (same as LU solver lines 410-412)
			double inv_chi_x = (KsiTensor.Str0.x > 1.0e-10) ? (1.0 / KsiTensor.Str0.x) : 1.0e10;
			double inv_chi_y = (KsiTensor.Str1.y > 1.0e-10) ? (1.0 / KsiTensor.Str1.y) : 1.0e10;
			double inv_chi_z = (KsiTensor.Str2.z > 1.0e-10) ? (1.0 / KsiTensor.Str2.z) : 1.0e10;

			inv_chi[offset + 0] = inv_chi_x;
			inv_chi[offset + 1] = inv_chi_y;
			inv_chi[offset + 2] = inv_chi_z;

			// RHS = H_ext + Mr/chi (same as LU solver lines 427-433)
			double Mr_over_chi_x = (KsiTensor.Str0.x > 1.0e-10) ? (MrVect.x / KsiTensor.Str0.x) : 0.0;
			double Mr_over_chi_y = (KsiTensor.Str1.y > 1.0e-10) ? (MrVect.y / KsiTensor.Str1.y) : 0.0;
			double Mr_over_chi_z = (KsiTensor.Str2.z > 1.0e-10) ? (MrVect.z / KsiTensor.Str2.z) : 0.0;

			rhs[offset + 0] = FlatExtern[offset + 0] + Mr_over_chi_x;
			rhs[offset + 1] = FlatExtern[offset + 1] + Mr_over_chi_y;
			rhs[offset + 2] = FlatExtern[offset + 2] + Mr_over_chi_z;
		}
		else
		{
			// For 6DOF: use isotropic chi (original behavior)
			double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
			if(chi < 1.0e-6) chi = 1.0e-6;
			double inv_chi_val = 1.0 / chi;

			for(int k = 0; k < dof; k++)
			{
				inv_chi[offset + k] = inv_chi_val;
				rhs[offset + k] = FlatExtern[offset + k];
			}
		}
	}

	// Update H-matrix diagonal with current inv_chi
	m_hacapk->UpdateDiagonal(inv_chi);

	// Initial guess
	for(int i = 0; i < totalDOF; i++)
	{
		sol[i] = FlatMagn[i];
	}

	// Build Jacobi preconditioner using H-matrix diagonal
	GetDiagonalElements_HMatrix_VariableDOF(diag_inv, inv_chi, totalDOF);
	for(int i = 0; i < totalDOF; i++)
	{
		diag_inv[i] = (std::abs(diag_inv[i]) > 1.0e-15) ? (1.0 / diag_inv[i]) : 1.0;
	}

	// Initialize: r0 = b - A*x0
	m_hacapk->MatVec(sol, v);  // v = A*x0 using H-matrix
	Copy(rhs, r, totalDOF);
	Axpy(-1.0, v, r, totalDOF);
	Copy(r, r0, totalDOF);

	double rho = 1.0, alpha_bicg = 1.0, omega = 1.0;
	std::fill(p.begin(), p.end(), 0.0);
	std::fill(v.begin(), v.end(), 0.0);

	double rhs_norm = Norm2(rhs, totalDOF);
	if(rhs_norm < 1.0e-30) rhs_norm = 1.0;

	int iter;
	for(iter = 1; iter <= max_iter; iter++)
	{
		double rho_old = rho;
		rho = Dot(r0, r, totalDOF);

		if(std::abs(rho) < 1.0e-30)
		{
			residual = Norm2(r, totalDOF) / rhs_norm;
			break;
		}

		if(iter == 1)
		{
			Copy(r, p, totalDOF);
		}
		else
		{
			if(std::abs(rho_old * omega) < 1.0e-30)
			{
				residual = Norm2(r, totalDOF) / rhs_norm;
				break;
			}
			double beta = (rho / rho_old) * (alpha_bicg / omega);
			Axpy(-omega, v, p, totalDOF);
			Scale(beta, p, totalDOF);
			Axpy(1.0, r, p, totalDOF);
		}

		// Apply preconditioner
		#pragma omp parallel for if(totalDOF > 100)
		for(int i = 0; i < totalDOF; i++)
		{
			p_hat[i] = diag_inv[i] * p[i];
		}

		// v = A * p_hat using H-matrix
		m_hacapk->MatVec(p_hat, v);

		double r0_dot_v = Dot(r0, v, totalDOF);
		if(std::abs(r0_dot_v) < 1.0e-30)
		{
			residual = Norm2(r, totalDOF) / rhs_norm;
			break;
		}
		alpha_bicg = rho / r0_dot_v;

		Copy(r, s, totalDOF);
		Axpy(-alpha_bicg, v, s, totalDOF);

		double s_norm = Norm2(s, totalDOF);
		if(s_norm / rhs_norm < tol)
		{
			Axpy(alpha_bicg, p_hat, sol, totalDOF);
			residual = s_norm / rhs_norm;
			break;
		}

		#pragma omp parallel for if(totalDOF > 100)
		for(int i = 0; i < totalDOF; i++)
		{
			s_hat[i] = diag_inv[i] * s[i];
		}

		// t = A * s_hat using H-matrix
		m_hacapk->MatVec(s_hat, t);

		double t_dot_s = Dot(t, s, totalDOF);
		double t_dot_t = Dot(t, t, totalDOF);
		if(std::abs(t_dot_t) < 1.0e-30)
		{
			Axpy(alpha_bicg, p_hat, sol, totalDOF);
			residual = s_norm / rhs_norm;
			break;
		}
		omega = t_dot_s / t_dot_t;

		Axpy(alpha_bicg, p_hat, sol, totalDOF);
		Axpy(omega, s_hat, sol, totalDOF);

		Copy(s, r, totalDOF);
		Axpy(-omega, t, r, totalDOF);

		double r_norm = Norm2(r, totalDOF);
		residual = r_norm / rhs_norm;
		if(residual < tol) break;

		if(std::abs(omega) < 1.0e-30) break;
	}

	// Copy solution back to flat array
	for(int i = 0; i < totalDOF; i++)
	{
		FlatMagn[i] = sol[i];
	}

	return iter;
}

//-------------------------------------------------------------------------
// GetDiagonalElements_HMatrix_VariableDOF
//-------------------------------------------------------------------------

void radTRelaxationMethNo_2::GetDiagonalElements_HMatrix_VariableDOF(std::vector<double>& diag,
                                                                      const std::vector<double>& inv_chi,
                                                                      int totalDOF)
{
	// Get diagonal elements A_ii = N_ii + 1/chi_i
	// Use cached N_ii values for efficiency (computed once during H-matrix build)
	if(m_hacapk->IsDiagonalCached())
	{
		const std::vector<double>& diag_N = m_hacapk->GetDiagonalN();
		#pragma omp parallel for if(totalDOF > 100)
		for(int i = 0; i < totalDOF; i++)
		{
			diag[i] = diag_N[i] + inv_chi[i];
		}
	}
	else
	{
		// Fallback: compute on-demand (slow)
		for(int i = 0; i < totalDOF; i++)
		{
			double N_ii = m_hacapk->GetInteractionMatrixElement(i, i);
			diag[i] = N_ii + inv_chi[i];
		}
	}
}

//-------------------------------------------------------------------------
// AutoRelax_VariableDOF (H-matrix version)
//-------------------------------------------------------------------------

int radTRelaxationMethNo_2::AutoRelax_VariableDOF(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	int AmOfMainElem = IntrctPtr->AmOfMainElem;
	if(AmOfMainElem <= 0) return 0;

	// HACApK supports both 3DOF tetrahedra and 6DOF hexahedra
	// HasVariableDOF() returns true only if DOF != 3, but HACApK works with uniform 3DOF too
	// For uniform 3DOF case, we need to set up the DOF tracking arrays
	if(!IntrctPtr->HasVariableDOF())
	{
		// Uniform 3DOF: set up DOF tracking manually
		// ComputeDOFOffsets sets up m_elemDOF, m_elemDOFOffset, m_totalDOF
		IntrctPtr->ComputeDOFOffsets();
		// Also need to set up flat arrays for HACApK
		IntrctPtr->SetupVariableDOFArrays();
	}

	int totalDOF = IntrctPtr->GetTotalDOF();
	if(totalDOF <= 0) return 0;

	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	double* FlatMagn = IntrctPtr->GetFlatMagnArray();
	double* FlatField = IntrctPtr->GetFlatFieldArray();
	double* FlatExtern = IntrctPtr->GetFlatExternFieldArray();

	if(FlatMagn == nullptr || FlatField == nullptr || FlatExtern == nullptr) return 0;

	// For 3DOF tetrahedra, we need the pre-computed InteractMatrix for element callbacks
	// This is because 3DOF uses component-to-component interaction which is expensive to compute on-demand
	// 6DOF hexahedra use on-demand computation (Yano-Sugahara MSC method)
	//
	// NOTE: This is currently O(N^2) and dominates the solve time for large 3DOF problems.
	// TODO: Implement on-demand 3x3 block computation to avoid full matrix pre-computation.
	bool need_precompute_matrix = !IntrctPtr->HasVariableDOF();  // uniform 3DOF
	if(need_precompute_matrix)
	{
		// First allocate InteractMatrix memory (may have been skipped when skipDenseMatrix=1)
		IntrctPtr->AllocateInteractMatrix();
		// Then compute the interaction matrix values
		IntrctPtr->SetupInteractMatrix();
	}

	// For 3DOF tetrahedra, initialize FlatExtern from ExternFieldArray
	if(need_precompute_matrix && IntrctPtr->ExternFieldArray != nullptr)
	{
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int offset = IntrctPtr->GetElementDOFOffset(elem);
			FlatExtern[offset + 0] = IntrctPtr->ExternFieldArray[elem].x;
			FlatExtern[offset + 1] = IntrctPtr->ExternFieldArray[elem].y;
			FlatExtern[offset + 2] = IntrctPtr->ExternFieldArray[elem].z;
		}
	}

	// Build H-matrix - must match current geometry
	// Delete old H-matrix if it was built for a different geometry
	if(m_hacapk && m_hacapk->GetInteraction() != IntrctPtr)
	{
		delete m_hacapk;
		m_hacapk = nullptr;
	}

	if(!m_hacapk)
	{
		m_hacapk = new RadHACApKManager(IntrctPtr);
	}

	// Reset timing statistics at start of solve
	rad.m_timing_hmatrix_build = 0.0;
	rad.m_timing_linear_solve = 0.0;
	rad.m_linear_iterations = 0;

	if(!m_hacapk->IsValid())
	{
		// Time H-matrix construction
		auto t_hmat_start = std::chrono::high_resolution_clock::now();

		if(!m_hacapk->BuildHMatrix(m_hacapk_params))
		{
			delete m_hacapk;
			m_hacapk = nullptr;
			return 0;
		}

		auto t_hmat_end = std::chrono::high_resolution_clock::now();
		rad.m_timing_hmatrix_build = std::chrono::duration<double>(t_hmat_end - t_hmat_start).count();
	}

	// For 3DOF tetrahedra: ensure flat N storage is ready after InteractMatrix is computed
	// This must happen AFTER InteractMatrix is set up (above) and AFTER BuildHMatrix
	// because BuildHMatrix may have returned early if InteractMatrix was NULL at that time
	if(need_precompute_matrix && !m_hacapk->IsFlatNReady())
	{
		m_hacapk->PrecomputeFlatInteractMatrix();
	}

	std::vector<double> OldMagn(totalDOF);
	double MisfitE2 = 1.0e30;
	int totalIterCount = 0;
	int outerIter = 0;

	// Initialize H field in NewFieldArray (used for chi(H) computation in nonlinear iteration)
	// Also initialize FlatField for compatibility
	const double H_init_mag = 1000.0;
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		if(dof == 3)
		{
			double H_ext_mag = 0.0;
			for(int k = 0; k < 3; k++)
			{
				H_ext_mag += FlatExtern[offset + k] * FlatExtern[offset + k];
			}
			H_ext_mag = std::sqrt(H_ext_mag);
			double scale = (H_ext_mag > 1.0e-10) ? std::min(1.0, H_init_mag / H_ext_mag) : 1.0;
			for(int k = 0; k < 3; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k] * scale;
			}
			// Initialize NewFieldArray for chi(H) computation
			IntrctPtr->NewFieldArray[elem].x = FlatExtern[offset + 0] * scale;
			IntrctPtr->NewFieldArray[elem].y = FlatExtern[offset + 1] * scale;
			IntrctPtr->NewFieldArray[elem].z = FlatExtern[offset + 2] * scale;
		}
		else if(dof == 6)
		{
			// 6DOF MSC hexahedra: initialize FlatField and estimate H from external field
			for(int k = 0; k < dof; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k];
			}
			// For 6DOF, estimate H from external field
			// The external field in RHS is H_ext (uniform applied field)
			// Estimate H direction from face normals weighted by sigma
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
			if(poly && poly->Use6DOF_MSC)
			{
				// Estimate H from sigma pattern: H ~ sum(sigma_i * n_i) / sum(n_i dot n_i)
				double Hx = 0.0, Hy = 0.0, Hz = 0.0;
				double wx = 0.0, wy = 0.0, wz = 0.0;
				for(int face = 0; face < 6; face++)
				{
					double sigma = FlatExtern[offset + face];
					TVector3d& n = poly->FaceNormal[face];
					Hx += sigma * n.x;
					Hy += sigma * n.y;
					Hz += sigma * n.z;
					wx += n.x * n.x;
					wy += n.y * n.y;
					wz += n.z * n.z;
				}
				IntrctPtr->NewFieldArray[elem].x = (wx > 1.0e-10) ? Hx / wx : 0.0;
				IntrctPtr->NewFieldArray[elem].y = (wy > 1.0e-10) ? Hy / wy : 0.0;
				IntrctPtr->NewFieldArray[elem].z = (wz > 1.0e-10) ? Hz / wz : 0.0;
			}
			else
			{
				// Fallback: use H_init_mag in z-direction
				IntrctPtr->NewFieldArray[elem].x = 0.0;
				IntrctPtr->NewFieldArray[elem].y = 0.0;
				IntrctPtr->NewFieldArray[elem].z = H_init_mag;
			}
		}
		else
		{
			for(int k = 0; k < dof; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k];
			}
			// Initialize NewFieldArray with H_init_mag in z-direction (fallback)
			IntrctPtr->NewFieldArray[elem].x = 0.0;
			IntrctPtr->NewFieldArray[elem].y = 0.0;
			IntrctPtr->NewFieldArray[elem].z = H_init_mag;
		}
	}

	// Cache polyhedron pointers for fast access (same as LU solver)
	std::vector<radTPolyhedron*> polyCache(AmOfMainElem, nullptr);
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
		polyCache[elem] = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
	}

	// Initialize CurrentChi with ELF-style initial value (same as LU solver lines 1494-1522)
	// This uses BH curve's 2nd point: chi = B2/(mu0*H2) - 1
	// Without this, CurrentChi starts at 1.0 causing slow convergence
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		if(dof == 6)
		{
			radTPolyhedron* poly = polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
				radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
				if(NonlinMater != nullptr)
				{
					double chi_init = NonlinMater->GetInitialChi_ELF_Style();
					if(chi_init > 0)
					{
						poly->CurrentChi = chi_init;
					}
				}
			}
		}
	}

	// B-field convergence tracking (same as LU solver, ELF mucal2)
	std::vector<double> OldBnorm(AmOfMainElem, 0.0);
	const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;

	// Outer nonlinear iteration (rewritten to match LU solver structure)
	for(outerIter = 0; outerIter < MaxIterNumber; outerIter++)
	{

		// Store old values
		for(int i = 0; i < totalDOF; i++)
		{
			OldMagn[i] = FlatMagn[i];
		}

		// Store old B norm for convergence check (same as LU solver)
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			if(dof == 6)
			{
				radTPolyhedron* poly = polyCache[elem];
				if(poly && poly->Use6DOF_MSC)
				{
					double chi = poly->CurrentChi;
					if(chi < 1.0e-6) chi = 1.0e-6;
					TVector3d& M = poly->Magn;
					TVector3d H(M.x / chi, M.y / chi, M.z / chi);
					TVector3d B(MU_0 * (H.x + M.x), MU_0 * (H.y + M.y), MU_0 * (H.z + M.z));
					OldBnorm[elem] = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
				}
			}
		}

		// Solve with BiCGSTAB using H-matrix
		// NOTE: max_iter for BiCGSTAB inner loop should be FIXED (not user's MaxIterNumber)
		// User's MaxIterNumber controls OUTER nonlinear iterations, not inner BiCGSTAB
		// bicg_tol is set via rad.SetBiCGSTABTol() Python API (default: 1e-4, ELF-compatible)
		double residual = 0.0;
		const double bicg_tol = rad.m_bicg_tol;
		const int bicg_max_iter = 10000;  // Inner BiCGSTAB max iterations (fixed)

		// Time BiCGSTAB solve
		auto t_bicg_start = std::chrono::high_resolution_clock::now();
		int n_iter = SolveBiCGSTAB_HMatrix_VariableDOF(totalDOF, bicg_tol, bicg_max_iter, residual);
		auto t_bicg_end = std::chrono::high_resolution_clock::now();
		rad.m_timing_linear_solve += std::chrono::duration<double>(t_bicg_end - t_bicg_start).count();
		rad.m_linear_iterations += n_iter;

		totalIterCount += n_iter;

		// Update element magnetization from flat array
		double M_diff_sq = 0.0;
		double M_norm_sq = 0.0;
		for(int i = 0; i < totalDOF; i++)
		{
			double diff = FlatMagn[i] - OldMagn[i];
			M_diff_sq += diff * diff;
			M_norm_sq += FlatMagn[i] * FlatMagn[i];
		}

		// Sync magnetization to element objects and compute H_new
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			int offset = IntrctPtr->GetElementDOFOffset(elem);
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];

			if(dof == 3)
			{
				g3dRelaxPtr->Magn.x = FlatMagn[offset];
				g3dRelaxPtr->Magn.y = FlatMagn[offset + 1];
				g3dRelaxPtr->Magn.z = FlatMagn[offset + 2];

				// Compute H_new = M / chi_current for 3DOF (same as LU solver)
				// This is needed for chi(H) update in nonlinear iteration
				if(IntrctPtr->NewFieldArray != nullptr)
				{
					radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
					if(MaterPtr != nullptr)
					{
						TVector3d& H_old = IntrctPtr->NewFieldArray[elem];
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_old, KsiTensor, MrVect);
						double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						if(chi < 1.0e-6) chi = 1.0e-6;
						IntrctPtr->NewFieldArray[elem].x = g3dRelaxPtr->Magn.x / chi;
						IntrctPtr->NewFieldArray[elem].y = g3dRelaxPtr->Magn.y / chi;
						IntrctPtr->NewFieldArray[elem].z = g3dRelaxPtr->Magn.z / chi;
					}
				}
			}
			else if(dof == 6)
			{
				radTPolyhedron* poly = polyCache[elem];
				if(poly && poly->Use6DOF_MSC)
				{
					// Store sigma values
					for(int k = 0; k < 6; k++)
					{
						poly->Sigma[k] = FlatMagn[offset + k];
					}

					// Compute effective magnetization from sigma (same as LU solver)
					double Mx = 0.0, My = 0.0, Mz = 0.0;
					double wx = 0.0, wy = 0.0, wz = 0.0;
					for(int face = 0; face < 6; face++)
					{
						double sigma = poly->Sigma[face];
						TVector3d& n = poly->FaceNormal[face];
						double nx2 = n.x * n.x;
						double ny2 = n.y * n.y;
						double nz2 = n.z * n.z;
						Mx += sigma * n.x;
						My += sigma * n.y;
						Mz += sigma * n.z;
						wx += nx2;
						wy += ny2;
						wz += nz2;
					}
					if(wx > 1.0e-10) poly->Magn.x = Mx / wx;
					if(wy > 1.0e-10) poly->Magn.y = My / wy;
					if(wz > 1.0e-10) poly->Magn.z = Mz / wz;

					// Compute H_new = M / chi_current (same as LU solver lines 1657-1668)
					if(IntrctPtr->NewFieldArray != nullptr)
					{
						double chi_used = poly->CurrentChi;
						if(chi_used < 1.0e-6) chi_used = 1.0e-6;
						IntrctPtr->NewFieldArray[elem].x = poly->Magn.x / chi_used;
						IntrctPtr->NewFieldArray[elem].y = poly->Magn.y / chi_used;
						IntrctPtr->NewFieldArray[elem].z = poly->Magn.z / chi_used;
					}
				}
			}
		}

		// Compute convergence and update chi (same structure as LU solver lines 1686-1793)
		double max_B_rel_change = 0.0;
		bool has_6dof_elements = false;

		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			if(dof == 6)
			{
				has_6dof_elements = true;
				radTPolyhedron* poly = polyCache[elem];
				if(poly && poly->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
				{
					TVector3d H_new = IntrctPtr->NewFieldArray[elem];
					radTMaterial* MaterPtr = (radTMaterial*)(IntrctPtr->g3dRelaxPtrVect[elem]->MaterHandle.rep);

					// Get chi used for this iteration's matrix (same as LU line 1699)
					double chi_matrix = poly->CurrentChi;
					double mu_old = chi_matrix + 1.0;

					// Compute H magnitude for chi update
					double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);

					// Use ELF-style dual-method chi update (same as LU lines 1708-1736)
					radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
					double chi_new;
					if(NonlinMater != nullptr)
					{
						chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old);
					}
					else
					{
						// Fallback for linear materials
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_new, KsiTensor, MrVect);
						chi_new = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
					}
					if(chi_new < 1.0e-6) chi_new = 1.0e-6;

					// Update chi for next iteration with optional under-relaxation
					// relax=0: full step, relax>0: chi = chi_new*(1-relax) + chi_old*relax
					double relax_hacapk = rad.m_relax;
					if(relax_hacapk > 0.0 && relax_hacapk <= 1.0)
					{
						chi_new = chi_new * (1.0 - relax_hacapk) + chi_matrix * relax_hacapk;
					}
					poly->CurrentChi = chi_new;

					// B-field convergence (same as LU lines 1741-1765)
					TVector3d& M_new = poly->Magn;
					double chi_for_B = chi_new;
					if(chi_for_B < 1.0e-6) chi_for_B = 1.0e-6;
					TVector3d H_for_B(M_new.x / chi_for_B, M_new.y / chi_for_B, M_new.z / chi_for_B);
					TVector3d B_new_vec(MU_0 * (H_for_B.x + M_new.x),
					                    MU_0 * (H_for_B.y + M_new.y),
					                    MU_0 * (H_for_B.z + M_new.z));
					double B_new_norm = std::sqrt(B_new_vec.x*B_new_vec.x + B_new_vec.y*B_new_vec.y + B_new_vec.z*B_new_vec.z);

					// Get B_sat from BH curve (same as LU lines 1753-1759)
					double B_sat = 1.0;
					if(NonlinMater != nullptr)
					{
						B_sat = NonlinMater->GetBsaturation();
						if(B_sat < 1.0e-10) B_sat = 1.0;
					}

					// B-field convergence: |B_new - B_old| / B_sat (same as LU lines 1761-1765)
					double B_old_norm = OldBnorm[elem];
					double B_rel_change = std::fabs(B_new_norm - B_old_norm) / B_sat;
					if(B_rel_change > max_B_rel_change)
						max_B_rel_change = B_rel_change;
				}
			}
		}

		// Convergence criterion (same as LU lines 1779-1794)
		double rel_change;
		if(has_6dof_elements)
		{
			// For 6DOF MSC: use ELF-style B-field change (mucal2)
			rel_change = max_B_rel_change;
		}
		else
		{
			// For 3DOF MMM: use M change
			rel_change = (M_norm_sq > 1.0e-30) ? std::sqrt(M_diff_sq / M_norm_sq) : std::sqrt(M_diff_sq);
		}
		MisfitE2 = rel_change * rel_change;

		if(rel_change <= PrecOnMagnetiz)
		{
			outerIter++;
			break;
		}

		if(radYield.Check() == 0) return outerIter;
	}

	IntrctPtr->RelaxStatusParam.MisfitM = std::sqrt(MisfitE2);
	return outerIter;
}

#endif // RADIA_USE_HACAPK
