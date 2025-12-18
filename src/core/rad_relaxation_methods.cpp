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

#include <time.h>

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
			// Row 0 of block (i,j): -Nij.Str0
			BaseMatrix_flat[(row_base + 0) + (col_base + 0)*ndof] = -Nij.Str0.x;
			BaseMatrix_flat[(row_base + 0) + (col_base + 1)*ndof] = -Nij.Str0.y;
			BaseMatrix_flat[(row_base + 0) + (col_base + 2)*ndof] = -Nij.Str0.z;

			// Row 1 of block (i,j): -Nij.Str1
			BaseMatrix_flat[(row_base + 1) + (col_base + 0)*ndof] = -Nij.Str1.x;
			BaseMatrix_flat[(row_base + 1) + (col_base + 1)*ndof] = -Nij.Str1.y;
			BaseMatrix_flat[(row_base + 1) + (col_base + 2)*ndof] = -Nij.Str1.z;

			// Row 2 of block (i,j): -Nij.Str2
			BaseMatrix_flat[(row_base + 2) + (col_base + 0)*ndof] = -Nij.Str2.x;
			BaseMatrix_flat[(row_base + 2) + (col_base + 1)*ndof] = -Nij.Str2.y;
			BaseMatrix_flat[(row_base + 2) + (col_base + 2)*ndof] = -Nij.Str2.z;
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
				double xj0 = x[3*j + 0];
				double xj1 = x[3*j + 1];
				double xj2 = x[3*j + 2];

				y0 -= Nij.Str0.x*xj0 + Nij.Str0.y*xj1 + Nij.Str0.z*xj2;
				y1 -= Nij.Str1.x*xj0 + Nij.Str1.y*xj1 + Nij.Str1.z*xj2;
				y2 -= Nij.Str2.x*xj0 + Nij.Str2.y*xj1 + Nij.Str2.z*xj2;
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
				// A[3i+k, 3j+l] = -N[i][j].Str_k.{x,y,z}[l]
				int row_base = 3*i;
				int col_base = 3*j;

				// Row 0 of block (i,j): -Nij.Str0
				A_flat[(row_base + 0) + (col_base + 0)*ndof] = -Nij.Str0.x;
				A_flat[(row_base + 0) + (col_base + 1)*ndof] = -Nij.Str0.y;
				A_flat[(row_base + 0) + (col_base + 2)*ndof] = -Nij.Str0.z;

				// Row 1 of block (i,j): -Nij.Str1
				A_flat[(row_base + 1) + (col_base + 0)*ndof] = -Nij.Str1.x;
				A_flat[(row_base + 1) + (col_base + 1)*ndof] = -Nij.Str1.y;
				A_flat[(row_base + 1) + (col_base + 2)*ndof] = -Nij.Str1.z;

				// Row 2 of block (i,j): -Nij.Str2
				A_flat[(row_base + 2) + (col_base + 0)*ndof] = -Nij.Str2.x;
				A_flat[(row_base + 2) + (col_base + 1)*ndof] = -Nij.Str2.y;
				A_flat[(row_base + 2) + (col_base + 2)*ndof] = -Nij.Str2.z;
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
		// IMPORTANT: Use stricter tolerance for inner BiCGSTAB
		// Separate tolerances: outer (e.g., 0.01) and inner bicg_tol (1e-6)
		// Using loose outer tolerance for BiCGSTAB causes premature convergence
		const double bicg_tol = 1.0e-6;  // Inner BiCGSTAB tolerance
		double residual = 0.0;
		int n_iter = SolveBiCGSTAB(ndof, bicg_tol, MaxIterNumber - totalIterCount, residual);
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

	// Negate only the 3 DOF blocks (MMM method needs -N)
	// MSC 6 DOF blocks already have correct sign (-K/(4pi))
	for(int row_elem = 0; row_elem < AmOfMainElem; row_elem++)
	{
		int dof_row = IntrctPtr->GetElementDOF(row_elem);
		int offset_row = IntrctPtr->GetElementDOFOffset(row_elem);

		for(int col_elem = 0; col_elem < AmOfMainElem; col_elem++)
		{
			int dof_col = IntrctPtr->GetElementDOF(col_elem);
			int offset_col = IntrctPtr->GetElementDOFOffset(col_elem);

			// Only negate blocks where at least one element is 3 DOF (MMM)
			// 6x6 blocks (both 6 DOF MSC) keep their sign
			if(dof_row == 3 || dof_col == 3)
			{
				for(int i = 0; i < dof_row; i++)
				{
					for(int j = 0; j < dof_col; j++)
					{
						int idx = (offset_row + i) * totalDOF + (offset_col + j);
						BaseMatrix[idx] = -BaseMatrix[idx];
					}
				}
			}
		}
	}

	// Store old values for convergence check
	std::vector<double> OldMagn(totalDOF);

	double PrecOnMagnetizE2 = PrecOnMagnetiz * PrecOnMagnetiz;
	double MisfitE2 = 1.0e30;
	int iterCount = 0;

	// Initialize H field
	const double H_init_mag = 1000.0;
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
		else
		{
			// MSC elements: initialize sigma to small values
			for(int k = 0; k < dof; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k];
			}
		}
	}

	// Work arrays
	std::vector<double> SystemMatrix(totalDOF * totalDOF);
	std::vector<double> RHS(totalDOF);

	for(iterCount = 0; iterCount < MaxIterNumber; iterCount++)
	{
		// Store old values
		for(int i = 0; i < totalDOF; i++)
		{
			OldMagn[i] = FlatMagn[i];
		}

		// Update H from M if not first iteration
		if(iterCount > 0)
		{
			for(int elem = 0; elem < AmOfMainElem; elem++)
			{
				int dof = IntrctPtr->GetElementDOF(elem);
				int offset = IntrctPtr->GetElementDOFOffset(elem);

				if(dof == 3)
				{
					// Standard element: H = M / chi
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
				// MSC elements: TODO - implement sigma -> H conversion
			}
		}

		// Copy base matrix and update diagonal with 1/chi
		// Base matrix already contains -K/(4pi) for 6 DOF MSC elements
		// Equation: (-K/(4pi) + 1/chi * I) * sigma = H_ext_n
		for(int i = 0; i < totalDOF * totalDOF; i++)
		{
			SystemMatrix[i] = BaseMatrix[i];
		}

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

				// Get average chi from material (assume linear isotropic for now)
				TVector3d H_est(0., 0., FlatExtern[offset]);  // Estimate H for material query
				TMatrix3d KsiTensor;
				TVector3d MrVect;
				MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);
				double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
				if(chi < 1.0e-6) chi = 1.0e-6;
				double inv_chi = 1.0 / chi;

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

		// Solve
		int ierr = SolveLU_Flat(SystemMatrix, RHS, totalDOF);
		if(ierr != 0) return iterCount;

		// Extract solution
		for(int i = 0; i < totalDOF; i++)
		{
			FlatMagn[i] = RHS[i];
		}

		// Compute convergence
		double M_diff_sq = 0.0;
		double M_norm_sq = 0.0;
		for(int i = 0; i < totalDOF; i++)
		{
			double dM = FlatMagn[i] - OldMagn[i];
			M_diff_sq += dM * dM;
			M_norm_sq += FlatMagn[i] * FlatMagn[i];
		}

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
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
				if(poly && poly->Use6DOF_MSC)
				{
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
				}
			}
		}

		double rel_change = (M_norm_sq > 1.0e-30) ? std::sqrt(M_diff_sq / M_norm_sq) : std::sqrt(M_diff_sq);
		MisfitE2 = rel_change * rel_change;

		if(rel_change <= PrecOnMagnetiz)
		{
			iterCount++;
			break;
		}

		if(radYield.Check() == 0) return iterCount;
	}

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
	// For 3 DOF elements (MMM): base matrix stores N, equation is (-N + 1/chi) * M = H_ext
	// For 6 DOF elements (MSC): base matrix stores -K/(4pi), equation is (-K/(4pi) + 1/chi) * sigma = H_ext_n
	//
	// Sign convention:
	// - 3 DOF blocks: need to negate (y -= N*x)
	// - 6 DOF blocks: already negative (y += (-K/(4pi))*x)

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

			// Get block from flat matrix
			const double* block = &FlatInteract[offset_row * totalDOF + offset_col];

			// Determine sign based on DOF types:
			// - 3x3 block (both 3 DOF): stores N, need to negate → sign = -1
			// - 6x6 block (both 6 DOF): stores -K/(4pi), use as-is → sign = +1
			// - 3x6 or 6x3 blocks: stores N or -K, need to negate → sign = -1
			double sign = (dof_row == 6 && dof_col == 6) ? 1.0 : -1.0;

			for(int i = 0; i < dof_row; i++)
			{
				double sum = 0.0;
				for(int j = 0; j < dof_col; j++)
				{
					sum += block[i * totalDOF + j] * x[offset_col + j];
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
	// For 3 DOF: interaction matrix stores N, diagonal = -N_ii + 1/chi
	// For 6 DOF: interaction matrix stores -K/(4pi), diagonal = -K/(4pi) + 1/chi (use as-is)
	const double* FlatInteract = IntrctPtr->GetFlatInteractMatrix();
	if(FlatInteract == nullptr) return;

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		// Get diagonal block
		const double* diag_block = &FlatInteract[offset * totalDOF + offset];

		// Sign convention: negate for 3 DOF (MMM), use as-is for 6 DOF (MSC)
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
			radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

			// Get average chi from material (assume linear isotropic for now)
			TVector3d H_est(0., 0., FlatExtern[offset]);  // Estimate H for material query
			TMatrix3d KsiTensor;
			TVector3d MrVect;
			MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);
			double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
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

	// Initial guess
	for(int i = 0; i < totalDOF; i++)
	{
		sol[i] = FlatMagn[i];
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
	double MisfitE2 = 1.0e30;
	int totalIterCount = 0;
	int outerIter = 0;

	// Initialize H field
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
		else
		{
			for(int k = 0; k < dof; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k];
			}
		}
	}

	// Outer nonlinear iteration
	for(outerIter = 0; outerIter < MaxIterNumber; outerIter++)
	{
		// Store old values
		for(int i = 0; i < totalDOF; i++)
		{
			OldMagn[i] = FlatMagn[i];
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
		double residual = 0.0;
		int n_iter = SolveBiCGSTAB_VariableDOF(totalDOF, PrecOnMagnetiz, MaxIterNumber - totalIterCount, residual);
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
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
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
					if(wx > 1.0e-10) g3dRelaxPtr->Magn.x = Mx / wx;
					if(wy > 1.0e-10) g3dRelaxPtr->Magn.y = My / wy;
					if(wz > 1.0e-10) g3dRelaxPtr->Magn.z = Mz / wz;
				}
			}
		}

		double rel_change = (M_norm_sq > 1.0e-30) ? std::sqrt(M_diff_sq / M_norm_sq) : std::sqrt(M_diff_sq);
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
