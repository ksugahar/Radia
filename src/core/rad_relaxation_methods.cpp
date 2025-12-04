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

#include <time.h>

//-------------------------------------------------------------------------

#ifdef _OPENMP
#include <omp.h>
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
	// Gaussian elimination with partial pivoting
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
}

int radTRelaxationMethNo_0::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

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

	// Store old magnetization for convergence check
	std::vector<TVector3d> OldMagnArray(AmOfMainElem);

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

	// Outer nonlinear iteration loop
	// For linear materials, this converges in 1 iteration
	// For nonlinear materials, chi(H) is updated each iteration
	for(iterCount = 0; iterCount < MaxIterNumber; iterCount++)
	{
		// Store old magnetization
		for(int i = 0; i < AmOfMainElem; i++)
		{
			OldMagnArray[i] = MagnAr[i];
		}

		// Update H field from current M: H = H_ext - N*M
		// Skip on first iteration (H already initialized above)
		if(iterCount > 0)
		{
			for(int i = 0; i < AmOfMainElem; i++)
			{
				TVector3d H_total = ExternFieldAr[i];
				for(int j = 0; j < AmOfMainElem; j++)
				{
					TMatrix3df& Nij = IntrcMat[i][j];
					H_total.x -= Nij.Str0.x*MagnAr[j].x + Nij.Str0.y*MagnAr[j].y + Nij.Str0.z*MagnAr[j].z;
					H_total.y -= Nij.Str1.x*MagnAr[j].x + Nij.Str1.y*MagnAr[j].y + Nij.Str1.z*MagnAr[j].z;
					H_total.z -= Nij.Str2.x*MagnAr[j].x + Nij.Str2.y*MagnAr[j].y + Nij.Str2.z*MagnAr[j].z;
				}
				NewFieldAr[i] = H_total;
			}
		}

		// Build system matrix using current chi(H)
		// A = -N + 1/chi(H), where chi depends on current H field
		std::vector<std::vector<double>> SystemMatrix(ndof, std::vector<double>(ndof, 0.0));
		std::vector<double> RHS(ndof, 0.0);

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

			// Build matrix rows for this element (3 rows: Mx, My, Mz)
			for(int comp_i = 0; comp_i < 3; comp_i++)
			{
				int row = 3 * i + comp_i;

				// Get chi for this component
				double chi_val = 0.0;
				if(comp_i == 0)      chi_val = KsiTensor.Str0.x;
				else if(comp_i == 1) chi_val = KsiTensor.Str1.y;
				else                 chi_val = KsiTensor.Str2.z;

				// Compute 1/chi (or large value if chi is very small)
				double inv_chi = (chi_val > 1.0e-10) ? (1.0 / chi_val) : 1.0e10;

				// Add N contributions: A = -N + 1/chi
				for(int j = 0; j < AmOfMainElem; j++)
				{
					TMatrix3df& Nij = IntrcMat[i][j];

					for(int comp_j = 0; comp_j < 3; comp_j++)
					{
						int col = 3 * j + comp_j;

						double Nij_val = 0.0;
						if(comp_i == 0)      Nij_val = (comp_j == 0) ? Nij.Str0.x : ((comp_j == 1) ? Nij.Str0.y : Nij.Str0.z);
						else if(comp_i == 1) Nij_val = (comp_j == 0) ? Nij.Str1.x : ((comp_j == 1) ? Nij.Str1.y : Nij.Str1.z);
						else                 Nij_val = (comp_j == 0) ? Nij.Str2.x : ((comp_j == 1) ? Nij.Str2.y : Nij.Str2.z);

						SystemMatrix[row][col] -= Nij_val;
					}
				}

				// Add 1/chi to diagonal
				SystemMatrix[row][row] += inv_chi;

				// RHS = H_ext + Mr/chi
				double Hext_comp = 0.0;
				if(comp_i == 0)      Hext_comp = ExternFieldAr[i].x;
				else if(comp_i == 1) Hext_comp = ExternFieldAr[i].y;
				else                 Hext_comp = ExternFieldAr[i].z;

				double Mr_comp = 0.0;
				if(comp_i == 0)      Mr_comp = MrVect.x;
				else if(comp_i == 1) Mr_comp = MrVect.y;
				else                 Mr_comp = MrVect.z;

				double Mr_over_chi = (chi_val > 1.0e-10) ? (Mr_comp / chi_val) : 0.0;
				RHS[row] = Hext_comp + Mr_over_chi;
			}
		}

		// Solve the linear system using LU decomposition
		int ierr = SolveLU(SystemMatrix, RHS, ndof);

		if(ierr != 0)
		{
			// Solver failed - singular matrix
			return iterCount;
		}

		// Extract LU solution (M values)
		for(int i = 0; i < AmOfMainElem; i++)
		{
			MagnAr[i].x = RHS[3 * i + 0];
			MagnAr[i].y = RHS[3 * i + 1];
			MagnAr[i].z = RHS[3 * i + 2];
		}

		// Newton-style M(H) update (Gauss-Seidel style):
		// After LU solves the linearized system, apply Newton-style correction
		// at each element. This uses the material's M(H) function directly.
		TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.);
		TMatrix3d E(E_Str0, E_Str1, E_Str2);

		for(int i = 0; i < AmOfMainElem; i++)
		{
			// Compute quasi-external field (field from all OTHER elements + external)
			TVector3d QuasiExtField = ExternFieldAr[i];
			TMatrix3df* MatrArrayPtr = IntrcMat[i];
			for(int j = 0; j < AmOfMainElem; j++)
			{
				if(j != i)
				{
					TMatrix3df& Nij = MatrArrayPtr[j];
					QuasiExtField.x += Nij.Str0.x*MagnAr[j].x + Nij.Str0.y*MagnAr[j].y + Nij.Str0.z*MagnAr[j].z;
					QuasiExtField.y += Nij.Str1.x*MagnAr[j].x + Nij.Str1.y*MagnAr[j].y + Nij.Str1.z*MagnAr[j].z;
					QuasiExtField.z += Nij.Str2.x*MagnAr[j].x + Nij.Str2.y*MagnAr[j].y + Nij.Str2.z*MagnAr[j].z;
				}
			}

			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
			radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

			// Get chi*Nii and Mr for current H estimate
			TMatrix3d MatrElemByInstKsi;
			TVector3d MatrElemByInstMr;
			MaterPtr->MultMatrByInstKsiAndMr(NewFieldAr[i], MatrArrayPtr[i], MatrElemByInstKsi, MatrElemByInstMr);

			// Solve local equation: H = (I - chi*Nii)^{-1} * (QuasiExtField + Mr)
			TMatrix3d BufMatr = E - MatrElemByInstKsi;
			TMatrix3d InvBufMatr;
			Matrix3d_inv(BufMatr, InvBufMatr);
			NewFieldAr[i] = InvBufMatr * (QuasiExtField + MatrElemByInstMr);

			// Use material's M(H) function directly
			MagnAr[i] = MaterPtr->M(NewFieldAr[i]);
		}

		// Compute convergence (MisfitE2) and update object magnetizations
		MisfitE2 = 0.0;
		for(int i = 0; i < AmOfMainElem; i++)
		{
			// Compute convergence (change in M)
			TVector3d dM = MagnAr[i] - OldMagnArray[i];
			MisfitE2 += dM.AmpE2();

			// Update the object's magnetization
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
			g3dRelaxPtr->Magn = MagnAr[i];
		}

		MisfitE2 /= AmOfMainElem;

		// Check convergence
		if(MisfitE2 <= PrecOnMagnetizE2)
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
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum) if(n > 100)
	for(int i = 0; i < n; i++)
	{
		sum += a[i] * b[i];
	}
	return sum;
}

double radTRelaxationMethNo_1::Norm2(const std::vector<double>& a, int n)
{
	return std::sqrt(Dot(a, a, n));
}

void radTRelaxationMethNo_1::Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n)
{
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		y[i] += alpha * x[i];
	}
}

void radTRelaxationMethNo_1::Copy(const std::vector<double>& src, std::vector<double>& dst, int n)
{
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		dst[i] = src[i];
	}
}

void radTRelaxationMethNo_1::Scale(double alpha, std::vector<double>& x, int n)
{
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		x[i] *= alpha;
	}
}

void radTRelaxationMethNo_1::GetDiagonalElements(std::vector<double>& diag, int n_elem)
{
	// Extract diagonal elements from interaction matrix for Jacobi preconditioner
	// Diagonal block [i][i] is a 3x3 matrix, we extract the diagonal of that
	// For nonlinear materials, chi is computed from current H stored in NewFieldArray
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;

	for(int i = 0; i < n_elem; i++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
		radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		// Use current H for chi(H) computation (important for nonlinear materials)
		TVector3d InstH = (NewFieldAr != nullptr) ? NewFieldAr[i] : TVector3d(0., 0., 0.);
		TMatrix3d KsiTensor;
		TVector3d MrVect;
		MaterPtr->DefineInstantKsiTensor(InstH, KsiTensor, MrVect);

		// For each component
		double chi_x = (KsiTensor.Str0.x > 1.0e-10) ? KsiTensor.Str0.x : 1.0e10;
		double chi_y = (KsiTensor.Str1.y > 1.0e-10) ? KsiTensor.Str1.y : 1.0e10;
		double chi_z = (KsiTensor.Str2.z > 1.0e-10) ? KsiTensor.Str2.z : 1.0e10;

		// Diagonal of system matrix: A = -N + 1/chi
		if(IntrcMat != nullptr)
		{
			// Nii is from InteractMatrix[i][i]
			TMatrix3df& Nii = IntrcMat[i][i];
			diag[3*i + 0] = -Nii.Str0.x + 1.0/chi_x;
			diag[3*i + 1] = -Nii.Str1.y + 1.0/chi_y;
			diag[3*i + 2] = -Nii.Str2.z + 1.0/chi_z;
		}
		else
		{
			// Fallback: just use 1/chi as diagonal (no N contribution)
			diag[3*i + 0] = 1.0/chi_x;
			diag[3*i + 1] = 1.0/chi_y;
			diag[3*i + 2] = 1.0/chi_z;
		}
	}
}

void radTRelaxationMethNo_1::DenseMatVec(const std::vector<double>& x, std::vector<double>& y, int ndof)
{
	// Computes y = A * x where A = -N + 1/chi
	// Uses H-matrix if available, otherwise dense matrix
	// For nonlinear materials, chi is computed from current H stored in NewFieldArray
	int n_elem = ndof / 3;
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;

	// Initialize y to zero
	std::fill(y.begin(), y.end(), 0.0);

	if(IntrcMat != nullptr)
	{
		// Dense matrix-vector product
		#pragma omp parallel for if(n_elem > 50)
		for(int i = 0; i < n_elem; i++)
		{
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
			radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

			// Use current H for chi(H) computation (important for nonlinear materials)
			TVector3d InstH = (NewFieldAr != nullptr) ? NewFieldAr[i] : TVector3d(0., 0., 0.);
			TMatrix3d KsiTensor;
			TVector3d MrVect;
			MaterPtr->DefineInstantKsiTensor(InstH, KsiTensor, MrVect);

			double inv_chi_x = (KsiTensor.Str0.x > 1.0e-10) ? 1.0/KsiTensor.Str0.x : 1.0e10;
			double inv_chi_y = (KsiTensor.Str1.y > 1.0e-10) ? 1.0/KsiTensor.Str1.y : 1.0e10;
			double inv_chi_z = (KsiTensor.Str2.z > 1.0e-10) ? 1.0/KsiTensor.Str2.z : 1.0e10;

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

int radTRelaxationMethNo_1::SolveBiCGSTAB(int ndof, double tol, int max_iter, double& residual)
{
	// BiCGSTAB with Jacobi preconditioner
	// Reference: van der Vorst, SIAM J. Sci. Stat. Comput. 13 (1992)
	// For nonlinear materials, chi is computed from current H stored in NewFieldArray

	int n_elem = ndof / 3;

	// Allocate work vectors
	std::vector<double> r(ndof), r0(ndof), p(ndof), v(ndof), s(ndof), t(ndof);
	std::vector<double> p_hat(ndof), s_hat(ndof), diag_inv(ndof);

	// Get RHS vector: b = H_external + Mr/chi(H)
	// For linear materials, Mr = 0 so b = H_external
	// For nonlinear materials with remanence, Mr/chi term is needed
	std::vector<double> rhs(ndof);
	if(IntrctPtr->ExternFieldArray == nullptr) return 0;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	for(int i = 0; i < n_elem; i++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
		radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		// Use current H for chi(H) computation
		TVector3d InstH = (NewFieldAr != nullptr) ? NewFieldAr[i] : TVector3d(0., 0., 0.);
		TMatrix3d KsiTensor;
		TVector3d MrVect;
		MaterPtr->DefineInstantKsiTensor(InstH, KsiTensor, MrVect);

		// b = H_ext + Mr/chi
		double inv_chi_x = (KsiTensor.Str0.x > 1.0e-10) ? 1.0/KsiTensor.Str0.x : 0.0;
		double inv_chi_y = (KsiTensor.Str1.y > 1.0e-10) ? 1.0/KsiTensor.Str1.y : 0.0;
		double inv_chi_z = (KsiTensor.Str2.z > 1.0e-10) ? 1.0/KsiTensor.Str2.z : 0.0;

		rhs[3*i + 0] = IntrctPtr->ExternFieldArray[i].x + MrVect.x * inv_chi_x;
		rhs[3*i + 1] = IntrctPtr->ExternFieldArray[i].y + MrVect.y * inv_chi_y;
		rhs[3*i + 2] = IntrctPtr->ExternFieldArray[i].z + MrVect.z * inv_chi_z;
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
	GetDiagonalElements(diag_inv, n_elem);
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

	// Initialize: r0 = b - A*x0
	DenseMatVec(sol, v, ndof);  // v = A*x0
	Copy(rhs, r, ndof);         // r = rhs
	Axpy(-1.0, v, r, ndof);     // r = r - v

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
		for(int i = 0; i < ndof; i++)
		{
			p_hat[i] = diag_inv[i] * p[i];
		}

		// v = A * p_hat
		DenseMatVec(p_hat, v, ndof);

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
		for(int i = 0; i < ndof; i++)
		{
			s_hat[i] = diag_inv[i] * s[i];
		}

		// t = A * s_hat
		DenseMatVec(s_hat, t, ndof);

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

	// Store old magnetization for convergence checking
	std::vector<TVector3d> OldMagnArray(AmOfMainElem);

	double PrecOnMagnetizE2 = PrecOnMagnetiz * PrecOnMagnetiz;
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
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;

	for(outerIter = 0; outerIter < MaxIterNumber; outerIter++)
	{
		// Store old magnetization for convergence check
		for(int i = 0; i < AmOfMainElem; i++)
		{
			OldMagnArray[i] = MagnAr[i];
		}

		// Update H field from current M: H = H_ext - N*M
		// This H is used to compute chi(H) for the linear system
		// Skip on first iteration (H already initialized above)
		if(outerIter > 0)
		{
			for(int i = 0; i < AmOfMainElem; i++)
			{
				TVector3d H_total = ExternFieldAr[i];
				for(int j = 0; j < AmOfMainElem; j++)
				{
					TMatrix3df& Nij = IntrcMat[i][j];
					H_total.x -= Nij.Str0.x*MagnAr[j].x + Nij.Str0.y*MagnAr[j].y + Nij.Str0.z*MagnAr[j].z;
					H_total.y -= Nij.Str1.x*MagnAr[j].x + Nij.Str1.y*MagnAr[j].y + Nij.Str1.z*MagnAr[j].z;
					H_total.z -= Nij.Str2.x*MagnAr[j].x + Nij.Str2.y*MagnAr[j].y + Nij.Str2.z*MagnAr[j].z;
				}
				NewFieldAr[i] = H_total;
			}
		}

		// Solve linear system using BiCGSTAB with current chi(H)
		// System: A*M = b where A = -N + 1/chi(H), b = H_ext + Mr/chi
		double residual = 0.0;
		int n_iter = SolveBiCGSTAB(ndof, PrecOnMagnetiz * 0.1, MaxIterNumber - totalIterCount, residual);
		totalIterCount += n_iter;

		// Newton-style M(H) update (Gauss-Seidel style):
		// After BiCGSTAB solves the linearized system, apply Newton-style correction
		// at each element. This uses the material's M(H) function directly.
		// Key: compute quasi-external field (excluding self-interaction), solve for H
		// using the local (I - chi*Nii)^{-1} inversion, then M = M(H).
		TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.);
		TMatrix3d E(E_Str0, E_Str1, E_Str2);

		for(int i = 0; i < AmOfMainElem; i++)
		{
			// Compute quasi-external field (field from all OTHER elements + external)
			TVector3d QuasiExtField = ExternFieldAr[i];
			TMatrix3df* MatrArrayPtr = IntrcMat[i];
			for(int j = 0; j < AmOfMainElem; j++)
			{
				if(j != i)
				{
					TMatrix3df& Nij = MatrArrayPtr[j];
					QuasiExtField.x += Nij.Str0.x*MagnAr[j].x + Nij.Str0.y*MagnAr[j].y + Nij.Str0.z*MagnAr[j].z;
					QuasiExtField.y += Nij.Str1.x*MagnAr[j].x + Nij.Str1.y*MagnAr[j].y + Nij.Str1.z*MagnAr[j].z;
					QuasiExtField.z += Nij.Str2.x*MagnAr[j].x + Nij.Str2.y*MagnAr[j].y + Nij.Str2.z*MagnAr[j].z;
				}
			}

			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
			radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

			// Get chi*Nii and Mr for current H estimate
			TMatrix3d MatrElemByInstKsi, KsiTensor;
			TVector3d MatrElemByInstMr, MrVect;
			MaterPtr->MultMatrByInstKsiAndMr(NewFieldAr[i], MatrArrayPtr[i], MatrElemByInstKsi, MatrElemByInstMr);

			// Solve local equation: H = (I - chi*Nii)^{-1} * (QuasiExtField + Mr)
			TMatrix3d BufMatr = E - MatrElemByInstKsi;
			TMatrix3d InvBufMatr;
			Matrix3d_inv(BufMatr, InvBufMatr);
			NewFieldAr[i] = InvBufMatr * (QuasiExtField + MatrElemByInstMr);

			// Use material's M(H) function directly
			MagnAr[i] = MaterPtr->M(NewFieldAr[i]);
		}

		// Compute magnetization change (MisfitE2) and update object magnetizations
		MisfitE2 = 0.0;
		for(int i = 0; i < AmOfMainElem; i++)
		{
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];

			// Compute difference
			TVector3d dM;
			dM.x = MagnAr[i].x - OldMagnArray[i].x;
			dM.y = MagnAr[i].y - OldMagnArray[i].y;
			dM.z = MagnAr[i].z - OldMagnArray[i].z;
			MisfitE2 += dM.x*dM.x + dM.y*dM.y + dM.z*dM.z;

			// Update the object's magnetization
			g3dRelaxPtr->Magn = MagnAr[i];
		}

		MisfitE2 /= AmOfMainElem;

		// Check convergence
		if(MisfitE2 <= PrecOnMagnetizE2)
		{
			outerIter++;
			break;
		}

		// Allow multitasking
		if(radYield.Check() == 0) return totalIterCount;
	}

	// Update relaxation status
	IntrctPtr->RelaxStatusParam.MisfitM = std::sqrt(MisfitE2);
	ComputeRelaxStatusParam(MagnAr, OldMagnArray.data(), NewFieldAr);

	return totalIterCount;
}
