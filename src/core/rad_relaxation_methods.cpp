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
// Intel MKL headers (includes dgesv_ declaration)
#include "mkl_cblas.h"
#include "mkl_lapack.h"
#include "mkl_trans.h"   // For mkl_domatcopy (matrix transpose)
#include "mkl_service.h" // For mkl_set_num_threads, mkl_get_max_threads
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
// Newton-style M(H) update is now integrated into unified VariableDOF solvers
// that handle both 3DOF (tetra) and 6DOF (hex) elements
//-------------------------------------------------------------------------

//=========================================================================
// Unified Nonlinear Iteration Helper Functions
// These functions encapsulate the common logic shared across all solvers
//=========================================================================

bool InitializeNonlinearContext(NonlinearContext& ctx, radTInteraction* IntrctPtr, bool MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return false;

	// Reset if needed
	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	ctx.AmOfMainElem = IntrctPtr->AmOfMainElem;
	if(ctx.AmOfMainElem <= 0) return false;

	ctx.totalDOF = IntrctPtr->GetTotalDOF();
	if(ctx.totalDOF <= 0) return false;

	// Get flat arrays (owned by IntrctPtr)
	ctx.FlatMagn = IntrctPtr->GetFlatMagnArray();
	ctx.FlatField = IntrctPtr->GetFlatFieldArray();
	ctx.FlatExtern = IntrctPtr->GetFlatExternFieldArray();

	if(ctx.FlatMagn == nullptr || ctx.FlatField == nullptr || ctx.FlatExtern == nullptr)
		return false;

	// Allocate state vectors
	ctx.OldMagn.resize(ctx.totalDOF);
	ctx.OldChi.resize(ctx.AmOfMainElem, 0.0);
	ctx.OldBnorm.resize(ctx.AmOfMainElem, 0.0);
	ctx.CurrentChiArray.resize(ctx.AmOfMainElem, 1.0);
	ctx.NewFieldArray.resize(ctx.totalDOF);
	ctx.polyCache.resize(ctx.AmOfMainElem, nullptr);

	// Initialize flags
	ctx.all_materials_linear = true;
	ctx.relax_param = rad.m_relax;

	// Cache polyhedron pointers
	for(int elem = 0; elem < ctx.AmOfMainElem; elem++)
	{
		ctx.polyCache[elem] = dynamic_cast<radTPolyhedron*>(IntrctPtr->g3dRelaxPtrVect[elem]);
	}

	// Initialize chi and H field (ELF mucal0 style)
	const double H_init_mag = 100.0;
	for(int elem = 0; elem < ctx.AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
		radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		// Get initial chi (ELF style - from BH curve 2nd point)
		double chi_init = 1.0;
		radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
		if(NonlinMater != nullptr)
		{
			chi_init = NonlinMater->GetInitialChi_ELF_Style();
			if(chi_init <= 0) chi_init = 1.0;
			ctx.all_materials_linear = false;
			ctx.B_sat = NonlinMater->GetBsaturation();
			if(ctx.B_sat < 1.0e-10) ctx.B_sat = 1.0;
		}
		else
		{
			// Linear material: use DefineInstantKsiTensor
			TVector3d H_est(0., 0., H_init_mag);
			TMatrix3d KsiTensor;
			TVector3d MrVect;
			MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);
			chi_init = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
			if(chi_init < 1.0e-6) chi_init = 1.0e-6;
		}
		ctx.CurrentChiArray[elem] = chi_init;

		// Store in poly->CurrentChi for 6DOF elements
		if(dof == 6)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				poly->CurrentChi = chi_init;
			}
		}

		// Initialize H field
		if(dof == 3)
		{
			double H_ext_mag = 0.0;
			for(int k = 0; k < 3; k++)
			{
				H_ext_mag += ctx.FlatExtern[offset + k] * ctx.FlatExtern[offset + k];
			}
			H_ext_mag = std::sqrt(H_ext_mag);
			double scale = (H_ext_mag > 1.0e-10) ? std::min(1.0, H_init_mag / H_ext_mag) : 1.0;
			for(int k = 0; k < 3; k++)
			{
				ctx.FlatField[offset + k] = ctx.FlatExtern[offset + k] * scale;
			}
		}
		else if(dof == 6)
		{
			for(int k = 0; k < dof; k++)
			{
				ctx.FlatField[offset + k] = ctx.FlatExtern[offset + k];
			}
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
				ctx.FlatField[offset + k] = ctx.FlatExtern[offset + k];
			}
		}
	}

	return true;
}

//-------------------------------------------------------------------------

void BuildBaseMatrix(NonlinearContext& ctx, radTInteraction* IntrctPtr)
{
	double* FlatInteract = IntrctPtr->GetFlatInteractMatrix();
	if(FlatInteract == nullptr) return;

	ctx.BaseMatrix.resize(ctx.totalDOF * ctx.totalDOF);

	// Copy interaction matrix
	std::memcpy(ctx.BaseMatrix.data(), FlatInteract, ctx.totalDOF * ctx.totalDOF * sizeof(double));

	// Sign convention:
	// MMM (3DOF): FlatInteract stores N, need -N -> negate
	// MSC (6DOF): FlatInteract stores -K/(4pi) -> use as-is
	for(int row_elem = 0; row_elem < ctx.AmOfMainElem; row_elem++)
	{
		int dof_row = IntrctPtr->GetElementDOF(row_elem);
		int offset_row = IntrctPtr->GetElementDOFOffset(row_elem);

		for(int col_elem = 0; col_elem < ctx.AmOfMainElem; col_elem++)
		{
			int dof_col = IntrctPtr->GetElementDOF(col_elem);
			int offset_col = IntrctPtr->GetElementDOFOffset(col_elem);

			// Negate blocks involving 3DOF (MMM) elements
			if(dof_row == 3 || dof_col == 3)
			{
				for(int i = 0; i < dof_row; i++)
				{
					for(int j = 0; j < dof_col; j++)
					{
						// Column-major: (row, col) = [col * totalDOF + row]
						int idx = (offset_col + j) * ctx.totalDOF + (offset_row + i);
						ctx.BaseMatrix[idx] = -ctx.BaseMatrix[idx];
					}
				}
			}
		}
	}
}

//-------------------------------------------------------------------------

void StoreOldValuesAndComputeBnorm(NonlinearContext& ctx, radTInteraction* IntrctPtr)
{
	const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;

	// Store old magnetization
	for(int i = 0; i < ctx.totalDOF; i++)
	{
		ctx.OldMagn[i] = ctx.FlatMagn[i];
	}

	// Store old chi and B-norm for each element
	for(int elem = 0; elem < ctx.AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		if(dof == 3)
		{
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			TVector3d M = g3dRelaxPtr->Magn;
			TVector3d H(ctx.FlatField[offset], ctx.FlatField[offset+1], ctx.FlatField[offset+2]);
			TVector3d B(MU_0 * (H.x + M.x), MU_0 * (H.y + M.y), MU_0 * (H.z + M.z));
			ctx.OldBnorm[elem] = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
		}
		else if(dof == 6)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				ctx.OldChi[elem] = poly->CurrentChi;
				TVector3d M = poly->Magn;
				double chi = poly->CurrentChi;
				if(chi < 1.0e-6) chi = 1.0e-6;
				TVector3d H(M.x / chi, M.y / chi, M.z / chi);
				TVector3d B(MU_0 * (H.x + M.x), MU_0 * (H.y + M.y), MU_0 * (H.z + M.z));
				ctx.OldBnorm[elem] = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
			}
		}
	}
}

//-------------------------------------------------------------------------

void UpdateMagnAndComputeH(NonlinearContext& ctx, radTInteraction* IntrctPtr)
{
	for(int elem = 0; elem < ctx.AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		if(dof == 3)
		{
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			g3dRelaxPtr->Magn.x = ctx.FlatMagn[offset + 0];
			g3dRelaxPtr->Magn.y = ctx.FlatMagn[offset + 1];
			g3dRelaxPtr->Magn.z = ctx.FlatMagn[offset + 2];

			// H = M / chi
			double chi = ctx.CurrentChiArray[elem];
			if(chi < 1.0e-6) chi = 1.0e-6;
			for(int k = 0; k < 3; k++)
			{
				ctx.FlatField[offset + k] = ctx.FlatMagn[offset + k] / chi;
			}
		}
		else if(dof == 6)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				// Store sigma values
				for(int k = 0; k < 6; k++)
				{
					poly->Sigma[k] = ctx.FlatMagn[offset + k];
				}

				// Compute effective magnetization from sigma
				double Mx = 0.0, My = 0.0, Mz = 0.0;
				double wx = 0.0, wy = 0.0, wz = 0.0;
				for(int face = 0; face < 6; face++)
				{
					double sigma = poly->Sigma[face];
					TVector3d& n = poly->FaceNormal[face];
					Mx += sigma * n.x;
					My += sigma * n.y;
					Mz += sigma * n.z;
					wx += n.x * n.x;
					wy += n.y * n.y;
					wz += n.z * n.z;
				}
				if(wx > 1.0e-10) poly->Magn.x = Mx / wx;
				if(wy > 1.0e-10) poly->Magn.y = My / wy;
				if(wz > 1.0e-10) poly->Magn.z = Mz / wz;

				// H = M / chi
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
}

//-------------------------------------------------------------------------

double UpdateChiAndCheckConvergence(NonlinearContext& ctx, radTInteraction* IntrctPtr)
{
	const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;
	double max_B_rel_change = 0.0;

	for(int elem = 0; elem < ctx.AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
		radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
		radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);

		double chi_matrix = ctx.CurrentChiArray[elem];
		double mu_old = chi_matrix + 1.0;
		double chi_new;

		if(dof == 3)
		{
			TVector3d H_new(ctx.FlatField[offset], ctx.FlatField[offset+1], ctx.FlatField[offset+2]);
			double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);

			if(NonlinMater != nullptr)
			{
				chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old, ctx.relax_param);
			}
			else
			{
				chi_new = chi_matrix;  // Linear: keep constant
			}
			ctx.CurrentChiArray[elem] = chi_new;

			// B-field convergence
			TVector3d M_new = g3dRelaxPtr->Magn;
			TVector3d B_new(MU_0 * (H_new.x + M_new.x), MU_0 * (H_new.y + M_new.y), MU_0 * (H_new.z + M_new.z));
			double B_new_norm = std::sqrt(B_new.x*B_new.x + B_new.y*B_new.y + B_new.z*B_new.z);

			double B_sat = ctx.B_sat;
			if(NonlinMater != nullptr)
			{
				B_sat = NonlinMater->GetBsaturation();
				if(B_sat < 1.0e-10) B_sat = 1.0;
			}

			double B_rel_change = std::fabs(B_new_norm - ctx.OldBnorm[elem]) / B_sat;
			if(B_rel_change > max_B_rel_change)
				max_B_rel_change = B_rel_change;
		}
		else if(dof == 6)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
			{
				TVector3d H_new = IntrctPtr->NewFieldArray[elem];
				double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);

				if(NonlinMater != nullptr)
				{
					chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old, ctx.relax_param);
				}
				else
				{
					chi_new = chi_matrix;
				}
				poly->CurrentChi = chi_new;
				ctx.CurrentChiArray[elem] = chi_new;

				// B-field convergence
				TVector3d M_new = poly->Magn;
				double chi_for_B = chi_new;
				if(chi_for_B < 1.0e-6) chi_for_B = 1.0e-6;
				TVector3d H_for_B(M_new.x / chi_for_B, M_new.y / chi_for_B, M_new.z / chi_for_B);
				TVector3d B_new(MU_0 * (H_for_B.x + M_new.x), MU_0 * (H_for_B.y + M_new.y), MU_0 * (H_for_B.z + M_new.z));
				double B_new_norm = std::sqrt(B_new.x*B_new.x + B_new.y*B_new.y + B_new.z*B_new.z);

				double B_sat = ctx.B_sat;
				if(NonlinMater != nullptr)
				{
					B_sat = NonlinMater->GetBsaturation();
					if(B_sat < 1.0e-10) B_sat = 1.0;
				}

				double B_rel_change = std::fabs(B_new_norm - ctx.OldBnorm[elem]) / B_sat;
				if(B_rel_change > max_B_rel_change)
					max_B_rel_change = B_rel_change;
			}
		}
	}

	ctx.max_B_rel_change = max_B_rel_change;
	return max_B_rel_change;
}

//=========================================================================
// Unified Nonlinear Iteration (base class implementation)
// Calls virtual SolveLinearStep which is overridden by each solver
//=========================================================================

int radTIterativeRelaxMeth::AutoRelax_Unified(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	// Initialize context
	NonlinearContext ctx;
	if(!InitializeNonlinearContext(ctx, IntrctPtr, MagnResetIsNotNeeded))
		return 0;

	// Build base matrix (geometric part without chi)
	BuildBaseMatrix(ctx, IntrctPtr);

	int iterCount = 0;
	double MisfitE2 = 1.0e30;

	// Nonlinear iteration loop
	for(iterCount = 0; iterCount < MaxIterNumber; iterCount++)
	{
		// Store old values for convergence check
		StoreOldValuesAndComputeBnorm(ctx, IntrctPtr);

		// Solve linear system (virtual - overridden by LU, BiCGSTAB, HACApK)
		int linearIter = SolveLinearStep(ctx, iterCount);
		(void)linearIter;  // May be used for statistics

		// Update element magnetization from solution
		UpdateMagnAndComputeH(ctx, IntrctPtr);

		// Update chi and check convergence
		double rel_change = UpdateChiAndCheckConvergence(ctx, IntrctPtr);
		MisfitE2 = rel_change * rel_change;

		// Linear materials: converge in exactly 1 iteration
		if(ctx.all_materials_linear)
		{
			iterCount++;
			break;
		}

		// Check convergence
		if(rel_change <= PrecOnMagnetiz)
		{
			iterCount++;
			break;
		}

		// Check for user abort
		if(radYield.Check() == 0)
		{
			return iterCount;
		}
	}

	IntrctPtr->RelaxStatusParam.MisfitM = std::sqrt(MisfitE2);

	return iterCount;
}

//=========================================================================
// Method 0: LU Direct Solver (Entry Point)
//=========================================================================

int radTRelaxationMethNo_0::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	// Use unified nonlinear iteration with LU linear solve
	return AutoRelax_Unified(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
}

//=========================================================================
// Method 1: BiCGSTAB Iterative Solver (default)
//=========================================================================

double radTRelaxationMethNo_1::Dot(const std::vector<double>& a, const std::vector<double>& b, int n)
{
#ifdef HAVE_LAPACK
	// Use Intel MKL CBLAS cblas_ddot for optimized dot product
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
	// Use Intel MKL CBLAS cblas_dnrm2 for optimized norm
	return cblas_dnrm2(n, a.data(), 1);
#else
	return std::sqrt(Dot(a, a, n));
#endif
}

void radTRelaxationMethNo_1::Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n)
{
#ifdef HAVE_LAPACK
	// Use Intel MKL CBLAS cblas_daxpy: y = alpha*x + y
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
	// Use Intel MKL CBLAS cblas_dcopy for optimized copy
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
	// Use Intel MKL CBLAS cblas_dscal for optimized scale
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
	// Use Intel MKL CBLAS cblas_dgemv: y = alpha*A*x + beta*y
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

int radTRelaxationMethNo_1::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	// Use unified nonlinear iteration with BiCGSTAB linear solve
	return AutoRelax_Unified(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
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

//-------------------------------------------------------------------------
// LU Solver: SolveLinearStep override
// Builds system matrix with current chi, solves with LU, stores result
//-------------------------------------------------------------------------

int radTRelaxationMethNo_0::SolveLinearStep(NonlinearContext& ctx, int iterCount)
{
	int totalDOF = ctx.totalDOF;
	int AmOfMainElem = ctx.AmOfMainElem;

	// Build system matrix: copy base matrix and add diagonal terms
	// LAPACK dgesv destroys the matrix, so we need a working copy
	// Only one copy is needed - dgesv works directly on this
	std::vector<double> SystemMatrix(totalDOF * totalDOF);
	std::memcpy(SystemMatrix.data(), ctx.BaseMatrix.data(), totalDOF * totalDOF * sizeof(double));

	// Build RHS vector (will be overwritten with solution by dgesv)
	std::vector<double> RHS(totalDOF);

	// Update diagonal and RHS based on current chi
	// Matrix is COLUMN-MAJOR: A(i,j) at [j * totalDOF + i]
	// Diagonal element A(k,k) is at [k * totalDOF + k] = [k * (totalDOF + 1)]
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		double chi = ctx.CurrentChiArray[elem];
		if(chi < 1.0e-6) chi = 1.0e-6;
		double inv_chi = 1.0 / chi;

		// Add 1/chi to diagonal and set RHS
		for(int k = 0; k < dof; k++)
		{
			int row = offset + k;
			// Column-major diagonal: A(row,row) at index [row * (totalDOF + 1)]
			SystemMatrix[row * (totalDOF + 1)] += inv_chi;
			RHS[row] = ctx.FlatExtern[row];
		}

		// Update poly->CurrentChi for 6DOF elements
		if(dof == 6)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				poly->CurrentChi = chi;
			}
		}
	}

	// Solve using LAPACK LU (dgesv solves A*x = b in-place)
	auto t_lu_start = std::chrono::high_resolution_clock::now();
#ifdef HAVE_LAPACK
	std::vector<int> ipiv(totalDOF);
	int nrhs = 1;
	int info = 0;

	// dgesv overwrites SystemMatrix with LU factors and RHS with solution
	// No need for extra copy - pass SystemMatrix directly
	dgesv_(&totalDOF, &nrhs, SystemMatrix.data(), &totalDOF, ipiv.data(), RHS.data(), &totalDOF, &info);
	if(info != 0) return -1;  // Singular matrix
#else
	// Fallback: transpose to row-major and use SolveLU_Flat
	for(int i = 0; i < totalDOF; i++)
	{
		for(int j = i + 1; j < totalDOF; j++)
		{
			std::swap(SystemMatrix[i * totalDOF + j], SystemMatrix[j * totalDOF + i]);
		}
	}
	int ierr = SolveLU_Flat(SystemMatrix, RHS, totalDOF);
	if(ierr != 0) return -1;
#endif
	auto t_lu_end = std::chrono::high_resolution_clock::now();
	double t_lu = std::chrono::duration<double>(t_lu_end - t_lu_start).count();
	rad.m_solve_t_lu_decomp += t_lu;
	rad.m_solve_t_linear_solve += t_lu;

	// Copy solution to FlatMagn
	for(int i = 0; i < totalDOF; i++)
	{
		ctx.FlatMagn[i] = RHS[i];
	}

	return 0;  // LU is direct solver, no iterations
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

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

#ifdef HAVE_LAPACK
	// FAST PATH: Use Intel MKL cblas_dgemv for uniform DOF systems
	// This is O(N^2) with highly optimized BLAS instead of manual loops

	if(!IntrctPtr->HasVariableDOF())
	{
		// Pure 3 DOF system (tetrahedra only) - use single BLAS call with alpha=-1.0
		// y = -1.0 * A * x + 0.0 * y
		// Then add diagonal: y[i] += inv_chi[i] * x[i]

		// Intel MKL cblas_dgemv:
		// y := alpha*A*x + beta*y (for CblasNoTrans)
		// A is m x n, x is n, y is m
		// For column-major (CblasColMajor): lda = leading dimension = m = totalDOF
		cblas_dgemv(CblasColMajor, CblasNoTrans,
		            totalDOF, totalDOF,      // m, n (matrix dimensions)
		            -1.0,                     // alpha = -1.0 (negate for 3DOF)
		            FlatInteract, totalDOF,   // A, lda
		            x.data(), 1,              // x, incx
		            0.0,                      // beta
		            y.data(), 1);             // y, incy

		// Add diagonal contribution: y[i] += inv_chi[i] * x[i]
		// This is element-wise multiplication and addition
		#pragma omp parallel for if(totalDOF > 1000)
		for(int i = 0; i < totalDOF; i++)
		{
			y[i] += inv_chi[i] * x[i];
		}
		return;
	}

	// Check if pure 6 DOF system (all hexahedra MSC)
	bool allMSC = true;
	for(int elem = 0; elem < AmOfMainElem && allMSC; elem++)
	{
		if(IntrctPtr->GetElementDOF(elem) != 6) allMSC = false;
	}

	if(allMSC)
	{
		// Pure 6 DOF MSC system - use single BLAS call with alpha=+1.0
		cblas_dgemv(CblasColMajor, CblasNoTrans,
		            totalDOF, totalDOF,
		            1.0,                      // alpha = +1.0 (no negation for MSC)
		            FlatInteract, totalDOF,
		            x.data(), 1,
		            0.0,
		            y.data(), 1);

		// Add diagonal contribution
		#pragma omp parallel for if(totalDOF > 1000)
		for(int i = 0; i < totalDOF; i++)
		{
			y[i] += inv_chi[i] * x[i];
		}
		return;
	}
#endif

	// SLOW PATH: Mixed DOF system (rare) - use block-wise loops
	// This handles the case where 3DOF and 6DOF elements are mixed
	std::fill(y.begin(), y.end(), 0.0);

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

int radTRelaxationMethNo_1::SolveBiCGSTAB_VariableDOF(int totalDOF, double tol, int max_iter, double& residual,
                                                       const std::vector<double>& elemChiArray)
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
			// 3DOF MMM element: use SAME equation as LU solver
			// Equation: (-N + I/chi) * M = H_ext
			// Use isotropic chi from elemChiArray (computed in outer loop)
			double chi = elemChiArray[elem];
			if(chi < 1.0e-6) chi = 1.0e-6;
			double inv_chi_val = 1.0 / chi;

			for(int k = 0; k < 3; k++)
			{
				inv_chi[offset + k] = inv_chi_val;  // Same chi for all 3 components
				// RHS = H_ext (same as LU solver, no Mr term for nonlinear isotropic materials)
				rhs[offset + k] = FlatExtern[offset + k];
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

	// Initial guess: use current magnetization (ELF-compatible)
	// Using the previous solution as initial guess significantly speeds up convergence
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

//-------------------------------------------------------------------------
// BiCGSTAB Solver: SolveLinearStep override
// Uses BiCGSTAB iterative solver with Jacobi preconditioner
//-------------------------------------------------------------------------

int radTRelaxationMethNo_1::SolveLinearStep(NonlinearContext& ctx, int iterCount)
{
	int totalDOF = ctx.totalDOF;
	int AmOfMainElem = ctx.AmOfMainElem;

	// Update poly->CurrentChi for 6DOF elements before BiCGSTAB
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		if(dof == 6)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				poly->CurrentChi = ctx.CurrentChiArray[elem];
			}
		}
	}

	// Call BiCGSTAB solver
	double residual = 0.0;
	const double bicg_tol = rad.m_bicg_tol;
	const int bicg_max_iter = 10000;

	auto t_bicg_start = std::chrono::high_resolution_clock::now();
	int n_iter = SolveBiCGSTAB_VariableDOF(totalDOF, bicg_tol, bicg_max_iter, residual, ctx.CurrentChiArray);
	auto t_bicg_end = std::chrono::high_resolution_clock::now();
	rad.m_solve_t_linear_solve += std::chrono::duration<double>(t_bicg_end - t_bicg_start).count();

	return n_iter;
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

int radTRelaxationMethNo_2::SolveBiCGSTAB_HMatrix_VariableDOF(int totalDOF, double tol, int max_iter, double& residual,
                                                               const std::vector<double>& elemChiArray)
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
	// FIX (2025-12-26): Use isotropic chi from elemChiArray for 3DOF elements (same as BiCGSTAB)
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

		if(dof == 3)
		{
			// 3DOF MMM element: use SAME equation as LU solver
			// Equation: (-N + I/chi) * M = H_ext
			// Use isotropic chi from elemChiArray (computed in outer loop)
			double chi = elemChiArray[elem];
			if(chi < 1.0e-6) chi = 1.0e-6;
			double inv_chi_val = 1.0 / chi;

			for(int k = 0; k < 3; k++)
			{
				inv_chi[offset + k] = inv_chi_val;  // Same chi for all 3 components
				// RHS = H_ext (same as LU solver, no Mr term for nonlinear isotropic materials)
				rhs[offset + k] = FlatExtern[offset + k];
			}
		}
		else if(dof == 6)
		{
			// For 6DOF MSC: use isotropic chi from poly->CurrentChi (same as BiCGSTAB)
			radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
			double chi = (poly && poly->Use6DOF_MSC && poly->CurrentChi > 1.0e-6)
			           ? poly->CurrentChi : 1.0;
			if(chi < 1.0e-6) chi = 1.0e-6;
			double inv_chi_val = 1.0 / chi;

			for(int k = 0; k < 6; k++)
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
	// FIX (2025-12-27): Recompute diagonal EVERY iteration (not cached)
	// Reason: ELF extracts diagonal from updated H-matrix blocks via HACApK_extract_diagonal_omp,
	// which includes the current inv_chi values. Radia was caching the diagonal computed with
	// INITIAL inv_chi, causing slower convergence (17 iter vs 14 iter).
	//
	// The fix is to recompute diag_inv = 1/(N_ii + inv_chi[i]) each iteration.
	// This is O(N) and uses cached N_ii values, so it's cheap.
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

	// OPTIMIZATION (2025-12-26): Removed O(N^2) SetupInteractMatrix() for 3DOF tetrahedra
	// PrecomputeGeometry3DOF() + Compute3x3BlockFast() provides on-demand matrix computation
	// using pre-computed face vertices. This is much faster than SetupInteractMatrix() because:
	// - HACApK ACA+ only needs a subset of matrix elements (not all N^2)
	// - Face geometry is pre-computed once, then reused for all matrix element computations
	// - Thread-local LRU cache avoids redundant computation for nearby elements
	bool is_uniform_3dof = !IntrctPtr->HasVariableDOF();  // uniform 3DOF

	// For 3DOF tetrahedra, initialize FlatExtern from ExternFieldArray
	if(is_uniform_3dof && IntrctPtr->ExternFieldArray != nullptr)
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

		// NOTE: No longer caching Jacobi preconditioner (FIX 2025-12-27)
		// We recompute diag_inv = 1/(N_ii + inv_chi[i]) each iteration
	}

	// NOTE: 3DOF tetrahedra use PrecomputeFlatInteractMatrix() for fast O(1) access
	// 6DOF hexahedra use PrecomputeGeometry() + Compute6x6BlockFast()

	std::vector<double> OldMagn(totalDOF);
	// Store current isotropic chi for ALL elements (unified 3DOF/6DOF handling, same as LU/BiCGSTAB)
	std::vector<double> CurrentChiArray_hacapk(AmOfMainElem, 1.0);
	double MisfitE2 = 1.0e30;
	int totalIterCount = 0;
	int outerIter = 0;

	// Linear material detection: if all materials are linear, converge in 1 iteration
	bool all_materials_linear = true;

	// Initialize H field in NewFieldArray (used for chi(H) computation in nonlinear iteration)
	// Also initialize FlatField for compatibility
	const double H_init_mag = 100.0;  // Same as LU/BiCGSTAB (100 A/m)
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

	// Initialize CurrentChi with ELF-style initial value (same as LU/BiCGSTAB)
	// FIX (2025-12-26): Initialize for BOTH 3DOF and 6DOF elements
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
		radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		// Get initial chi (ELF mucal0 style - from BH curve 2nd point)
		double chi_init = 1.0;
		radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
		if(NonlinMater != nullptr)
		{
			chi_init = NonlinMater->GetInitialChi_ELF_Style();
			if(chi_init <= 0) chi_init = 1.0;
			all_materials_linear = false;  // At least one nonlinear material found
		}
		else
		{
			// Fallback for linear materials: use DefineInstantKsiTensor
			TVector3d H_est(0., 0., H_init_mag);
			TMatrix3d KsiTensor;
			TVector3d MrVect;
			MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);
			chi_init = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
			if(chi_init < 1.0e-6) chi_init = 1.0e-6;
		}
		CurrentChiArray_hacapk[elem] = chi_init;

		// For 6DOF elements, also store in poly->CurrentChi
		if(dof == 6)
		{
			radTPolyhedron* poly = polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				poly->CurrentChi = chi_init;
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

		// Store old B norm for convergence check (same as LU/BiCGSTAB)
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			int offset = IntrctPtr->GetElementDOFOffset(elem);

			if(dof == 3)
			{
				// 3DOF element: store chi and compute B-field norm
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				TVector3d M = g3dRelaxPtr->Magn;
				double chi = CurrentChiArray_hacapk[elem];
				if(chi < 1.0e-6) chi = 1.0e-6;
				TVector3d H(M.x / chi, M.y / chi, M.z / chi);
				TVector3d B(MU_0 * (H.x + M.x), MU_0 * (H.y + M.y), MU_0 * (H.z + M.z));
				OldBnorm[elem] = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
			}
			else if(dof == 6)
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
		int n_iter = SolveBiCGSTAB_HMatrix_VariableDOF(totalDOF, bicg_tol, bicg_max_iter, residual, CurrentChiArray_hacapk);
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

				// Update H field: H = M / chi (same as LU/BiCGSTAB)
				double chi = CurrentChiArray_hacapk[elem];
				if(chi < 1.0e-6) chi = 1.0e-6;
				for(int k = 0; k < 3; k++)
				{
					FlatField[offset + k] = FlatMagn[offset + k] / chi;
				}
				// Also update NewFieldArray for chi(H) computation
				if(IntrctPtr->NewFieldArray != nullptr)
				{
					IntrctPtr->NewFieldArray[elem].x = FlatField[offset];
					IntrctPtr->NewFieldArray[elem].y = FlatField[offset + 1];
					IntrctPtr->NewFieldArray[elem].z = FlatField[offset + 2];
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

		// Compute convergence and update chi (same structure as LU/BiCGSTAB)
		double max_B_rel_change = 0.0;
		bool has_6dof_elements = false;

		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			int offset = IntrctPtr->GetElementDOFOffset(elem);

			if(dof == 3)
			{
				// 3DOF element: update chi using ELF-style dual method (same as LU/BiCGSTAB)
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

				// Get H from updated FlatField (H = M / chi)
				TVector3d H_new(FlatField[offset], FlatField[offset+1], FlatField[offset+2]);
				double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);

				// Chi update using ELF-style dual-method
				double chi_matrix = CurrentChiArray_hacapk[elem];
				double mu_old = chi_matrix + 1.0;

				radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
				double chi_new;
				double relax_hacapk = rad.m_relax;
				if(NonlinMater != nullptr)
				{
					// ELF-style dual-method with relax parameter
					chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old, relax_hacapk);
				}
				else
				{
					// Fallback for linear materials
					TMatrix3d KsiTensor;
					TVector3d MrVect;
					MaterPtr->DefineInstantKsiTensor(H_new, KsiTensor, MrVect);
					chi_new = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
					if(chi_new < 1.0e-6) chi_new = 1.0e-6;
					// Apply under-relaxation for linear materials
					if(relax_hacapk > 0.0 && relax_hacapk <= 1.0)
					{
						chi_new = chi_new * (1.0 - relax_hacapk) + chi_matrix * relax_hacapk;
					}
				}
				CurrentChiArray_hacapk[elem] = chi_new;

				// B-field convergence (ELF mucal2 style) for 3DOF elements
				TVector3d M_new = g3dRelaxPtr->Magn;
				TVector3d B_new_vec(MU_0 * (H_new.x + M_new.x),
				                    MU_0 * (H_new.y + M_new.y),
				                    MU_0 * (H_new.z + M_new.z));
				double B_new_norm = std::sqrt(B_new_vec.x*B_new_vec.x + B_new_vec.y*B_new_vec.y + B_new_vec.z*B_new_vec.z);

				double B_sat = 1.0;
				if(NonlinMater != nullptr)
				{
					B_sat = NonlinMater->GetBsaturation();
					if(B_sat < 1.0e-10) B_sat = 1.0;
				}

				double B_old_norm = OldBnorm[elem];
				double B_rel_change = std::fabs(B_new_norm - B_old_norm) / B_sat;
				if(B_rel_change > max_B_rel_change)
					max_B_rel_change = B_rel_change;
			}
			else if(dof == 6)
			{
				has_6dof_elements = true;
				radTPolyhedron* poly = polyCache[elem];
				if(poly && poly->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
				{
					TVector3d H_new = IntrctPtr->NewFieldArray[elem];
					radTMaterial* MaterPtr = (radTMaterial*)(IntrctPtr->g3dRelaxPtrVect[elem]->MaterHandle.rep);

					// Get chi used for this iteration's matrix
					double chi_matrix = poly->CurrentChi;
					double mu_old = chi_matrix + 1.0;

					// Compute H magnitude for chi update
					double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);

					// Use ELF-style dual-method chi update
					radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
					double chi_new;
					double relax_hacapk = rad.m_relax;
					if(NonlinMater != nullptr)
					{
						// ELF-style dual-method with relax parameter
						chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old, relax_hacapk);
					}
					else
					{
						// Fallback for linear materials
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_new, KsiTensor, MrVect);
						chi_new = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						if(chi_new < 1.0e-6) chi_new = 1.0e-6;
						// Apply under-relaxation for linear materials
						if(relax_hacapk > 0.0 && relax_hacapk <= 1.0)
						{
							chi_new = chi_new * (1.0 - relax_hacapk) + chi_matrix * relax_hacapk;
						}
					}
					poly->CurrentChi = chi_new;

					// B-field convergence
					TVector3d& M_new = poly->Magn;
					double chi_for_B = chi_new;
					if(chi_for_B < 1.0e-6) chi_for_B = 1.0e-6;
					TVector3d H_for_B(M_new.x / chi_for_B, M_new.y / chi_for_B, M_new.z / chi_for_B);
					TVector3d B_new_vec(MU_0 * (H_for_B.x + M_new.x),
					                    MU_0 * (H_for_B.y + M_new.y),
					                    MU_0 * (H_for_B.z + M_new.z));
					double B_new_norm = std::sqrt(B_new_vec.x*B_new_vec.x + B_new_vec.y*B_new_vec.y + B_new_vec.z*B_new_vec.z);

					// Get B_sat from BH curve
					double B_sat = 1.0;
					if(NonlinMater != nullptr)
					{
						B_sat = NonlinMater->GetBsaturation();
						if(B_sat < 1.0e-10) B_sat = 1.0;
					}

					// B-field convergence: |B_new - B_old| / B_sat
					double B_old_norm = OldBnorm[elem];
					double B_rel_change = std::fabs(B_new_norm - B_old_norm) / B_sat;
					if(B_rel_change > max_B_rel_change)
						max_B_rel_change = B_rel_change;
				}
			}
		}

		// Convergence criterion: use B-field change for both 3DOF and 6DOF (same as LU/BiCGSTAB)
		double rel_change = max_B_rel_change;
		MisfitE2 = rel_change * rel_change;

		// Linear materials: converge in exactly 1 iteration (chi is constant, no Newton needed)
		// ELF pattern: linear materials don't need permeability update - solution is exact
		// Check BEFORE normal convergence to ensure single iteration
		if(all_materials_linear)
		{
			outerIter++;
			break;
		}

		if(rel_change <= PrecOnMagnetiz)
		{
			outerIter++;
			break;
		}

		if(radYield.Check() == 0) return outerIter;
	}

	IntrctPtr->RelaxStatusParam.MisfitM = std::sqrt(MisfitE2);

	// Set statistics for GetSolveStats (ELF-compatible)
	rad.m_solve_linear_iterations = rad.m_linear_iterations;  // Copy from HACApK-specific counter
	rad.m_solve_t_linear_solve = rad.m_timing_linear_solve;   // Copy linear solve time
	rad.m_solve_t_matrix_build = rad.m_timing_hmatrix_build;  // H-matrix build time as matrix build

	return outerIter;
}

#endif // RADIA_USE_HACAPK
