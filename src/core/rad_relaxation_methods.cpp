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
#include "rad_constants.h"    // For RadConst::INV_FOUR_PI
#include "rad_bicgstab.h"     // For templated BiCGSTAB

#include <time.h>
#include <chrono>   // For timing instrumentation
#include <cstring>  // For std::memcpy
#include <stdexcept> // For std::runtime_error (loop-star antisym-IMA fail-loud)
#include <cstdio>   // For fprintf in debug logging
#include <cstdlib>  // For getenv
#include <array>    // For std::array in IMA mirror computation
#include <atomic>   // For std::atomic in parallel early-exit patterns
#include <algorithm> // For std::min/std::max
#include <cmath>    // For std::isfinite
#include <map>      // For the sparse coarse-op (E) CSR assembly, Step 3b

// Uncomment to enable chi value debugging
// #define RADIA_DEBUG_CHI

// External access to radTApplication for NonlinearMethod setting
extern radTApplication rad;

//-------------------------------------------------------------------------

#include "rad_parallel.h"

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

	radTRelaxStatusParam& RelStatParR = IntrctPtr->RelaxStatusParam;

	ngcore::ParallelFor(ngcore::IntRange(IntrctPtr->AmOfMainElem), [&](size_t i) {
		double LocalTestBufMaxModM = 0., LocalTestBufMaxModH = 0.;
		if(RelStatParR.MisfitM >= 0. && OldMagnArray != nullptr)
		{
			TVector3d Mnew_mi_MoldVect = NewMagnArray[i] - OldMagnArray[i];
			double local_misfit = Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y
						+ Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;
			ngcore::AtomicAdd(BufMisfitM, local_misfit);
		}
		if(RelStatParR.MaxModM >= 0.)
		{
			LocalTestBufMaxModM = sqrt(NewMagnArray[i].x*NewMagnArray[i].x
								+ NewMagnArray[i].y*NewMagnArray[i].y
								+ NewMagnArray[i].z*NewMagnArray[i].z);
			ngcore::AtomicMax(BufMaxModM, LocalTestBufMaxModM);
		}
		if(RelStatParR.MaxModH >= 0.)
		{
			LocalTestBufMaxModH = sqrt(NewFieldArray[i].x*NewFieldArray[i].x
								+ NewFieldArray[i].y*NewFieldArray[i].y
								+ NewFieldArray[i].z*NewFieldArray[i].z);
			ngcore::AtomicMax(BufMaxModH, LocalTestBufMaxModH);
		}
	});
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
// that handle compact 3-component elements. Mesh-backed magnetic-material solves route through HDiv-VIM.
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
		radTHysteresisMaterial* HystMater = dynamic_cast<radTHysteresisMaterial*>(MaterPtr);
		if(HystMater != nullptr)
		{
			chi_init = HystMater->GetInitialChi_ELF_Style();
			if(chi_init <= 0) chi_init = 1.0;
			ctx.all_materials_linear = false;
			ctx.B_sat = HystMater->GetBsaturation();
			if(ctx.B_sat < 1.0e-10) ctx.B_sat = 1.0;
		}
		else if(NonlinMater != nullptr)
		{
			chi_init = NonlinMater->GetInitialChi_ELF_Style();
			if(chi_init <= 0) chi_init = 1.0;
			double H_ext_mag = 0.0;
			if(dof == 3)
			{
				for(int k = 0; k < 3; k++) H_ext_mag += ctx.FlatExtern[offset + k] * ctx.FlatExtern[offset + k];
				H_ext_mag = std::sqrt(H_ext_mag);
			}
			else if(dof >= 4)
			{
				const TVector3d& Hext = IntrctPtr->ExternFieldArray[elem];
				H_ext_mag = std::sqrt(Hext.x*Hext.x + Hext.y*Hext.y + Hext.z*Hext.z);
			}
			if(H_ext_mag > 1.0e-10)
			{
				double chi_ext = NonlinMater->ComputeChiDualMethod(H_ext_mag, chi_init + 1.0, 0.0);
				if(chi_ext > 1.0e-6 && chi_ext < 1.0e10) chi_init = chi_ext;
			}
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

		// Store in poly->CurrentChi for face-coefficient elements.
		if(dof >= 4)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->UseFaceChargeDOF)
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
		else if(dof >= 4)
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

bool BuildBaseMatrix(NonlinearContext& ctx, radTInteraction* IntrctPtr)
{
	double* FlatInteract = IntrctPtr->GetFlatInteractMatrix();
	if(FlatInteract == nullptr) return false;

	// CRITICAL: Use size_t to avoid int32 overflow for DOF > 46340
	size_t matrix_size = (size_t)ctx.totalDOF * (size_t)ctx.totalDOF;
	try {
		ctx.BaseMatrix.resize(matrix_size);
	} catch (const std::bad_alloc&) {
		double required_gb = (double)matrix_size * 8 / (1024.0 * 1024.0 * 1024.0);
		fprintf(stderr, "Radia::Solve> Matrix requires %.1f GB memory for DOF=%d. Problem too large.\n", required_gb, ctx.totalDOF);
		return false;
	}

	// Copy interaction matrix and NEGATE all entries
	// Physical equation: A = -K/(4pi) + diag(1/chi) for all element types
	// FlatInteract stores the physical interaction kernel.
	// Negate it to match the relaxation equation convention.
	std::memcpy(ctx.BaseMatrix.data(), FlatInteract, matrix_size * sizeof(double));

	// Negate entire matrix.
	for(size_t i = 0; i < matrix_size; i++)
	{
		ctx.BaseMatrix[i] = -ctx.BaseMatrix[i];
	}

	// NOTE: Do NOT symmetrize the K matrix. The asymmetric K[i,j] != K[j,i] is
	// physically correct for non-orthogonal hexahedra where face normals and
	// eval points differ between elements. BiCGSTAB uses FlatInteract directly
	// and gets correct results; LU should use the same unsymmetrized matrix.

	return true;
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
		ctx.OldChi[elem] = ctx.CurrentChiArray[elem];

		if(dof == 3)
		{
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			TVector3d M = g3dRelaxPtr->Magn;
			TVector3d H(ctx.FlatField[offset], ctx.FlatField[offset+1], ctx.FlatField[offset+2]);
			TVector3d B(MU_0 * (H.x + M.x), MU_0 * (H.y + M.y), MU_0 * (H.z + M.z));
			ctx.OldBnorm[elem] = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
		}
		else if(dof >= 4)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->UseFaceChargeDOF)
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
// Helper function to compute actual H field at element centers from sigma values
// This is the ELF-compatible approach: H = H_ext + H_demag (from all sigma contributions)
//
// FIX (2026-02-05): For IMA (Image Method Analysis), add demagnetizing field from
// mirror elements. Without this, nonlinear IMA produces ~22% error because chi
// update uses H = M/chi which doesn't include mirror demagnetizing contributions.
//-------------------------------------------------------------------------

void ComputeActualHFieldFromSigma(NonlinearContext& ctx, radTInteraction* IntrctPtr)
{
	// First pass: store magnetization/sigma values for all elements
	for(int elem = 0; elem < ctx.AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		if(dof == 3)
		{
			// Compact 3-component elements: FlatMagn contains M directly.
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			g3dRelaxPtr->Magn.x = ctx.FlatMagn[offset + 0];
			g3dRelaxPtr->Magn.y = ctx.FlatMagn[offset + 1];
			g3dRelaxPtr->Magn.z = ctx.FlatMagn[offset + 2];
		}
		else if(dof >= 4)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->UseFaceChargeDOF)
			{
				// Store sigma values in polyhedron (5 for wedge, 6 for hex)
				for(int k = 0; k < dof; k++)
				{
					poly->Sigma[k] = ctx.FlatMagn[offset + k];
				}

				// Compute effective magnetization from sigma using weighted least-squares
				// M_x = sum(sigma_i * n_x_i) / sum(n_x_i^2) (ELF-compatible)
				double Mx = 0.0, My = 0.0, Mz = 0.0;
				double wx = 0.0, wy = 0.0, wz = 0.0;
				for(int face = 0; face < dof; face++)
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
			}
		}
	}

	// Second pass: compute H field for chi update
	// For 6DOF hexahedral elements: use H = M/chi (constitutive relation)
	// This is the same approach used by HACApK and matches ELF behavior.
	//
	// NOTE (FIX 2026-02-02): The previous approach computed H = H_ext + H_demag via B_comp,
	// but this caused ~2% error vs ELF because of numerical issues in self-field computation.
	// HACApK uses H = M/chi and matches ELF perfectly, so we adopt the same approach here.
	//
	// KNOWN ISSUE (2026-02-05): For IMA (Image Method Analysis), the H = M/chi approach
	// does NOT properly account for nonlinear material behavior with mirror elements.
	// This causes ~22% error for nonlinear materials with IMA.
	//
	// Investigation notes:
	// - Adding mirror demagnetizing field to H made results WORSE (74% error vs 22%)
	// - The issue may be in the nonlinear iteration formulation, not just chi update
	// - For now, nonlinear IMA is NOT supported
	//
	// Workaround: Use full model (no IMA) for nonlinear problems.

	for(int elem_j = 0; elem_j < ctx.AmOfMainElem; elem_j++)
	{
		int dof_j = IntrctPtr->GetElementDOF(elem_j);
		int offset_j = IntrctPtr->GetElementDOFOffset(elem_j);

		if(dof_j == 3)
		{
			// 3DOF elements: compute H = M / chi (constitutive relation)
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem_j];
			double chi = ctx.CurrentChiArray[elem_j];
			if(chi < 1.0e-6) chi = 1.0e-6;

			if(IntrctPtr->NewFieldArray != nullptr)
			{
				IntrctPtr->NewFieldArray[elem_j].x = g3dRelaxPtr->Magn.x / chi;
				IntrctPtr->NewFieldArray[elem_j].y = g3dRelaxPtr->Magn.y / chi;
				IntrctPtr->NewFieldArray[elem_j].z = g3dRelaxPtr->Magn.z / chi;
			}

			// Also update FlatField for convergence check
			ctx.FlatField[offset_j + 0] = g3dRelaxPtr->Magn.x / chi;
			ctx.FlatField[offset_j + 1] = g3dRelaxPtr->Magn.y / chi;
			ctx.FlatField[offset_j + 2] = g3dRelaxPtr->Magn.z / chi;
		}
		else if(dof_j >= 4)
		{
			radTPolyhedron* poly_j = ctx.polyCache[elem_j];
			if(poly_j && poly_j->UseFaceChargeDOF && IntrctPtr->NewFieldArray != nullptr)
			{
				// Compute H = M / chi (constitutive relation, same as HACApK)
				double chi = ctx.CurrentChiArray[elem_j];
				if(chi < 1.0e-6) chi = 1.0e-6;

				IntrctPtr->NewFieldArray[elem_j].x = poly_j->Magn.x / chi;
				IntrctPtr->NewFieldArray[elem_j].y = poly_j->Magn.y / chi;
				IntrctPtr->NewFieldArray[elem_j].z = poly_j->Magn.z / chi;

			}
		}
	}

}

//-------------------------------------------------------------------------

void UpdateMagnAndComputeH(NonlinearContext& ctx, radTInteraction* IntrctPtr)
{
	// For 3DOF elements: use H = M / chi (constitutive relation)
	// For 6DOF elements: compute actual H from sigma (ELF-compatible)

	// First, update magnetization for all elements
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

			// H = M / chi (for 3DOF elements)
			double chi = ctx.CurrentChiArray[elem];
			if(chi < 1.0e-6) chi = 1.0e-6;
			for(int k = 0; k < 3; k++)
			{
				ctx.FlatField[offset + k] = ctx.FlatMagn[offset + k] / chi;
			}
		}
		else if(dof >= 4)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->UseFaceChargeDOF)
			{
				// Store sigma values (5 for wedge, 6 for hex)
				for(int k = 0; k < dof; k++)
				{
					poly->Sigma[k] = ctx.FlatMagn[offset + k];
				}

				// Compute effective magnetization from sigma (ELF-compatible weighted least-squares)
				double Mx = 0.0, My = 0.0, Mz = 0.0;
				double wx = 0.0, wy = 0.0, wz = 0.0;
				for(int face = 0; face < dof; face++)
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
			}
		}
	}

	// For 6DOF elements: compute H = M / chi (constitutive relation)
	// This approach is stable and gives convergent solutions, though with ~18% error vs ELF.
	// The error is likely due to different matrix formulations between Radia and ELF.
	//
	// Note: Using H = |sum(sigma*n)| (ELF mucal1 style) was tested but gave 85% error,
	// suggesting that Radia's matrix formulation requires the H = M/chi approach.
	for(int elem = 0; elem < ctx.AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);

		if(dof >= 4)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->UseFaceChargeDOF && IntrctPtr->NewFieldArray != nullptr)
			{
				// Compute H = M / chi (constitutive relation)
				double chi = ctx.CurrentChiArray[elem];
				if(chi < 1.0e-6) chi = 1.0e-6;

				TVector3d H_from_M;
				H_from_M.x = poly->Magn.x / chi;
				H_from_M.y = poly->Magn.y / chi;
				H_from_M.z = poly->Magn.z / chi;

				IntrctPtr->NewFieldArray[elem] = H_from_M;
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
		radTHysteresisMaterial* HystMater = dynamic_cast<radTHysteresisMaterial*>(MaterPtr);

		double chi_matrix = ctx.CurrentChiArray[elem];
		double mu_old = chi_matrix + 1.0;
		double chi_new;

		if(dof == 3)
		{
			TVector3d H_new(ctx.FlatField[offset], ctx.FlatField[offset+1], ctx.FlatField[offset+2]);
			double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);

			if(HystMater != nullptr)
			{
				// Energy-based hysteresis: pass full H vector for correct direction
				chi_new = HystMater->ComputeChiFromH(H_new);
				if(ctx.use_newton && !ctx.DifferentialChiArray.empty())
				{
					ctx.DifferentialChiArray[elem] = HystMater->ComputeDifferentialChi(H_mag);
				}
			}
			else if(NonlinMater != nullptr)
			{
				chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old, ctx.relax_param);
				if(ctx.use_newton && !ctx.DifferentialChiArray.empty())
				{
					ctx.DifferentialChiArray[elem] = NonlinMater->ComputeDifferentialChi(H_mag);
				}
			}
			else
			{
				chi_new = chi_matrix;  // Linear: keep constant
				if(ctx.use_newton && !ctx.DifferentialChiArray.empty())
				{
					ctx.DifferentialChiArray[elem] = chi_matrix;
				}
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
		else if(dof >= 4)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->UseFaceChargeDOF && IntrctPtr->NewFieldArray != nullptr)
			{
				TVector3d H_new = IntrctPtr->NewFieldArray[elem];
				double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);
				TVector3d M_poly = poly->Magn;
				double M_mag = std::sqrt(M_poly.x*M_poly.x + M_poly.y*M_poly.y + M_poly.z*M_poly.z);


				if(HystMater != nullptr)
				{
					chi_new = HystMater->ComputeChiFromH(H_new);
					if(ctx.use_newton && !ctx.DifferentialChiArray.empty())
					{
						ctx.DifferentialChiArray[elem] = HystMater->ComputeDifferentialChi(H_mag);
					}
				}
				else if(NonlinMater != nullptr)
				{
					chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old, ctx.relax_param);
					if(ctx.use_newton && !ctx.DifferentialChiArray.empty())
					{
						ctx.DifferentialChiArray[elem] = NonlinMater->ComputeDifferentialChi(H_mag);
					}
				}
				else
				{
					chi_new = chi_matrix;
					if(ctx.use_newton && !ctx.DifferentialChiArray.empty())
					{
						ctx.DifferentialChiArray[elem] = chi_matrix;
					}
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

//-------------------------------------------------------------------------

double ApplyLineSearchDamping(NonlinearContext& ctx, radTInteraction* IntrctPtr,
                              const std::vector<double>& sigma_trial)
{
	// Only apply if Newton is active and damping is enabled
	if(!ctx.use_newton || !ctx.newton_damping_enabled) {
		return 1.0;  // Full step (no damping)
	}

	const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;

	// Store original solution for backtracking
	std::vector<double> sigma_old = ctx.OldSigma;

	double omega = 1.0;  // Start with full Newton step
	double best_residual = 1.0e30;

	for(int ls_iter = 0; ls_iter < ctx.newton_ls_max_iter; ls_iter++)
	{
		// Apply damped update: sigma = omega * trial + (1-omega) * old
		for(int i = 0; i < ctx.totalDOF; i++)
		{
			ctx.FlatMagn[i] = omega * sigma_trial[i] + (1.0 - omega) * sigma_old[i];
		}

		// Update magnetization from FlatMagn and compute H field
		// This syncs poly->Magn and computes IntrctPtr->NewFieldArray
		ComputeActualHFieldFromSigma(ctx, IntrctPtr);

		// Compute B-field change metric (same as convergence check)
		// Returns max_elem |B_new - B_old| / B_sat
		double residual = 0.0;

		for(int elem = 0; elem < ctx.AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			int offset = IntrctPtr->GetElementDOFOffset(elem);

			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
			radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);

			// Get current chi (not updated yet - use old chi for line search)
			double chi = ctx.CurrentChiArray[elem];
			if(chi < 1.0e-6) chi = 1.0e-6;

			double B_sat = ctx.B_sat;
			if(NonlinMater != nullptr) {
				B_sat = NonlinMater->GetBsaturation();
				if(B_sat < 1.0e-10) B_sat = 1.0;
			}

			if(dof == 3)
			{
				// 3DOF: compute B from current FlatMagn.
				// FIX (line-search metric): use the FRESH secant chi(H) -- the same chi the
				// outer convergence test recomputes -- NOT the stale ctx.CurrentChiArray.
				// With the stale chi the damping minimized a different quantity than the
				// convergence metric and green-lit Newton steps that convergence then judged
				// WORSE -> Newton diverged after activation.
				TVector3d M(ctx.FlatMagn[offset], ctx.FlatMagn[offset+1], ctx.FlatMagn[offset+2]);
				double chi_b = chi;
				if(NonlinMater != nullptr) {
					double H_mag_ls = std::sqrt(M.x*M.x + M.y*M.y + M.z*M.z) / chi;
					chi_b = NonlinMater->ComputeChiDualMethod(H_mag_ls, chi + 1.0, rad.m_relax);
					if(chi_b < 1.0e-6) chi_b = 1.0e-6;
				}
				TVector3d H(M.x / chi_b, M.y / chi_b, M.z / chi_b);
				TVector3d B(MU_0 * (H.x + M.x), MU_0 * (H.y + M.y), MU_0 * (H.z + M.z));
				double B_new_norm = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);

				double B_rel_change = std::fabs(B_new_norm - ctx.OldBnorm[elem]) / B_sat;
				if(B_rel_change > residual)
					residual = B_rel_change;
			}
			else if(dof >= 4)
			{
				// 6DOF: compute B from poly->Magn (already updated by ComputeActualHFieldFromSigma)
				radTPolyhedron* poly = ctx.polyCache[elem];
				if(poly && poly->UseFaceChargeDOF && IntrctPtr->NewFieldArray != nullptr)
				{
					TVector3d M = poly->Magn;
					// FIX (line-search metric): fresh secant chi(H), mirroring convergence (see 3DOF note).
					double chi_b = chi;
					if(NonlinMater != nullptr) {
						double H_mag_ls = std::sqrt(M.x*M.x + M.y*M.y + M.z*M.z) / chi;
						chi_b = NonlinMater->ComputeChiDualMethod(H_mag_ls, chi + 1.0, rad.m_relax);
						if(chi_b < 1.0e-6) chi_b = 1.0e-6;
					}
					TVector3d H(M.x / chi_b, M.y / chi_b, M.z / chi_b);
					TVector3d B(MU_0 * (H.x + M.x), MU_0 * (H.y + M.y), MU_0 * (H.z + M.z));
					double B_new_norm = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);

					double B_rel_change = std::fabs(B_new_norm - ctx.OldBnorm[elem]) / B_sat;
					if(B_rel_change > residual)
						residual = B_rel_change;
				}
			}
		}

		// Accept step if residual decreased (or first iteration)
		if(ls_iter == 0)
		{
			best_residual = residual;
			// Check if full step is acceptable (at least 1% improvement from previous iteration)
			// Note: ctx.max_B_rel_change contains residual from PREVIOUS iteration
			if(ctx.max_B_rel_change > 1.0e-12 && residual < ctx.max_B_rel_change * 0.99)
			{
				ctx.accepted_omegas.push_back(omega);
				return omega;  // Accept full Newton step
			}
			// If first nonlinear iteration, always accept full step
			if(ctx.max_B_rel_change < 1.0e-12)
			{
				ctx.accepted_omegas.push_back(omega);
				return omega;
			}
		}
		else if(ctx.max_B_rel_change > 1.0e-12 && residual < ctx.max_B_rel_change * 0.99)
		{
			// FIX (acceptance reference): accept a damped step only when it beats the
			// PREVIOUS ITERATE (ctx.max_B_rel_change), not merely the undamped full Newton
			// step (best_residual).  The old reference accepted steps that beat the bad full
			// step yet still WORSENED the residual vs where we were -> monotonic divergence.
			ctx.accepted_omegas.push_back(omega);
			return omega;
		}

		// Reject step: reduce omega and retry
		omega *= 0.5;
		ctx.total_ls_backtracks++;

		// Stop if omega too small
		if(omega < ctx.newton_ls_min_omega)
		{
			omega = ctx.newton_ls_min_omega;
			ctx.accepted_omegas.push_back(omega);
			// Accept minimal step to avoid stall
			for(int i = 0; i < ctx.totalDOF; i++)
				ctx.FlatMagn[i] = omega * sigma_trial[i] + (1.0 - omega) * sigma_old[i];
			ComputeActualHFieldFromSigma(ctx, IntrctPtr);
			return omega;
		}
	}

	// Max iterations reached: accept last omega
	ctx.accepted_omegas.push_back(omega);
	return omega;
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
	ctx.nonlinear_tol = PrecOnMagnetiz;

	// Newton-Raphson initialization
	ctx.use_newton = rad.m_use_newton;
	if(ctx.use_newton)
	{
		ctx.DifferentialChiArray.resize(ctx.AmOfMainElem, 1.0);
		ctx.OldSigma.resize(ctx.totalDOF, 0.0);

		// Line search damping parameters
		ctx.newton_damping_enabled = rad.m_newton_damping_enabled;
		ctx.newton_ls_max_iter = rad.m_newton_ls_max_iter;
		ctx.newton_ls_min_omega = rad.m_newton_ls_min_omega;
		ctx.total_ls_backtracks = 0;
		ctx.accepted_omegas.reserve(MaxIterNumber);
	}

	// Build base matrix (geometric part without chi)
	if(NeedsDenseMatrix())
	{
		if(!BuildBaseMatrix(ctx, IntrctPtr))
			return 0;  // Memory allocation failed
	}

	int iterCount = 0;
	double MisfitE2 = 1.0e30;

	// Nonlinear iteration loop
	for(iterCount = 0; iterCount < MaxIterNumber; iterCount++)
	{
		// Store old sigma for Newton RHS correction
		if(ctx.use_newton)
		{
			for(int i = 0; i < ctx.totalDOF; i++)
				ctx.OldSigma[i] = ctx.FlatMagn[i];
		}

		// Store old values for convergence check
		StoreOldValuesAndComputeBnorm(ctx, IntrctPtr);
		// Solve linear system (virtual - overridden by LU, BiCGSTAB, HACApK)
		int linearIter = SolveLinearStep(ctx, iterCount);
		if(linearIter < 0)
		{
			IntrctPtr->RelaxStatusParam.MisfitM = 1.0e30;
			return MaxIterNumber;
		}

		// Update element magnetization and compute actual H field from the RAW solver iterate (FlatMagn)
		ComputeActualHFieldFromSigma(ctx, IntrctPtr);

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
// B-input Newton-Raphson Solver for Energy-Based Hysteresis
//
// Solves: F(M) = M - Inverse(mu_0*(H_ext + N*M + M)) / mu_0 = 0
// Jacobian: dF/dM = I - (dJ/dB) * (N + I)
// where dJ/dB is block-diagonal (3x3 per element), computed analytically
// from ComputeJacobian(dBdH, chi_d) -> dJ/dB = I - mu_0 * inv(dB/dH)
//
// Converges in 2-4 Newton iterations (vs hundreds for Hantila/Picard).
// Requires all elements to have radTEnergyHysteresisMaterial.
//=========================================================================

int radTIterativeRelaxMeth::AutoRelax_BInput_Newton(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	static const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;
	static const double NU_0 = 1.0 / MU_0;

	// Initialize context
	// For B-input Newton, NEVER reset magnetization.
	// Element M and hysteresis material state (m_Jk_prev, m_Jk_current) must
	// persist between Solve() calls for proper hysteresis stepping.
	// Virgin state has M=0 naturally (no explicit reset needed).
	NonlinearContext ctx;
	if(!InitializeNonlinearContext(ctx, IntrctPtr, /*MagnResetIsNotNeeded=*/1))
		return 0;

	ctx.use_b_input = true;

	int n_elem = ctx.AmOfMainElem;
	int dof = ctx.totalDOF;

	// Verify all elements are hysteresis materials (energy or play)
	std::vector<radTHysteresisMaterial*> hys_mats(n_elem, nullptr);
	for(int i = 0; i < n_elem; i++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
		radTMaterial* mat = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
		hys_mats[i] = dynamic_cast<radTHysteresisMaterial*>(mat);
		if(hys_mats[i] == nullptr) return 0;  // Not all hysteresis - abort
	}

	// Build base matrix (geometric interaction matrix N)
	if(NeedsDenseMatrix())
	{
		if(!BuildBaseMatrix(ctx, IntrctPtr))
			return 0;
	}

	// Build NpI = N + I (column-major, dof x dof)
	// BaseMatrix stores the geometric part as -N.
	// So NpI = -BaseMatrix + I  (to get N + I)
	// N = -BaseMatrix, so NpI = -BaseMatrix + I.
	std::vector<double> NpI(static_cast<size_t>(dof) * dof, 0.0);
	for(size_t j = 0; j < (size_t)dof; j++)
	{
		for(size_t i = 0; i < (size_t)dof; i++)
		{
			// BaseMatrix is column-major: A(i,j) at index [j*dof + i]
			// BaseMatrix = -N, so N = -BaseMatrix.
			double N_ij = -ctx.BaseMatrix[j * dof + i];
			NpI[j * dof + i] = N_ij + (i == j ? 1.0 : 0.0);
		}
	}

	// Save hysteresis states for all elements (reference point for play operators)
	ctx.saved_hys_states.resize(n_elem);
	for(int i = 0; i < n_elem; i++)
		hys_mats[i]->SaveState(ctx.saved_hys_states[i]);

	// Initialize M: use current element magnetization as initial guess
	// This is much better than Forward(H_ext_only) for hysteresis stepping,
	// because M changes gradually between field steps.
	// On first solve (virgin state), element M is already zero (correct start).
	{
		bool all_zero = true;
		for(int i = 0; i < n_elem; i++)
		{
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
			int offset = IntrctPtr->GetElementDOFOffset(i);
			ctx.FlatMagn[offset+0] = g3dRelaxPtr->Magn.x;
			ctx.FlatMagn[offset+1] = g3dRelaxPtr->Magn.y;
			ctx.FlatMagn[offset+2] = g3dRelaxPtr->Magn.z;
			if(g3dRelaxPtr->Magn.x != 0 || g3dRelaxPtr->Magn.y != 0 || g3dRelaxPtr->Magn.z != 0)
				all_zero = false;
		}

		// If all M = 0 (virgin state), initialize from Forward(H_ext)
		if(all_zero)
		{
			for(int i = 0; i < n_elem; i++)
			{
				int offset = IntrctPtr->GetElementDOFOffset(i);
				TVector3d H_ext_i(ctx.FlatExtern[offset], ctx.FlatExtern[offset+1], ctx.FlatExtern[offset+2]);

				hys_mats[i]->RestoreState(ctx.saved_hys_states[i]);
				TVector3d B = hys_mats[i]->Inverse(H_ext_i);  // H -> B
				TVector3d M_init = (1.0 / MU_0) * B - H_ext_i;

				ctx.FlatMagn[offset+0] = M_init.x;
				ctx.FlatMagn[offset+1] = M_init.y;
				ctx.FlatMagn[offset+2] = M_init.z;
			}
		}
	}

	// Working arrays
	std::vector<double> M_vec(dof);
	std::vector<double> H_vec(dof);
	std::vector<double> B_vec(dof);
	std::vector<double> M_model(dof);
	std::vector<double> F_vec(dof);
	std::vector<double> dJdB_blocks(n_elem * 9);  // n_elem 3x3 blocks, column-major
	std::vector<double> dM(dof);

	// Copy initial M
	for(int i = 0; i < dof; i++)
		M_vec[i] = ctx.FlatMagn[i];


	int iterCount = 0;
	double final_residual = 1.0;

	// Interaction matrix N (column-major): stored as -BaseMatrix.
	// For MatVec: H = H_ext + N*M, where N = -BaseMatrix
	// So H = H_ext - BaseMatrix*M

	for(iterCount = 0; iterCount < MaxIterNumber; iterCount++)
	{
		// Restore all hysteresis states to beginning-of-step reference
		ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
			hys_mats[i]->RestoreState(ctx.saved_hys_states[i]);
		});

		// Compute H = H_ext + N*M
		// N*M = -BaseMatrix*M.
#ifdef HAVE_LAPACK
		{
			// H = H_ext (copy)
			cblas_dcopy(dof, ctx.FlatExtern, 1, H_vec.data(), 1);
			// H += (-1) * BaseMatrix * M  (column-major MatVec)
			cblas_dgemv(CblasColMajor, CblasNoTrans, dof, dof,
			            -1.0, ctx.BaseMatrix.data(), dof,
			            M_vec.data(), 1,
			            1.0, H_vec.data(), 1);
		}
#else
		{
			for(int i = 0; i < dof; i++) H_vec[i] = ctx.FlatExtern[i];
			for(int j = 0; j < dof; j++)
			{
				double Mj = M_vec[j];
				for(int i = 0; i < dof; i++)
					H_vec[i] -= ctx.BaseMatrix[(size_t)j * dof + i] * Mj;
			}
		}
#endif

		// Compute B = mu_0 * (H + M)
		for(int i = 0; i < dof; i++)
			B_vec[i] = MU_0 * (H_vec[i] + M_vec[i]);

		// Evaluate Forward(B) per element -> M_model, and compute dJ/dB analytically
		ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
			int offset = IntrctPtr->GetElementDOFOffset((int)i);
			TVector3d B_i(B_vec[offset], B_vec[offset+1], B_vec[offset+2]);

			// Restore to start-of-step reference state before each Forward.
			// Forward() overwrites m_Jk_prev, so without this the pinning
			// reference drifts across Newton iterations.
			hys_mats[i]->RestoreState(ctx.saved_hys_states[i]);

			// Forward(B) -> H (natural for B-input Play), sets m_Jk_current internally
			TVector3d H_inv = hys_mats[i]->Forward(B_i);

			// M_model = J_total / mu_0 = (B - mu_0*H) / mu_0
			TVector3d M_i = (1.0 / MU_0) * B_i - H_inv;
			M_model[offset+0] = M_i.x;
			M_model[offset+1] = M_i.y;
			M_model[offset+2] = M_i.z;

			// Analytical Jacobian: ComputeJacobian gives dB/dH (3x3)
			// dJ/dB = I - mu_0 * inv(dB/dH)
			TMatrix3d dBdH;
			double chi_d;
			hys_mats[i]->ComputeJacobian(dBdH, chi_d);

			TMatrix3d dBdH_inv = Matrix3d_inv(dBdH);
			// dJ/dB = I - mu_0 * dBdH_inv
			// Store as column-major 3x3 block
			size_t blk = i * 9;
			// Column 0: dJ/dB[:, 0]
			dJdB_blocks[blk + 0] = 1.0 - MU_0 * dBdH_inv.Str0.x;
			dJdB_blocks[blk + 1] =      - MU_0 * dBdH_inv.Str1.x;
			dJdB_blocks[blk + 2] =      - MU_0 * dBdH_inv.Str2.x;
			// Column 1: dJ/dB[:, 1]
			dJdB_blocks[blk + 3] =      - MU_0 * dBdH_inv.Str0.y;
			dJdB_blocks[blk + 4] = 1.0 - MU_0 * dBdH_inv.Str1.y;
			dJdB_blocks[blk + 5] =      - MU_0 * dBdH_inv.Str2.y;
			// Column 2: dJ/dB[:, 2]
			dJdB_blocks[blk + 6] =      - MU_0 * dBdH_inv.Str0.z;
			dJdB_blocks[blk + 7] =      - MU_0 * dBdH_inv.Str1.z;
			dJdB_blocks[blk + 8] = 1.0 - MU_0 * dBdH_inv.Str2.z;
		});

		// Compute residual F = M - M_model
		double F_norm = 0.0;
		double M_norm = 0.0;
		for(int i = 0; i < dof; i++)
		{
			F_vec[i] = M_vec[i] - M_model[i];
			F_norm += F_vec[i] * F_vec[i];
			M_norm += M_model[i] * M_model[i];
		}
		F_norm = std::sqrt(F_norm);
		M_norm = std::sqrt(M_norm);
		if(M_norm < 1.0) M_norm = 1.0;

		double rel_residual = F_norm / M_norm;
		final_residual = rel_residual;


		// Check convergence
		if(rel_residual < PrecOnMagnetiz)
		{
			// Use model M for consistency
			for(int i = 0; i < dof; i++)
				M_vec[i] = M_model[i];
			iterCount++;
			break;
		}

		// Solve J_F * dM = -F where J_F = I - dJ/dB * (N+I)
		int solve_err = SolveBInputLinearStep(ctx, NpI, dJdB_blocks, F_vec, dM);
		if(solve_err != 0) break;

		// Line search: find step that reduces ||F||
		double step = 1.0;
		bool ls_success = false;
		for(int ls = 0; ls < 10; ls++)
		{
			// M_trial = M + step * dM
			std::vector<double> M_trial(dof);
			for(int i = 0; i < dof; i++)
				M_trial[i] = M_vec[i] + step * dM[i];

			// Compute H_trial = H_ext + N*M_trial
			std::vector<double> H_trial(dof);
#ifdef HAVE_LAPACK
			cblas_dcopy(dof, ctx.FlatExtern, 1, H_trial.data(), 1);
			cblas_dgemv(CblasColMajor, CblasNoTrans, dof, dof,
			            -1.0, ctx.BaseMatrix.data(), dof,
			            M_trial.data(), 1,
			            1.0, H_trial.data(), 1);
#else
			for(int i = 0; i < dof; i++) H_trial[i] = ctx.FlatExtern[i];
			for(int j = 0; j < dof; j++)
			{
				double Mj = M_trial[j];
				for(int i = 0; i < dof; i++)
					H_trial[i] -= ctx.BaseMatrix[(size_t)j * dof + i] * Mj;
			}
#endif

			// B_trial = mu_0*(H_trial + M_trial)
			std::vector<double> B_trial(dof);
			for(int i = 0; i < dof; i++)
				B_trial[i] = MU_0 * (H_trial[i] + M_trial[i]);

			// Evaluate F_trial: restore states, call Forward (B->H)
			double F_trial_norm2 = 0.0;
			for(int i = 0; i < n_elem; i++)
			{
				hys_mats[i]->RestoreState(ctx.saved_hys_states[i]);
				int offset = IntrctPtr->GetElementDOFOffset(i);
				TVector3d B_i(B_trial[offset], B_trial[offset+1], B_trial[offset+2]);
				TVector3d H_inv = hys_mats[i]->Forward(B_i);
				TVector3d M_mod = (1.0 / MU_0) * B_i - H_inv;
				double fx = M_trial[offset+0] - M_mod.x;
				double fy = M_trial[offset+1] - M_mod.y;
				double fz = M_trial[offset+2] - M_mod.z;
				F_trial_norm2 += fx*fx + fy*fy + fz*fz;
			}

			if(std::sqrt(F_trial_norm2) < F_norm)
			{
				ls_success = true;
				break;
			}
			step *= 0.5;
		}


		// Update M
		for(int i = 0; i < dof; i++)
			M_vec[i] += step * dM[i];

		// Check for user abort
		if(radYield.Check() == 0)
			break;
	}

	// Store converged M back to elements
	for(int i = 0; i < dof; i++)
		ctx.FlatMagn[i] = M_vec[i];

	// Final Inverse + CommitState: set definitive play operator state for each element
	// The last Newton iteration did RestoreState + Inverse at converged B,
	// so m_Jk_current is correct. But we must ensure RestoreState + Inverse
	// is done with the final converged B (not line search trial B).
	ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
		int offset = IntrctPtr->GetElementDOFOffset((int)i);
		TVector3d M_final(M_vec[offset], M_vec[offset+1], M_vec[offset+2]);

		// Update element magnetization
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[(int)i];
		g3dRelaxPtr->SetM(M_final);

		// Restore state to beginning-of-step, then Forward(B_converged)
		// to set m_Jk_current to the definitive values
		hys_mats[i]->RestoreState(ctx.saved_hys_states[i]);
		TVector3d B_final = MU_0 * TVector3d(
			H_vec[offset] + M_vec[offset],
			H_vec[offset+1] + M_vec[offset+1],
			H_vec[offset+2] + M_vec[offset+2]);
		hys_mats[i]->Forward(B_final);  // B -> H (natural)

		// CommitState: m_Jk_prev = m_Jk_current
		// This ensures the next Solve() step starts from the correct reference
		hys_mats[i]->CommitState();
	});

	// Recompute final H for all elements
	{
#ifdef HAVE_LAPACK
		cblas_dcopy(dof, ctx.FlatExtern, 1, H_vec.data(), 1);
		cblas_dgemv(CblasColMajor, CblasNoTrans, dof, dof,
		            -1.0, ctx.BaseMatrix.data(), dof,
		            M_vec.data(), 1,
		            1.0, H_vec.data(), 1);
#else
		for(int i = 0; i < dof; i++) H_vec[i] = ctx.FlatExtern[i];
		for(int j = 0; j < dof; j++)
		{
			double Mj = M_vec[j];
			for(int i = 0; i < dof; i++)
				H_vec[i] -= ctx.BaseMatrix[(size_t)j * dof + i] * Mj;
		}
#endif
		ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
			int offset = IntrctPtr->GetElementDOFOffset((int)i);
			ctx.FlatField[offset+0] = H_vec[offset+0];
			ctx.FlatField[offset+1] = H_vec[offset+1];
			ctx.FlatField[offset+2] = H_vec[offset+2];
		});
	}

	// Update RelaxStatusParam
	IntrctPtr->RelaxStatusParam.MisfitM = final_residual;

	return iterCount;
}

//=========================================================================
// B-input Hantila Solver for Energy-Based Hysteresis
//
// Uses constant LHS (I - alpha*N), LU-factored ONCE.
// Each iteration: O(N^2) back-substitution + Forward(B) per element.
// Hantila splits M = alpha*H + R, so LHS = (I - alpha*N) is constant.
//
// Algorithm:
//   LHS = I + alpha * BaseMatrix  (where BaseMatrix = -N)
//   LU factor LHS once
//   Each iteration:
//     1. H = LU_solve(H_ext - BaseMatrix * R)   [= (I-alpha*N)^{-1} * (H_ext + N*R)]
//     2. M = R + alpha * H
//     3. B = mu_0 * (H + M)
//     4. Forward(B_i) -> Jk, M_new_i = J/mu_0
//     5. R_new = M_new - alpha * H
//     6. Convergence: max|dB|/B_sat < tol
//=========================================================================

int radTIterativeRelaxMeth::AutoRelax_BInput_Hantila(
	double PrecOnMagnetiz, int MaxIterNumber,
	double alpha, double relax, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	static const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;

	// Initialize context (preserve element M for hysteresis state continuity)
	NonlinearContext ctx;
	if(!InitializeNonlinearContext(ctx, IntrctPtr, /*MagnResetIsNotNeeded=*/1))
		return 0;

	ctx.use_b_input = true;

	int dof = ctx.totalDOF;
	int n_elem = ctx.AmOfMainElem;
	if(n_elem <= 0 || dof <= 0) return 0;

	// Collect hysteresis material pointers (energy or play).
	// Parallel with atomic early-exit flag: each element's dynamic_cast is
	// independent, but if ANY element fails we abort and return 0.
	std::vector<radTHysteresisMaterial*> hys_mats(n_elem);
	std::atomic<bool> all_hys{true};
	ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
		radTg3dRelax* g3d = IntrctPtr->g3dRelaxPtrVect[(int)i];
		radTMaterial* mat = (radTMaterial*)(g3d->MaterHandle.rep);
		hys_mats[i] = dynamic_cast<radTHysteresisMaterial*>(mat);
		if(!hys_mats[i]) all_hys.store(false, std::memory_order_relaxed);
	});
	if(!all_hys.load()) return 0;

	// Build base matrix (geometric interaction matrix N)
	if(NeedsDenseMatrix())
	{
		if(!BuildBaseMatrix(ctx, IntrctPtr))
			return 0;
	}

	// Auto-compute alpha if not provided: alpha >= max(dM/dH)
	// Probe multiple H values to find max susceptibility across full range
	// (descending hysteresis branch can have much higher chi than ascending)
	if(alpha <= 0.0)
	{
		double max_chi = 0.0;
		// Probe H magnitudes: logarithmically spaced from 10 to 100000 A/m
		double H_probes[] = {10.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0, 50000.0};
		int n_probes = 7;
		ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
			// Also check at actual H_ext
			int offset = IntrctPtr->GetElementDOFOffset((int)i);
			TVector3d H_ext_i(ctx.FlatExtern[offset], ctx.FlatExtern[offset+1], ctx.FlatExtern[offset+2]);

			// Save state before probing (probes are temporary)
			std::vector<TVector3d> probe_save;
			int ss = hys_mats[i]->GetStateSize();
			probe_save.resize(ss);
			hys_mats[i]->SaveState(probe_save);

			double local_max_chi = 0.0;
			for(int p = -1; p < n_probes; p++)
			{
				TVector3d H_probe;
				if(p < 0)
					H_probe = H_ext_i;
				else
					H_probe = TVector3d(0.0, 0.0, H_probes[p]);

				hys_mats[i]->RestoreState(probe_save);
				hys_mats[i]->Inverse(H_probe);  // H -> B (probing for max chi)
				TMatrix3d dBdH; double chi_d;
				hys_mats[i]->ComputeJacobian(dBdH, chi_d);
				if(chi_d > local_max_chi) local_max_chi = chi_d;
			}
			// Restore original state
			hys_mats[i]->RestoreState(probe_save);

			// Reduction across threads
			ngcore::AtomicMax(max_chi, local_max_chi);
		});
		alpha = max_chi * 1.5;  // 50% safety margin for hysteresis state variations
		if(alpha < 10.0) alpha = 10.0;
	}

	// Save hysteresis states (start-of-step reference)
	ctx.saved_hys_states.resize(n_elem);
	ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
		int state_size = hys_mats[i]->GetStateSize();
		ctx.saved_hys_states[i].resize(state_size);
		// SaveState returns m_Jk_prev as flat TVector3d array
		hys_mats[i]->SaveState(ctx.saved_hys_states[i]);
	});

	// Build LHS = I + alpha * BaseMatrix  (BaseMatrix = -N, so LHS = I - alpha*N)
	size_t mat_size = (size_t)dof * dof;
	std::vector<double> LHS(mat_size);
	for(size_t k = 0; k < mat_size; k++)
		LHS[k] = alpha * ctx.BaseMatrix[k];
	for(int i = 0; i < dof; i++)
		LHS[(size_t)i * dof + i] += 1.0;

	// LU factorize LHS ONCE
#ifdef HAVE_LAPACK
	std::vector<int> ipiv(dof);
	{
		int info = 0;
		int n = dof;
		dgetrf_(&n, &n, LHS.data(), &n, ipiv.data(), &info);
		if(info != 0) return 0;
	}
#else
	return 0;  // Hantila requires LAPACK for LU
#endif

	// Initialize: M from element current magnetization, R = M - alpha * H
	std::vector<double> M_vec(dof);
	std::vector<double> H_vec(dof);
	std::vector<double> R_vec(dof, 0.0);

	// Get current element M
	for(int i = 0; i < dof; i++)
		M_vec[i] = ctx.FlatMagn[i];

	// If M is all zero (virgin), initialize via Forward(H_ext)
	{
		double M_norm2 = 0.0;
		for(int i = 0; i < dof; i++) M_norm2 += M_vec[i] * M_vec[i];
		if(M_norm2 < 1e-30)
		{
			ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
				int offset = IntrctPtr->GetElementDOFOffset((int)i);
				TVector3d H_ext_i(ctx.FlatExtern[offset], ctx.FlatExtern[offset+1], ctx.FlatExtern[offset+2]);

				hys_mats[i]->RestoreState(ctx.saved_hys_states[i]);
				TVector3d B_fwd = hys_mats[i]->Inverse(H_ext_i);  // H -> B
				TVector3d M_fwd = (1.0 / MU_0) * B_fwd - H_ext_i;
				M_vec[offset+0] = M_fwd.x;
				M_vec[offset+1] = M_fwd.y;
				M_vec[offset+2] = M_fwd.z;
			});
		}
	}

	// Compute initial H = H_ext + N*M, then R = M - alpha*H
#ifdef HAVE_LAPACK
	cblas_dcopy(dof, ctx.FlatExtern, 1, H_vec.data(), 1);
	cblas_dgemv(CblasColMajor, CblasNoTrans, dof, dof,
	            -1.0, ctx.BaseMatrix.data(), dof,
	            M_vec.data(), 1, 1.0, H_vec.data(), 1);
#endif
	for(int i = 0; i < dof; i++)
		R_vec[i] = M_vec[i] - alpha * H_vec[i];

	// Estimate B_sat for convergence check
	double B_sat = 2.0;
	{
		double Js_sum = 0.0;
		ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
			ngcore::AtomicMax(Js_sum, hys_mats[i]->GetBsaturation());
		});
		if(Js_sum > B_sat) B_sat = Js_sum;
	}

	double final_residual = 1.0;
	int iterCount = 0;

	for(int iter = 0; iter < MaxIterNumber; iter++)
	{
		iterCount = iter + 1;

		// Step 1: Solve (I - alpha*N)*H = H_ext + N*R
		// RHS = H_ext + N*R = H_ext - BaseMatrix*R
		std::vector<double> rhs(dof);
#ifdef HAVE_LAPACK
		cblas_dcopy(dof, ctx.FlatExtern, 1, rhs.data(), 1);
		cblas_dgemv(CblasColMajor, CblasNoTrans, dof, dof,
		            -1.0, ctx.BaseMatrix.data(), dof,
		            R_vec.data(), 1, 1.0, rhs.data(), 1);

		// Back-substitution with pre-factored LHS
		{
			int n = dof;
			int nrhs = 1;
			int info = 0;
			cblas_dcopy(dof, rhs.data(), 1, H_vec.data(), 1);
			char trans = 'N';
			dgetrs_(&trans, &n, &nrhs, LHS.data(), &n, ipiv.data(),
			        H_vec.data(), &n, &info);
			if(info != 0) break;
		}
#endif

		// Step 2: M = R + alpha * H
		for(int i = 0; i < dof; i++)
			M_vec[i] = R_vec[i] + alpha * H_vec[i];

		// Step 3: B = mu_0 * (H + M), then Forward(B) per element
		std::vector<double> M_new(dof);

		ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
			int offset = IntrctPtr->GetElementDOFOffset((int)i);
			TVector3d H_i(H_vec[offset], H_vec[offset+1], H_vec[offset+2]);
			TVector3d M_i(M_vec[offset], M_vec[offset+1], M_vec[offset+2]);
			TVector3d B_i = MU_0 * (H_i + M_i);

			// Forward(B) -> H with restored start-of-step state
			hys_mats[i]->RestoreState(ctx.saved_hys_states[i]);
			TVector3d H_inv = hys_mats[i]->Forward(B_i);
			TVector3d M_model = (1.0 / MU_0) * B_i - H_inv;

			M_new[offset+0] = M_model.x;
			M_new[offset+1] = M_model.y;
			M_new[offset+2] = M_model.z;
		});

		// Convergence: ||M_model - M_linear||^2 / ||M_model||^2
		double sum_dM2 = 0.0, sum_M2 = 0.0;
		for(int i = 0; i < dof; i++)
		{
			double dM = M_new[i] - M_vec[i];
			sum_dM2 += dM * dM;
			sum_M2 += M_new[i] * M_new[i];
		}
		double rel_residual = (sum_M2 > 1e-30) ? std::sqrt(sum_dM2 / sum_M2) : std::sqrt(sum_dM2);

		// Step 4: Update R = M_new - alpha * H  (with optional under-relaxation)
		double omega = 1.0 - relax;
		for(int i = 0; i < dof; i++)
		{
			double R_target = M_new[i] - alpha * H_vec[i];
			if(relax > 0.0 && iter > 0)
				R_vec[i] = (1.0 - omega) * R_vec[i] + omega * R_target;
			else
				R_vec[i] = R_target;
		}

		// Convergence check (skip first iteration to let B_old settle)
		if(iter > 0)
		{
			final_residual = rel_residual;
			if(rel_residual < PrecOnMagnetiz)
				break;
		}

		if(radYield.Check() == 0) break;
	}

	// Store converged M back to elements and commit hysteresis state
	ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
		int offset = IntrctPtr->GetElementDOFOffset((int)i);
		TVector3d M_final(M_vec[offset], M_vec[offset+1], M_vec[offset+2]);
		ctx.FlatMagn[offset+0] = M_final.x;
		ctx.FlatMagn[offset+1] = M_final.y;
		ctx.FlatMagn[offset+2] = M_final.z;

		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[(int)i];
		g3dRelaxPtr->SetM(M_final);

		// Final Forward at converged B + CommitState
		hys_mats[i]->RestoreState(ctx.saved_hys_states[i]);
		TVector3d B_final = MU_0 * TVector3d(
			H_vec[offset] + M_final.x,
			H_vec[offset+1] + M_final.y,
			H_vec[offset+2] + M_final.z);
		hys_mats[i]->Forward(B_final);
		hys_mats[i]->CommitState();
	});

	// Update field arrays
#ifdef HAVE_LAPACK
	{
		cblas_dcopy(dof, ctx.FlatExtern, 1, H_vec.data(), 1);
		cblas_dgemv(CblasColMajor, CblasNoTrans, dof, dof,
		            -1.0, ctx.BaseMatrix.data(), dof,
		            M_vec.data(), 1, 1.0, H_vec.data(), 1);
		ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
			int offset = IntrctPtr->GetElementDOFOffset((int)i);
			ctx.FlatField[offset+0] = H_vec[offset+0];
			ctx.FlatField[offset+1] = H_vec[offset+1];
			ctx.FlatField[offset+2] = H_vec[offset+2];
		});
	}
#endif

	IntrctPtr->RelaxStatusParam.MisfitM = final_residual;
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

bool radTRelaxationMethNo_0::NeedsDenseMatrix() const
{
	return true;
}

//=========================================================================
// Method 0: B-input Newton Linear Step (dense Jacobian + LAPACK dgesv_)
//
// Assembles J_F = I - block_diag(dJ/dB) * (N+I)  (dof x dof)
// Solves J_F * dM = -F via LAPACK dgesv_ (LU with partial pivoting)
//=========================================================================

int radTRelaxationMethNo_0::SolveBInputLinearStep(
	NonlinearContext& ctx,
	const std::vector<double>& NpI,
	const std::vector<double>& dJdB_blocks,
	const std::vector<double>& F,
	std::vector<double>& dM)
{
	int dof = ctx.totalDOF;
	int n_elem = ctx.AmOfMainElem;
	size_t matrix_size = (size_t)dof * dof;

	// Allocate Jacobian matrix J_F (column-major for LAPACK)
	std::vector<double> J_F;
	try {
		J_F.resize(matrix_size, 0.0);
	} catch (const std::bad_alloc&) {
		return -2;
	}

	// Initialize J_F = I
	for(int i = 0; i < dof; i++)
		J_F[(size_t)i * dof + i] = 1.0;

	// Subtract dJ/dB * (N+I) block by block
	// J_F[o_i:o_i+3, :] -= dJdB_i @ NpI[o_i:o_i+3, :]
	for(int elem = 0; elem < n_elem; elem++)
	{
		int o_i = IntrctPtr->GetElementDOFOffset(elem);
		const double* dJdB_i = &dJdB_blocks[elem * 9];  // 3x3, column-major

#ifdef HAVE_LAPACK
		// BLAS dgemm: C(3,dof) -= A(3,3) * B(3,dof)
		// In column-major:
		//   A = dJdB_i at pointer, lda=3
		//   B = NpI rows o_i..o_i+2 = NpI + o_i, ldb=dof
		//   C = J_F rows o_i..o_i+2 = J_F + o_i, ldc=dof
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		            3, dof, 3,
		            -1.0,
		            dJdB_i, 3,
		            &NpI[o_i], dof,
		            1.0,
		            &J_F[o_i], dof);
#else
		// Manual: J_F[o_i+r, c] -= sum_d dJdB_i[r,d] * NpI[o_i+d, c]
		for(int c = 0; c < dof; c++)
		{
			for(int r = 0; r < 3; r++)
			{
				double sum = 0.0;
				for(int d = 0; d < 3; d++)
					sum += dJdB_i[d * 3 + r] * NpI[(size_t)c * dof + o_i + d];
				J_F[(size_t)c * dof + o_i + r] -= sum;
			}
		}
#endif
	}

	// RHS = -F (will be overwritten with solution by dgesv)
	dM.resize(dof);
	for(int i = 0; i < dof; i++)
		dM[i] = -F[i];

	// Solve J_F * dM = -F via LAPACK dgesv_
#ifdef HAVE_LAPACK
	std::vector<int> ipiv(dof);
	int nrhs = 1;
	int info = 0;

	{
		ngcore::SuspendTaskManager stm;
		radia::MKLThreadGuard mkl_guard(radia::GetNumThreads());
		dgesv_(&dof, &nrhs, J_F.data(), &dof, ipiv.data(), dM.data(), &dof, &info);
	}

	return (info == 0) ? 0 : -1;
#else
	// Fallback: use SolveLU_Flat
	return SolveLU_Flat(J_F, dM, dof);
#endif
}

//=========================================================================
// Variable DOF Solver Methods for hybrid surface-charge + standard element analysis.
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

	{
		ngcore::SuspendTaskManager stm;
		radia::MKLThreadGuard mkl_guard(radia::GetNumThreads());
		dgesv_(&n, &nrhs, A_col.data(), &n, ipiv.data(), b.data(), &n, &info);
	}

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
	// CRITICAL: Use size_t to avoid int32 overflow for DOF > 46340
	size_t matrix_size = (size_t)totalDOF * (size_t)totalDOF;
	// SystemMatrix (the dense O(N^2) working copy for dgesv) is allocated lazily;
	// ensureSystemMatrix() returns false on OOM (-> caller returns -2).
	std::vector<double> SystemMatrix;
	bool systemMatrixReady = false;
	auto ensureSystemMatrix = [&]() -> bool {
		if(systemMatrixReady) return true;
		try { SystemMatrix.resize(matrix_size); }
		catch (const std::bad_alloc&) {
			double required_gb = (double)matrix_size * 8 / (1024.0 * 1024.0 * 1024.0);
			fprintf(stderr, "Radia::Solve> LU solver requires %.1f GB memory for DOF=%d. Use BiCGSTAB (method 1) or HACApK (method 2) for large problems.\n", required_gb, totalDOF);
			return false;
		}
		systemMatrixReady = true;
		return true;
	};
	// System equation: (-K/(4pi) + I/chi) * sigma = H_ext_n (ELF-compatible); BaseMatrix already negated
	// in SetupBaseMatrix_VariableDOF.  Build RHS vector (will be overwritten with solution by dgesv).
	std::vector<double> RHS(totalDOF);

	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		if(IntrctPtr->GetElementDOF(elem) != 3)
		{
			fprintf(stderr, "Radia::Solve> mesh-backed magnetic material solve was removed from the legacy C++ relaxation path; use HDiv-VIM.\n");
			return -4;
		}
	}
	if(!ensureSystemMatrix()) return -2;
	std::memcpy(SystemMatrix.data(), ctx.BaseMatrix.data(), matrix_size * sizeof(double));

	// Update diagonal and RHS based on current chi
	// Matrix is ROW-MAJOR: A(i,j) at [i * totalDOF + j]
	// Diagonal element A(k,k) is at [k * totalDOF + k] = [k * (totalDOF + 1)]
	// Newton: use chi_d for system matrix (start after 10 Picard iterations)
	const int newton_start_iter = 10;
	bool newton_active = ctx.use_newton && iterCount >= newton_start_iter && !ctx.DifferentialChiArray.empty();

	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		double chi_abs = ctx.CurrentChiArray[elem];
		if(chi_abs < 1.0e-6) chi_abs = 1.0e-6;

		// Newton: use chi_d for matrix diagonal; Picard: use chi_abs
		// Newton: use chi_d for diagonal, add RHS correction
		double chi_matrix = chi_abs;
		if(newton_active)
		{
			chi_matrix = ctx.DifferentialChiArray[elem];
			if(chi_matrix < 1.0e-6) chi_matrix = 1.0e-6;
		}
		double inv_chi = 1.0 / chi_matrix;

		double newton_correction = 0.0;
		if(newton_active)
		{
			newton_correction = inv_chi - 1.0 / chi_abs;
		}


		// Add 1/chi to diagonal and set RHS
		// Newton: RHS = H_ext + (1/chi_d - 1/chi_abs) * sigma_old
		for(int k = 0; k < dof; k++)
		{
			int row = offset + k;
			SystemMatrix[row * (totalDOF + 1)] += inv_chi;
			RHS[row] = ctx.FlatExtern[row];
			if(newton_active)
			{
				RHS[row] += newton_correction * ctx.OldSigma[row];
			}
		}
	}

	// Solve using LAPACK LU (dgesv solves A*x = b in-place)
	auto t_lu_start = std::chrono::high_resolution_clock::now();
#ifdef HAVE_LAPACK
	std::vector<int> ipiv(totalDOF);
	int nrhs = 1;
	int info = 0;

	// CRITICAL: dgesv expects COLUMN-MAJOR format, but BaseMatrix is ROW-MAJOR
	// Transpose in-place: swap A[i,j] with A[j,i] for i < j
	// Each (i,j) pair is touched exactly once (j>i restriction), so different
	// outer-i rows never write the same cell -> safe to parallelize over i.
	ngcore::ParallelForRange(ngcore::IntRange(totalDOF), [&](ngcore::IntRange r) {
		for (auto i : r)
		{
			for(int j = (int)i + 1; j < totalDOF; j++)
			{
				// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
				std::swap(SystemMatrix[(size_t)i * totalDOF + j], SystemMatrix[(size_t)j * totalDOF + i]);
			}
		}
	});

		{
			// dgesv overwrites SystemMatrix with LU factors and RHS with the solution x0 = A^-1 b
			ngcore::SuspendTaskManager stm;
			radia::MKLThreadGuard mkl_guard(radia::GetNumThreads());
			dgesv_(&totalDOF, &nrhs, SystemMatrix.data(), &totalDOF, ipiv.data(), RHS.data(), &totalDOF, &info);
		}
		if(info != 0) return -1;  // Singular matrix
#else
	// Fallback: transpose to row-major and use SolveLU_Flat
	// Each (i,j) pair touched once (j>i); safe to parallelize over outer i.
	ngcore::ParallelForRange(ngcore::IntRange(totalDOF), [&](ngcore::IntRange r) {
		for (auto i : r)
		{
			for(int j = (int)i + 1; j < totalDOF; j++)
			{
				// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
				std::swap(SystemMatrix[(size_t)i * totalDOF + j], SystemMatrix[(size_t)j * totalDOF + i]);
			}
		}
	});
	int ierr = SolveLU_Flat(SystemMatrix, RHS, totalDOF);
	if(ierr != 0) return -1;
#endif
	auto t_lu_end = std::chrono::high_resolution_clock::now();
	double t_lu = std::chrono::duration<double>(t_lu_end - t_lu_start).count();
	rad.m_solve_t_lu_decomp += t_lu;
	rad.m_solve_t_linear_solve += t_lu;

	// Apply line search damping if Newton is active
	std::vector<double> sigma_trial = RHS;  // RHS contains solution after dgesv
	double omega = ApplyLineSearchDamping(ctx, IntrctPtr, sigma_trial);

	// If line search already updated FlatMagn (omega < 0.999), we're done
	// Otherwise copy trial solution (omega=1.0 case, full step)
	if(omega >= 0.999)
	{
		for(int i = 0; i < totalDOF; i++)
		{
			ctx.FlatMagn[i] = sigma_trial[i];
		}
	}
	// else: ApplyLineSearchDamping already updated FlatMagn with damped solution


	return 0;  // LU is direct solver, no iterations
}
