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
#include "rad_gmres.h"        // For templated restarted GMRES (KKT saddle, non-normal/indefinite)

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
// that handle both 3-DOF dipole and 4-6 DOF face-charge MSC elements
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

	// B-input PLAY-model moment Picard: cache hysteresis material pointers and
	// snapshot per-element start-of-step play state.  When set, the dof>=4
	// branch of UpdateChiAndCheckConvergence uses ComputeChiFromB(B) (B-input)
	// and restores the snapshot before each material evaluation so each element
	// keeps its own play trajectory; states are committed at end of the solve.
	ctx.b_input_play = rad.m_b_input_moment;
	if(ctx.b_input_play)
	{
		ctx.hys_mat_cache.assign(ctx.AmOfMainElem, nullptr);
		ctx.hys_play_state.assign(ctx.AmOfMainElem, std::vector<double>());
	}

	// Cache polyhedron pointers
	for(int elem = 0; elem < ctx.AmOfMainElem; elem++)
	{
		ctx.polyCache[elem] = dynamic_cast<radTPolyhedron*>(IntrctPtr->g3dRelaxPtrVect[elem]);
		if(ctx.b_input_play)
		{
			radTMaterial* MaterPtr0 = (radTMaterial*)(IntrctPtr->g3dRelaxPtrVect[elem]->MaterHandle.rep);
			radTHysteresisMaterial* HystMater0 = dynamic_cast<radTHysteresisMaterial*>(MaterPtr0);
			ctx.hys_mat_cache[elem] = HystMater0;
			if(HystMater0 != nullptr)
			{
				int ssz = HystMater0->GetStateSize();
				ctx.hys_play_state[elem].assign((size_t)(ssz > 0 ? ssz : 0), 0.0);
				if(ssz > 0) HystMater0->SaveStateToArray(ctx.hys_play_state[elem].data());
			}
		}
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

		// Store in poly->CurrentChi for face-charge MSC elements
		if(dof >= 4)
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
	// FlatInteract stores K/(4pi) for MSC and N for MMM
	// Both need negation to get -K/(4pi) and -N respectively
	std::memcpy(ctx.BaseMatrix.data(), FlatInteract, matrix_size * sizeof(double));

	// Negate entire matrix (physically correct for both MMM and MSC)
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
			// 3DOF MMM dipole elements (RecMag): FlatMagn contains M directly
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			g3dRelaxPtr->Magn.x = ctx.FlatMagn[offset + 0];
			g3dRelaxPtr->Magn.y = ctx.FlatMagn[offset + 1];
			g3dRelaxPtr->Magn.z = ctx.FlatMagn[offset + 2];
		}
		else if(dof >= 4)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
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
			if(poly_j && poly_j->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
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
			if(poly && poly->Use6DOF_MSC)
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
			if(poly && poly->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
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
			if(poly && poly->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
			{
				TVector3d H_new = IntrctPtr->NewFieldArray[elem];
				double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);
				TVector3d M_poly = poly->Magn;
				double M_mag = std::sqrt(M_poly.x*M_poly.x + M_poly.y*M_poly.y + M_poly.z*M_poly.z);


				if(HystMater != nullptr)
				{
					if(ctx.b_input_play)
					{
						// B-INPUT play update: build the element flux density from the
						// current solved state, B = mu0*(H + M), then chi = |M|/|H| with
						// H = Forward(B).  Restore the start-of-step play state before the
						// evaluation (Forward advances state) and save it back afterwards
						// so each element keeps its own play trajectory across iterations.
						TVector3d B_elem(MU_0 * (H_new.x + M_poly.x),
						                 MU_0 * (H_new.y + M_poly.y),
						                 MU_0 * (H_new.z + M_poly.z));
						if(elem < (int)ctx.hys_play_state.size() && !ctx.hys_play_state[elem].empty())
							HystMater->RestoreStateFromArray(ctx.hys_play_state[elem].data());
						chi_new = HystMater->ComputeChiFromB(B_elem);
						if(elem < (int)ctx.hys_play_state.size() && !ctx.hys_play_state[elem].empty())
							HystMater->SaveStateToArray(ctx.hys_play_state[elem].data());
					}
					else
					{
						chi_new = HystMater->ComputeChiFromH(H_new);
					}
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

static void CopyVectorToFlatMagn(NonlinearContext& ctx, const std::vector<double>& src)
{
	const int n = std::min(ctx.totalDOF, (int)src.size());
	for(int i = 0; i < n; i++) ctx.FlatMagn[i] = src[(size_t)i];
}

static void RestoreChiArray(NonlinearContext& ctx, radTInteraction* IntrctPtr, const std::vector<double>& chi)
{
	ctx.CurrentChiArray = chi;
	for(int elem = 0; elem < ctx.AmOfMainElem && elem < (int)chi.size(); elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		if(dof >= 4)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC) poly->CurrentChi = chi[(size_t)elem];
		}
	}
}

static double TryMomentAndersonAcceleration(NonlinearContext& ctx, radTInteraction* IntrctPtr, double baseRelChange)
{
	if(rad.m_moment_anderson_depth <= 0 || !ctx.last_solve_was_moment_hacapk || ctx.use_newton ||
	   ctx.totalDOF <= 0 || ctx.OldMagn.size() != (size_t)ctx.totalDOF)
		return baseRelChange;

	std::vector<double> baseMagn((size_t)ctx.totalDOF);
	std::vector<double> residual((size_t)ctx.totalDOF);
	for(int i = 0; i < ctx.totalDOF; i++)
	{
		baseMagn[(size_t)i] = ctx.FlatMagn[i];
		residual[(size_t)i] = baseMagn[(size_t)i] - ctx.OldMagn[(size_t)i];
	}

	if(!ctx.moment_anderson_have_prev ||
	   ctx.MomentAndersonPrevResidual.size() != (size_t)ctx.totalDOF ||
	   ctx.MomentAndersonPrevImage.size() != (size_t)ctx.totalDOF)
	{
		ctx.MomentAndersonPrevResidual = residual;
		ctx.MomentAndersonPrevImage = baseMagn;
		ctx.moment_anderson_have_prev = true;
		return baseRelChange;
	}

	double num = 0.0, den = 0.0;
	for(int i = 0; i < ctx.totalDOF; i++)
	{
		const double df = residual[(size_t)i] - ctx.MomentAndersonPrevResidual[(size_t)i];
		num += residual[(size_t)i] * df;
		den += df * df;
	}
	if(den <= 1.0e-300)
	{
		ctx.MomentAndersonPrevResidual = residual;
		ctx.MomentAndersonPrevImage = baseMagn;
		ctx.moment_anderson_rejected++;
		return baseRelChange;
	}

	const double gamma = num / den;
	if(!std::isfinite(gamma) || gamma < -1.0 || gamma > 1.0)
	{
		ctx.MomentAndersonPrevResidual = residual;
		ctx.MomentAndersonPrevImage = baseMagn;
		ctx.moment_anderson_rejected++;
		return baseRelChange;
	}

	std::vector<double> accelerated((size_t)ctx.totalDOF);
	for(int i = 0; i < ctx.totalDOF; i++)
	{
		const double gi = baseMagn[(size_t)i];
		const double gip = ctx.MomentAndersonPrevImage[(size_t)i];
		accelerated[(size_t)i] = gi - gamma * (gi - gip);
	}

	const std::vector<double> baseChi = ctx.CurrentChiArray;
	RestoreChiArray(ctx, IntrctPtr, ctx.OldChi);
	CopyVectorToFlatMagn(ctx, accelerated);
	ComputeActualHFieldFromSigma(ctx, IntrctPtr);
	const double accelRelChange = UpdateChiAndCheckConvergence(ctx, IntrctPtr);

	if(std::isfinite(accelRelChange) && accelRelChange < baseRelChange)
	{
		std::vector<double> acceptedResidual((size_t)ctx.totalDOF);
		for(int i = 0; i < ctx.totalDOF; i++)
			acceptedResidual[(size_t)i] = accelerated[(size_t)i] - ctx.OldMagn[(size_t)i];
		ctx.MomentAndersonPrevResidual.swap(acceptedResidual);
		ctx.MomentAndersonPrevImage.swap(accelerated);
		ctx.moment_anderson_accepted++;
		return accelRelChange;
	}

	CopyVectorToFlatMagn(ctx, baseMagn);
	RestoreChiArray(ctx, IntrctPtr, baseChi);
	ComputeActualHFieldFromSigma(ctx, IntrctPtr);
	ctx.max_B_rel_change = baseRelChange;
	ctx.MomentAndersonPrevResidual = residual;
	ctx.MomentAndersonPrevImage = baseMagn;
	ctx.moment_anderson_rejected++;
	return baseRelChange;
}

//-------------------------------------------------------------------------
/**
 * Apply adaptive line search damping to Newton-Raphson update.
 *
 * Finds optimal damping factor omega in [min_omega, 1.0] such that:
 *   sigma_new = omega * sigma_trial + (1 - omega) * sigma_old
 *
 * Uses backtracking line search with residual = max_B_rel_change as merit function.
 *
 * @param ctx Nonlinear context with trial solution to be damped
 * @param IntrctPtr Interaction data
 * @param sigma_trial Trial solution from linear solve [totalDOF]
 * @return Accepted omega value (1.0 = full Newton step, <1.0 = damped)
 */
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
				if(poly && poly->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
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
		ctx.last_solve_was_moment_hacapk = false;
		ctx.last_moment_linear_tol = 0.0;
		ctx.last_moment_krylov_solver = rad.m_moment_krylov_solver;

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
		rel_change = TryMomentAndersonAcceleration(ctx, IntrctPtr, rel_change);
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

	// B-input PLAY moment solve: commit each element's converged play state so the
	// next quasi-static step starts from the correct reference.  The per-element
	// hys_play_state holds (start-of-step pinning, converged m_pk_current) from the
	// final UpdateChiAndCheckConvergence evaluation; RestoreStateFromArray then
	// CommitState sets m_pk_prev = m_pk_pinning = converged m_pk_current.
	if(ctx.b_input_play)
	{
		for(int elem = 0; elem < ctx.AmOfMainElem; elem++)
		{
			radTHysteresisMaterial* HystMater = (elem < (int)ctx.hys_mat_cache.size()) ? ctx.hys_mat_cache[elem] : nullptr;
			if(HystMater != nullptr && elem < (int)ctx.hys_play_state.size() && !ctx.hys_play_state[elem].empty())
			{
				HystMater->RestoreStateFromArray(ctx.hys_play_state[elem].data());
				HystMater->CommitState();
			}
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
	// BaseMatrix stores the geometric part: -N for MMM (negated during BuildBaseMatrix)
	// So NpI = -BaseMatrix + I  (to get N + I)
	// Actually BaseMatrix = -N for MMM, so N = -BaseMatrix, and NpI = -BaseMatrix + I
	std::vector<double> NpI(static_cast<size_t>(dof) * dof, 0.0);
	for(size_t j = 0; j < (size_t)dof; j++)
	{
		for(size_t i = 0; i < (size_t)dof; i++)
		{
			// BaseMatrix is column-major: A(i,j) at index [j*dof + i]
			// For MMM: BaseMatrix = -N, so N = -BaseMatrix
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

	// Interaction matrix N (column-major): stored as -BaseMatrix for MMM
	// For MatVec: H = H_ext + N*M, where N = -BaseMatrix
	// So H = H_ext - BaseMatrix*M

	for(iterCount = 0; iterCount < MaxIterNumber; iterCount++)
	{
		// Restore all hysteresis states to beginning-of-step reference
		ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
			hys_mats[i]->RestoreState(ctx.saved_hys_states[i]);
		});

		// Compute H = H_ext + N*M
		// N*M = -BaseMatrix*M (since BaseMatrix = -N for MMM)
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

// The multipole-moment MMM path (BuildMomentSystemCore) assembles its own per-element moment system and never
// reads the dense interaction/base matrix, so the dense BaseMatrix can be skipped (no O(N^2)).  EXCEPTION: the
// B-input Newton/Hantila path (rad.m_b_input_newton / _hantila) builds NpI from the dense BaseMatrix, so it
// still needs it -> keep dense there.
bool radTRelaxationMethNo_0::NeedsDenseMatrix() const
{
	// The multipole-moment MMM path (BuildMomentSystemCore) assembles its own per-element moment system and
	// NEVER reads the dense interaction/base matrix.  When EVERY element is a surface-charge moment element
	// (tet 4 / wedge,pyramid 5 / hex 6 DOF) and we are not on the B-input hysteresis path (which builds NpI
	// from the dense BaseMatrix), the dense base matrix is not needed.  This lets the method-2 non-hex moment
	// reroute-to-LU (built with skipDenseMatrix) run -- without it BuildBaseMatrix returns false (no dense
	// InteractMatrix) and AutoRelax silently returns 0 iterations.  Mirrors radTRelaxationMethNo_1.
	// The dense BaseMatrix is needed ONLY by the dense 3-DOF B-input Newton/Hantila
	// path (which builds NpI from it).  The all-moment B-input PLAY route
	// (m_b_input_moment) uses BuildMomentSystemCore and does NOT read the dense
	// matrix, so treat it like the ordinary moment path (no dense build).
	bool b_input_dense = (rad.m_b_input_newton || rad.m_b_input_hantila) && !rad.m_b_input_moment;
	if(!b_input_dense && IntrctPtr != nullptr)
	{
		int nElem = IntrctPtr->AmOfMainElem;
		if(nElem > 0)
		{
			bool allMoment = true;
			for(int elem = 0; elem < nElem; elem++)
			{
				int dof = IntrctPtr->GetElementDOF(elem);
				if(dof != 6 && dof != 5 && dof != 4) { allMoment = false; break; }
			}
			if(allMoment) return false;
		}
	}
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
// Method 1: BiCGSTAB Iterative Solver (default)
//=========================================================================

double radTRelaxationMethNo_1::Dot(const std::vector<double>& a, const std::vector<double>& b, int n)
{
	double sum = 0.0;
	ngcore::ParallelForRange(ngcore::IntRange(n), [&](ngcore::IntRange r) {
		double local_sum = 0.0;
		for (auto i : r) {
			local_sum += a[i] * b[i];
		}
		ngcore::AtomicAdd(sum, local_sum);
	});
	return sum;
}

double radTRelaxationMethNo_1::Norm2(const std::vector<double>& a, int n)
{
	return std::sqrt(Dot(a, a, n));
}

void radTRelaxationMethNo_1::Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n)
{
	ngcore::ParallelFor(ngcore::IntRange(n), [&](size_t i) {
		y[i] += alpha * x[i];
	});
}

void radTRelaxationMethNo_1::Copy(const std::vector<double>& src, std::vector<double>& dst, int n)
{
	ngcore::ParallelFor(ngcore::IntRange(n), [&](size_t i) {
		dst[i] = src[i];
	});
}

void radTRelaxationMethNo_1::Scale(double alpha, std::vector<double>& x, int n)
{
	ngcore::ParallelFor(ngcore::IntRange(n), [&](size_t i) {
		x[i] *= alpha;
	});
}

void radTRelaxationMethNo_1::GetDiagonalElements(std::vector<double>& diag, const std::vector<double>& inv_chi, int n_elem)
{
	// Extract diagonal elements from interaction matrix for legacy scalar-Jacobi diagnostics.
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

		// Diagonal of system matrix: A = -N + 1/chi (physically correct)
		// M = chi * H_total => (I/chi - N) * M = H_ext => (-N + 1/chi) * M = H_ext
		if(IntrcMat != nullptr)
		{
			// Nii is from InteractMatrix[i][i]
			TMatrix3df& Nii = IntrcMat[i][i];
			diag[3*i + 0] = -Nii.Str0.x + inv_chi_x;  // Physically correct: +1/chi
			diag[3*i + 1] = -Nii.Str1.y + inv_chi_y;  // Physically correct: +1/chi
			diag[3*i + 2] = -Nii.Str2.z + inv_chi_z;  // Physically correct: +1/chi
		}
		else
		{
			// Fallback: just use +1/chi as diagonal (no N contribution)
			diag[3*i + 0] = inv_chi_x;
			diag[3*i + 1] = inv_chi_y;
			diag[3*i + 2] = inv_chi_z;
		}
	}
}

void radTRelaxationMethNo_1::DenseMatVec(const std::vector<double>& x, std::vector<double>& y,
                                         const std::vector<double>& inv_chi, int ndof)
{
	// Computes y = A * x where A = -N + 1/chi (physically correct)
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
		ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
			// Use pre-computed 1/chi values
			double inv_chi_x = inv_chi[3*i + 0];
			double inv_chi_y = inv_chi[3*i + 1];
			double inv_chi_z = inv_chi[3*i + 2];

			// y[i] = -sum(N[i][j] * x[j]) + (1/chi) * x[i] (physically correct)
			double y0 = inv_chi_x * x[3*i + 0];  // +1/chi (physically correct)
			double y1 = inv_chi_y * x[3*i + 1];  // +1/chi (physically correct)
			double y2 = inv_chi_z * x[3*i + 2];  // +1/chi (physically correct)

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
		});
	}
}

void radTRelaxationMethNo_1::BuildFlatMatrix(std::vector<double>& A_flat, const std::vector<double>& inv_chi, int ndof)
{
	// Build flat matrix A = -N - diag(1/chi) for BLAS dgemv (ELF-compatible)
	// Stored in ROW-MAJOR order: A[target][source]
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

	// CRITICAL: Use size_t to avoid int32 overflow for DOF > 46340
	A_flat.resize((size_t)ndof * (size_t)ndof, 0.0);

	if(IntrcMat != nullptr)
	{
		ngcore::ParallelFor(ngcore::IntRange(n_elem), [&](size_t i) {
			for(int j = 0; j < n_elem; j++)
			{
				TMatrix3df& Nij = IntrcMat[i][j];

				// ROW-MAJOR: A[row*ndof + col] for A[row][col]
				int row_base = 3*i;
				int col_base = 3*j;

				// Row 0 (Hx response): A[3i+0, 3j+l] = -dHx/dM_l = -Str_l.x
				A_flat[(row_base + 0)*ndof + (col_base + 0)] = -Nij.Str0.x;  // -dHx/dMx
				A_flat[(row_base + 0)*ndof + (col_base + 1)] = -Nij.Str1.x;  // -dHx/dMy
				A_flat[(row_base + 0)*ndof + (col_base + 2)] = -Nij.Str2.x;  // -dHx/dMz

				// Row 1 (Hy response): A[3i+1, 3j+l] = -dHy/dM_l = -Str_l.y
				A_flat[(row_base + 1)*ndof + (col_base + 0)] = -Nij.Str0.y;  // -dHy/dMx
				A_flat[(row_base + 1)*ndof + (col_base + 1)] = -Nij.Str1.y;  // -dHy/dMy
				A_flat[(row_base + 1)*ndof + (col_base + 2)] = -Nij.Str2.y;  // -dHy/dMz

				// Row 2 (Hz response): A[3i+2, 3j+l] = -dHz/dM_l = -Str_l.z
				A_flat[(row_base + 2)*ndof + (col_base + 0)] = -Nij.Str0.z;  // -dHz/dMx
				A_flat[(row_base + 2)*ndof + (col_base + 1)] = -Nij.Str1.z;  // -dHz/dMy
				A_flat[(row_base + 2)*ndof + (col_base + 2)] = -Nij.Str2.z;  // -dHz/dMz
			}

			// Add diagonal 1/chi terms (physically correct: +1/chi)
			A_flat[(3*i + 0)*ndof + (3*i + 0)] += inv_chi[3*i + 0];
			A_flat[(3*i + 1)*ndof + (3*i + 1)] += inv_chi[3*i + 1];
			A_flat[(3*i + 2)*ndof + (3*i + 2)] += inv_chi[3*i + 2];
		});
	}
}

void radTRelaxationMethNo_1::DenseMatVec_BLAS(const std::vector<double>& A_flat, const std::vector<double>& x,
                                              std::vector<double>& y, int ndof)
{
#ifdef HAVE_LAPACK
	// Use Intel MKL CBLAS cblas_dgemv: y = alpha*A*x + beta*y
	// A is stored in ROW-MAJOR order: A[target][source]
	cblas_dgemv(CblasRowMajor, CblasNoTrans, ndof, ndof,
	            1.0, A_flat.data(), ndof, x.data(), 1,
	            0.0, y.data(), 1);
#else
	// Fallback: manual matrix-vector multiply (ROW-MAJOR)
	std::fill(y.begin(), y.end(), 0.0);
	for(int i = 0; i < ndof; i++)
	{
		for(int j = 0; j < ndof; j++)
		{
			y[i] += A_flat[i * ndof + j] * x[j];  // row-major: A[i][j] at i*n+j
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

bool radTRelaxationMethNo_1::NeedsDenseMatrix() const
{
	if(IntrctPtr == nullptr) return true;
	int nElem = IntrctPtr->AmOfMainElem;
	if(nElem <= 0) return true;
	for(int elem = 0; elem < nElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		if(dof != 6 && dof != 5 && dof != 4) return true;
	}
	return false;
}

//=========================================================================
// Variable DOF Solver Methods for hybrid collocation MMMM + standard element analysis.
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
	// SystemMatrix (the dense O(N^2) working copy for dgesv) is allocated LAZILY (Phase 2 Increment 4):
	// the moment+method2 H-BiCGSTAB path returns before any dense solve, so on that scalable path we never
	// pay the O(N^2).  Dense paths (moment dense-LU method-0 path and MMM dense LU) call
	// ensureSystemMatrix() first; it returns false on OOM (-> caller returns -2).
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

	// multipole-moment MMM: assemble the parameter-free MOMENT system (BuildMomentSystemCore) for surface-charge
	// polyhedra (tet 4-DOF + wedge/pyramid 5-DOF + hex 6-DOF).  A's COLUMNS are the face DOF, so dgesv's solution is
	// sigma in DOF order -- a drop-in for the retired EIEM2 transpose + dgesv + write-back.  moment is now
	// UNCONDITIONAL (the old EIEM2/moment opt-out was removed in Phase 3b-1).  Pure tet (3 DOF = MMM)
	// uses the dense MMM path below; mixed tet+MSC is rejected fail-loud in MakeAutoRelax.
	bool useMoment = true;
	if(useMoment) { for(int e = 0; e < AmOfMainElem; e++) { int dd = IntrctPtr->GetElementDOF(e); if(dd != 6 && dd != 5 && dd != 4) { useMoment = false; break; } } }

	if(useMoment)
	{
		if(ctx.use_newton)
		{
			fprintf(stderr, "Radia::Solve> newton_method=True is not implemented for the multipole-moment surface-charge path; refusing to run Picard silently.\n");
			return -4;
		}
		std::vector<int> momElem; IntrctPtr->CollectMomentElems(momElem);   // hex(6)+wedge/pyramid(5), matches BuildMomentSystemCore
		int nMom = (int)momElem.size();
		std::vector<double> chiPerHex((size_t)nMom), HextPerHex((size_t)nMom*3);
		for(int h = 0; h < nMom; h++)
		{
			int e = momElem[h];
			double chi_abs = ctx.CurrentChiArray[e]; if(chi_abs < 1.0e-6) chi_abs = 1.0e-6;
			chiPerHex[h] = chi_abs;
			const TVector3d& He = IntrctPtr->ExternFieldArray[e];      // external field at element centroid
			HextPerHex[(size_t)h*3+0] = He.x; HextPerHex[(size_t)h*3+1] = He.y; HextPerHex[(size_t)h*3+2] = He.z;
			radTPolyhedron* poly = ctx.polyCache[e]; if(poly && poly->Use6DOF_MSC) poly->CurrentChi = chi_abs;
		}
		// MMMM multipole-moment surface-charge solve: assemble the dense per-element moment system
		// (BuildMomentSystemCore) and solve it with the dense LU below.  MMMM does NOT connect to HACApK;
		// large-scale moment routing through the H-matrix was removed.
		if(!ensureSystemMatrix()) return -2;   // dense moment method-0 path
		IntrctPtr->BuildMomentSystemCore(chiPerHex.data(), HextPerHex.data(), SystemMatrix, RHS);
		rad.m_solve_t_moment_fieldgrad += IntrctPtr->LastMomentFieldGradTime();
		rad.m_solve_t_moment_system_build += IntrctPtr->LastMomentSystemBuildTime();
	}
	else
	{
	if(!ensureSystemMatrix()) return -2;   // pure-MMM dense path needs the dense BaseMatrix copy
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

		// Update poly->CurrentChi for face-charge MSC elements
		if(dof >= 4)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				poly->CurrentChi = chi_abs;  // Always store absolute chi
			}
		}
	}
	}  // end else (pure-MMM dense path; the moment path above already filled SystemMatrix + RHS)

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

//-------------------------------------------------------------------------
// BiCGSTAB solver variable DOF methods
//-------------------------------------------------------------------------

void radTRelaxationMethNo_1::MatVec_VariableDOF(const std::vector<double>& x, std::vector<double>& y,
                                                 const std::vector<double>& inv_chi, int totalDOF)
{
	// Computes y = A * x where A = (base matrix) - diag(1/chi) (ELF-compatible)
	// Sign convention:
	// - MMM (3 DOF): FlatInteract stores N, equation is (-N - I/chi) -> negate N, subtract I/chi
	// - MSC (6 DOF): FlatInteract stores -K/(4pi), equation is (-K/(4pi) - I/chi) -> use as-is, subtract I/chi
	//
	// IMPORTANT: FlatInteract is stored in ROW-MAJOR format (C/NumPy style)
	// Element at (row, col) is at index [row * totalDOF + col]

	const double* FlatInteract = IntrctPtr->GetFlatInteractMatrix();
	if(FlatInteract == nullptr) return;

	// TaskManager path for all dense BiCGSTAB systems (pure 3DOF, pure 6DOF, and mixed DOF).
	// FlatInteract stores the demag/base matrix in row-major target/source order.  The BiCGSTAB
	// equation uses -FlatInteract plus the positive diagonal 1/chi contribution.
	ngcore::ParallelFor(ngcore::IntRange(totalDOF), [&](size_t row) {
		const double* Arow = &FlatInteract[(size_t)row * totalDOF];
		double sum = 0.0;
		for(int col = 0; col < totalDOF; col++) sum += Arow[col] * x[col];
		y[row] = -sum + inv_chi[row] * x[row];
	});
}

void radTRelaxationMethNo_1::GetDiagonalElements_VariableDOF(std::vector<double>& diag,
                                                              const std::vector<double>& inv_chi, int totalDOF)
{
	// Extract diagonal elements for legacy scalar-Jacobi diagnostics.
	// Physical equation: A = -K/(4pi) + diag(1/chi)
	// FlatInteract stores K/(4pi) for MSC and N for MMM - negate both
	const double* FlatInteract = IntrctPtr->GetFlatInteractMatrix();
	if(FlatInteract == nullptr) return;

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		// Get diagonal block
		// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
		const double* diag_block = &FlatInteract[(size_t)offset * totalDOF + offset];

		// Negate for all element types (physical equation: A = -K + diag(1/chi))
		double sign = -1.0;

		for(int k = 0; k < dof; k++)
		{
			// Diagonal element: sign*matrix_ii + 1/chi (physically correct)
			// CRITICAL: Use size_t cast for indexing with totalDOF
			diag[offset + k] = sign * diag_block[(size_t)k * totalDOF + k] + inv_chi[offset + k];
		}
	}
}

#ifdef HAVE_LAPACK
// dgetrf_ and dgetri_ are provided by Intel MKL headers (mkl_lapack.h)
// included at the top of this file

bool radTRelaxationMethNo_1::BuildBlockJacobiPreconditioner_VariableDOF(
	std::vector<double>& blockInverse, std::vector<int>& blockOffsets,
	const std::vector<double>& inv_chi, int totalDOF)
{
	// Build element-block Jacobi preconditioner by inverting each element's natural diagonal block.
	// This is much better than scalar Jacobi for poorly conditioned MSC matrices
	const double* FlatInteract = IntrctPtr->GetFlatInteractMatrix();
	if(FlatInteract == nullptr) return false;

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	// Calculate total storage needed for block inverses
	int total_block_storage = 0;
	blockOffsets.resize(AmOfMainElem + 1);
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		blockOffsets[elem] = total_block_storage;
		total_block_storage += dof * dof;
	}
	blockOffsets[AmOfMainElem] = total_block_storage;
	blockInverse.resize(total_block_storage);

	// Process each element's diagonal block
	int max_dof = 6;  // Maximum DOF per element (hexahedra)
	std::vector<double> block_copy(max_dof * max_dof);
	std::vector<int> ipiv(max_dof);
	std::vector<double> work(max_dof * max_dof);
	int lwork = max_dof * max_dof;

	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int mat_offset = IntrctPtr->GetElementDOFOffset(elem);
		int block_offset = blockOffsets[elem];

		// Extract diagonal block: A_block = -K_block + inv_chi * I
		// FlatInteract stores K/(4pi) in ROW-MAJOR [target][source] format
		for(int i = 0; i < dof; i++)
		{
			for(int j = 0; j < dof; j++)
			{
				// Row-major access: A[row][col] at FlatInteract[(row)*totalDOF + (col)]
				double K_ij = FlatInteract[(size_t)(mat_offset + i) * totalDOF + (mat_offset + j)];
				// Store in COLUMN-MAJOR for LAPACK (block_copy[i + j*dof] = element at row i, col j)
				block_copy[i + j * dof] = -K_ij;
				if(i == j)
				{
					block_copy[i + j * dof] += inv_chi[mat_offset + i];
				}
			}
		}

		// Invert the block using LAPACK: LU factorization then inverse
		int info = 0;
		dgetrf_(&dof, &dof, block_copy.data(), &dof, ipiv.data(), &info);
		if(info != 0)
		{
			fprintf(stderr, "[Block Jacobi] Element %d: singular diagonal block (info=%d)\n", elem, info);
			return false;
		}
		dgetri_(&dof, block_copy.data(), &dof, ipiv.data(), work.data(), &lwork, &info);
		if(info != 0)
		{
			fprintf(stderr, "[Block Jacobi] Element %d: inversion failed (info=%d)\n", elem, info);
			return false;
		}

		// Store inverse block in ROW-MAJOR format for efficient application
		for(int i = 0; i < dof; i++)
		{
			for(int j = 0; j < dof; j++)
			{
				// Convert from LAPACK column-major to row-major
				blockInverse[block_offset + i * dof + j] = block_copy[i + j * dof];
			}
		}
	}

	return true;
}

void radTRelaxationMethNo_1::ApplyBlockJacobiPreconditioner_VariableDOF(
	const std::vector<double>& x, std::vector<double>& y,
	const std::vector<double>& blockInverse, const std::vector<int>& blockOffsets)
{
	// Apply element-block Jacobi preconditioner: y = M^{-1} * x.
	// where M is the block-diagonal of A
	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	ngcore::ParallelFor(ngcore::IntRange(AmOfMainElem), [&](size_t elem) {
		int dof = IntrctPtr->GetElementDOF(elem);
		int mat_offset = IntrctPtr->GetElementDOFOffset(elem);
		int block_offset = blockOffsets[elem];

		// y_elem = inv_block * x_elem (block is stored row-major)
		for(int i = 0; i < dof; i++)
		{
			double sum = 0.0;
			for(int j = 0; j < dof; j++)
			{
				sum += blockInverse[block_offset + i * dof + j] * x[mat_offset + j];
			}
			y[mat_offset + i] = sum;
		}
	});
}
#endif

int radTRelaxationMethNo_1::SolveBiCGSTAB_VariableDOF(NonlinearContext& ctx,
                                                       int totalDOF, double tol, int max_iter, double& residual,
                                                       const std::vector<double>& elemChiArray,
                                                       bool use_newton,
                                                       const std::vector<double>* absChiArray,
                                                       const double* oldSigma)
{
	// BiCGSTAB with element-block Jacobi preconditioner for variable DOF systems.
	// Keep a TaskManager region active across the whole method-1 solve so every ParallelFor
	// matvec/vector/preconditioner operation is multi-threaded even from a bare rad.Solve(..., 1).
	ngcore::RegionTaskManager rtm(radia::GetMaxThreads());

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	// Allocate work vectors
	std::vector<double> r(totalDOF), r0(totalDOF), p(totalDOF), v(totalDOF), s(totalDOF), t(totalDOF);
	std::vector<double> p_hat(totalDOF), s_hat(totalDOF);
	std::vector<double> inv_chi(totalDOF);
	std::vector<double> rhs(totalDOF);
	std::vector<double> sol(totalDOF);

	double* FlatMagn = IntrctPtr->GetFlatMagnArray();
	double* FlatField = IntrctPtr->GetFlatFieldArray();
	double* FlatExtern = IntrctPtr->GetFlatExternFieldArray();

	if(FlatMagn == nullptr || FlatField == nullptr || FlatExtern == nullptr) return 0;

	// Pre-compute 1/chi and RHS for all elements
	// Picard: (D(1/chi_abs) + G) sigma = H_ext
	// Newton: (D(1/chi_d) + G) sigma_new = H_ext + D(1/chi_d - 1/chi_abs) * sigma_old
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		// elemChiArray: chi_d (Newton) or chi_abs (Picard)
		double chi_matrix = elemChiArray[elem];
		if(chi_matrix < 1.0e-6) chi_matrix = 1.0e-6;
		double inv_chi_val = 1.0 / chi_matrix;

		double newton_correction = 0.0;
		if(use_newton && absChiArray && oldSigma)
		{
			double chi_abs = (*absChiArray)[elem];
			if(chi_abs < 1.0e-6) chi_abs = 1.0e-6;
			newton_correction = inv_chi_val - 1.0 / chi_abs;
		}

		for(int k = 0; k < dof; k++)
		{
			inv_chi[offset + k] = inv_chi_val;
			rhs[offset + k] = FlatExtern[offset + k];
			if(use_newton && oldSigma)
			{
				rhs[offset + k] += newton_correction * oldSigma[offset + k];
			}
		}
	}

	// Initial guess: use current magnetization (ELF-compatible)
	// Using the previous solution as initial guess significantly speeds up convergence
	for(int i = 0; i < totalDOF; i++)
	{
		sol[i] = FlatMagn[i];
	}

	// Build preconditioner
	std::vector<double> blockInverse;
	std::vector<int> blockOffsets;

#ifdef HAVE_LAPACK
	if(!BuildBlockJacobiPreconditioner_VariableDOF(blockInverse, blockOffsets, inv_chi, totalDOF))
	{
		fprintf(stderr, "[BiCG] Error: block-Jacobi build failed; scalar Jacobi fallback is disabled.\n");
		return -11;
	}
#else
	fprintf(stderr, "[BiCG] Error: block-Jacobi requires LAPACK; scalar Jacobi fallback is disabled.\n");
	return -11;
#endif

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
		rho = this->Dot(r0, r, totalDOF);

		if(std::abs(rho) < 1.0e-30)
		{
			residual = this->Norm2(r, totalDOF) / rhs_norm;
			break;
		}

		if(iter == 1)
		{
			this->Copy(r, p, totalDOF);
		}
		else
		{
			if(std::abs(rho_old * omega) < 1.0e-30)
			{
				residual = this->Norm2(r, totalDOF) / rhs_norm;
				break;
			}
			double beta = (rho / rho_old) * (alpha_bicg / omega);
			this->Axpy(-omega, v, p, totalDOF);
			this->Scale(beta, p, totalDOF);
			this->Axpy(1.0, r, p, totalDOF);
		}

		// Apply preconditioner: p_hat = M^{-1} * p
		ApplyBlockJacobiPreconditioner_VariableDOF(p, p_hat, blockInverse, blockOffsets);

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

		// Apply preconditioner: s_hat = M^{-1} * s
		ApplyBlockJacobiPreconditioner_VariableDOF(s, s_hat, blockInverse, blockOffsets);

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

		// Only detect true numerical blowup (NaN/Inf or extreme divergence)
		// BiCGSTAB naturally has non-monotonic convergence - do NOT stop on temporary spikes
		if(std::isnan(residual) || std::isinf(residual) || residual > 1.0e15)
		{
			fprintf(stderr, "[BiCG] Exit: numerical blowup at iter %d (residual=%.4e)\n", iter, residual);
			break;
		}

		if(std::abs(omega) < 1.0e-30) {
			fprintf(stderr, "[BiCG] Exit: omega=%.4e at iter %d\n", omega, iter);
			break;
		}
	}

	// Apply line search damping if Newton is active
	std::vector<double> sigma_trial = sol;
	double omega_ls = ApplyLineSearchDamping(ctx, this->IntrctPtr, sigma_trial);

	// If line search already updated FlatMagn (omega_ls < 0.999), we're done
	// Otherwise copy trial solution (omega_ls=1.0 case, full step)
	if(omega_ls >= 0.999)
	{
		for(int i = 0; i < totalDOF; i++)
		{
			FlatMagn[i] = sigma_trial[i];
		}
	}
	// else: ApplyLineSearchDamping already updated FlatMagn with damped solution

	return iter;
}

//-------------------------------------------------------------------------
// BiCGSTAB Solver: SolveLinearStep override
// Uses BiCGSTAB iterative solver with element-block Jacobi preconditioner.
//-------------------------------------------------------------------------

int radTRelaxationMethNo_1::SolveLinearStep(NonlinearContext& ctx, int iterCount)
{
	int totalDOF = ctx.totalDOF;
	int AmOfMainElem = ctx.AmOfMainElem;

	bool useMoment = (AmOfMainElem > 0);
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		if(dof != 6 && dof != 5 && dof != 4)
		{
			useMoment = false;
			break;
		}
	}

	if(useMoment)
	{
		if(ctx.use_newton)
		{
			fprintf(stderr, "Radia::Solve> newton_method=True is not implemented for the multipole-moment surface-charge path; refusing to run Picard silently.\n");
			return -4;
		}
		// Keep the TaskManager active for dense moment assembly, matvec, and block-Jacobi application.
		// Without this region the method-1 moment branch falls back to effectively serial ngcore::ParallelFor
		// execution, which makes MDX scaling measurements misleading.
		ngcore::RegionTaskManager rtm(radia::GetMaxThreads());

		std::vector<int> momElem;
		IntrctPtr->CollectMomentElems(momElem);
		int nMom = (int)momElem.size();
		if(nMom <= 0) return 0;

		std::vector<double> chiPerHex((size_t)nMom), HextPerHex((size_t)nMom * 3);
		for(int h = 0; h < nMom; h++)
		{
			int elem = momElem[h];
			double chi_abs = ctx.CurrentChiArray[elem];
			if(chi_abs < 1.0e-6) chi_abs = 1.0e-6;
			chiPerHex[h] = chi_abs;

			const TVector3d& He = IntrctPtr->ExternFieldArray[elem];
			HextPerHex[(size_t)h * 3 + 0] = He.x;
			HextPerHex[(size_t)h * 3 + 1] = He.y;
			HextPerHex[(size_t)h * 3 + 2] = He.z;

			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC) poly->CurrentChi = chi_abs;
		}

		const int nHex = IntrctPtr->GetNumHexElements();
		const std::vector<int>& hexElem = IntrctPtr->GetHexaElemIndices();
		bool pureHexMatrixFree = (nMom == nHex && totalDOF == 6 * nHex);
		for(int h = 0; pureHexMatrixFree && h < nHex; h++)
		{
			const int elem = hexElem[h];
			if(h >= (int)momElem.size() || momElem[h] != elem ||
			   IntrctPtr->GetElementDOF(elem) != 6 ||
			   IntrctPtr->GetElementDOFOffset(elem) != 6 * h)
			{
				pureHexMatrixFree = false;
			}
		}

		if(pureHexMatrixFree)
		{
			const bool serialElementOps = (nHex <= 8);
			std::vector<double> RHS((size_t)totalDOF, 0.0);
			if(serialElementOps)
			{
				for(int h = 0; h < nHex; h++)
					for(int t = 0; t < 3; t++) RHS[(size_t)6 * h + t] = chiPerHex[h] * HextPerHex[(size_t)3 * h + t];
			}
			else
			{
				ngcore::ParallelFor(ngcore::IntRange(nHex), [&](size_t h) {
					for(int t = 0; t < 3; t++) RHS[(size_t)6 * h + t] = chiPerHex[h] * HextPerHex[(size_t)3 * h + t];
				});
			}

			auto t_setup_start = std::chrono::high_resolution_clock::now();
			IntrctPtr->PrecomputeMomentGeometry();

			std::vector<double> localLBlock((size_t)nHex * 36, 0.0);
			std::vector<double> diagKBlock((size_t)nHex * 36, 0.0);
			std::vector<double> zeroChi((size_t)nHex, 0.0);
			auto buildLocalBlocks = [&](size_t hh) {
				const int h = (int)hh;
				IntrctPtr->MomentSystemBlock6x6(h, h, zeroChi.data(), &localLBlock[(size_t)h * 36]);
				IntrctPtr->MomentSystemBlock6x6(h, h, nullptr, &diagKBlock[(size_t)h * 36], true);
			};
			if(serialElementOps) { for(int h = 0; h < nHex; h++) buildLocalBlocks((size_t)h); }
			else { ngcore::ParallelFor(ngcore::IntRange(nHex), buildLocalBlocks); }

			std::vector<double> blockInverse((size_t)nHex * 36, 0.0);
#ifdef HAVE_LAPACK
			std::atomic<int> bj_bad_elem{-1};
			std::atomic<int> bj_bad_info{0};
			auto buildPrecondBlock = [&](size_t hh) {
				const int h = (int)hh;
				double rawBlock[36];
				double block_copy[36];
				int ipiv[6];
				double work[36];
				int six = 6, lwork = 36, info = 0;
				const double* Lh = &localLBlock[(size_t)h * 36];
				const double* Kh = &diagKBlock[(size_t)h * 36];
				for(int k = 0; k < 36; k++) rawBlock[k] = Lh[k] + chiPerHex[h] * Kh[k];
				for(int i = 0; i < 6; i++)
					for(int j = 0; j < 6; j++)
						block_copy[i + j * 6] = rawBlock[i * 6 + j];   // LAPACK column-major

				dgetrf_(&six, &six, block_copy, &six, ipiv, &info);
				if(info == 0) dgetri_(&six, block_copy, &six, ipiv, work, &lwork, &info);
				if(info != 0)
				{
					int expected = -1;
					if(bj_bad_elem.compare_exchange_strong(expected, h)) bj_bad_info.store(info);
					return;
				}

				double* Binv = &blockInverse[(size_t)h * 36];
				for(int i = 0; i < 6; i++)
					for(int j = 0; j < 6; j++)
						Binv[i * 6 + j] = block_copy[i + j * 6];
			};
			if(serialElementOps) { for(int h = 0; h < nHex; h++) buildPrecondBlock((size_t)h); }
			else { ngcore::ParallelFor(ngcore::IntRange(nHex), buildPrecondBlock); }
			if(bj_bad_elem.load() >= 0)
			{
				fprintf(stderr, "[BiCG] Error: matrix-free moment block-Jacobi failed at hex %d (info=%d); scalar/identity substitute is disabled.\n",
				        bj_bad_elem.load(), bj_bad_info.load());
				return -11;
			}
#else
			fprintf(stderr, "[BiCG] Error: matrix-free moment block-Jacobi requires LAPACK; scalar/identity substitute is disabled.\n");
			return -11;
#endif

			// BUILD-ONCE geometry coupling K (Sugahara 2026-06-27): the matrix-free MomentKernelMatVec6x6
			// RE-evaluates the full O(nHex^2) kernel on EVERY matvec, so a BiCGSTAB solve (~tens of matvecs)
			// repeats ~N matrix builds.  Instead build the chi-INDEPENDENT geometry block matrix K ONCE
			// (K[h][g] = MomentSystemBlock6x6(h,g,kernelOnly=true), incl. INV4PI + IMA images); the matvec
			// applies y = diag(chi)(K x) + L_local x, and a nonlinear Picard loop reuses K (only chi, the
			// block diagonal, changes).
			//   method 0/1 -> DENSE K (momentKdense, O(N^2) storage, cheap dense GEMV).
			//   method 2   -> H-MATRIX K (RadHACApKMomentSystem kernel-only, O(N log N) matvec + storage) --
			//                 the collocation-MMMM COARSE tier (Sugahara 2026-07-02, loop-free abandoned;
			//                 field-correct but loop-polluted internal M, acceptable for coarse/optimization).
			const bool useHMatrixK = rad.m_moment_use_hmatrix;
			std::vector<double> momentKdense;              // dense K   (built when !useHMatrixK)
			std::vector<double> momentHKx;                 // matvec input copy (the H-matvec uses the vector API)
			if(useHMatrixK)
			{
				// BUILD-ONCE across the Picard loop: m_momentHK is a member, so the H-matrix K survives
				// between SolveLinearStep calls of one solve (K is chi-free; only chi changes per iteration).
				if(m_momentHK && m_momentHK->GetInteraction() != IntrctPtr)
				{
					delete m_momentHK;
					m_momentHK = nullptr;
				}
				if(!m_momentHK)
				{
					m_momentHK = new RadHACApKMomentSystem(IntrctPtr);
					RadHACApKParams hkParams;
					hkParams.aca_eps   = rad.m_hacapk_eps;
					hkParams.leaf_size = rad.m_hacapk_leaf_size;
					hkParams.eta       = rad.m_hacapk_eta;
					hkParams.print_level = 0;
					if(!m_momentHK->BuildHMatrix(hkParams))   // self-wraps a RegionTaskManager; nests under rtm harmlessly
					{
						delete m_momentHK;
						m_momentHK = nullptr;
						fprintf(stderr, "[BiCG] Error: collocation-MMMM moment H-matrix (method 2) build failed; no dense fallback.\n");
						return -12;
					}
				}
				momentHKx.assign((size_t)totalDOF, 0.0);
			}
			else
			{
				// O(N^2) dense; research code -- NO artificial size cap.  A too-large allocation fails loud
				// (std::bad_alloc); use method 2 (HACApK, O(N log N)) for very large N.
				const size_t momentDenseElems = (size_t)totalDOF * (size_t)totalDOF;
				momentKdense.assign(momentDenseElems, 0.0);
				auto buildKdenseRow = [&](size_t hh) {
					const int h = (int)hh;
					double Kblk[36];
					for(int g = 0; g < nHex; g++)
					{
						IntrctPtr->MomentSystemBlock6x6(h, g, nullptr, Kblk, true);   // chi-independent geometry block (incl. images)
						for(int i = 0; i < 6; i++)
							for(int j = 0; j < 6; j++)
								momentKdense[(size_t)(6*h+i) * (size_t)totalDOF + (size_t)(6*g+j)] = Kblk[i*6+j];
					}
				};
				if(serialElementOps) { for(int h = 0; h < nHex; h++) buildKdenseRow((size_t)h); }
				else { ngcore::ParallelFor(ngcore::IntRange(nHex), buildKdenseRow); }
			}
			std::vector<double> momentKx((size_t)totalDOF, 0.0);

			auto t_setup_end = std::chrono::high_resolution_clock::now();
			rad.m_solve_t_moment_system_build += std::chrono::duration<double>(t_setup_end - t_setup_start).count();

			std::vector<double> sigma(totalDOF);
			for(int i = 0; i < totalDOF; i++) sigma[i] = ctx.FlatMagn[i];

			auto matvec = [&](const double* x, double* y) {
				// build-once: y = diag(chi)(K x) + L_local x.  K x via dense GEMV (method 0/1) or the
				// RadHACApKMomentSystem H-matvec (method 2, O(N log N)); NO kernel recompute either way.
				if(useHMatrixK)
				{
					std::memcpy(momentHKx.data(), x, (size_t)totalDOF * sizeof(double));
					m_momentHK->MatVec(momentHKx, momentKx);   // momentKx = K x  (y overwritten, no accumulate)
				}
				else
				{
					cblas_dgemv(CblasRowMajor, CblasNoTrans, totalDOF, totalDOF, 1.0,
					            momentKdense.data(), totalDOF, x, 1, 0.0, momentKx.data(), 1);
				}
				auto applyRow = [&](size_t hh) {
					const int h = (int)hh;
					const double chih = chiPerHex[h];
					const double* Lh = &localLBlock[(size_t)h * 36];
					const double* xh = &x[(size_t)6 * h];
					const double* Kxh = &momentKx[(size_t)6 * h];
					double* yh = &y[(size_t)6 * h];
					for(int i = 0; i < 6; i++)
					{
						double s = chih * Kxh[i];
						const double* Li = &Lh[(size_t)i * 6];
						for(int j = 0; j < 6; j++) s += Li[j] * xh[j];
						yh[i] = s;
					}
				};
				if(serialElementOps) { for(int h = 0; h < nHex; h++) applyRow((size_t)h); }
				else { ngcore::ParallelFor(ngcore::IntRange(nHex), applyRow); }
			};
			auto precond = [&](const double* x, double* y) {
				auto applyBlock = [&](size_t hh) {
					const int h = (int)hh;
					const double* Binv = &blockInverse[(size_t)h * 36];
					const double* xh = &x[(size_t)6 * h];
					double* yh = &y[(size_t)6 * h];
					for(int i = 0; i < 6; i++)
					{
						double sum = 0.0;
						for(int j = 0; j < 6; j++) sum += Binv[i * 6 + j] * xh[j];
						yh[i] = sum;
					}
				};
				if(serialElementOps) { for(int h = 0; h < nHex; h++) applyBlock((size_t)h); }
				else { ngcore::ParallelFor(ngcore::IntRange(nHex), applyBlock); }
			};

			const double bicg_tol = rad.m_bicg_tol;
			const int bicg_max_iter = 10000;

			auto t_bicg_start = std::chrono::high_resolution_clock::now();
			radia::bicgstab::Result result;
			if(rad.m_moment_krylov_solver == 1)
			{
				// GMRES for the plain (non-loop-free) solve -- makes the m_moment_krylov_solver=1 flag LIVE
				// (was a dead stub).  GMRES is robust where BiCGSTAB stalls on the non-normal operator.
				radia::gmres::Result gres = radia::gmres::Solve(totalDOF, matvec, precond, RHS.data(), sigma.data(),
				                                                bicg_tol, bicg_max_iter, rad.m_moment_gmres_restart);
				result.iterations = gres.iterations; result.residual = gres.residual; result.converged = gres.converged;
			}
			else
			{
				result = radia::bicgstab::Solve<double>(totalDOF, matvec, precond, RHS.data(), sigma.data(),
				                                        bicg_tol, bicg_max_iter);
			}
			auto t_bicg_end = std::chrono::high_resolution_clock::now();
			rad.m_solve_t_linear_solve += std::chrono::duration<double>(t_bicg_end - t_bicg_start).count();
			rad.m_solve_linear_iterations += result.iterations;
			ctx.last_solve_was_moment_hacapk = useHMatrixK;
			ctx.last_moment_linear_tol = bicg_tol;
			ctx.last_moment_krylov_solver = rad.m_moment_krylov_solver;

			if(!result.converged)
			{
				fprintf(stderr, "[BiCG] Warning: matrix-free moment BiCGSTAB did not reach tol %.3e after %d iterations (residual=%.4e)\n",
				        bicg_tol, result.iterations, result.residual);
				return -3;
			}

			std::vector<double> sigma_trial = sigma;
			double omega_ls = ApplyLineSearchDamping(ctx, this->IntrctPtr, sigma_trial);
			if(omega_ls >= 0.999)
			{
				for(int i = 0; i < totalDOF; i++) ctx.FlatMagn[i] = sigma_trial[i];
			}
			return result.iterations;
		}

		const bool serialMomentOps = (nMom <= 8);
		std::vector<int> momDof((size_t)nMom), momOff((size_t)nMom);
		for(int h = 0; h < nMom; h++)
		{
			const int elem = momElem[h];
			momDof[h] = IntrctPtr->GetElementDOF(elem);
			momOff[h] = IntrctPtr->GetElementDOFOffset(elem);
		}

		std::vector<double> RHS((size_t)totalDOF, 0.0);
		for(int h = 0; h < nMom; h++)
		{
			const int off = momOff[h];
			for(int t = 0; t < 3; t++) RHS[(size_t)off + t] = chiPerHex[h] * HextPerHex[(size_t)3 * h + t];
		}

		auto t_setup_start = std::chrono::high_resolution_clock::now();
		IntrctPtr->PrecomputeMomentAnyGeometry();
		std::vector<double> localLBlock((size_t)nMom * 36, 0.0);
		std::vector<double> diagKBlock((size_t)nMom * 36, 0.0);
		std::vector<double> zeroChi((size_t)nMom, 0.0);
		auto buildMomentLocal = [&](size_t hh) {
			const int h = (int)hh;
			IntrctPtr->MomentSystemBlockAny(h, h, zeroChi.data(), &localLBlock[(size_t)h * 36]);
			IntrctPtr->MomentSystemBlockAny(h, h, nullptr, &diagKBlock[(size_t)h * 36], true);
		};
		if(serialMomentOps) { for(int h = 0; h < nMom; h++) buildMomentLocal((size_t)h); }
		else { ngcore::ParallelFor(ngcore::IntRange(nMom), buildMomentLocal); }

		std::vector<double> blockInverse;
		std::vector<int> blockOffsets;

#ifdef HAVE_LAPACK
		// Moment rows strongly couple the local face DOF (dipole/monopole/quadrupole constraints).
		// A scalar diagonal Jacobi preconditioner is too weak even for small nonlinear cubes.
		// The production local choice is the natural per-element 5x5/6x6 moment block.
		blockOffsets.resize((size_t)nMom + 1);
		int total_block_storage = 0;
		for(int h = 0; h < nMom; h++)
		{
			int dof = momDof[h];
			blockOffsets[h] = total_block_storage;
			total_block_storage += dof * dof;
		}
		blockOffsets[nMom] = total_block_storage;
		blockInverse.assign((size_t)total_block_storage, 0.0);

		const int max_dof = 6;
		std::vector<double> block_copy((size_t)max_dof * max_dof);
		std::vector<int> ipiv(max_dof);
		std::vector<double> work((size_t)max_dof * max_dof);
		int lwork = max_dof * max_dof;
		for(int h = 0; h < nMom; h++)
		{
			int elem = momElem[h];
			int dof = momDof[h];
			int boff = blockOffsets[h];
			const double* Lh = &localLBlock[(size_t)h * 36];
			const double* Kh = &diagKBlock[(size_t)h * 36];

			for(int i = 0; i < dof; i++)
			{
				for(int j = 0; j < dof; j++)
				{
					// LAPACK wants column-major storage: block_copy[row + col*dof].
					block_copy[i + j * dof] = Lh[i * 6 + j] + chiPerHex[h] * Kh[i * 6 + j];
				}
			}

			int info = 0;
			dgetrf_(&dof, &dof, block_copy.data(), &dof, ipiv.data(), &info);
			if(info == 0) dgetri_(&dof, block_copy.data(), &dof, ipiv.data(), work.data(), &lwork, &info);
			if(info != 0)
			{
				fprintf(stderr, "[BiCG] Error: moment block-Jacobi failed at element %d (info=%d); scalar Jacobi fallback is disabled.\n",
				        elem, info);
				return -11;
			}

			for(int i = 0; i < dof; i++)
			{
				for(int j = 0; j < dof; j++)
				{
					blockInverse[boff + i * dof + j] = block_copy[i + j * dof];
				}
			}
		}
#else
		fprintf(stderr, "[BiCG] Error: moment block-Jacobi requires LAPACK; scalar Jacobi fallback is disabled.\n");
		return -11;
#endif
		// BUILD-ONCE dense geometry kernel for the variable-DOF (tet/wedge/mixed) moment path.
		// Mirrors the pure-hex momentKdense optimization (commit 1e5c9b50) and the classic MMM/MSC
		// cached-matrix matvec (MatVec_VariableDOF + FlatInteract, commit ceb8f9ea).  Previously the
		// matvec re-evaluated MomentSystemBlockAny for EVERY (h,g) pair on EVERY BiCGSTAB iteration
		// (~one full O(nMom^2) kernel build per matvec) -- the dominant cost that made tet/wedge
		// method-1 ~100x slower per iteration than the pure-hex path (regular tet 648 DOF: 161 iters
		// but ~58 ms/iter).  Build the chi-INDEPENDENT kernel block K[h][g] ONCE and scatter it into
		// a dense totalDOF x totalDOF row-major matrix (variable-DOF offsets via momOff/momDof); the
		// matvec then applies y = diag(chi)(K x) + L_local x as a cheap dense GEMV -- bit-identical to
		// the old per-iteration recompute, just computed once.  O(N^2) storage (= method 0 / hex);
		// for very large N use HDiv-VIM (loop-free).
		std::vector<double> momentKdenseAny((size_t)totalDOF * (size_t)totalDOF, 0.0);
		{
			auto buildKdenseRowAny = [&](size_t hh) {
				const int h = (int)hh;
				const int rdof = momDof[h];
				const int roff = momOff[h];
				double Kblk[36];
				for(int g = 0; g < nMom; g++)
				{
					const int cdof = momDof[g];
					const int coff = momOff[g];
					IntrctPtr->MomentSystemBlockAny(h, g, nullptr, Kblk, true);   // chi-independent geometry block (incl. images)
					for(int i = 0; i < rdof; i++)
						for(int j = 0; j < cdof; j++)
							momentKdenseAny[(size_t)(roff + i) * (size_t)totalDOF + (size_t)(coff + j)] = Kblk[i * 6 + j];
				}
			};
			if(serialMomentOps) { for(int h = 0; h < nMom; h++) buildKdenseRowAny((size_t)h); }
			else { ngcore::ParallelFor(ngcore::IntRange(nMom), buildKdenseRowAny); }
		}
		std::vector<double> momentKxAny((size_t)totalDOF, 0.0);

		auto t_setup_end = std::chrono::high_resolution_clock::now();
		rad.m_solve_t_moment_system_build += std::chrono::duration<double>(t_setup_end - t_setup_start).count();

		std::vector<double> sigma(totalDOF);
		for(int i = 0; i < totalDOF; i++) sigma[i] = ctx.FlatMagn[i];

		auto matvec = [&](const double* x, double* y) {
			// build-once: y = diag(chi)(Kdense x) + L_local x  -- cheap dense GEMV, NO kernel recompute
			cblas_dgemv(CblasRowMajor, CblasNoTrans, totalDOF, totalDOF, 1.0,
			            momentKdenseAny.data(), totalDOF, x, 1, 0.0, momentKxAny.data(), 1);
			auto applyRow = [&](size_t hh) {
				const int h = (int)hh;
				const int rdof = momDof[h];
				const int roff = momOff[h];
				const double chih = chiPerHex[h];
				const double* Lh = &localLBlock[(size_t)h * 36];
				for(int i = 0; i < rdof; i++)
				{
					double s = chih * momentKxAny[(size_t)roff + i];
					const double* Li = &Lh[(size_t)i * 6];
					for(int j = 0; j < rdof; j++) s += Li[j] * x[(size_t)roff + j];
					y[(size_t)roff + i] = s;
				}
			};
			if(serialMomentOps) { for(int h = 0; h < nMom; h++) applyRow((size_t)h); }
			else { ngcore::ParallelFor(ngcore::IntRange(nMom), applyRow); }
		};
		auto precond = [&](const double* x, double* y) {
			ngcore::ParallelFor(ngcore::IntRange(nMom), [&](size_t hh) {
				int h = (int)hh;
				int dof = momDof[h];
				int off = momOff[h];
				int boff = blockOffsets[h];
				for(int i = 0; i < dof; i++)
				{
					double sum = 0.0;
					for(int j = 0; j < dof; j++) sum += blockInverse[boff + i * dof + j] * x[off + j];
					y[off + i] = sum;
				}
			});
		};

		const double bicg_tol = rad.m_bicg_tol;
		const int bicg_max_iter = 10000;
		auto t_bicg_start = std::chrono::high_resolution_clock::now();
		radia::bicgstab::Result result =
			radia::bicgstab::Solve<double>(totalDOF, matvec, precond, RHS.data(), sigma.data(),
			                               bicg_tol, bicg_max_iter);
		auto t_bicg_end = std::chrono::high_resolution_clock::now();
		rad.m_solve_t_linear_solve += std::chrono::duration<double>(t_bicg_end - t_bicg_start).count();
		rad.m_solve_linear_iterations += result.iterations;

		if(!result.converged)
		{
			fprintf(stderr, "[BiCG] Warning: variable-DOF matrix-free moment BiCGSTAB did not reach tol %.3e after %d iterations (residual=%.4e)\n",
			        bicg_tol, result.iterations, result.residual);
			return -3;
		}

		std::vector<double> sigma_trial = sigma;
		double omega_ls = ApplyLineSearchDamping(ctx, this->IntrctPtr, sigma_trial);
		if(omega_ls >= 0.999)
		{
			for(int i = 0; i < totalDOF; i++) ctx.FlatMagn[i] = sigma_trial[i];
		}
		return result.iterations;
	}

	// Update poly->CurrentChi for 6DOF elements before BiCGSTAB
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		if(dof >= 4)
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

	// Newton: use chi_d in system matrix with RHS correction (start after 10 Picard iters)
	const int newton_start_iter_bicg = 10;
	bool newton_active = ctx.use_newton && iterCount >= newton_start_iter_bicg && !ctx.DifferentialChiArray.empty();

	auto t_bicg_start = std::chrono::high_resolution_clock::now();
	int n_iter;
	if(newton_active)
	{
		// Newton: pass chi_d for system matrix, chi_abs for RHS correction
		n_iter = SolveBiCGSTAB_VariableDOF(ctx, totalDOF, bicg_tol, bicg_max_iter, residual,
		                                    ctx.DifferentialChiArray, true,
		                                    &ctx.CurrentChiArray, ctx.OldSigma.data());
	}
	else
	{
		n_iter = SolveBiCGSTAB_VariableDOF(ctx, totalDOF, bicg_tol, bicg_max_iter, residual,
		                                    ctx.CurrentChiArray);
	}
	auto t_bicg_end = std::chrono::high_resolution_clock::now();
	rad.m_solve_t_linear_solve += std::chrono::duration<double>(t_bicg_end - t_bicg_start).count();
	rad.m_solve_linear_iterations += n_iter;

	return n_iter;
}

//=========================================================================
// Method 2: BiCGSTAB with H-matrix (HACApK ACA+)
//=========================================================================


double radTRelaxationMethNo_2::Dot(const std::vector<double>& a, const std::vector<double>& b, int n)
{
	double sum = 0.0;
	ngcore::ParallelForRange(ngcore::IntRange(n), [&](ngcore::IntRange r) {
		double local_sum = 0.0;
		for (auto i : r) {
			local_sum += a[i] * b[i];
		}
		ngcore::AtomicAdd(sum, local_sum);
	});
	return sum;
}

double radTRelaxationMethNo_2::Norm2(const std::vector<double>& a, int n)
{
	return std::sqrt(Dot(a, a, n));
}

void radTRelaxationMethNo_2::Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n)
{
	ngcore::ParallelFor(ngcore::IntRange(n), [&](size_t i) {
		y[i] += alpha * x[i];
	});
}

void radTRelaxationMethNo_2::Copy(const std::vector<double>& src, std::vector<double>& dst, int n)
{
	ngcore::ParallelFor(ngcore::IntRange(n), [&](size_t i) {
		dst[i] = src[i];
	});
}

void radTRelaxationMethNo_2::Scale(double alpha, std::vector<double>& x, int n)
{
	ngcore::ParallelFor(ngcore::IntRange(n), [&](size_t i) {
		x[i] *= alpha;
	});
}

int radTRelaxationMethNo_2::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	// HACApK supports 3-DOF dipole and 4-6 DOF face-charge MSC elements.
	return AutoRelax_VariableDOF(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
}

//-------------------------------------------------------------------------
// SolveBiCGSTAB_HMatrix_VariableDOF
// BiCGSTAB with H-matrix for 3-DOF dipole and 4-6 DOF face-charge MSC elements
//-------------------------------------------------------------------------

int radTRelaxationMethNo_2::SolveBiCGSTAB_HMatrix_VariableDOF(NonlinearContext& ctx,
                                                               int totalDOF, double tol, int max_iter, double& residual,
                                                               const std::vector<double>& elemChiArray,
                                                               bool use_newton,
                                                               const std::vector<double>* absChiArray,
                                                               const double* oldSigma)
{
	if (!m_hacapk || !m_hacapk->IsValid()) return 0;

	// TaskManager self-wrap (AGENTS.md "Parallelization: NGSolve TaskManager"): one region around the
	// whole MMM/MSC method-2 BiCGSTAB (init + loop matvecs) so it is multi-threaded even when
	// driven by a bare rad.Solve(...,2) without `with TaskManager()`.  Nested -> no-op.
	ngcore::RegionTaskManager rtm(radia::GetMaxThreads());

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	// Allocate work vectors
	std::vector<double> r(totalDOF), r0(totalDOF), p(totalDOF), v(totalDOF), s(totalDOF), t(totalDOF);
	std::vector<double> p_hat(totalDOF), s_hat(totalDOF);
	std::vector<double> inv_chi(totalDOF);
	std::vector<double> rhs(totalDOF);
	std::vector<double> sol(totalDOF);

	double* FlatMagn = IntrctPtr->GetFlatMagnArray();
	double* FlatField = IntrctPtr->GetFlatFieldArray();
	double* FlatExtern = IntrctPtr->GetFlatExternFieldArray();

	if(FlatMagn == nullptr || FlatField == nullptr || FlatExtern == nullptr) return 0;

	// Pre-compute 1/chi and RHS for all elements
	// Picard: (D(1/chi_abs) + G) sigma = H_ext
	// Newton: (D(1/chi_d) + G) sigma_new = H_ext + D(1/chi_d - 1/chi_abs) * sigma_old
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		if(dof != 3 && dof != 4 && dof != 5 && dof != 6)
		{
			std::cerr << "[HACApK] Error: Element " << elem << " has " << dof
			          << " DOF, expected 3 (dipole) or 4-6 (face-charge MSC)" << std::endl;
			return 0;
		}

		// elemChiArray contains chi_d (Newton) or chi_abs (Picard) for system matrix
		double chi_matrix = elemChiArray[elem];
		if(chi_matrix < 1.0e-6) chi_matrix = 1.0e-6;
		double inv_chi_val = 1.0 / chi_matrix;

		// Newton correction: (1/chi_d - 1/chi_abs) per element
		double newton_correction = 0.0;
		if(use_newton && absChiArray && oldSigma)
		{
			double chi_abs = (*absChiArray)[elem];
			if(chi_abs < 1.0e-6) chi_abs = 1.0e-6;
			newton_correction = inv_chi_val - 1.0 / chi_abs;
		}

		for(int k = 0; k < dof; k++)
		{
			inv_chi[offset + k] = inv_chi_val;
			rhs[offset + k] = FlatExtern[offset + k];
			if(use_newton && oldSigma)
			{
				rhs[offset + k] += newton_correction * oldSigma[offset + k];
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

	// Build preconditioner
	std::vector<double> hmat_blockInverse;
	std::vector<int> hmat_blockOffsets;

#ifdef HAVE_LAPACK
	// Element-wise block Jacobi is the minimum meaningful local preconditioner for mixed DOF systems.
	if(!BuildBlockJacobiPreconditioner_HMatrix(hmat_blockInverse, hmat_blockOffsets, inv_chi, totalDOF))
	{
		fprintf(stderr, "[HACApK BiCG] Error: block-Jacobi build failed; scalar Jacobi fallback is disabled.\n");
		return -11;
	}
#else
	fprintf(stderr, "[HACApK BiCG] Error: block-Jacobi requires LAPACK; scalar Jacobi fallback is disabled.\n");
	return -11;
#endif

	// Initialize: r0 = b - A*x0
	m_hacapk->MatVec(sol, v);  // v = A*x0 using H-matrix
	this->Copy(rhs, r, totalDOF);
	this->Axpy(-1.0, v, r, totalDOF);
	this->Copy(r, r0, totalDOF);

	double rho = 1.0, alpha_bicg = 1.0, omega = 1.0;
	std::fill(p.begin(), p.end(), 0.0);
	std::fill(v.begin(), v.end(), 0.0);

	double rhs_norm = this->Norm2(rhs, totalDOF);
	if(rhs_norm < 1.0e-30) rhs_norm = 1.0;

	int iter;
	for(iter = 1; iter <= max_iter; iter++)
	{
		double rho_old = rho;
		rho = this->Dot(r0, r, totalDOF);

		if(std::abs(rho) < 1.0e-30)
		{
			residual = this->Norm2(r, totalDOF) / rhs_norm;
			break;
		}

		if(iter == 1)
		{
			this->Copy(r, p, totalDOF);
		}
		else
		{
			if(std::abs(rho_old * omega) < 1.0e-30)
			{
				residual = this->Norm2(r, totalDOF) / rhs_norm;
				break;
			}
			double beta = (rho / rho_old) * (alpha_bicg / omega);
			this->Axpy(-omega, v, p, totalDOF);
			this->Scale(beta, p, totalDOF);
			this->Axpy(1.0, r, p, totalDOF);
		}

		// Apply preconditioner
		this->ApplyBlockJacobiPreconditioner_HMatrix(p, p_hat, hmat_blockInverse, hmat_blockOffsets);

		// v = A * p_hat using H-matrix
		m_hacapk->MatVec(p_hat, v);

		double r0_dot_v = this->Dot(r0, v, totalDOF);
		if(std::abs(r0_dot_v) < 1.0e-30)
		{
			residual = this->Norm2(r, totalDOF) / rhs_norm;
			break;
		}
		alpha_bicg = rho / r0_dot_v;

		this->Copy(r, s, totalDOF);
		this->Axpy(-alpha_bicg, v, s, totalDOF);

		double s_norm = this->Norm2(s, totalDOF);
		if(s_norm / rhs_norm < tol)
		{
			this->Axpy(alpha_bicg, p_hat, sol, totalDOF);
			residual = s_norm / rhs_norm;
			break;
		}

		ApplyBlockJacobiPreconditioner_HMatrix(s, s_hat, hmat_blockInverse, hmat_blockOffsets);

		// t = A * s_hat using H-matrix
		m_hacapk->MatVec(s_hat, t);

		double t_dot_s = this->Dot(t, s, totalDOF);
		double t_dot_t = this->Dot(t, t, totalDOF);
		if(std::abs(t_dot_t) < 1.0e-30)
		{
			this->Axpy(alpha_bicg, p_hat, sol, totalDOF);
			residual = s_norm / rhs_norm;
			break;
		}
		omega = t_dot_s / t_dot_t;

		this->Axpy(alpha_bicg, p_hat, sol, totalDOF);
		this->Axpy(omega, s_hat, sol, totalDOF);

		this->Copy(s, r, totalDOF);
		this->Axpy(-omega, t, r, totalDOF);

		double r_norm = this->Norm2(r, totalDOF);
		residual = r_norm / rhs_norm;
		if(residual < tol) break;

		if(std::abs(omega) < 1.0e-30) break;
	}

	// Apply line search damping if Newton is active
	std::vector<double> sigma_trial = sol;
	double omega_ls = ApplyLineSearchDamping(ctx, this->IntrctPtr, sigma_trial);

	// If line search already updated FlatMagn (omega_ls < 0.999), we're done
	// Otherwise copy trial solution (omega_ls=1.0 case, full step)
	if(omega_ls >= 0.999)
	{
		for(int i = 0; i < totalDOF; i++)
		{
			FlatMagn[i] = sigma_trial[i];
		}
	}
	// else: ApplyLineSearchDamping already updated FlatMagn with damped solution

	return iter;
}

//-------------------------------------------------------------------------
// GetDiagonalElements_HMatrix_VariableDOF
//-------------------------------------------------------------------------

void radTRelaxationMethNo_2::GetDiagonalElements_HMatrix_VariableDOF(std::vector<double>& diag,
                                                                      const std::vector<double>& inv_chi,
                                                                      int totalDOF)
{
	// Get diagonal elements A_ii = -K_ii/(4pi) + 1/chi_i (physically correct)
	// GetDiagonalN() returns raw K/(4pi), so we negate it
	// Use cached N_ii values for efficiency (computed once during H-matrix build)
	if(m_hacapk->IsDiagonalCached())
	{
		const std::vector<double>& diag_N = m_hacapk->GetDiagonalN();
		ngcore::ParallelFor(ngcore::IntRange(totalDOF), [&](size_t i) {
			diag[i] = -diag_N[i] + inv_chi[i];  // Physical: -K/(4pi) + 1/chi
		});
	}
	else
	{
		// Fallback: compute on-demand (slow)
		for(int i = 0; i < totalDOF; i++)
		{
			double N_ii = m_hacapk->GetInteractionMatrixElement(i, i);
			diag[i] = -N_ii + inv_chi[i];  // Physical: -K/(4pi) + 1/chi
		}
	}
}

//-------------------------------------------------------------------------
// Element-block Jacobi preconditioner for H-matrix BiCGSTAB.
// Extracts each element's natural diagonal block DOF-generically through GetInteractionMatrixElement.
//-------------------------------------------------------------------------

#ifdef HAVE_LAPACK
bool radTRelaxationMethNo_2::BuildBlockJacobiPreconditioner_HMatrix(
	std::vector<double>& blockInverse, std::vector<int>& blockOffsets,
	const std::vector<double>& inv_chi, int totalDOF)
{
	if(!m_hacapk || !IntrctPtr) return false;

	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	// Calculate total storage needed for block inverses
	int total_block_storage = 0;
	blockOffsets.resize(AmOfMainElem + 1);
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		blockOffsets[elem] = total_block_storage;
		total_block_storage += dof * dof;
	}
	blockOffsets[AmOfMainElem] = total_block_storage;
	blockInverse.resize(total_block_storage);

	int max_dof = 6;
	std::vector<double> block_copy(max_dof * max_dof);
	std::vector<int> ipiv(max_dof);
	std::vector<double> work(max_dof * max_dof);
	int lwork = max_dof * max_dof;

	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int mat_offset = IntrctPtr->GetElementDOFOffset(elem);
		int block_offset = blockOffsets[elem];

		// Diagonal block A_block = -N_block + (1/chi) I, extracted DOF-generically from the H-matrix
		// kernel via GetInteractionMatrixElement (returns +N).  Works for any element DOF; the multipole-moment MMM
		// MSC path never reaches method 2 here (pure hex/wedge/pyramid is rerouted to the LU/Picard moment driver
		// in SolveGen), so this block-Jacobi only ever sees 3-DOF tetrahedra (MMM).
		for(int i = 0; i < dof; i++)
		{
			for(int j = 0; j < dof; j++)
			{
				// column-major for LAPACK
				block_copy[i + j * dof] = -m_hacapk->GetInteractionMatrixElement(mat_offset + i, mat_offset + j);
				if(i == j)
				{
					block_copy[i + j * dof] += inv_chi[mat_offset + i];
				}
			}
		}

		// Invert with LAPACK
		int info = 0;
		dgetrf_(&dof, &dof, block_copy.data(), &dof, ipiv.data(), &info);
		if(info != 0)
		{
			fprintf(stderr, "[HACApK Block Jacobi] Element %d: singular diagonal block (info=%d)\n", elem, info);
			return false;
		}
		dgetri_(&dof, block_copy.data(), &dof, ipiv.data(), work.data(), &lwork, &info);
		if(info != 0)
		{
			fprintf(stderr, "[HACApK Block Jacobi] Element %d: inversion failed (info=%d)\n", elem, info);
			return false;
		}

		// Store inverse in row-major format
		for(int i = 0; i < dof; i++)
		{
			for(int j = 0; j < dof; j++)
			{
				blockInverse[block_offset + i * dof + j] = block_copy[i + j * dof];
			}
		}
	}

	return true;
}
#endif

void radTRelaxationMethNo_2::ApplyBlockJacobiPreconditioner_HMatrix(
	const std::vector<double>& x, std::vector<double>& y,
	const std::vector<double>& blockInverse, const std::vector<int>& blockOffsets)
{
	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	ngcore::ParallelFor(ngcore::IntRange(AmOfMainElem), [&](size_t elem) {
		int dof = IntrctPtr->GetElementDOF(elem);
		int mat_offset = IntrctPtr->GetElementDOFOffset(elem);
		int block_offset = blockOffsets[elem];

		for(int i = 0; i < dof; i++)
		{
			double sum = 0.0;
			for(int j = 0; j < dof; j++)
			{
				sum += blockInverse[block_offset + i * dof + j] * x[mat_offset + j];
			}
			y[mat_offset + i] = sum;
		}
	});
}

//-------------------------------------------------------------------------
// AutoRelax_VariableDOF (H-matrix version)
//-------------------------------------------------------------------------

int radTRelaxationMethNo_2::AutoRelax_VariableDOF(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	int AmOfMainElem = IntrctPtr->AmOfMainElem;
	if(AmOfMainElem <= 0) return 0;

	// HACApK supports 3-DOF dipole and 4-6 DOF face-charge MSC elements.
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
		m_hacapk = new RadHACApKMMMManager(IntrctPtr);
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

		// NOTE: The element-block Jacobi preconditioner is rebuilt each nonlinear iteration
		// because the diagonal blocks depend on the current per-element susceptibility.
	}


	// NOTE: 3DOF tetrahedra use precomputed geometry for fast O(1) access.
	// Surface-charge MMM multipole-moment solves use the dense LU / matrix-free moment path (no HACApK).

	std::vector<double> OldMagn(totalDOF);
	// Store current isotropic chi for ALL elements (unified 3DOF/6DOF handling, same as LU/BiCGSTAB)
	std::vector<double> CurrentChiArray_hacapk(AmOfMainElem, 1.0);
	// Newton-Raphson: differential chi (chi_d = (dB/dH)/mu_0 - 1) and old sigma
	std::vector<double> DifferentialChiArray(AmOfMainElem, 1.0);
	std::vector<double> OldSigma(totalDOF, 0.0);
	bool use_newton = rad.m_use_newton;
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
		else if(dof >= 4)
		{
			// Face-charge MSC elements: initialize FlatField and estimate H from external field
			for(int k = 0; k < dof; k++)
			{
				FlatField[offset + k] = FlatExtern[offset + k];
			}
			// For face-charge MSC, estimate H from external field
			// The external field in RHS is H_ext (uniform applied field)
			// Estimate H direction from face normals weighted by sigma
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
			if(poly && poly->Use6DOF_MSC)
			{
				// Estimate H from sigma pattern: H ~ sum(sigma_i * n_i) / sum(n_i dot n_i)
				double Hx = 0.0, Hy = 0.0, Hz = 0.0;
				double wx = 0.0, wy = 0.0, wz = 0.0;
				for(int face = 0; face < dof; face++)
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
			// Initialize NewFieldArray with H_init_mag in z-direction when no element direction is available.
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
	// FIX (2025-12-26): Initialize for BOTH 3-DOF dipole and face-charge MSC elements
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

		// For face-charge MSC elements, also store in poly->CurrentChi
		if(dof >= 4)
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
	double max_B_rel_change = 1.0e30;  // Initialize for first iteration

	// Outer nonlinear iteration (rewritten to match LU solver structure)
	for(outerIter = 0; outerIter < MaxIterNumber; outerIter++)
	{

		// Store old values
		for(int i = 0; i < totalDOF; i++)
		{
			OldMagn[i] = FlatMagn[i];
			OldSigma[i] = FlatMagn[i];  // Store sigma_old for Newton RHS correction
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
			else if(dof >= 4)
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

		// Build temporary NonlinearContext for line search
		NonlinearContext ctx_temp;
		ctx_temp.totalDOF = totalDOF;
		ctx_temp.AmOfMainElem = AmOfMainElem;
		ctx_temp.FlatMagn = FlatMagn;
		ctx_temp.FlatField = FlatField;  // Use FlatField (double*) not NewFieldArray (TVector3d*)
		ctx_temp.FlatExtern = IntrctPtr->GetFlatExternFieldArray();
		ctx_temp.OldSigma = OldSigma;
		ctx_temp.OldBnorm = OldBnorm;
		ctx_temp.CurrentChiArray = CurrentChiArray_hacapk;
		ctx_temp.DifferentialChiArray = DifferentialChiArray;
		ctx_temp.polyCache = polyCache;
		ctx_temp.use_newton = use_newton;
		ctx_temp.newton_damping_enabled = rad.m_newton_damping_enabled;
		ctx_temp.newton_ls_max_iter = rad.m_newton_ls_max_iter;
		ctx_temp.newton_ls_min_omega = rad.m_newton_ls_min_omega;
		ctx_temp.max_B_rel_change = max_B_rel_change;  // Use previous iteration's value
		ctx_temp.B_sat = 1.0;  // Default, will be overridden by material

		// Solve with BiCGSTAB using H-matrix
		// NOTE: max_iter for BiCGSTAB inner loop should be FIXED (not user's MaxIterNumber)
		// User's MaxIterNumber controls OUTER nonlinear iterations, not inner BiCGSTAB
		// bicg_tol is set via rad.SetBiCGSTABTol() Python API (default: 1e-4, ELF-compatible)
		double residual = 0.0;
		const double bicg_tol = rad.m_bicg_tol;
		const int bicg_max_iter = 10000;  // Inner BiCGSTAB max iterations (fixed)

		// Time BiCGSTAB solve
		// Newton-Raphson: use chi_d for system matrix with RHS correction
		// Start with Picard for first 10 iterations to approach solution, then switch to Newton
		const int newton_start_iter = 10;
		auto t_bicg_start = std::chrono::high_resolution_clock::now();
		int n_iter;
		if(use_newton && outerIter >= newton_start_iter)
		{
			// Newton: system matrix uses chi_d, RHS adds correction with chi_abs
			n_iter = SolveBiCGSTAB_HMatrix_VariableDOF(ctx_temp, totalDOF, bicg_tol, bicg_max_iter, residual,
			                                            DifferentialChiArray, true,
			                                            &CurrentChiArray_hacapk, OldSigma.data());
		}
		else
		{
			// Picard: use chi_abs
			n_iter = SolveBiCGSTAB_HMatrix_VariableDOF(ctx_temp, totalDOF, bicg_tol, bicg_max_iter, residual,
			                                            CurrentChiArray_hacapk);
		}
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
			else if(dof >= 4)
			{
				radTPolyhedron* poly = polyCache[elem];
				if(poly && poly->Use6DOF_MSC)
				{
					// Store sigma values (5 for wedge, 6 for hex)
					for(int k = 0; k < dof; k++)
					{
						poly->Sigma[k] = FlatMagn[offset + k];
					}

					// Compute effective magnetization from sigma (same as LU solver)
					double Mx = 0.0, My = 0.0, Mz = 0.0;
					double wx = 0.0, wy = 0.0, wz = 0.0;
					for(int face = 0; face < dof; face++)
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
		max_B_rel_change = 0.0;  // Reset for this iteration (declared outside loop)
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
					// Newton: also compute differential chi
					if(use_newton)
					{
						DifferentialChiArray[elem] = NonlinMater->ComputeDifferentialChi(H_mag);
					}
				}
				else
				{
					// Fallback for linear materials
					TMatrix3d KsiTensor;
					TVector3d MrVect;
					MaterPtr->DefineInstantKsiTensor(H_new, KsiTensor, MrVect);
					chi_new = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
					if(chi_new < 1.0e-6) chi_new = 1.0e-6;
					if(relax_hacapk > 0.0 && relax_hacapk <= 1.0)
					{
						chi_new = chi_new * (1.0 - relax_hacapk) + chi_matrix * relax_hacapk;
					}
					if(use_newton) DifferentialChiArray[elem] = chi_new;
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
			else if(dof >= 4)
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
						chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old, relax_hacapk);
						if(use_newton)
						{
							DifferentialChiArray[elem] = NonlinMater->ComputeDifferentialChi(H_mag);
						}
					}
					else
					{
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_new, KsiTensor, MrVect);
						chi_new = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						if(chi_new < 1.0e-6) chi_new = 1.0e-6;
						if(relax_hacapk > 0.0 && relax_hacapk <= 1.0)
						{
							chi_new = chi_new * (1.0 - relax_hacapk) + chi_matrix * relax_hacapk;
						}
						if(use_newton) DifferentialChiArray[elem] = chi_new;
					}
					poly->CurrentChi = chi_new;
					CurrentChiArray_hacapk[elem] = chi_new;

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

