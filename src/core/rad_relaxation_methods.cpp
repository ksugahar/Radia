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

#include <time.h>
#include <chrono>   // For timing instrumentation
#include <cstring>  // For std::memcpy
#include <cstdio>   // For fprintf in debug logging
#include <cstdlib>  // For getenv
#include <array>    // For std::array in IMA mirror computation

// Uncomment to enable chi value debugging
// #define RADIA_DEBUG_CHI

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
		if(dof >= 5)
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
		else if(dof >= 5)
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

		if(dof == 3)
		{
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			TVector3d M = g3dRelaxPtr->Magn;
			TVector3d H(ctx.FlatField[offset], ctx.FlatField[offset+1], ctx.FlatField[offset+2]);
			TVector3d B(MU_0 * (H.x + M.x), MU_0 * (H.y + M.y), MU_0 * (H.z + M.z));
			ctx.OldBnorm[elem] = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
		}
		else if(dof >= 5)
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
			// 3DOF elements (tetrahedra, wedges): FlatMagn contains M directly
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			g3dRelaxPtr->Magn.x = ctx.FlatMagn[offset + 0];
			g3dRelaxPtr->Magn.y = ctx.FlatMagn[offset + 1];
			g3dRelaxPtr->Magn.z = ctx.FlatMagn[offset + 2];
		}
		else if(dof >= 5)
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
		else if(dof_j >= 5)
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

#ifdef RADIA_DEBUG_CHI
				TVector3d H_new = IntrctPtr->NewFieldArray[elem_j];
				double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);
				fprintf(stderr, "ComputeActualH elem %d (6DOF): H = M/chi = [%.1f, %.1f, %.1f], |H|=%.1f, chi=%.2f\n",
				        elem_j, H_new.x, H_new.y, H_new.z, H_mag, chi);
#endif
			}
		}
	}

	// Third pass: Add IMA mirror contributions to H field if IMA is enabled
	// NOTE (2026-02-05): This fix attempt made results WORSE (74% error vs 22% before).
	// The IMA demagnetizing contribution is being computed but has wrong effect.
	// DISABLED pending further investigation.
	// TODO: Investigate why adding mirror demagnetizing makes things worse.
	if(false && IntrctPtr->IsIMAEnabled() && IntrctPtr->NewFieldArray != nullptr)
	{
		int imaSym = IntrctPtr->GetIMASymmetry();
		int signX = IntrctPtr->GetIMASignX();
		int signY = IntrctPtr->GetIMASignY();
		int signZ = IntrctPtr->GetIMASignZ();

		// Compute IMA mirror contributions for each element
		for(int elem_i = 0; elem_i < ctx.AmOfMainElem; elem_i++)
		{
			radTPolyhedron* poly_i = ctx.polyCache[elem_i];
			if(!poly_i || !poly_i->Use6DOF_MSC) continue;

			// Get observation point (element i's center)
			TVector3d obs = poly_i->CentrPoint;

			// Accumulate field from all elements' IMA mirrors
			TVector3d H_ima(0.0, 0.0, 0.0);

			for(int elem_j = 0; elem_j < ctx.AmOfMainElem; elem_j++)
			{
				radTPolyhedron* poly_j = ctx.polyCache[elem_j];
				if(!poly_j || !poly_j->Use6DOF_MSC) continue;

				// Get source element's face data
				int nFaces = poly_j->AmOfFaces;
				if(nFaces != 6) continue;  // Only hexahedra

				// Get face vertices from VectHandlePgnAndTrans
				// Transform 2D polygon vertices to 3D global coordinates
				std::array<std::array<TVector3d, 4>, 8> faceVertices;
				int faceIdx = 0;
				for(auto& hpt : poly_j->VectHandlePgnAndTrans)
				{
					if(faceIdx >= 8) break;
					radTPolygon* pgn = hpt.PgnHndl.rep;
					radTrans* tr = hpt.TransHndl.rep;
					if(!pgn || pgn->AmOfEdgePoints < 4) {
						faceIdx++;
						continue;
					}
					// Transform 2D local (x, y, CoordZ) to 3D global via polygon's transform
					const radTVect2dVect& verts2d = pgn->EdgePointsVector;
					faceVertices[faceIdx][0] = tr->TrPoint(TVector3d(verts2d[0].x, verts2d[0].y, pgn->CoordZ));
					faceVertices[faceIdx][1] = tr->TrPoint(TVector3d(verts2d[1].x, verts2d[1].y, pgn->CoordZ));
					faceVertices[faceIdx][2] = tr->TrPoint(TVector3d(verts2d[2].x, verts2d[2].y, pgn->CoordZ));
					faceVertices[faceIdx][3] = tr->TrPoint(TVector3d(verts2d[3].x, verts2d[3].y, pgn->CoordZ));
					faceIdx++;
				}

				// Lambda to compute field from a mirrored element
				auto computeMirroredFieldContribution = [&](int mirrorAxis, int sign) -> TVector3d {
					TVector3d H_mirror(0.0, 0.0, 0.0);

					// Create mirrored face vertices
					std::array<std::array<TVector3d, 4>, 8> mirrorVerts;
					for(int i = 0; i < nFaces; i++)
					{
						for(int j = 0; j < 4; j++)
						{
							mirrorVerts[i][j] = faceVertices[i][j];
							if(mirrorAxis & radTInteraction::IMA_X) mirrorVerts[i][j].x = -mirrorVerts[i][j].x;
							if(mirrorAxis & radTInteraction::IMA_Y) mirrorVerts[i][j].y = -mirrorVerts[i][j].y;
							if(mirrorAxis & radTInteraction::IMA_Z) mirrorVerts[i][j].z = -mirrorVerts[i][j].z;
						}
					}

					// Count number of mirror axes (odd number = need winding reversal)
					int numAxes = 0;
					if(mirrorAxis & radTInteraction::IMA_X) numAxes++;
					if(mirrorAxis & radTInteraction::IMA_Y) numAxes++;
					if(mirrorAxis & radTInteraction::IMA_Z) numAxes++;
					bool reverseWinding = (numAxes % 2 == 1);

					// Compute mirrored center point
					TVector3d mirrorCenter = poly_j->CentrPoint;
					if(mirrorAxis & radTInteraction::IMA_X) mirrorCenter.x = -mirrorCenter.x;
					if(mirrorAxis & radTInteraction::IMA_Y) mirrorCenter.y = -mirrorCenter.y;
					if(mirrorAxis & radTInteraction::IMA_Z) mirrorCenter.z = -mirrorCenter.z;

					// Compute field from mirrored faces using solid angle integration
					for(int i = 0; i < nFaces; i++)
					{
						const TVector3d& MV0 = mirrorVerts[i][0];
						const TVector3d& MV1 = mirrorVerts[i][1];
						const TVector3d& MV2 = mirrorVerts[i][2];
						const TVector3d& MV3 = mirrorVerts[i][3];

						// Mirror sigma = sign * original sigma
						double mirrorSigma = sign * poly_j->Sigma[i];

						// Use FieldFromQuadFaceMirrored for accurate solid angle integration
						TVector3d H_face = poly_j->FieldFromQuadFaceMirrored(
							obs, MV0, MV1, MV2, MV3, mirrorSigma, reverseWinding, mirrorCenter);

						H_mirror.x += H_face.x;
						H_mirror.y += H_face.y;
						H_mirror.z += H_face.z;
					}

					// Apply 1/(4*pi) factor
					H_mirror.x *= RadConst::INV_FOUR_PI;
					H_mirror.y *= RadConst::INV_FOUR_PI;
					H_mirror.z *= RadConst::INV_FOUR_PI;

					return H_mirror;
				};

				// Single axis contributions
				if(imaSym & radTInteraction::IMA_X)
					H_ima += computeMirroredFieldContribution(radTInteraction::IMA_X, signX);
				if(imaSym & radTInteraction::IMA_Y)
					H_ima += computeMirroredFieldContribution(radTInteraction::IMA_Y, signY);
				if(imaSym & radTInteraction::IMA_Z)
					H_ima += computeMirroredFieldContribution(radTInteraction::IMA_Z, signZ);

				// Dual axis contributions
				if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Y))
					H_ima += computeMirroredFieldContribution(radTInteraction::IMA_XY, signX * signY);
				if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Z))
					H_ima += computeMirroredFieldContribution(radTInteraction::IMA_XZ, signX * signZ);
				if((imaSym & radTInteraction::IMA_Y) && (imaSym & radTInteraction::IMA_Z))
					H_ima += computeMirroredFieldContribution(radTInteraction::IMA_YZ, signY * signZ);

				// Triple axis contribution
				if((imaSym & radTInteraction::IMA_X) && (imaSym & radTInteraction::IMA_Y) && (imaSym & radTInteraction::IMA_Z))
					H_ima += computeMirroredFieldContribution(radTInteraction::IMA_XYZ, signX * signY * signZ);
			}

			// Add IMA contributions to H field
			IntrctPtr->NewFieldArray[elem_i].x += H_ima.x;
			IntrctPtr->NewFieldArray[elem_i].y += H_ima.y;
			IntrctPtr->NewFieldArray[elem_i].z += H_ima.z;

#ifdef RADIA_DEBUG_CHI
			if(elem_i < 3) {  // Only print first few elements
				fprintf(stderr, "ComputeActualH elem %d: H_ima = [%.1f, %.1f, %.1f]\n",
				        elem_i, H_ima.x, H_ima.y, H_ima.z);
			}
#endif
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
		else if(dof >= 5)
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

		if(dof >= 5)
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
		else if(dof >= 5)
		{
			radTPolyhedron* poly = ctx.polyCache[elem];
			if(poly && poly->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
			{
				TVector3d H_new = IntrctPtr->NewFieldArray[elem];
				double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);
				TVector3d M_poly = poly->Magn;
				double M_mag = std::sqrt(M_poly.x*M_poly.x + M_poly.y*M_poly.y + M_poly.z*M_poly.z);

#ifdef RADIA_DEBUG_CHI
				fprintf(stderr, "  Element %d (6DOF): H_mag = %.2f, M_mag = %.2f, M/H = %.2f, chi_old = %.2f\n",
				        elem, H_mag, M_mag, H_mag > 1e-6 ? M_mag/H_mag : 0.0, mu_old - 1.0);
#endif

				if(NonlinMater != nullptr)
				{
					chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old, ctx.relax_param);
#ifdef RADIA_DEBUG_CHI
					fprintf(stderr, "    -> chi_new from B-H = %.2f\n", chi_new);
#endif
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
	if(!BuildBaseMatrix(ctx, IntrctPtr))
		return 0;  // Memory allocation failed

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

		// Update element magnetization and compute actual H field from sigma
		// Uses H = H_ext + H_demag (ELF-compatible) instead of circular H = M/chi
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
		#pragma omp parallel for if(n_elem > 50)
		for(int i = 0; i < n_elem; i++)
		{
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
		}
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
		#pragma omp parallel for if(n_elem > 50)
		for(int i = 0; i < n_elem; i++)
		{
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
		}
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
	// CRITICAL: Use size_t to avoid int32 overflow for DOF > 46340
	size_t matrix_size = (size_t)totalDOF * (size_t)totalDOF;
	std::vector<double> SystemMatrix;
	try {
		SystemMatrix.resize(matrix_size);
	} catch (const std::bad_alloc&) {
		// Memory allocation failed - likely DOF is too large for LU
		double required_gb = (double)matrix_size * 8 / (1024.0 * 1024.0 * 1024.0);
		fprintf(stderr, "Radia::Solve> LU solver requires %.1f GB memory for DOF=%d. Use BiCGSTAB (method 1) or HACApK (method 2) for large problems.\n", required_gb, totalDOF);
		return -2;  // Memory allocation failure
	}
	// Copy base matrix (already contains -K/(4pi), see SetupBaseMatrix_VariableDOF)
	// System equation: (-K/(4pi) + I/chi) * sigma = H_ext_n (ELF-compatible)
	// BaseMatrix is already negated in SetupBaseMatrix_VariableDOF (line 282-285)
	std::memcpy(SystemMatrix.data(), ctx.BaseMatrix.data(), matrix_size * sizeof(double));

	// Build RHS vector (will be overwritten with solution by dgesv)
	std::vector<double> RHS(totalDOF);

	// Update diagonal and RHS based on current chi
	// Matrix is COLUMN-MAJOR: A(i,j) at [j * totalDOF + i]
	// Diagonal element A(k,k) is at [k * totalDOF + k] = [k * (totalDOF + 1)]
#ifdef RADIA_DEBUG_CHI
	fprintf(stderr, "=== LU Solver Debug: chi and matrix values ===\n");
#endif
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);

		double chi = ctx.CurrentChiArray[elem];
		if(chi < 1.0e-6) chi = 1.0e-6;
		double inv_chi = 1.0 / chi;

#ifdef RADIA_DEBUG_CHI
		fprintf(stderr, "Element %d: chi = %.6f, inv_chi = %.6e, dof = %d\n", elem, chi, inv_chi, dof);
#endif

		// Add 1/chi to diagonal and set RHS = H_ext_n
		// MSC equation: (-K/(4pi) + 1/chi * I) * sigma = H_ext_n (ELF-compatible)
		for(int k = 0; k < dof; k++)
		{
			int row = offset + k;
#ifdef RADIA_DEBUG_CHI
			double diag_before = SystemMatrix[row * (totalDOF + 1)];
#endif
			SystemMatrix[row * (totalDOF + 1)] += inv_chi;
#ifdef RADIA_DEBUG_CHI
			double diag_after = SystemMatrix[row * (totalDOF + 1)];
			fprintf(stderr, "  row %d: diag_before = %.6e, diag_after = %.6e, RHS = %.6e\n",
			        row, diag_before, diag_after, ctx.FlatExtern[row]);
#endif
			RHS[row] = ctx.FlatExtern[row];
		}

		// Update poly->CurrentChi for 6DOF elements
		if(dof >= 5)
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

	// CRITICAL: dgesv expects COLUMN-MAJOR format, but BaseMatrix is ROW-MAJOR
	// Transpose in-place: swap A[i,j] with A[j,i] for i < j
	for(int i = 0; i < totalDOF; i++)
	{
		for(int j = i + 1; j < totalDOF; j++)
		{
			// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
			std::swap(SystemMatrix[(size_t)i * totalDOF + j], SystemMatrix[(size_t)j * totalDOF + i]);
		}
	}

	// dgesv overwrites SystemMatrix with LU factors and RHS with solution
	dgesv_(&totalDOF, &nrhs, SystemMatrix.data(), &totalDOF, ipiv.data(), RHS.data(), &totalDOF, &info);

	if(info != 0) return -1;  // Singular matrix
#else
	// Fallback: transpose to row-major and use SolveLU_Flat
	for(int i = 0; i < totalDOF; i++)
	{
		for(int j = i + 1; j < totalDOF; j++)
		{
			// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
			std::swap(SystemMatrix[(size_t)i * totalDOF + j], SystemMatrix[(size_t)j * totalDOF + i]);
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

#ifdef RADIA_DEBUG_CHI
	fprintf(stderr, "=== LU Solver Debug: Solution (sigma/M) values ===\n");
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		int offset = IntrctPtr->GetElementDOFOffset(elem);
		fprintf(stderr, "Element %d (dof=%d): ", elem, dof);
		for(int k = 0; k < dof; k++)
		{
			fprintf(stderr, "%.3e ", RHS[offset + k]);
		}
		fprintf(stderr, "\n");
	}
#endif

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
		// Then subtract diagonal: y[i] -= inv_chi[i] * x[i] (ELF-compatible)

		// Intel MKL cblas_dgemv:
		// y := alpha*A*x + beta*y (for CblasNoTrans)
		// A is m x n, x is n, y is m
		// ROW-MAJOR (CblasRowMajor): lda = leading dimension = n = totalDOF (columns)
		// A[target][source] format: cache-efficient for matvec
		cblas_dgemv(CblasRowMajor, CblasNoTrans,
		            totalDOF, totalDOF,      // m, n (matrix dimensions)
		            -1.0,                     // alpha = -1.0 (negate for 3DOF)
		            FlatInteract, totalDOF,   // A, lda (row-major: lda = n)
		            x.data(), 1,              // x, incx
		            0.0,                      // beta
		            y.data(), 1);             // y, incy

		// Add diagonal contribution: y[i] += inv_chi[i] * x[i] (physically correct)
		// This is element-wise multiplication and addition
		#pragma omp parallel for if(totalDOF > 1000)
		for(int i = 0; i < totalDOF; i++)
		{
			y[i] += inv_chi[i] * x[i];  // Physically correct: +1/chi
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
		// Pure 6 DOF MSC system - use single BLAS call with alpha=-1.0
		// Physical equation: A = -K/(4pi) + diag(1/chi)
		// FlatInteract stores K/(4pi), so we need to negate
		// ROW-MAJOR (CblasRowMajor): A[target][source] format
		cblas_dgemv(CblasRowMajor, CblasNoTrans,
		            totalDOF, totalDOF,
		            -1.0,                     // alpha = -1.0 (negate to get -K/(4pi))
		            FlatInteract, totalDOF,   // A, lda (row-major: lda = n)
		            x.data(), 1,
		            0.0,
		            y.data(), 1);

		// Add diagonal contribution (physically correct: +1/chi)
		#pragma omp parallel for if(totalDOF > 1000)
		for(int i = 0; i < totalDOF; i++)
		{
			y[i] += inv_chi[i] * x[i];  // Physically correct: +1/chi
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

		// Diagonal contribution: +(1/chi) * x (physically correct)
		for(int k = 0; k < dof_row; k++)
		{
			y[offset_row + k] = inv_chi[offset_row + k] * x[offset_row + k];  // Physically correct: +1/chi
		}

		// Matrix-vector product
		for(int col_elem = 0; col_elem < AmOfMainElem; col_elem++)
		{
			int dof_col = IntrctPtr->GetElementDOF(col_elem);
			int offset_col = IntrctPtr->GetElementDOFOffset(col_elem);

			// Get block from flat matrix - COLUMN-MAJOR: block starts at [col * totalDOF + row]
			// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
			const double* block = &FlatInteract[(size_t)offset_col * totalDOF + offset_row];

			// Physical equation: A = -K/(4pi) + diag(1/chi)
			// FlatInteract stores K/(4pi) for MSC and N for MMM - negate both
			double sign = -1.0;

			// Column-major block access: element (i, j) within block is at [j * totalDOF + i]
			// CRITICAL: Use size_t cast for indexing with totalDOF
			for(int i = 0; i < dof_row; i++)
			{
				double sum = 0.0;
				for(int j = 0; j < dof_col; j++)
				{
					sum += block[(size_t)j * totalDOF + i] * x[offset_col + j];
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
	// Build block-Jacobi preconditioner by inverting each element's diagonal block
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
			// Singular block - use identity as fallback
			fprintf(stderr, "[Block Jacobi] Element %d: singular diagonal block (info=%d), using identity\n", elem, info);
			for(int i = 0; i < dof * dof; i++) block_copy[i] = 0;
			for(int i = 0; i < dof; i++) block_copy[i + i * dof] = 1.0;
		}
		else
		{
			dgetri_(&dof, block_copy.data(), &dof, ipiv.data(), work.data(), &lwork, &info);
			if(info != 0)
			{
				fprintf(stderr, "[Block Jacobi] Element %d: inversion failed (info=%d), using identity\n", elem, info);
				for(int i = 0; i < dof * dof; i++) block_copy[i] = 0;
				for(int i = 0; i < dof; i++) block_copy[i + i * dof] = 1.0;
			}
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
	// Apply block-Jacobi preconditioner: y = M^{-1} * x
	// where M is the block-diagonal of A
	int AmOfMainElem = IntrctPtr->AmOfMainElem;

	#pragma omp parallel for if(AmOfMainElem > 20)
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
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
	}
}
#endif

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
		else if(dof >= 5)
		{
			// MSC hexahedron (6 DOF): equation (K/(4pi) + 1/chi * I) * sigma = H_ext_n
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
			radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);

			// Use CurrentChi which was initialized in AutoRelax_VariableDOF
			double chi = (poly && poly->Use6DOF_MSC && poly->CurrentChi > 1.0e-6)
			           ? poly->CurrentChi : 1.0;
			if(chi < 1.0e-6) chi = 1.0e-6;

			double inv_chi_val = 1.0 / chi;

			for(int k = 0; k < dof; k++)
			{
				inv_chi[offset + k] = inv_chi_val;
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

	// Analyze matrix conditioning to choose preconditioner
	bool use_block_jacobi = false;
	double diag_ratio = 1.0;
	double min_dominance_val = 1.0;
	{
		const double* FlatInteract = IntrctPtr->GetFlatInteractMatrix();
		double min_diag = 1e30, max_diag = 0;
		double min_rowsum = 1e30, max_rowsum = 0;
		double min_dominance = 1e30;
		int weak_rows = 0;

		for(int i = 0; i < totalDOF; i++)
		{
			// Diagonal: -K[i,i] + 1/chi[i]
			double diag = -FlatInteract[(size_t)i * totalDOF + i] + inv_chi[i];
			double abs_diag = std::abs(diag);
			if(abs_diag < min_diag) min_diag = abs_diag;
			if(abs_diag > max_diag) max_diag = abs_diag;

			// Row sum (off-diagonal): sum |K[i,j]| for j != i
			double rowsum = 0;
			for(int j = 0; j < totalDOF; j++)
			{
				if(j != i) rowsum += std::abs(FlatInteract[(size_t)i * totalDOF + j]);
			}
			if(rowsum < min_rowsum) min_rowsum = rowsum;
			if(rowsum > max_rowsum) max_rowsum = rowsum;

			// Diagonal dominance ratio
			double dominance = abs_diag / (rowsum + 1e-15);
			if(dominance < min_dominance) min_dominance = dominance;
			if(dominance < 0.1) weak_rows++;
		}

		diag_ratio = max_diag / (min_diag + 1e-15);
		min_dominance_val = min_dominance;

		fprintf(stderr, "[BiCG Conditioning] DOF=%d\n", totalDOF);
		fprintf(stderr, "[BiCG Conditioning] Diagonal: min=%.3e, max=%.3e, ratio=%.1f\n",
		        min_diag, max_diag, diag_ratio);
		fprintf(stderr, "[BiCG Conditioning] Row sum: min=%.3e, max=%.3e\n", min_rowsum, max_rowsum);
		fprintf(stderr, "[BiCG Conditioning] Min dominance: %.4f, weak rows: %d/%d (%.1f%%)\n",
		        min_dominance, weak_rows, totalDOF, 100.0 * weak_rows / totalDOF);

		// Use block Jacobi if matrix is poorly conditioned
		// Criteria: diagonal ratio > 10 OR min dominance < 0.1
		if(diag_ratio > 10.0 || min_dominance < 0.1)
		{
			use_block_jacobi = true;
			fprintf(stderr, "[BiCG] Using BLOCK Jacobi preconditioner (better for ill-conditioned matrix)\n");
		}
	}

	// Build preconditioner
	std::vector<double> blockInverse;
	std::vector<int> blockOffsets;

#ifdef HAVE_LAPACK
	if(use_block_jacobi)
	{
		// Build block Jacobi preconditioner
		if(!BuildBlockJacobiPreconditioner_VariableDOF(blockInverse, blockOffsets, inv_chi, totalDOF))
		{
			fprintf(stderr, "[BiCG] Warning: Block Jacobi build failed, falling back to scalar Jacobi\n");
			use_block_jacobi = false;
		}
	}
#else
	use_block_jacobi = false;  // Block Jacobi requires LAPACK
#endif

	if(!use_block_jacobi)
	{
		// Build scalar Jacobi preconditioner
		GetDiagonalElements_VariableDOF(diag_inv, inv_chi, totalDOF);
		for(int i = 0; i < totalDOF; i++)
		{
			diag_inv[i] = (std::abs(diag_inv[i]) > 1.0e-15) ? (1.0 / diag_inv[i]) : 1.0;
		}
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

		// Apply preconditioner: p_hat = M^{-1} * p
#ifdef HAVE_LAPACK
		if(use_block_jacobi)
		{
			ApplyBlockJacobiPreconditioner_VariableDOF(p, p_hat, blockInverse, blockOffsets);
		}
		else
#endif
		{
			#pragma omp parallel for if(totalDOF > 100)
			for(int i = 0; i < totalDOF; i++)
			{
				p_hat[i] = diag_inv[i] * p[i];
			}
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

		// Apply preconditioner: s_hat = M^{-1} * s
#ifdef HAVE_LAPACK
		if(use_block_jacobi)
		{
			ApplyBlockJacobiPreconditioner_VariableDOF(s, s_hat, blockInverse, blockOffsets);
		}
		else
#endif
		{
			#pragma omp parallel for if(totalDOF > 100)
			for(int i = 0; i < totalDOF; i++)
			{
				s_hat[i] = diag_inv[i] * s[i];
			}
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
			fprintf(stderr, "[BiCG] Exit: divergence at iter %d\n", iter);
			break;
		}

		if(std::abs(omega) < 1.0e-30) {
			fprintf(stderr, "[BiCG] Exit: omega=%.4e at iter %d\n", omega, iter);
			break;
		}
	}

	// DEBUG: Always print solution info
	fprintf(stderr, "[BiCG] %d iters, residual=%.4e, IMA=%d\n", iter, residual, IntrctPtr->IsIMAEnabled() ? 1 : 0);
	if(residual > 1e-3)
	{
		fprintf(stderr, "[BiCG WARNING] High residual! Solution may be incorrect.\n");
		// Verify MatVec: compute A*sol and compare to rhs
		std::vector<double> Asol(totalDOF);
		MatVec_VariableDOF(sol, Asol, inv_chi, totalDOF);
		double err_norm = 0.0;
		for(int i = 0; i < totalDOF; i++)
		{
			double err = Asol[i] - rhs[i];
			err_norm += err * err;
		}
		err_norm = std::sqrt(err_norm);
		fprintf(stderr, "[BiCG] ||A*sol - rhs|| = %.4e, ||rhs|| = %.4e\n", err_norm, rhs_norm);
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
		if(dof >= 5)
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

		if(dof != 3 && dof < 5)
		{
			// HACApK supports 3DOF tetrahedra and 5/6DOF MSC elements (wedges/hexahedra)
			std::cerr << "[HACApK] Error: Element " << elem << " has " << dof
			          << " DOF, expected 3 (tetrahedra), 5 (wedges), or 6 (hexahedra)" << std::endl;
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
		else if(dof >= 5)
		{
			// For 6DOF MSC: use isotropic chi from poly->CurrentChi (same as BiCGSTAB)
			radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
			double chi = (poly && poly->Use6DOF_MSC && poly->CurrentChi > 1.0e-6)
			           ? poly->CurrentChi : 1.0;
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
	// FIX (2025-12-27): Recompute diagonal EVERY iteration (not cached)
	// Reason: ELF extracts diagonal from updated H-matrix blocks via HACApK_extract_diagonal_omp,
	// which includes the current inv_chi values. Radia was caching the diagonal computed with
	// INITIAL inv_chi, causing slower convergence (17 iter vs 14 iter).
	//
	// The fix is to recompute diag_inv = 1/(N_ii - inv_chi[i]) each iteration (ELF-compatible).
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
	// Get diagonal elements A_ii = -K_ii/(4pi) + 1/chi_i (physically correct)
	// GetDiagonalN() returns raw K/(4pi), so we negate it
	// Use cached N_ii values for efficiency (computed once during H-matrix build)
	if(m_hacapk->IsDiagonalCached())
	{
		const std::vector<double>& diag_N = m_hacapk->GetDiagonalN();
		#pragma omp parallel for if(totalDOF > 100)
		for(int i = 0; i < totalDOF; i++)
		{
			diag[i] = -diag_N[i] + inv_chi[i];  // Physical: -K/(4pi) + 1/chi
		}
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
		// We recompute diag_inv = 1/(N_ii - inv_chi[i]) each iteration (ELF-compatible)
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
		else if(dof >= 5)
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
		if(dof >= 5)
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
			else if(dof >= 5)
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
			else if(dof >= 5)
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
			else if(dof >= 5)
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
