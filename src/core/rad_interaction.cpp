/*-------------------------------------------------------------------------
*
* File name:      radintrc.cpp
*
* Project:        RADIA
*
* Description:    Magnetic interaction between "relaxable" field source objects
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#include "rad_interaction.h"
// #include "rad_subdivided_rectangle.h" REMOVED (Phase C, 2026-04-16)
#include "rad_polyhedron.h"  // For IsTetrahedron() check in N_self fix
#include "rad_constants.h"   // For RadConst::INV_FOUR_PI
#include <array>
#include <unordered_map>
#include <map>
#include <memory>
#include <mutex>
#include <cmath>
#include <algorithm>
#include <utility>
#include <chrono>
#include <atomic>

#include "rad_parallel.h"

//-------------------------------------------------------------------------
// Static member definitions for RadIMAFieldContext
//-------------------------------------------------------------------------
std::atomic<bool> RadIMAFieldContext::s_active{false};
std::atomic<int> RadIMAFieldContext::s_symmetry{0};
std::atomic<int> RadIMAFieldContext::s_signX{1};
std::atomic<int> RadIMAFieldContext::s_signY{1};
std::atomic<int> RadIMAFieldContext::s_signZ{1};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTInteraction::radTInteraction(const radThg& In_hg, const radThg& In_hgMoreExtSrc, const radTCompCriterium& InCompCriterium, short InMemAllocTotAtOnce, char ExtraExternFieldArrayIsNeeded, char KeepTransData, int rankMPI, int nProcMPI, char skipDenseMatrix) //OC08012020 + skipDenseMatrix
//radTInteraction::radTInteraction(const radThg& In_hg, const radThg& In_hgMoreExtSrc, const radTCompCriterium& InCompCriterium, short InMemAllocTotAtOnce, char ExtraExternFieldArrayIsNeeded, char KeepTransData)
{
	if(!Setup(In_hg, In_hgMoreExtSrc, InCompCriterium, InMemAllocTotAtOnce, ExtraExternFieldArrayIsNeeded, KeepTransData, rankMPI, nProcMPI, skipDenseMatrix)) //OC08012020 + skipDenseMatrix
	//if(!Setup(In_hg, In_hgMoreExtSrc, InCompCriterium, InMemAllocTotAtOnce, ExtraExternFieldArrayIsNeeded, KeepTransData))
	{
		SomethingIsWrong = 1;
		Send.ErrorMessage("Radia::Error118");
		throw 0;
	}
}

//-------------------------------------------------------------------------

radTInteraction::radTInteraction()
{
	AmOfMainElem = 0;
	AmOfExtElem = 0;
	InteractMatrix = nullptr;
	ExternFieldArray = nullptr;
	AuxOldMagnArray = nullptr;
	AuxOldFieldArray = nullptr;

	NewMagnArray = nullptr;
	NewFieldArray = nullptr;
	IdentTransPtr = nullptr;

	RelaxSubIntervArray = nullptr; // New
	mKeepTransData = 0;

	// Tetrahedron/Wedge/Hexahedron geometry cache
	m_tetraGeomReady = false;
	m_wedgeGeomReady = false;
	m_hexaGeomReady = false;
	m_hexaTriDataReady = false;

	// PEEC element flag
	m_hasPEECElements = false;

	// IMA flags - initialize to disabled
	m_imaEnabled = false;
	m_imaSymmetry = 0;
	m_imaSignX = 1;
	m_imaSignY = 1;
	m_imaSignZ = 1;
	m_imaNumElements = 0;
}

//-------------------------------------------------------------------------

int radTInteraction::Setup(const radThg& In_hg, const radThg& In_hgMoreExtSrc, const radTCompCriterium& InCompCriterium, short InMemAllocTotAtOnce, char AuxOldMagnArrayIsNeeded, char KeepTransData, int rankMPI, int nProcMPI, char skipDenseMatrix) //OC08012020 + skipDenseMatrix
//int radTInteraction::Setup(const radThg& In_hg, const radThg& In_hgMoreExtSrc, const radTCompCriterium& InCompCriterium, short InMemAllocTotAtOnce, char AuxOldMagnArrayIsNeeded, char KeepTransData)
{
	SomethingIsWrong = 0;

	AmOfMainElem = 0;
	AmOfExtElem = 0;
	InteractMatrix = nullptr;
	ExternFieldArray = nullptr;
	AuxOldMagnArray = nullptr;
	AuxOldFieldArray = nullptr;

	NewMagnArray = nullptr;
	NewFieldArray = nullptr;
	IdentTransPtr = nullptr;

	RelaxSubIntervArray = nullptr; // New
	AmOfRelaxSubInterv = 0; // New

	// Tetrahedron/Wedge/Hexahedron geometry cache
	m_tetraGeomReady = false;
	m_wedgeGeomReady = false;
	m_hexaGeomReady = false;
	m_hexaTriDataReady = false;

	// IMA flags - initialize to disabled
	m_imaEnabled = false;
	m_imaSymmetry = 0;
	m_imaSignX = 1;
	m_imaSignY = 1;
	m_imaSignZ = 1;
	m_imaNumElements = 0;

	SourceHandle = In_hg;
	CompCriterium = InCompCriterium;
	FillInMainTransOnly = 0;
	RelaxationStarted = 0;

	MoreExtSourceHandle = In_hgMoreExtSrc;

	MemAllocTotAtOnce = InMemAllocTotAtOnce;

	IdentTransPtr = new radIdentTrans();

	radTlphgPtr NewListOfTransPtr;
	CountMainRelaxElems(static_cast<radTg3d*>(SourceHandle.rep), &NewListOfTransPtr);

	if(!NotEmpty()) return 0;

	//m_rankMPI = -1; //OC20122019 (to set from Application?) 
	//m_nProcMPI = 0;
	m_rankMPI = rankMPI; //OC08012019 (to set from Application?) 
	m_nProcMPI = nProcMPI; 

	bool IntrctMatrMemAllocShouldBeDone = true;
	if(m_rankMPI > 0) IntrctMatrMemAllocShouldBeDone = false;

//	if(MPI_Comm_size(MPI_COMM_WORLD, &m_nProcMPI) != MPI_SUCCESS) { Send.ErrorMessage("Radia::Error601"); return 0;}
//	if(MPI_Comm_rank(MPI_COMM_WORLD, &m_rankMPI) != MPI_SUCCESS) { Send.ErrorMessage("Radia::Error601"); return 0;} //Get the rank of the process
//	if(m_rankMPI > 0) IntrctMatrMemAllocShouldBeDone = false;
//#endif

	if(IntrctMatrMemAllocShouldBeDone) //OC20122019
	{
		// For HACApK solver (skipDenseMatrix=1), skip InteractMatrix allocation
		// because HACApK builds its own H-matrix, dense matrix is unnecessary overhead
		// and would consume O(N^2) memory which exceeds available memory for large N
		AllocateMemory(AuxOldMagnArrayIsNeeded, skipDenseMatrix); //In case of MPI-parallelization, this has to be executed by master only

		if(SomethingIsWrong)
		{
			EmptyVectOfPtrToListsOfTrans(); return 0;
		}
		FillInRelaxSubIntervArray(); //New
	}
	FillInMainTransPtrArray();

	// Check if any element has variable DOF (e.g., 4-6 independent face coefficients).
	// If so, use the variable DOF interaction matrix setup
	ComputeDOFOffsets();

	int setupResult = 1;
	if(skipDenseMatrix)
	{
		// Skip dense matrix construction for H-matrix solver (HACApK)
		// H-matrix builds its own compressed matrix, no need for dense matrix
		// Just allocate the auxiliary arrays needed for nonlinear iteration
		// Allocate flat arrays (always needed now for unified solver)
		// Note: ComputeDOFOffsets() already called above at line 134
		m_flatExternFieldArray.resize(m_totalDOF, 0.0);
		m_flatMagnArray.resize(m_totalDOF, 0.0);
		m_flatFieldArray.resize(m_totalDOF, 0.0);
		// Allocate standard arrays (needed for chi update in nonlinear iteration)
		vNewMagnArray.resize(AmOfMainElem);
		vNewFieldArray.resize(AmOfMainElem);
		vExternFieldArray.resize(AmOfMainElem);
		NewMagnArray = vNewMagnArray.data();
		NewFieldArray = vNewFieldArray.data();
		ExternFieldArray = vExternFieldArray.data();
		// Keep face geometry caches coherent when dense matrix construction is skipped.
		// Mesh-backed magnetic-material solves are handled by HDiv-VIM; this cache is
		// retained only for fixed-magnetization field and diagnostic paths.
		PrecomputeHexaGeometry();
	}
	else
	{
		// Build full dense interaction matrix (for LU/BiCGSTAB solvers)
		// Always use VariableDOF matrix format for unified nonlinear solver
		// This works for both 3DOF (tetrahedra) and 6DOF (hexahedra) elements
		setupResult = SetupInteractMatrix_VariableDOF();
	}
	if(!setupResult) { DeallocateMemory(); return 0;} //OC26122019 //Most CPU-intensive

	if(IntrctMatrMemAllocShouldBeDone) //OC29122019
	{
		SetupExternFieldArray();
		AddExternFieldFromMoreExtSource();
		//ZeroAuxOldMagnArray();
		ZeroAuxOldArrays();

		InitAuxArrays();
	}

	mKeepTransData = KeepTransData;
	if(!KeepTransData) //OC021103
	{
		DestroyMainTransPtrArray();
		EmptyVectOfPtrToListsOfTrans();
	}

	////ResetM();
	//InitAuxArrays(); //OC30122019 (moved up)

	return 1;
}

//-------------------------------------------------------------------------

radTInteraction::~radTInteraction()
{
	DeallocateMemory(); //OC27122019
}

//-------------------------------------------------------------------------

bool radTInteraction::HasSurfaceChargeElements() const
{
	// True if any relaxable element carries independent face coefficients.
	// Production mesh-backed magnetic-material solves are routed to HDiv-VIM.
	for(int i = 0; i < AmOfMainElem; i++)
	{
		radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[i]);
		if(poly && poly->UseFaceChargeDOF) return true;
	}
	return false;
}

//-------------------------------------------------------------------------

void radTInteraction::DeallocateMemory() //OC27122019
{
	// RAII: automatic cleanup via vInteractMatrix, vInteractMatrixPtrs, and vGenMatrStorage

	g3dExternPtrVect.erase(g3dExternPtrVect.begin(), g3dExternPtrVect.end()); //OC240408, to enable current scaling/update

	// Automatic cleanup via RAII for vector arrays
	// RelaxSubIntervArray: automatic cleanup via vRelaxSubIntervArray

	if(mKeepTransData) //OC021103
	{
		DestroyMainTransPtrArray();
		EmptyVectOfPtrToListsOfTrans();
	}
	if(IdentTransPtr != nullptr) delete IdentTransPtr; //required by EmptyVectOfPtrToListsOfTrans();
}

//-------------------------------------------------------------------------

void radTInteraction::CountMainRelaxElems(radTg3d* g3dPtr, radTlphgPtr* CurrListOfTransPtrPtr)
{
	radTGroup* GroupPtr = Cast.GroupCast(g3dPtr);
	if(GroupPtr == 0)
	{
		radTg3dRelax* g3dRelaxPtr = Cast.g3dRelaxCast(g3dPtr);
		if((g3dRelaxPtr != 0) && (g3dRelaxPtr->MaterHandle.rep != 0))
		{
			g3dRelaxPtrVect.push_back(g3dRelaxPtr);
			AmOfMainElem++;

			radTlphgPtr* TotalListOfElemTransPtrPtr = new radTlphgPtr(*CurrListOfTransPtrPtr);
			PushFrontNativeElemTransList(g3dRelaxPtr, TotalListOfElemTransPtrPtr);
			IntVectOfPtrToListsOfTransPtr.push_back(TotalListOfElemTransPtrPtr);
		}
		else 
		{
			g3dExternPtrVect.push_back(g3dPtr);
			AmOfExtElem++;

			radTlphgPtr* TotalListOfElemTransPtrPtr	= new radTlphgPtr(*CurrListOfTransPtrPtr);
			PushFrontNativeElemTransList(g3dPtr, TotalListOfElemTransPtrPtr);
			ExtVectOfPtrToListsOfTransPtr.push_back(TotalListOfElemTransPtrPtr);
		}
	}
	else
	{
		// radTSubdividedRecMag handling REMOVED (Phase C, 2026-04-16). Generic Group traversal only.
		radTlphgPtr* LocListOfTransPtrPtr = CurrListOfTransPtrPtr;

		short GroupListOfTransIsNotEmpty = 1;
		if(GroupPtr->g3dListOfTransform.empty()) GroupListOfTransIsNotEmpty = 0;

		if(GroupListOfTransIsNotEmpty)
		{
			LocListOfTransPtrPtr = new radTlphgPtr(*CurrListOfTransPtrPtr);
			PushFrontNativeElemTransList(GroupPtr, LocListOfTransPtrPtr);
		}

		for(radTmhg::iterator iter = GroupPtr->GroupMapOfHandlers.begin();
			iter != GroupPtr->GroupMapOfHandlers.end(); ++iter)
			CountMainRelaxElems(static_cast<radTg3d*>(iter->second.rep), LocListOfTransPtrPtr);

		if(GroupListOfTransIsNotEmpty) delete LocListOfTransPtrPtr;
	}
}

//-------------------------------------------------------------------------

void radTInteraction::FillInRelaxSubIntervArray() // New
{
	// Reset and rebuild from scratch
	AmOfRelaxSubInterv = 0;

	if(RelaxSubIntervConstrVect.size() == 0) return;

	int CurrentStartNo = 0;
	int PlainCount = -1;

	for(auto Iter = RelaxSubIntervConstrVect.begin(); Iter != RelaxSubIntervConstrVect.end(); ++Iter)
	{
		int LocStartNo = (*Iter).StartNo;
		if(LocStartNo != CurrentStartNo)
		{
			RelaxSubIntervArray[++PlainCount] = radTRelaxSubInterval(CurrentStartNo, LocStartNo-1, TRelaxSubIntervalID::RelaxApart);
		}
		RelaxSubIntervArray[++PlainCount] = *Iter;
		CurrentStartNo = (*Iter).FinNo + 1;
	}
	if(CurrentStartNo != AmOfMainElem)
		RelaxSubIntervArray[++PlainCount] = radTRelaxSubInterval(CurrentStartNo, AmOfMainElem-1, TRelaxSubIntervalID::RelaxApart);

	AmOfRelaxSubInterv = ++PlainCount;

	// Do NOT erase RelaxSubIntervConstrVect - keep it for future rebuilds
	// RelaxSubIntervConstrVect.erase(RelaxSubIntervConstrVect.begin(), RelaxSubIntervConstrVect.end());
}

//-------------------------------------------------------------------------

void radTInteraction::AddRelaxSubInterval(int StartNo, int FinNo, TRelaxSubIntervalID SubIntervalID)
{
	if(StartNo < 0 || FinNo < 0) return;
	if(StartNo > FinNo) return;
	if(FinNo >= AmOfMainElem) return;

	radTRelaxSubInterval SubInterval(StartNo, FinNo, SubIntervalID);
	RelaxSubIntervConstrVect.push_back(SubInterval);

	// Reallocate RelaxSubIntervArray with sufficient size
	int MaxSubIntervArraySize = 2 * ((int)(RelaxSubIntervConstrVect.size())) + 1;
	if(MaxSubIntervArraySize > (int)vRelaxSubIntervArray.size())
	{
		vRelaxSubIntervArray.resize(MaxSubIntervArraySize);
		RelaxSubIntervArray = vRelaxSubIntervArray.data();
	}

	// Rebuild RelaxSubIntervArray after adding new interval
	FillInRelaxSubIntervArray();
}

//-------------------------------------------------------------------------

void radTInteraction::AllocateMemory(char AuxOldMagnArrayIsNeeded, char skipInteractMatrix)
{
	vExternFieldArray.resize(AmOfMainElem);
	ExternFieldArray = vExternFieldArray.data();

	if(AuxOldMagnArrayIsNeeded)
	{
		vAuxOldMagnArray.resize(AmOfMainElem);
		vAuxOldFieldArray.resize(AmOfMainElem);
		AuxOldMagnArray = vAuxOldMagnArray.data();
		AuxOldFieldArray = vAuxOldFieldArray.data();
	}

	vNewMagnArray.resize(AmOfMainElem);
	vNewFieldArray.resize(AmOfMainElem);
	NewMagnArray = vNewMagnArray.data();
	NewFieldArray = vNewFieldArray.data();

	// Skip InteractMatrix allocation for HACApK solver
	// HACApK builds its own H-matrix, dense matrix is unnecessary overhead
	// For 100k DOF, dense matrix would require ~80 GB of memory
	if(skipInteractMatrix)
	{
		// Just initialize pointers to nullptr
		vInteractMatrixPtrs.resize(AmOfMainElem, nullptr);
		InteractMatrix = vInteractMatrixPtrs.data();
	}
	else
	{
		// Check memory requirements before allocation (legacy LU needs a dense matrix)
		// Dense matrix requires N^2 * sizeof(TMatrix3df) = N^2 * 36 bytes
		size_t matrix_size = (size_t)AmOfMainElem * (size_t)AmOfMainElem;
		size_t required_bytes = matrix_size * sizeof(TMatrix3df);
		const size_t MAX_DENSE_MATRIX_BYTES = 100ULL * 1024 * 1024 * 1024;  // 100 GB limit

		if(required_bytes > MAX_DENSE_MATRIX_BYTES)
		{
			std::cerr << "[Radia] Error: Dense matrix too large for the legacy LU solver." << std::endl;
			std::cerr << "[Radia] Elements=" << AmOfMainElem << ", required memory="
			          << (required_bytes / (1024*1024*1024)) << " GB" << std::endl;
			std::cerr << "[Radia] Use a mesh-backed soft-iron model with HDiv-VIM for large problems." << std::endl;
			std::cerr.flush();
			SomethingIsWrong = 1;
			return;
		}

		vInteractMatrixPtrs.resize(AmOfMainElem, nullptr);
		InteractMatrix = vInteractMatrixPtrs.data();

		try {
			if(MemAllocTotAtOnce)
			{
				vGenMatrStorage.resize(AmOfMainElem * AmOfMainElem);
				TMatrix3df* GenMatrPtr = vGenMatrStorage.data();

				for(int i=0; i<AmOfMainElem; i++)
				{
					InteractMatrix[i] = &(GenMatrPtr[i*AmOfMainElem]);
					vInteractMatrixPtrs[i] = InteractMatrix[i];
				}
			}
			else
			{
				vInteractMatrix.resize(AmOfMainElem);
				for(int i=0; i<AmOfMainElem; i++)
				{
					vInteractMatrix[i].resize(AmOfMainElem);
					InteractMatrix[i] = vInteractMatrix[i].data();
					vInteractMatrixPtrs[i] = InteractMatrix[i];
				}
			}
		} catch(const std::bad_alloc&) {
			std::cerr << "[Radia] Error: Memory allocation failed for dense interaction matrix." << std::endl;
			std::cerr << "[Radia] Use a mesh-backed soft-iron model with HDiv-VIM for large problems." << std::endl;
			std::cerr.flush();
			SomethingIsWrong = 1;
			return;
		}
	}

	int MaxSubIntervArraySize = 2 * ((int)(RelaxSubIntervConstrVect.size())) + 1; // New
	//try
	//{
		if(MaxSubIntervArraySize > 1)
		{
			vRelaxSubIntervArray.resize(MaxSubIntervArraySize);
			RelaxSubIntervArray = vRelaxSubIntervArray.data();
		}
	//}
	//catch (radTException* radExceptionPtr)
	//{
	//	Send.ErrorMessage(radExceptionPtr->what());	return;
	//}
	//catch (...)
	//{
	//	Send.ErrorMessage("Radia::Error999"); return;
	//}
}

//-------------------------------------------------------------------------

void radTInteraction::AllocateInteractMatrix()
{
	// Allocate only InteractMatrix (for HACApK 3DOF case where Setup was called with skipDenseMatrix=1)
	// This is needed when HACApK needs to use pre-computed interaction matrix for 3DOF tetrahedra

	if(AmOfMainElem <= 0) return;
	if(InteractMatrix != nullptr && vInteractMatrixPtrs.size() > 0 && vInteractMatrixPtrs[0] != nullptr)
	{
		// Already allocated
		return;
	}

	vInteractMatrixPtrs.resize(AmOfMainElem, nullptr);
	InteractMatrix = vInteractMatrixPtrs.data();

	if(MemAllocTotAtOnce)
	{
		vGenMatrStorage.resize(AmOfMainElem * AmOfMainElem);
		TMatrix3df* GenMatrPtr = vGenMatrStorage.data();

		for(int i=0; i<AmOfMainElem; i++)
		{
			InteractMatrix[i] = &(GenMatrPtr[i*AmOfMainElem]);
			vInteractMatrixPtrs[i] = InteractMatrix[i];
		}
	}
	else
	{
		vInteractMatrix.resize(AmOfMainElem);
		for(int i=0; i<AmOfMainElem; i++)
		{
			vInteractMatrix[i].resize(AmOfMainElem);
			InteractMatrix[i] = vInteractMatrix[i].data();
			vInteractMatrixPtrs[i] = InteractMatrix[i];
		}
	}
}

//-------------------------------------------------------------------------

void radTInteraction::NestedFor_Trans(radTrans* BaseTransPtr, const radTlphgPtr::const_iterator& Iter, int ElemLocInd, char I_or_E)
{
	radTrans* TransPtr = (radTrans*)(((**Iter).Handler_g).rep);
	radTrans* LocTotTransPtr = BaseTransPtr;
	radTrans LocTotTrans;

	radTlphgPtr::const_iterator LocalNextIter = Iter;
	LocalNextIter++;
	int Mult = (**Iter).m;

	if(Mult == 1)
	{
		TrProduct(LocTotTransPtr, TransPtr, LocTotTrans);
		AddTransOrNestedFor(&LocTotTrans, LocalNextIter, ElemLocInd, I_or_E);
	}
	else
	{
		AddTransOrNestedFor(LocTotTransPtr, LocalNextIter, ElemLocInd, I_or_E);
		if(FillInMainTransOnly) return;
		for(int km = 1; km < Mult; km++)
		{
			TrProduct(LocTotTransPtr, TransPtr, LocTotTrans);
			LocTotTransPtr = &LocTotTrans;
			AddTransOrNestedFor(LocTotTransPtr, LocalNextIter, ElemLocInd, I_or_E);
		}
	}
}

//-------------------------------------------------------------------------

void radTInteraction::FillInMainTransPtrArray()
{
	vMainTransPtrArray.resize(AmOfMainElem);
	MainTransPtrArray = vMainTransPtrArray.data();
	FillInMainTransOnly = 1;

	for(int i=0; i<AmOfMainElem; i++)
	{
		FillInTransPtrVectForElem(i, 'I');
		if(Cast.IdentTransCast(TransPtrVect[0]) == 0) 
		{
			MainTransPtrArray[i] = new radTrans(*(TransPtrVect[0]));
		}
		else MainTransPtrArray[i] = IdentTransPtr;
		EmptyTransPtrVect();
	}
	FillInMainTransOnly = 0;
}

//-------------------------------------------------------------------------

int radTInteraction::CountRelaxElemsWithSym()
{
	int AmOfElemWithSym = 0;

	for(int i=0; i<AmOfMainElem; i++)
	{
		radTlphgPtr& Loc_lphgPtr = *(IntVectOfPtrToListsOfTransPtr[i]);
		int LocTotMult = 1;

		for(radTlphgPtr::iterator TrIter = Loc_lphgPtr.begin();	
			TrIter != Loc_lphgPtr.end(); ++TrIter)
		{
			LocTotMult *= (**TrIter).m;
		}
		AmOfElemWithSym += LocTotMult;
	}
	return AmOfElemWithSym;
}

//-------------------------------------------------------------------------

int radTInteraction::SetupInteractMatrix() //OC26122019
//void radTInteraction::SetupInteractMatrix()
{
	radTFieldKey FieldKeyInteract; FieldKeyInteract.B_=FieldKeyInteract.H_=FieldKeyInteract.PreRelax_=1;
	TVector3d ZeroVect(0.,0.,0.);

	//--New
	int AmOfElemWithSym = CountRelaxElemsWithSym();
	//--EndNew

	if(m_nProcMPI < 2) //OC01012020
	{
		// Simplified global coordinate version:
		// Compute interaction matrix directly in global coordinates
		// without intermediate coordinate transformations.
		// This matches the behavior of rad.Fld() which gives correct results.
		//
		// OPTIMIZATION (2025-12-11): OpenMP parallelization of O(N^2) matrix build
		// Parallelize outer loop only (MSVC OpenMP 2.0 doesn't support collapse)
		// Each ColNo iteration is independent, enabling parallel computation.

		ngcore::ParallelFor(ngcore::IntRange(AmOfMainElem), [&](size_t ColNo)
		{
			radTg3dRelax* g3dRelaxPtrColNo = g3dRelaxPtrVect[(int)ColNo];

			for(int StrNo=0; StrNo<AmOfMainElem; StrNo++)
			{
				// Get observation point (element center) directly in global coordinates
				TVector3d ObsPoiVect = (g3dRelaxPtrVect[StrNo])->ReturnCentrPoint();

				// Create thread-local Field object to avoid race conditions
				radTField Field(FieldKeyInteract, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
				Field.AmOfIntrctElemWithSym = AmOfElemWithSym;

				// Compute field contribution in global coordinates
				g3dRelaxPtrColNo->B_comp(&Field);

				// Store result directly (no transformation)
				TMatrix3d SubMatrix;
				SubMatrix.Str0 = Field.B;  // dH/dMx
				SubMatrix.Str1 = Field.H;  // dH/dMy
				SubMatrix.Str2 = Field.A;  // dH/dMz

				InteractMatrix[StrNo][(int)ColNo] = SubMatrix;
			}
		});

		//--New
		for(int ClNo=0; ClNo<AmOfMainElem; ClNo++)
		{
			radTg3dRelax* g3dRelaxPtrClNo = g3dRelaxPtrVect[ClNo];
			g3dRelaxPtrVect[ClNo] = g3dRelaxPtrClNo->FormalIntrctMemberPtr();
		}
		//--EndNew
	}
	return 1; //OC26122019
}

//-------------------------------------------------------------------------
//=========================================================================
// Variable DOF support for hybrid surface-charge + standard element analysis.
//=========================================================================
//-------------------------------------------------------------------------

void radTInteraction::ComputeDOFOffsets()
{
	// Compute DOF offsets for each element
	// This allows mixing elements with different DOF counts.

	m_elemDOF.resize(AmOfMainElem);
	m_elemDOFOffset.resize(AmOfMainElem + 1);  // +1 for end sentinel

	m_totalDOF = 0;
	m_hasVariableDOF = false;
	m_hasPEECElements = false;  // Will be set by coupled solver if PEEC elements are present

	for(int i = 0; i < AmOfMainElem; i++)
	{
		m_elemDOFOffset[i] = m_totalDOF;
		int dof = g3dRelaxPtrVect[i]->NumberOfDegOfFreedom();
		m_elemDOF[i] = dof;
		m_totalDOF += dof;

		// Check if any element has non-standard DOF
		if(dof != 3)
		{
			m_hasVariableDOF = true;
		}
	}
	m_elemDOFOffset[AmOfMainElem] = m_totalDOF;  // End sentinel
}

//-------------------------------------------------------------------------

double* radTInteraction::GetInteractBlock(int row_elem, int col_elem)
{
	// Return pointer to interaction block (row_elem, col_elem) in flattened matrix
	// Matrix is ROW-MAJOR: A(i,j) at index [i * m_totalDOF + j]
	if(m_flatInteractMatrix.empty()) return nullptr;

	int offset_row = m_elemDOFOffset[row_elem];
	int offset_col = m_elemDOFOffset[col_elem];

	// Block starts at row offset_row, column offset_col in the totalDOF x totalDOF matrix
	// ROW-MAJOR: element at (row, col) is at index [row * m_totalDOF + col]
	// A[target][source] format: ELF-compatible, BiCGSTAB/HACApK-optimal
	// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
	return &m_flatInteractMatrix[(size_t)offset_row * m_totalDOF + offset_col];
}

const double* radTInteraction::GetInteractBlock(int row_elem, int col_elem) const
{
	if(m_flatInteractMatrix.empty()) return nullptr;

	int offset_row = m_elemDOFOffset[row_elem];
	int offset_col = m_elemDOFOffset[col_elem];

	// ROW-MAJOR: element at (row, col) is at index [row * m_totalDOF + col]
	// A[target][source] format: ELF-compatible, BiCGSTAB/HACApK-optimal
	// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
	return &m_flatInteractMatrix[(size_t)offset_row * m_totalDOF + offset_col];
}

//-------------------------------------------------------------------------

void radTInteraction::SetupVariableDOFArrays()
{
	// Set up flat arrays for HACApK solver without building interaction matrix
	// This is used when all elements are 3DOF (tetrahedra) but HACApK is requested
	// HACApK computes matrix elements on-demand using callbacks

	if(m_totalDOF <= 0)
	{
		// Need to compute DOF offsets first
		ComputeDOFOffsets();
	}

	// Allocate flattened field arrays (NOT the interaction matrix)
	m_flatExternFieldArray.resize(m_totalDOF, 0.0);
	m_flatMagnArray.resize(m_totalDOF, 0.0);
	m_flatFieldArray.resize(m_totalDOF, 0.0);

	// Also allocate standard arrays (NewFieldArray, NewMagnArray, ExternFieldArray)
	// These are needed for chi update in nonlinear iteration
	vNewMagnArray.resize(AmOfMainElem);
	vNewFieldArray.resize(AmOfMainElem);
	vExternFieldArray.resize(AmOfMainElem);
	NewMagnArray = vNewMagnArray.data();
	NewFieldArray = vNewFieldArray.data();
	ExternFieldArray = vExternFieldArray.data();
}

//-------------------------------------------------------------------------

int radTInteraction::SetupInteractMatrix_VariableDOF()
{
	// Build interaction matrix with variable DOF blocks
	// This function now builds flat matrix format for BOTH 3DOF and 6DOF elements
	// This enables using the unified nonlinear solver for all element types

	// First compute DOF offsets
	ComputeDOFOffsets();

	// Note: We no longer redirect to SetupInteractMatrix() for 3DOF-only cases
	// The flat matrix format is needed for the unified VariableDOF solver

	// Check memory requirements before allocation
	// Dense matrix requires N^2 doubles = N^2 * 8 bytes
	// For 100k DOF, this is ~80 GB which exceeds typical system memory
	size_t matrix_size = (size_t)m_totalDOF * (size_t)m_totalDOF;
	size_t required_bytes = matrix_size * sizeof(double);
	const size_t MAX_DENSE_MATRIX_BYTES = 100ULL * 1024 * 1024 * 1024;  // 100 GB limit

	if(required_bytes > MAX_DENSE_MATRIX_BYTES)
	{
		std::cerr << "[Radia] Error: Dense matrix too large for the legacy LU solver." << std::endl;
		std::cerr << "[Radia] DOF=" << m_totalDOF << ", required memory="
		          << (required_bytes / (1024*1024*1024)) << " GB" << std::endl;
		std::cerr << "[Radia] Use a mesh-backed soft-iron model with HDiv-VIM for large problems." << std::endl;
		std::cerr.flush();
		return 0;  // Signal failure
	}

	// Allocate flattened interaction matrix
	// CRITICAL: Use size_t to avoid int32 overflow for DOF > 46340 (sqrt(INT32_MAX))
	try {
		m_flatInteractMatrix.resize(matrix_size, 0.0);
	} catch(const std::bad_alloc&) {
		std::cerr << "[Radia] Error: Memory allocation failed for dense interaction matrix." << std::endl;
		std::cerr << "[Radia] DOF=" << m_totalDOF << ", required memory="
		          << (required_bytes / (1024*1024*1024)) << " GB" << std::endl;
		std::cerr << "[Radia] Use a mesh-backed soft-iron model with HDiv-VIM for large problems." << std::endl;
		return 0;  // Signal failure
	}

	// Allocate flattened field arrays
	m_flatExternFieldArray.resize(m_totalDOF, 0.0);
	m_flatMagnArray.resize(m_totalDOF, 0.0);
	m_flatFieldArray.resize(m_totalDOF, 0.0);

	// Also allocate standard arrays (NewFieldArray, NewMagnArray, ExternFieldArray)
	// These are needed for chi update in nonlinear iteration
	vNewMagnArray.resize(AmOfMainElem);
	vNewFieldArray.resize(AmOfMainElem);
	vExternFieldArray.resize(AmOfMainElem);
	NewMagnArray = vNewMagnArray.data();
	NewFieldArray = vNewFieldArray.data();
	ExternFieldArray = vExternFieldArray.data();

	radTFieldKey FieldKeyInteract;
	FieldKeyInteract.B_ = FieldKeyInteract.H_ = FieldKeyInteract.PreRelax_ = 1;
	TVector3d ZeroVect(0., 0., 0.);

	int AmOfElemWithSym = CountRelaxElemsWithSym();

	// Check if we have symmetry transformations
	// If no symmetries, we can use simplified global coordinate computation with OpenMP
	bool hasSymmetry = (AmOfElemWithSym > AmOfMainElem);

	// Check if any element uses independent face coefficients.  Dense relaxation
	// matrix assembly is kept only for compact 3-component elements; mesh-backed
	// magnetic-material solves route through HDiv-VIM.
	bool hasFaceChargeElements = false;
	for(int i = 0; i < AmOfMainElem && !hasFaceChargeElements; i++)
	{
		if(m_elemDOF[i] >= 4) hasFaceChargeElements = true;
	}

	// Build interaction matrix with variable-size blocks
	// For each pair (row_elem, col_elem), compute the interaction block

	// FAST PATH: pure compact 3-component elements without symmetry.
	if(!hasSymmetry && !hasFaceChargeElements)
	{
		// Pre-compute tetrahedron geometry for fast block computation
		// This caches face vertices, normals, and areas to avoid B_comp overhead
		PrecomputeTetraGeometry();

		if(m_tetraGeomReady)
		{
			// ULTRA-FAST PATH: Use pre-computed geometry
			// Compute3x3BlockFast uses cached geometry arrays (no B_comp overhead)
			ngcore::ParallelFor(ngcore::IntRange(AmOfMainElem), [&](size_t col)
			{
				int offset_col = m_elemDOFOffset[(int)col];

				for(int row = 0; row < AmOfMainElem; row++)
				{
					int offset_row = m_elemDOFOffset[row];

					// Get pointer to this block in the flattened matrix
					// ROW-MAJOR: A(row, col) at index [row * m_totalDOF + col]
					// A[target][source] format: ELF-compatible, BiCGSTAB/HACApK-optimal
					// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
					double* block = &m_flatInteractMatrix[(size_t)offset_row * m_totalDOF + offset_col];

					// Compute 3x3 block using cached geometry
					double N_mat[9];
					Compute3x3BlockFast(row, (int)col, N_mat);

					// Copy to flat matrix (both row-major, direct copy)
					// N_mat is [target][source] row-major, matrix is [target][source] row-major
					// CRITICAL: Use size_t cast for all indexing with m_totalDOF
					block[(size_t)0 * m_totalDOF + 0] = N_mat[0];  // (0,0)
					block[(size_t)0 * m_totalDOF + 1] = N_mat[1];  // (0,1)
					block[(size_t)0 * m_totalDOF + 2] = N_mat[2];  // (0,2)
					block[(size_t)1 * m_totalDOF + 0] = N_mat[3];  // (1,0)
					block[(size_t)1 * m_totalDOF + 1] = N_mat[4];  // (1,1)
					block[(size_t)1 * m_totalDOF + 2] = N_mat[5];  // (1,2)
					block[(size_t)2 * m_totalDOF + 0] = N_mat[6];  // (2,0)
					block[(size_t)2 * m_totalDOF + 1] = N_mat[7];  // (2,1)
					block[(size_t)2 * m_totalDOF + 2] = N_mat[8];  // (2,2)
				}
			});
			return 1;
		}

		// Fallback: FAST PATH without geometry cache (original B_comp method)
		ngcore::ParallelFor(ngcore::IntRange(AmOfMainElem), [&](size_t col)
		{
			radTg3dRelax* elem_col = g3dRelaxPtrVect[(int)col];
			int dof_col = m_elemDOF[(int)col];
			int offset_col = m_elemDOFOffset[(int)col];

			for(int row = 0; row < AmOfMainElem; row++)
			{
				radTg3dRelax* elem_row = g3dRelaxPtrVect[row];
				int dof_row = m_elemDOF[row];
				int offset_row = m_elemDOFOffset[row];

				// Get pointer to this block in the flattened matrix
				// ROW-MAJOR: A(row, col) at index [row * m_totalDOF + col]
				// A[target][source] format: ELF-compatible, BiCGSTAB/HACApK-optimal
				// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
				double* block = &m_flatInteractMatrix[(size_t)offset_row * m_totalDOF + offset_col];

				// Compute the interaction block based on DOF types.
				// FAST PATH: only for compact 3x3 blocks.
				if(dof_row == 3 && dof_col == 3)
				{
					// 3x3 N-matrix computation: H_field at row center from col magnetization
					// Must compute N[:, k] = H-field from unit magnetization M_k
					// This is required for wedge and other 3DOF polyhedra
					TVector3d ObsPoiVect = elem_row->ReturnCentrPoint();

					// Save original magnetization
					radTPolyhedron* poly_col = dynamic_cast<radTPolyhedron*>(elem_col);
					TVector3d orig_magn(0., 0., 0.);
					if(poly_col) {
						orig_magn = poly_col->Magn;
					}

					// Compute N-matrix columns by setting unit magnetizations
					for(int m_dir = 0; m_dir < 3; m_dir++)
					{
						// Set unit magnetization in direction m_dir
						TVector3d unit_M(0., 0., 0.);
						if(m_dir == 0) unit_M.x = 1.0;
						else if(m_dir == 1) unit_M.y = 1.0;
						else unit_M.z = 1.0;

						if(poly_col) {
							poly_col->Magn = unit_M;
						}

						// Compute H field at observation point
						radTField Field(FieldKeyInteract, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
						Field.AmOfIntrctElemWithSym = AmOfElemWithSym;
						elem_col->B_comp(&Field);

						// Store H-field as column m_dir of N-matrix (row-major storage)
						// N[i][m_dir] = H_i from unit M in direction m_dir
						block[(size_t)0 * m_totalDOF + m_dir] = Field.H.x;
						block[(size_t)1 * m_totalDOF + m_dir] = Field.H.y;
						block[(size_t)2 * m_totalDOF + m_dir] = Field.H.z;
					}

					// Restore original magnetization
					if(poly_col) {
						poly_col->Magn = orig_magn;
					}
				}
				// Non-3x3 blocks are no longer assembled by this legacy dense path.
			}
		});
		return 1;
	}

	// Independent face-coefficient blocks are no longer assembled by the legacy
	// dense relaxation path.  Mesh-backed magnetic-material solves route through
	// HDiv-VIM; fixed-M field evaluation owns the face geometry below.
	if(hasFaceChargeElements)
	{
		PrecomputeHexaGeometry();
		PrecomputeWedgeGeometry();
		return 1;
	}

	// What remains below is the symmetry-aware compact 3-component path.

	// SLOW PATH: With symmetry transformations (original code)
	for(int col = 0; col < AmOfMainElem; col++)
	{
		FillInTransPtrVectForElem(col, 'I');

		radTg3dRelax* elem_col = g3dRelaxPtrVect[col];
		int dof_col = m_elemDOF[col];
		int offset_col = m_elemDOFOffset[col];

		for(int row = 0; row < AmOfMainElem; row++)
		{
			radTg3dRelax* elem_row = g3dRelaxPtrVect[row];
			int dof_row = m_elemDOF[row];
			int offset_row = m_elemDOFOffset[row];

			// Get pointer to this block in the flattened matrix
			// ROW-MAJOR: A(row, col) at index [row * m_totalDOF + col]
			// A[target][source] format: ELF-compatible, BiCGSTAB/HACApK-optimal
			// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
			double* block = &m_flatInteractMatrix[(size_t)offset_row * m_totalDOF + offset_col];

			// Compute the interaction block based on DOF types
			// 3x3 blocks use B_comp; 6x6/5x5/cross-DOF use ComputeMixedBlockFast

			if(dof_row == 3 && dof_col == 3)
			{
				// Standard 3x3 interaction: use existing B_comp method
				TVector3d InitObsPoiVect = MainTransPtrArray[row]->TrPoint(elem_row->ReturnCentrPoint());

				TMatrix3d SubMatrix(ZeroVect, ZeroVect, ZeroVect), BufSubMatrix;
				for(unsigned i = 0; i < TransPtrVect.size(); i++)
				{
					TVector3d ObsPoiVect = TransPtrVect[i]->TrPoint_inv(InitObsPoiVect);

					radTField Field(FieldKeyInteract, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
					Field.AmOfIntrctElemWithSym = AmOfElemWithSym;

					elem_col->B_comp(&Field);

					BufSubMatrix.Str0 = Field.B;
					BufSubMatrix.Str1 = Field.H;
					BufSubMatrix.Str2 = Field.A;

					TransPtrVect[i]->TrMatrix(BufSubMatrix);
					SubMatrix += BufSubMatrix;
				}
				MainTransPtrArray[row]->TrMatrix_inv(SubMatrix);

				// Copy 3x3 matrix to flattened block (ROW-MAJOR)
				// A[i][j] at [i * stride + j] where stride = m_totalDOF
				// CRITICAL: Use size_t cast for indexing with m_totalDOF
				block[(size_t)0 * m_totalDOF + 0] = SubMatrix.Str0.x;  // (0,0)
				block[(size_t)0 * m_totalDOF + 1] = SubMatrix.Str0.y;  // (0,1)
				block[(size_t)0 * m_totalDOF + 2] = SubMatrix.Str0.z;  // (0,2)
				block[(size_t)1 * m_totalDOF + 0] = SubMatrix.Str1.x;  // (1,0)
				block[(size_t)1 * m_totalDOF + 1] = SubMatrix.Str1.y;  // (1,1)
				block[(size_t)1 * m_totalDOF + 2] = SubMatrix.Str1.z;  // (1,2)
				block[(size_t)2 * m_totalDOF + 0] = SubMatrix.Str2.x;  // (2,0)
				block[(size_t)2 * m_totalDOF + 1] = SubMatrix.Str2.y;  // (2,1)
				block[(size_t)2 * m_totalDOF + 2] = SubMatrix.Str2.z;  // (2,2)
			}
			// Mesh-backed magnetic-material blocks are no longer assembled here.
			// In this compact path dof is always 3; the defensive else zeroes any
			// unexpected block.
			else
			{
				// Unknown DOF combination - zero out the block (ROW-MAJOR)
				for(int i = 0; i < dof_row; i++)
				{
					for(int j = 0; j < dof_col; j++)
					{
						// ROW-MAJOR: A[i][j] at [i * stride + j]
						block[(size_t)i * m_totalDOF + j] = 0.0;
					}
				}
			}
		}
		EmptyTransPtrVect();
	}

	// Update formal interaction member pointers
	for(int ClNo = 0; ClNo < AmOfMainElem; ClNo++)
	{
		radTg3dRelax* g3dRelaxPtrClNo = g3dRelaxPtrVect[ClNo];
		g3dRelaxPtrVect[ClNo] = g3dRelaxPtrClNo->FormalIntrctMemberPtr();
	}

	return 1;
}

//-------------------------------------------------------------------------

void radTInteraction::SetupExternFieldArray()
{
	radTFieldKey FieldKeyExtern; FieldKeyExtern.H_=1;
	TVector3d ZeroVect(0.,0.,0.), InitObsPoiVect(0.,0.,0.), ObsPoiVect(0.,0.,0.);

	for(int k=0; k<AmOfMainElem; k++) ExternFieldArray[k] = ZeroVect;

	// Also zero m_flatExternFieldArray (always used now with unified solver)
	if(!m_flatExternFieldArray.empty())
	{
		for(size_t i = 0; i < m_flatExternFieldArray.size(); i++)
			m_flatExternFieldArray[i] = 0.0;
	}

	for(int ExtElNo=0; ExtElNo<AmOfExtElem; ExtElNo++)
	{
		FillInTransPtrVectForElem(ExtElNo, 'E');
		radTg3d* ExtElPtr = g3dExternPtrVect[ExtElNo];

		for(int StrNo=0; StrNo<AmOfMainElem; StrNo++)
		{
			InitObsPoiVect = MainTransPtrArray[StrNo]->TrPoint((g3dRelaxPtrVect[StrNo])->CentrPoint);
			TVector3d BufVect(0.,0.,0.);
			for(unsigned i=0; i<TransPtrVect.size(); i++)
			{
				TVector3d ObsPoiVect = TransPtrVect[i]->TrPoint_inv(InitObsPoiVect);
				radTField Field(FieldKeyExtern, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.); // Improve
				ExtElPtr->B_comp(&Field);
				// H field is a pseudo-vector (axial vector): H' = det(T) * T * H
				BufVect += TransPtrVect[i]->TrAxialVect(Field.H);
			}
			// Inverse pseudo-vector transformation: H_local = det(T) * T_inv * H_world
			ExternFieldArray[StrNo] += MainTransPtrArray[StrNo]->TrAxialVect_inv(BufVect);
		}
		EmptyTransPtrVect();
	}

	// Populate m_flatExternFieldArray (always used now with unified solver)
	if(!m_flatExternFieldArray.empty() && AmOfExtElem > 0)
	{
		for(int StrNo = 0; StrNo < AmOfMainElem; StrNo++)
		{
			int dof = m_elemDOF[StrNo];
			int offset = m_elemDOFOffset[StrNo];
			radTg3dRelax* elem = g3dRelaxPtrVect[StrNo];

			if(dof == 3)
			{
				// Standard element: store H_ext components
				TVector3d& H_ext = ExternFieldArray[StrNo];
				m_flatExternFieldArray[offset + 0] += H_ext.x;
				m_flatExternFieldArray[offset + 1] += H_ext.y;
				m_flatExternFieldArray[offset + 2] += H_ext.z;
			}
			// Face-coefficient elements do not store a per-face external field here.
		}
	}
	//g3dExternPtrVect.erase(g3dExternPtrVect.begin(), g3dExternPtrVect.end()); //OC240408, to enable current scaling/update
}

//-------------------------------------------------------------------------

void radTInteraction::AddExternFieldFromMoreExtSource()
{
	if(MoreExtSourceHandle.rep != 0)
	{
		radTFieldKey FieldKeyExtern; FieldKeyExtern.H_=1;
		TVector3d ZeroVect(0.,0.,0.), InitObsPoiVect(0.,0.,0.);

		for(int StrNo=0; StrNo<AmOfMainElem; StrNo++)
		{
			radTrans* ATransPtr = MainTransPtrArray[StrNo];

			InitObsPoiVect = MainTransPtrArray[StrNo]->TrPoint((g3dRelaxPtrVect[StrNo])->CentrPoint);
			radTField Field(FieldKeyExtern, CompCriterium, InitObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.); // Improve

			(static_cast<radTg3d*>(MoreExtSourceHandle.rep))->B_genComp(&Field);

			//TVector3d BufVect = ExternFieldArray[StrNo];

			// H field is a pseudo-vector: inverse transformation uses TrAxialVect_inv
			ExternFieldArray[StrNo] += MainTransPtrArray[StrNo]->TrAxialVect_inv(Field.H);
		}

		// Also populate m_flatExternFieldArray (always used now with unified solver)
		if(!m_flatExternFieldArray.empty())
		{
			for(int StrNo = 0; StrNo < AmOfMainElem; StrNo++)
			{
				int dof = m_elemDOF[StrNo];
				int offset = m_elemDOFOffset[StrNo];
				radTg3dRelax* elem = g3dRelaxPtrVect[StrNo];

				if(dof == 3)
				{
					// Standard element: store H_ext components
					TVector3d& H_ext = ExternFieldArray[StrNo];
					m_flatExternFieldArray[offset + 0] = H_ext.x;
					m_flatExternFieldArray[offset + 1] = H_ext.y;
					m_flatExternFieldArray[offset + 2] = H_ext.z;
				}
				// Face-coefficient elements do not store a per-face external field here.
			}
		}
	}
}

//-------------------------------------------------------------------------

void radTInteraction::AddMoreExternField(const radThg& hExtraExtSrc)
{
	if(hExtraExtSrc.rep == 0) return;

	radTg3d* pExtraExtSrc = static_cast<radTg3d*>(hExtraExtSrc.rep);

	radTFieldKey FieldKeyExtern; FieldKeyExtern.H_=1;
	TVector3d ZeroVect(0.,0.,0.), InitObsPoiVect(0.,0.,0.);

	for(int StrNo=0; StrNo<AmOfMainElem; StrNo++) 
	{
		radTrans* aTransPtr = MainTransPtrArray[StrNo];
		InitObsPoiVect = MainTransPtrArray[StrNo]->TrPoint((g3dRelaxPtrVect[StrNo])->CentrPoint);

		radTField Field(FieldKeyExtern, CompCriterium, InitObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.); // Improve
		pExtraExtSrc->B_genComp(&Field);

		// H field is a pseudo-vector: inverse transformation uses TrAxialVect_inv
		ExternFieldArray[StrNo] += MainTransPtrArray[StrNo]->TrAxialVect_inv(Field.H);
	}
}

//-------------------------------------------------------------------------

void radTInteraction::ZeroAuxOldArrays()
{
	if(AmOfMainElem <= 0) return;

	if(AuxOldMagnArray != nullptr)
	{
		TVector3d *tAuxOldMagn = AuxOldMagnArray;
		for(int StrNo=0; StrNo<AmOfMainElem; StrNo++) 
		{
			tAuxOldMagn->x = 0;
			tAuxOldMagn->y = 0;
			(tAuxOldMagn++)->z = 0;
		}
	}
	if(AuxOldFieldArray != nullptr)
	{
		TVector3d *tAuxOldField = AuxOldFieldArray;
		for(int StrNo=0; StrNo<AmOfMainElem; StrNo++) 
		{
			tAuxOldField->x = 0;
			tAuxOldField->y = 0;
			(tAuxOldField++)->z = 0;
		}
	}
}

//-------------------------------------------------------------------------

void radTInteraction::SubstractOldMagn()
{
	if((AuxOldMagnArray == nullptr) || (AmOfMainElem <= 0)) return;

	TVector3d *tAuxOldMagn = AuxOldMagnArray;
	for(int StNo=0; StNo<AmOfMainElem; StNo++)
	{
		TVector3d &M = (g3dRelaxPtrVect[StNo])->Magn;
		M -= *(tAuxOldMagn++); 
	}
}

//-------------------------------------------------------------------------

void radTInteraction::AddOldMagn()
{
	if((AuxOldMagnArray == nullptr) || (AmOfMainElem <= 0)) return;

	TVector3d *tAuxOldMagn = AuxOldMagnArray;
	for(int StNo=0; StNo<AmOfMainElem; StNo++)
	{
		TVector3d &M = (g3dRelaxPtrVect[StNo])->Magn;
		M += *(tAuxOldMagn++); 
	}
}

//-------------------------------------------------------------------------

double radTInteraction::CalcQuadNewOldMagnDif()
{
	if((AuxOldMagnArray == nullptr) || (AmOfMainElem <= 0)) return 0;

	double SumE2 = 0;
	TVector3d *tAuxOldMagn = AuxOldMagnArray;
	for(int StNo=0; StNo<AmOfMainElem; StNo++)
	{
		TVector3d CurDifM = (g3dRelaxPtrVect[StNo])->Magn - *(tAuxOldMagn++); 
		SumE2 += CurDifM.AmpE2(); //CurDifM*CurDifM;
	}
	return SumE2;
}

//-------------------------------------------------------------------------

void radTInteraction::FindMaxModMandH(double& MaxModM, double& MaxModH)
{
	double BufMaxModMe2, BufMaxModHe2, TestBufMaxModMe2, TestBufMaxModHe2;
	BufMaxModMe2 = BufMaxModHe2 = TestBufMaxModMe2 = TestBufMaxModHe2 = 1.E-17;

	for(int i=0; i<AmOfMainElem; i++)
	{
		TVector3d &NewMagn = NewMagnArray[i];
		TestBufMaxModMe2 = NewMagn.x*NewMagn.x + NewMagn.y*NewMagn.y + NewMagn.z*NewMagn.z;
		if(BufMaxModMe2 < TestBufMaxModMe2) BufMaxModMe2 = TestBufMaxModMe2;

		TVector3d &NewField = NewFieldArray[i];
		TestBufMaxModHe2 = NewField.x*NewField.x + NewField.y*NewField.y + NewField.z*NewField.z;
		if(BufMaxModHe2 < TestBufMaxModHe2) BufMaxModHe2 = TestBufMaxModHe2;
	}
	MaxModM = sqrt(BufMaxModMe2);
	MaxModH = sqrt(BufMaxModHe2);
}

// DumpBin / DumpBinParse / CAuxBinStrVect ctor REMOVED (Phase B2c, 2026-04-15)

//=========================================================================
// PrecomputeTetraGeometry: Pre-compute tetrahedron face geometry
// Extracts vertices, normals, and areas for fast 3x3 block computation
// Shared compact-tetrahedron geometry cache for dense interaction assembly.
//=========================================================================

void radTInteraction::PrecomputeTetraGeometry()
{
	if(m_tetraGeomReady || AmOfMainElem == 0) return;

	// Count tetrahedra and build index map (like PrecomputeHexaGeometry)
	// Works for both pure-tet and mixed meshes
	int nTet = 0;
	m_tetraElemIndices.clear();
	m_globalToTetraIdx.assign(AmOfMainElem, -1);
	for(int e = 0; e < AmOfMainElem; e++)
	{
		if(m_elemDOF[e] == 3)
		{
			m_globalToTetraIdx[e] = nTet;
			m_tetraElemIndices.push_back(e);
			nTet++;
		}
	}
	if(nTet == 0) return;

	// Allocate arrays indexed by type-specific tet index (not global element index)
	m_tetraCenters.resize(nTet * 3);
	m_tetraFaceVertices.resize(nTet * 4 * 3 * 3);
	m_tetraFaceNormals.resize(nTet * 4 * 3);
	m_tetraFaceAreas.resize(nTet * 4);

	for(int t = 0; t < nTet; t++)
	{
		int elemIdx = m_tetraElemIndices[t];
		radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[elemIdx]);
		if(!poly || poly->AmOfFaces != 4) continue;

		// Store element center
		int cIdx = t * 3;
		m_tetraCenters[cIdx + 0] = poly->CentrPoint.x;
		m_tetraCenters[cIdx + 1] = poly->CentrPoint.y;
		m_tetraCenters[cIdx + 2] = poly->CentrPoint.z;

		// Store face data for each of the 4 triangular faces
		for(int f = 0; f < 4; f++)
		{
			const radTHandlePgnAndTrans& hpt = poly->VectHandlePgnAndTrans[f];
			radTPolygon* pgn = hpt.PgnHndl.rep;
			radTrans* tr = hpt.TransHndl.rep;

			// Get 3 vertices of this triangular face
			const radTVect2dVect& verts2d = pgn->EdgePointsVector;
			if(verts2d.size() < 3) continue;

			int fvIdx = (t * 4 + f) * 3 * 3;
			TVector3d V[3];
			for(int v = 0; v < 3; v++)
			{
				V[v] = tr->TrPoint(TVector3d(verts2d[v].x, verts2d[v].y, pgn->CoordZ));
				m_tetraFaceVertices[fvIdx + v * 3 + 0] = V[v].x;
				m_tetraFaceVertices[fvIdx + v * 3 + 1] = V[v].y;
				m_tetraFaceVertices[fvIdx + v * 3 + 2] = V[v].z;
			}

			// Compute face normal (outward pointing)
			TVector3d e1 = {V[1].x - V[0].x, V[1].y - V[0].y, V[1].z - V[0].z};
			TVector3d e2 = {V[2].x - V[0].x, V[2].y - V[0].y, V[2].z - V[0].z};
			TVector3d n = {e1.y*e2.z - e1.z*e2.y, e1.z*e2.x - e1.x*e2.z, e1.x*e2.y - e1.y*e2.x};
			double nLen = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);

			// Face area = 0.5 * |cross product|
			m_tetraFaceAreas[t * 4 + f] = 0.5 * nLen;

			// Normalize and check orientation (outward from centroid)
			if(nLen > 1e-20)
			{
				n.x /= nLen; n.y /= nLen; n.z /= nLen;

				// Face center
				TVector3d fc = {(V[0].x + V[1].x + V[2].x) / 3.0,
				                (V[0].y + V[1].y + V[2].y) / 3.0,
				                (V[0].z + V[1].z + V[2].z) / 3.0};
				// Vector from centroid to face center
				TVector3d toFace = {fc.x - poly->CentrPoint.x,
				                    fc.y - poly->CentrPoint.y,
				                    fc.z - poly->CentrPoint.z};
				// If normal points inward, flip it
				if(n.x*toFace.x + n.y*toFace.y + n.z*toFace.z < 0)
				{
					n.x = -n.x; n.y = -n.y; n.z = -n.z;
				}
			}

			// Store normalized outward normal
			int fnIdx = (t * 4 + f) * 3;
			m_tetraFaceNormals[fnIdx + 0] = n.x;
			m_tetraFaceNormals[fnIdx + 1] = n.y;
			m_tetraFaceNormals[fnIdx + 2] = n.z;
		}
	}

	m_tetraGeomReady = true;
}

//=========================================================================
// FieldFromChargedTriangleLocal: Compute H field from charged triangle
// Uses analytic formula (van Oosterom & Strackee, 1983)
// Returns field WITHOUT 4pi divisor
//=========================================================================

static void FieldFromChargedTriangleLocal(const double* obs,
                                          const double* v0, const double* v1, const double* v2,
                                          double sigma, double* H_out)
{
	const double EPS = 1.0e-20;

	// Triangle edges
	double e1[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
	double e2[3] = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};

	// Orthonormal basis: a = e1/|e1|, c = e1 x e2 / |e1 x e2|, b = c x a
	double e1_len = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
	if(e1_len < EPS)
	{
		H_out[0] = H_out[1] = H_out[2] = 0.0;
		return;
	}

	double basis_a[3] = {e1[0]/e1_len, e1[1]/e1_len, e1[2]/e1_len};

	// c = e1 x e2 (normal direction)
	double c[3] = {e1[1]*e2[2] - e1[2]*e2[1],
	               e1[2]*e2[0] - e1[0]*e2[2],
	               e1[0]*e2[1] - e1[1]*e2[0]};
	double c_len = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
	if(c_len < EPS)
	{
		H_out[0] = H_out[1] = H_out[2] = 0.0;
		return;
	}
	double basis_c[3] = {c[0]/c_len, c[1]/c_len, c[2]/c_len};

	// b = c x a
	double basis_b[3] = {basis_c[1]*basis_a[2] - basis_c[2]*basis_a[1],
	                     basis_c[2]*basis_a[0] - basis_c[0]*basis_a[2],
	                     basis_c[0]*basis_a[1] - basis_c[1]*basis_a[0]};

	// Local coordinates of triangle vertices (v0 at origin in local frame)
	double xy0_x = 0.0, xy0_y = 0.0;
	double xy1_x = e1_len, xy1_y = 0.0;
	double xy2_x = e2[0]*basis_a[0] + e2[1]*basis_a[1] + e2[2]*basis_a[2];
	double xy2_y = e2[0]*basis_b[0] + e2[1]*basis_b[1] + e2[2]*basis_b[2];

	double XY[3][2] = {{xy0_x, xy0_y}, {xy1_x, xy1_y}, {xy2_x, xy2_y}};
	double DS[3], AM[3], XD[3], YD[3];
	double EPSG = 0.0;

	for(int j = 0; j < 3; j++)
	{
		int l = (j + 1) % 3;
		double dx = XY[l][0] - XY[j][0];
		double dy = XY[l][1] - XY[j][1];
		if(fabs(dx) < EPS) dx = (dx >= 0) ? EPS : -EPS;

		DS[j] = sqrt(dx*dx + dy*dy);
		AM[j] = dy / dx;
		XD[j] = -dx / DS[j];
		YD[j] =  dy / DS[j];

		if(DS[j] > EPSG) EPSG = DS[j];
	}
	EPSG *= 1.0e-12;

	// Transform observation point to local coordinates
	double d[3] = {obs[0]-v0[0], obs[1]-v0[1], obs[2]-v0[2]};
	double EE1 = d[0]*basis_a[0] + d[1]*basis_a[1] + d[2]*basis_a[2];
	double EE2 = d[0]*basis_b[0] + d[1]*basis_b[1] + d[2]*basis_b[2];
	double EE3 = d[0]*basis_c[0] + d[1]*basis_c[1] + d[2]*basis_c[2];

	double X[3], Y[3], H[3], E[3], R[3];
	for(int j = 0; j < 3; j++)
	{
		X[j] = EE1 - XY[j][0];
		Y[j] = EE2 - XY[j][1];
		H[j] = Y[j] * X[j];
		E[j] = EE3*EE3 + X[j]*X[j];
		R[j] = sqrt(X[j]*X[j] + Y[j]*Y[j] + EE3*EE3);
	}

	double Z = EE3;

	// Edge contributions
	double RM[3], RP[3], RR[3], AL[3];
	for(int j = 0; j < 3; j++)
	{
		int jp1 = (j + 1) % 3;
		RM[j] = R[j] + R[jp1] - DS[j];
		RP[j] = R[j] + R[jp1] + DS[j];
		RR[j] = (RM[j] / RP[j] > EPS) ? (RM[j] / RP[j]) : EPS;
		AL[j] = log(RR[j]);
	}

	// Field components in local frame WITHOUT 4pi divisor
	double HH1 = sigma * (-YD[0]*AL[0] - YD[1]*AL[1] - YD[2]*AL[2]);
	double HH2 = sigma * (-XD[0]*AL[0] - XD[1]*AL[1] - XD[2]*AL[2]);
	double HH3 = 0.0;

	// Normal component (atan terms)
	if(fabs(Z) > EPSG)
	{
		double ZR[3];
		for(int j = 0; j < 3; j++)
		{
			ZR[j] = Z * R[j];
		}

		double AT[3], BT[3];
		for(int j = 0; j < 3; j++)
		{
			int jp1 = (j + 1) % 3;
			AT[j] = (AM[j]*E[j] - H[j]) / ZR[j];
			BT[j] = (AM[j]*E[jp1] - H[jp1]) / ZR[jp1];
		}

		HH3 = sigma * (-atan(AT[0]) - atan(AT[1]) - atan(AT[2])
		               +atan(BT[0]) + atan(BT[1]) + atan(BT[2]));
	}

	// Transform back to global coordinates
	H_out[0] = HH1*basis_a[0] + HH2*basis_b[0] + HH3*basis_c[0];
	H_out[1] = HH1*basis_a[1] + HH2*basis_b[1] + HH3*basis_c[1];
	H_out[2] = HH1*basis_a[2] + HH2*basis_b[2] + HH3*basis_c[2];
}

//=========================================================================
// FieldGradFromChargedTriangleLocal: analytic H field AND its gradient gH = grad_obs(H) (the moment
// formulation's quadrupole field-gradient, = the Hessian -grad grad I0 of the single-layer potential)
// from a UNIFORMLY charged FLAT triangle.  H reproduces FieldFromChargedTriangleLocal exactly; gH is the
// CLOSED FORM obtained by symbolic differentiation of that H (Mathematica-verified vs the 64pt-Gauss
// moment kernel: gH rel ~ Gauss accuracy, symmetric + traceless to machine eps.  gH is
// assembled from ONLY the tangential (log-term HH1/HH2) derivatives
// plus tracelessness (Gzz = -(Gxx+Gyy)) and symmetry (Gxz = dHH1/de3, Gyz = dHH2/de3) -> NO atan
// derivative -> well-conditioned near the source plane.  Returns field+grad WITHOUT the 1/4pi factor
// (caller applies it, as for FieldFromChargedTriangleLocal).  gH order: (xx,yy,zz,xy,xz,yz).
// A degenerate (zero-area) triangle contributes nothing (c_len guard) -- so a fan-triangulated quad
// with a collapsed edge integrates the triangle automatically.
//=========================================================================
static void FieldGradFromChargedTriangleLocal(const double* obs,
                                              const double* v0, const double* v1, const double* v2,
                                              double sigma, double* H_out, double* gH_out)
{
	H_out[0] = H_out[1] = H_out[2] = 0.0;
	for(int k = 0; k < 6; k++) gH_out[k] = 0.0;

	const double EPS = 1.0e-20;

	double e1[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
	double e2[3] = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};
	double e1_len = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
	if(e1_len < EPS) return;
	double basis_a[3] = {e1[0]/e1_len, e1[1]/e1_len, e1[2]/e1_len};
	double c[3] = {e1[1]*e2[2] - e1[2]*e2[1],
	               e1[2]*e2[0] - e1[0]*e2[2],
	               e1[0]*e2[1] - e1[1]*e2[0]};
	double c_len = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
	if(c_len < EPS) return;                            // degenerate (zero-area) triangle
	double basis_c[3] = {c[0]/c_len, c[1]/c_len, c[2]/c_len};
	double basis_b[3] = {basis_c[1]*basis_a[2] - basis_c[2]*basis_a[1],
	                     basis_c[2]*basis_a[0] - basis_c[0]*basis_a[2],
	                     basis_c[0]*basis_a[1] - basis_c[1]*basis_a[0]};

	double xy2_x = e2[0]*basis_a[0] + e2[1]*basis_a[1] + e2[2]*basis_a[2];
	double xy2_y = e2[0]*basis_b[0] + e2[1]*basis_b[1] + e2[2]*basis_b[2];
	double XY[3][2] = {{0.0, 0.0}, {e1_len, 0.0}, {xy2_x, xy2_y}};
	double DS[3], AM[3], XD[3], YD[3], EPSG = 0.0;
	for(int j = 0; j < 3; j++)
	{
		int l = (j + 1) % 3;
		double dx = XY[l][0] - XY[j][0];
		double dy = XY[l][1] - XY[j][1];
		if(fabs(dx) < EPS) dx = (dx >= 0) ? EPS : -EPS;
		DS[j] = sqrt(dx*dx + dy*dy);
		AM[j] = dy / dx;
		XD[j] = -dx / DS[j];
		YD[j] =  dy / DS[j];
		if(DS[j] > EPSG) EPSG = DS[j];
	}
	EPSG *= 1.0e-12;

	double d[3] = {obs[0]-v0[0], obs[1]-v0[1], obs[2]-v0[2]};
	double EE1 = d[0]*basis_a[0] + d[1]*basis_a[1] + d[2]*basis_a[2];
	double EE2 = d[0]*basis_b[0] + d[1]*basis_b[1] + d[2]*basis_b[2];
	double EE3 = d[0]*basis_c[0] + d[1]*basis_c[1] + d[2]*basis_c[2];

	double X[3], Y[3], Hh[3], E[3], R[3];
	for(int j = 0; j < 3; j++)
	{
		X[j] = EE1 - XY[j][0];
		Y[j] = EE2 - XY[j][1];
		Hh[j] = Y[j] * X[j];
		E[j] = EE3*EE3 + X[j]*X[j];
		R[j] = sqrt(X[j]*X[j] + Y[j]*Y[j] + EE3*EE3);
	}
	double Z = EE3;
	double RM[3], RP[3], RR[3], AL[3];
	for(int j = 0; j < 3; j++)
	{
		int jp1 = (j + 1) % 3;
		RM[j] = R[j] + R[jp1] - DS[j];
		RP[j] = R[j] + R[jp1] + DS[j];
		RR[j] = (RM[j] / RP[j] > EPS) ? (RM[j] / RP[j]) : EPS;
		AL[j] = log(RR[j]);
	}
	double HH1 = sigma * (-YD[0]*AL[0] - YD[1]*AL[1] - YD[2]*AL[2]);
	double HH2 = sigma * (-XD[0]*AL[0] - XD[1]*AL[1] - XD[2]*AL[2]);
	double HH3 = 0.0;
	if(fabs(Z) > EPSG)
	{
		double AT[3], BT[3];
		for(int j = 0; j < 3; j++)
		{
			int jp1 = (j + 1) % 3;
			AT[j] = (AM[j]*E[j]   - Hh[j])   / (Z*R[j]);
			BT[j] = (AM[j]*E[jp1] - Hh[jp1]) / (Z*R[jp1]);
		}
		HH3 = sigma * (-atan(AT[0]) - atan(AT[1]) - atan(AT[2])
		               +atan(BT[0]) + atan(BT[1]) + atan(BT[2]));
	}
	H_out[0] = HH1*basis_a[0] + HH2*basis_b[0] + HH3*basis_c[0];
	H_out[1] = HH1*basis_a[1] + HH2*basis_b[1] + HH3*basis_c[1];
	H_out[2] = HH1*basis_a[2] + HH2*basis_b[2] + HH3*basis_c[2];

	// ---- local-frame gradient G (Mathematica-verified closed form, log-term derivatives only) ----
	// dR[j]/dEE = (X[j]/R[j], Y[j]/R[j], EE3/R[j]);  AL[j] = log(RM[j]/RP[j]) with
	// dRM[j]/dEE_k = dRP[j]/dEE_k = (dR[j]+dR[jp1])_k  ->  dAL[j]/dEE_k = (dR[j]+dR[jp1])_k*(1/RM-1/RP).
	double Gxx = 0.0, Gyy = 0.0, Gxy1 = 0.0, Gxy2 = 0.0, Gxz = 0.0, Gyz = 0.0;
	for(int j = 0; j < 3; j++)
	{
		int jp1 = (j + 1) % 3;
		double invRj = 1.0/R[j], invRjp = 1.0/R[jp1];
		double dR1 = X[j]*invRj + X[jp1]*invRjp;        // d(R[j]+R[jp1])/dEE1
		double dR2 = Y[j]*invRj + Y[jp1]*invRjp;        // /dEE2
		double dR3 = EE3*invRj + EE3*invRjp;            // /dEE3
		double diff = 1.0/RM[j] - 1.0/RP[j];
		double dAL1 = dR1*diff, dAL2 = dR2*diff, dAL3 = dR3*diff;
		Gxx  += -YD[j]*dAL1;     // dHH1/dEE1
		Gyy  += -XD[j]*dAL2;     // dHH2/dEE2
		Gxy1 += -YD[j]*dAL2;     // dHH1/dEE2
		Gxy2 += -XD[j]*dAL1;     // dHH2/dEE1  (== Gxy1 analytically; averaged for exact symmetry)
		Gxz  += -YD[j]*dAL3;     // dHH1/dEE3  (== dHH3/dEE1 by symmetry)
		Gyz  += -XD[j]*dAL3;     // dHH2/dEE3  (== dHH3/dEE2 by symmetry)
	}
	Gxx *= sigma; Gyy *= sigma; Gxy1 *= sigma; Gxy2 *= sigma; Gxz *= sigma; Gyz *= sigma;
	double Gxy = 0.5*(Gxy1 + Gxy2);
	double Gzz = -(Gxx + Gyy);                          // tracelessness (div H = 0 off source)
	double gl[3][3] = {{Gxx, Gxy, Gxz}, {Gxy, Gyy, Gyz}, {Gxz, Gyz, Gzz}};
	// global tensor = B * gl * B^T,  B columns = (basis_a, basis_b, basis_c)
	double B[3][3] = {{basis_a[0], basis_b[0], basis_c[0]},
	                  {basis_a[1], basis_b[1], basis_c[1]},
	                  {basis_a[2], basis_b[2], basis_c[2]}};
	double Gg[3][3];
	for(int i = 0; i < 3; i++) for(int jj = 0; jj < 3; jj++)
	{
		double s = 0.0;
		for(int a = 0; a < 3; a++) for(int b = 0; b < 3; b++) s += B[i][a]*gl[a][b]*B[jj][b];
		Gg[i][jj] = s;
	}
	gH_out[0] = Gg[0][0]; gH_out[1] = Gg[1][1]; gH_out[2] = Gg[2][2];
	gH_out[3] = Gg[0][1]; gH_out[4] = Gg[0][2]; gH_out[5] = Gg[1][2];
}

//=========================================================================
// Compute3x3BlockFast: Fast 3x3 interaction block for tetrahedra
// Uses pre-computed geometry (no B_comp overhead)
// ELF-style compact-tetrahedron field kernel used by dense interaction assembly.
//=========================================================================

void radTInteraction::Compute3x3BlockFast(int elem_i, int elem_j, double* N_mat) const
{
	std::memset(N_mat, 0, 9 * sizeof(double));

	if(!m_tetraGeomReady || elem_i < 0 || elem_i >= AmOfMainElem ||
	   elem_j < 0 || elem_j >= AmOfMainElem)
	{
		return;
	}

	// Convert global element indices to type-specific tet indices
	// m_globalToTetraIdx provides O(1) lookup (built in PrecomputeTetraGeometry)
	if(m_globalToTetraIdx.empty()) return;
	int tet_i = m_globalToTetraIdx[elem_i];
	int tet_j = m_globalToTetraIdx[elem_j];
	if(tet_i < 0 || tet_j < 0) return;

	// Observation point: center of element i
	const double* obs = &m_tetraCenters[tet_i * 3];

	// Column element center (for point charge cancellation)
	const double* col_center = &m_tetraCenters[tet_j * 3];

	// Unit magnetization vectors
	const double M_x[3] = {1.0, 0.0, 0.0};
	const double M_y[3] = {0.0, 1.0, 0.0};
	const double M_z[3] = {0.0, 0.0, 1.0};

	// Accumulate H field for each unit M direction
	double H_from_Mx[3] = {0.0, 0.0, 0.0};
	double H_from_My[3] = {0.0, 0.0, 0.0};
	double H_from_Mz[3] = {0.0, 0.0, 0.0};

	// Track total magnetic charge for centroid cancellation
	double total_charge_Mx = 0.0;
	double total_charge_My = 0.0;
	double total_charge_Mz = 0.0;

	// Process each of the 4 triangular faces
	for(int f = 0; f < 4; f++)
	{
		int fnIdx = (tet_j * 4 + f) * 3;
		const double* n_f = &m_tetraFaceNormals[fnIdx];

		// Surface charge density sigma = M dot n for each unit M
		double sigma_Mx = M_x[0]*n_f[0] + M_x[1]*n_f[1] + M_x[2]*n_f[2];
		double sigma_My = M_y[0]*n_f[0] + M_y[1]*n_f[1] + M_y[2]*n_f[2];
		double sigma_Mz = M_z[0]*n_f[0] + M_z[1]*n_f[1] + M_z[2]*n_f[2];

		// Accumulate total charge for each M direction
		double area = m_tetraFaceAreas[tet_j * 4 + f];
		total_charge_Mx += sigma_Mx * area;
		total_charge_My += sigma_My * area;
		total_charge_Mz += sigma_Mz * area;

		// Get face vertices
		int fvIdx = (tet_j * 4 + f) * 3 * 3;
		const double* V0 = &m_tetraFaceVertices[fvIdx + 0];
		const double* V1 = &m_tetraFaceVertices[fvIdx + 3];
		const double* V2 = &m_tetraFaceVertices[fvIdx + 6];

		double H_f[3];

		// H from Mx contribution
		if(fabs(sigma_Mx) > 1e-20)
		{
			FieldFromChargedTriangleLocal(obs, V0, V1, V2, sigma_Mx, H_f);
			H_from_Mx[0] += H_f[0];
			H_from_Mx[1] += H_f[1];
			H_from_Mx[2] += H_f[2];
		}

		// H from My contribution
		if(fabs(sigma_My) > 1e-20)
		{
			FieldFromChargedTriangleLocal(obs, V0, V1, V2, sigma_My, H_f);
			H_from_My[0] += H_f[0];
			H_from_My[1] += H_f[1];
			H_from_My[2] += H_f[2];
		}

		// H from Mz contribution
		if(fabs(sigma_Mz) > 1e-20)
		{
			FieldFromChargedTriangleLocal(obs, V0, V1, V2, sigma_Mz, H_f);
			H_from_Mz[0] += H_f[0];
			H_from_Mz[1] += H_f[1];
			H_from_Mz[2] += H_f[2];
		}
	}

	// Add point charge cancellation at centroid
	double r[3] = {obs[0] - col_center[0], obs[1] - col_center[1], obs[2] - col_center[2]};
	double dist_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	double dist = sqrt(dist_sq);

	if(dist > 1e-15)
	{
		double inv_dist3 = 1.0 / (dist * dist_sq);

		// Point charge H = -Q * r / |r|^3 (Q = total surface charge)
		double H_point_Mx[3] = {-total_charge_Mx * r[0] * inv_dist3,
		                        -total_charge_Mx * r[1] * inv_dist3,
		                        -total_charge_Mx * r[2] * inv_dist3};
		double H_point_My[3] = {-total_charge_My * r[0] * inv_dist3,
		                        -total_charge_My * r[1] * inv_dist3,
		                        -total_charge_My * r[2] * inv_dist3};
		double H_point_Mz[3] = {-total_charge_Mz * r[0] * inv_dist3,
		                        -total_charge_Mz * r[1] * inv_dist3,
		                        -total_charge_Mz * r[2] * inv_dist3};

		H_from_Mx[0] += H_point_Mx[0]; H_from_Mx[1] += H_point_Mx[1]; H_from_Mx[2] += H_point_Mx[2];
		H_from_My[0] += H_point_My[0]; H_from_My[1] += H_point_My[1]; H_from_My[2] += H_point_My[2];
		H_from_Mz[0] += H_point_Mz[0]; H_from_Mz[1] += H_point_Mz[1]; H_from_Mz[2] += H_point_Mz[2];
	}

	// =========== IMA: Mirrored source contributions ===========
	// Mirror source geometry inline.  M is a pseudovector: sign matrix S[beta]
	// per component.
	if(m_imaEnabled)
	{
		auto addMirrorTet = [&](int mirrorAxis, int combinedSign) {
			// Mirror source center
			double mirCenter[3] = {col_center[0], col_center[1], col_center[2]};
			if(mirrorAxis & IMA_X) mirCenter[0] = -mirCenter[0];
			if(mirrorAxis & IMA_Y) mirCenter[1] = -mirCenter[1];
			if(mirrorAxis & IMA_Z) mirCenter[2] = -mirCenter[2];

			int numMirrors = 0;
			if(mirrorAxis & IMA_X) numMirrors++;
			if(mirrorAxis & IMA_Y) numMirrors++;
			if(mirrorAxis & IMA_Z) numMirrors++;
			bool flipWinding = (numMirrors % 2 == 1);

			// Sign matrix for pseudovector M: S[k] = combinedSign, flip mirrored axes
			double S[3] = {(double)combinedSign, (double)combinedSign, (double)combinedSign};
			if(mirrorAxis & IMA_X) S[0] = -S[0];
			if(mirrorAxis & IMA_Y) S[1] = -S[1];
			if(mirrorAxis & IMA_Z) S[2] = -S[2];

			// Accumulate mirror H for each unit M direction
			double mirH_Mx[3] = {0,0,0}, mirH_My[3] = {0,0,0}, mirH_Mz[3] = {0,0,0};
			double mirCharge_Mx = 0, mirCharge_My = 0, mirCharge_Mz = 0;

			for(int f = 0; f < 4; f++)
			{
				// Copy and mirror face vertices
				int fvIdx = (tet_j * 4 + f) * 3 * 3;
				double V0[3], V1[3], V2[3];
				for(int k = 0; k < 3; k++)
				{
					V0[k] = m_tetraFaceVertices[fvIdx + 0*3 + k];
					V1[k] = m_tetraFaceVertices[fvIdx + 1*3 + k];
					V2[k] = m_tetraFaceVertices[fvIdx + 2*3 + k];
				}
				if(mirrorAxis & IMA_X) { V0[0] = -V0[0]; V1[0] = -V1[0]; V2[0] = -V2[0]; }
				if(mirrorAxis & IMA_Y) { V0[1] = -V0[1]; V1[1] = -V1[1]; V2[1] = -V2[1]; }
				if(mirrorAxis & IMA_Z) { V0[2] = -V0[2]; V1[2] = -V1[2]; V2[2] = -V2[2]; }
				if(flipWinding) {
					for(int k = 0; k < 3; k++) std::swap(V1[k], V2[k]);
				}

				// Recompute normal from mirrored vertices
				double e1[3] = {V1[0]-V0[0], V1[1]-V0[1], V1[2]-V0[2]};
				double e2[3] = {V2[0]-V0[0], V2[1]-V0[1], V2[2]-V0[2]};
				double n_f[3] = {e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]};
				double nLen = sqrt(n_f[0]*n_f[0] + n_f[1]*n_f[1] + n_f[2]*n_f[2]);
				double area = 0.5 * nLen;
				if(nLen > 1e-20) { n_f[0] /= nLen; n_f[1] /= nLen; n_f[2] /= nLen; }

				// sigma = M_unit dot n for each unit M direction
				double sigma_Mx = n_f[0], sigma_My = n_f[1], sigma_Mz = n_f[2];
				mirCharge_Mx += sigma_Mx * area;
				mirCharge_My += sigma_My * area;
				mirCharge_Mz += sigma_Mz * area;

				double H_f[3];
				if(fabs(sigma_Mx) > 1e-20) {
					FieldFromChargedTriangleLocal(obs, V0, V1, V2, sigma_Mx, H_f);
					mirH_Mx[0] += H_f[0]; mirH_Mx[1] += H_f[1]; mirH_Mx[2] += H_f[2];
				}
				if(fabs(sigma_My) > 1e-20) {
					FieldFromChargedTriangleLocal(obs, V0, V1, V2, sigma_My, H_f);
					mirH_My[0] += H_f[0]; mirH_My[1] += H_f[1]; mirH_My[2] += H_f[2];
				}
				if(fabs(sigma_Mz) > 1e-20) {
					FieldFromChargedTriangleLocal(obs, V0, V1, V2, sigma_Mz, H_f);
					mirH_Mz[0] += H_f[0]; mirH_Mz[1] += H_f[1]; mirH_Mz[2] += H_f[2];
				}
			}

			// Point charge cancellation at mirrored center
			double r[3] = {obs[0]-mirCenter[0], obs[1]-mirCenter[1], obs[2]-mirCenter[2]};
			double dist_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
			if(dist_sq > 1e-30)
			{
				double dist = sqrt(dist_sq);
				double inv_dist3 = 1.0 / (dist * dist_sq);
				mirH_Mx[0] += -mirCharge_Mx*r[0]*inv_dist3; mirH_Mx[1] += -mirCharge_Mx*r[1]*inv_dist3; mirH_Mx[2] += -mirCharge_Mx*r[2]*inv_dist3;
				mirH_My[0] += -mirCharge_My*r[0]*inv_dist3; mirH_My[1] += -mirCharge_My*r[1]*inv_dist3; mirH_My[2] += -mirCharge_My*r[2]*inv_dist3;
				mirH_Mz[0] += -mirCharge_Mz*r[0]*inv_dist3; mirH_Mz[1] += -mirCharge_Mz*r[1]*inv_dist3; mirH_Mz[2] += -mirCharge_Mz*r[2]*inv_dist3;
			}

			// Add with sign matrix: H_from_Mx[alpha] += S[0] * mirH_Mx[alpha]
			for(int a = 0; a < 3; a++) {
				H_from_Mx[a] += S[0] * mirH_Mx[a];
				H_from_My[a] += S[1] * mirH_My[a];
				H_from_Mz[a] += S[2] * mirH_Mz[a];
			}
		};

		bool hasX = (m_imaSymmetry & IMA_X) != 0;
		bool hasY = (m_imaSymmetry & IMA_Y) != 0;
		bool hasZ = (m_imaSymmetry & IMA_Z) != 0;

		if(hasX) addMirrorTet(IMA_X, m_imaSignX);
		if(hasY) addMirrorTet(IMA_Y, m_imaSignY);
		if(hasZ) addMirrorTet(IMA_Z, m_imaSignZ);
		if(hasX && hasY) addMirrorTet(IMA_XY, m_imaSignX * m_imaSignY);
		if(hasX && hasZ) addMirrorTet(IMA_XZ, m_imaSignX * m_imaSignZ);
		if(hasY && hasZ) addMirrorTet(IMA_YZ, m_imaSignY * m_imaSignZ);
		if(hasX && hasY && hasZ) addMirrorTet(IMA_XYZ, m_imaSignX * m_imaSignY * m_imaSignZ);
	}

	// Store in ROW-MAJOR format (same as flat matrix)
	// N[i][j] = H_i due to unit M_j (i = Hx,Hy,Hz; j = Mx,My,Mz)
	N_mat[0] = H_from_Mx[0] * RadConst::INV_FOUR_PI;
	N_mat[1] = H_from_My[0] * RadConst::INV_FOUR_PI;
	N_mat[2] = H_from_Mz[0] * RadConst::INV_FOUR_PI;
	N_mat[3] = H_from_Mx[1] * RadConst::INV_FOUR_PI;
	N_mat[4] = H_from_My[1] * RadConst::INV_FOUR_PI;
	N_mat[5] = H_from_Mz[1] * RadConst::INV_FOUR_PI;
	N_mat[6] = H_from_Mx[2] * RadConst::INV_FOUR_PI;
	N_mat[7] = H_from_My[2] * RadConst::INV_FOUR_PI;
	N_mat[8] = H_from_Mz[2] * RadConst::INV_FOUR_PI;
}

//=========================================================================
// PrecomputeHexaGeometry: Pre-compute hexahedron face geometry
// Hexahedra have 6 quadrilateral faces, each split into 2 triangles
// Reference path for hexahedral surface-charge collocation elements.
//=========================================================================

void radTInteraction::PrecomputeHexaGeometry()
{
	if(m_hexaGeomReady || AmOfMainElem == 0) return;

	// Count hexahedra and build index map (forward + reverse)
	int nHex = 0;
	m_hexaElemIndices.clear();
	m_globalToHexIdx.assign(AmOfMainElem, -1);
	for(int e = 0; e < AmOfMainElem; e++)
	{
		if(m_elemDOF[e] == 6)
		{
			m_globalToHexIdx[e] = nHex;
			m_hexaElemIndices.push_back(e);
			nHex++;
		}
	}
	if(nHex == 0) return;

	// Allocate arrays
	m_hexaCenters.resize(nHex * 3);
	m_hexaFaceNormals.resize(nHex * 6 * 3); // 6 faces, xyz
	m_hexaFaceAreas.resize(nHex * 6);       // 6 faces
	m_hexaTriVertices.resize(nHex * 6 * 2 * 3 * 3);  // 6 faces, 2 tris, 3 verts, xyz
	m_hexaTriSigns.resize(nHex * 6 * 2);    // 6 faces, 2 tris


	for(int h = 0; h < nHex; h++)
	{
		int elemIdx = m_hexaElemIndices[h];
		radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[elemIdx]);
		if(!poly || poly->AmOfFaces != 6) continue;

		// Store element center
		int cIdx = h * 3;
		m_hexaCenters[cIdx + 0] = poly->CentrPoint.x;
		m_hexaCenters[cIdx + 1] = poly->CentrPoint.y;
		m_hexaCenters[cIdx + 2] = poly->CentrPoint.z;

		// Process each of the 6 faces
		for(int f = 0; f < 6; f++)
		{
			// Store face normal (already computed in poly->FaceNormal)
			int fnIdx = (h * 6 + f) * 3;
			m_hexaFaceNormals[fnIdx + 0] = poly->FaceNormal[f].x;
			m_hexaFaceNormals[fnIdx + 1] = poly->FaceNormal[f].y;
			m_hexaFaceNormals[fnIdx + 2] = poly->FaceNormal[f].z;

			// Store face area
			m_hexaFaceAreas[h * 6 + f] = poly->FaceArea[f];

			// Get face vertices and split into 2 triangles.  Hexahedron faces use the stored
			// REAL vertices (exact 2-triangle geometry for planar and non-planar quads alike);
			// faces without stored verts fall back to the flattened-polygon reconstruction.
			TVector3d V[4];
			TVector3d RV[4]; int rnv = 0;
			if(poly->GetRealFaceVerts(f, RV, rnv) && rnv >= 4)
			{
				V[0] = RV[0]; V[1] = RV[1]; V[2] = RV[2]; V[3] = RV[3];
			}
			else
			{
				const radTHandlePgnAndTrans& hpt = poly->VectHandlePgnAndTrans[f];
				radTPolygon* pgn = hpt.PgnHndl.rep;
				radTrans* tr = hpt.TransHndl.rep;

				const radTVect2dVect& verts2d = pgn->EdgePointsVector;
				if(verts2d.size() < 4) continue;

				for(int v = 0; v < 4; v++)
					V[v] = tr->TrPoint(TVector3d(verts2d[v].x, verts2d[v].y, pgn->CoordZ));
			}

			// Triangle 1: V0, V1, V2
			// Triangle 2: V0, V2, V3
			TVector3d tri_verts[2][3] = {
				{V[0], V[1], V[2]},
				{V[0], V[2], V[3]}
			};

			for(int t = 0; t < 2; t++)
			{
				const TVector3d& T0 = tri_verts[t][0];
				const TVector3d& T1 = tri_verts[t][1];
				const TVector3d& T2 = tri_verts[t][2];

				// Store triangle vertices
				int tvIdx = ((h * 6 + f) * 2 + t) * 3 * 3;
				m_hexaTriVertices[tvIdx + 0] = T0.x;
				m_hexaTriVertices[tvIdx + 1] = T0.y;
				m_hexaTriVertices[tvIdx + 2] = T0.z;
				m_hexaTriVertices[tvIdx + 3] = T1.x;
				m_hexaTriVertices[tvIdx + 4] = T1.y;
				m_hexaTriVertices[tvIdx + 5] = T1.z;
				m_hexaTriVertices[tvIdx + 6] = T2.x;
				m_hexaTriVertices[tvIdx + 7] = T2.y;
				m_hexaTriVertices[tvIdx + 8] = T2.z;

				// Compute triangle normal and sign correction
				TVector3d e1 = {T1.x - T0.x, T1.y - T0.y, T1.z - T0.z};
				TVector3d e2 = {T2.x - T0.x, T2.y - T0.y, T2.z - T0.z};
				TVector3d tri_n = {e1.y*e2.z - e1.z*e2.y,
				                   e1.z*e2.x - e1.x*e2.z,
				                   e1.x*e2.y - e1.y*e2.x};
				double nLen = sqrt(tri_n.x*tri_n.x + tri_n.y*tri_n.y + tri_n.z*tri_n.z);

				double sign = 1.0;
				if(nLen > 1e-20)
				{
					tri_n.x /= nLen; tri_n.y /= nLen; tri_n.z /= nLen;

					// Triangle center
					TVector3d tc = {(T0.x + T1.x + T2.x) / 3.0,
					                (T0.y + T1.y + T2.y) / 3.0,
					                (T0.z + T1.z + T2.z) / 3.0};
					// Vector from element center to triangle center
					TVector3d toTri = {tc.x - poly->CentrPoint.x,
					                   tc.y - poly->CentrPoint.y,
					                   tc.z - poly->CentrPoint.z};
					// Sign: +1 if normal points outward, -1 if inward
					double dot = tri_n.x*toTri.x + tri_n.y*toTri.y + tri_n.z*toTri.z;
					sign = (dot >= 0) ? 1.0 : -1.0;
				}
				m_hexaTriSigns[(h * 6 + f) * 2 + t] = sign;
			}
		}

	}

	m_hexaGeomReady = true;

	// Also pre-compute triangle local coordinate systems
	const_cast<radTInteraction*>(this)->PrecomputeHexaTriangleData();
}

//=========================================================================
// BuildFaceGeom: per-DOF hex face geometry in the matrix DOF order.
// Row-major (m_totalDOF x 11): [elem_local, area, cx,cy,cz, nx,ny,nz(outward), ecx,ecy,ecz].
// DOF<->face mapping (DOF = m_elemDOFOffset[elem] + f) aligns the rows 1:1 with
// GetInteractMatrix.  Lets Python form the div(B)=0 constraint (Sum_f area_f*sigma_f
// = 0 per element), the uniform-field RHS (n.H), and the dipole moment (Sum_f sigma_f*area_f*(c_f-c_e)).
//=========================================================================

void radTInteraction::BuildFaceGeom(std::vector<double>& Gflat) const
{
	const int STRIDE = 11;
	Gflat.assign((size_t)m_totalDOF * STRIDE, 0.0);
	for(int d = 0; d < m_totalDOF; d++) Gflat[(size_t)d * STRIDE + 0] = -1.0;  // elem_local default -1 (non-hex)

	int nHex = (int)m_hexaElemIndices.size();
	for(int h = 0; h < nHex; h++)
	{
		int elemIdx = m_hexaElemIndices[h];
		radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[elemIdx]);
		if(!poly || poly->AmOfFaces != 6) continue;
		int off = m_elemDOFOffset[elemIdx];
		const TVector3d ec = poly->CentrPoint;
		for(int f = 0; f < 6; f++)
		{
			const radTHandlePgnAndTrans& hpt = poly->VectHandlePgnAndTrans[f];
			radTPolygon* pgn = hpt.PgnHndl.rep;
			radTrans* tr = hpt.TransHndl.rep;
			if(!pgn || !tr) continue;
			const radTVect2dVect& v2d = pgn->EdgePointsVector;
			if(v2d.size() < 4) continue;
			TVector3d V4[4];
			for(int v = 0; v < 4; v++)
				V4[v] = tr->TrPoint(TVector3d(v2d[v].x, v2d[v].y, pgn->CoordZ));
			// quad area via two triangles (V0V1V2 + V0V2V3)
			double area = 0.0;
			for(int t = 0; t < 2; t++)
			{
				const TVector3d& A = V4[0]; const TVector3d& B = V4[t + 1]; const TVector3d& C = V4[t + 2];
				double ux = B.x - A.x, uy = B.y - A.y, uz = B.z - A.z;
				double wx = C.x - A.x, wy = C.y - A.y, wz = C.z - A.z;
				double rx = uy * wz - uz * wy, ry = uz * wx - ux * wz, rz = ux * wy - uy * wx;
				area += 0.5 * std::sqrt(rx * rx + ry * ry + rz * rz);
			}
			double cx = 0.25 * (V4[0].x + V4[1].x + V4[2].x + V4[3].x);
			double cy = 0.25 * (V4[0].y + V4[1].y + V4[2].y + V4[3].y);
			double cz = 0.25 * (V4[0].z + V4[1].z + V4[2].z + V4[3].z);
			// outward unit normal: stored FaceNormal, flipped to point away from the element center
			TVector3d nrm = poly->FaceNormal[f];
			double nlen = std::sqrt(nrm.x*nrm.x + nrm.y*nrm.y + nrm.z*nrm.z);
			if(nlen > 1e-20) { nrm.x /= nlen; nrm.y /= nlen; nrm.z /= nlen; }
			double outdot = nrm.x*(cx - ec.x) + nrm.y*(cy - ec.y) + nrm.z*(cz - ec.z);
			if(outdot < 0.0) { nrm.x = -nrm.x; nrm.y = -nrm.y; nrm.z = -nrm.z; }

			size_t b = (size_t)(off + f) * STRIDE;
			Gflat[b + 0] = (double)h;
			Gflat[b + 1] = area;
			Gflat[b + 2] = cx; Gflat[b + 3] = cy; Gflat[b + 4] = cz;
			Gflat[b + 5] = nrm.x; Gflat[b + 6] = nrm.y; Gflat[b + 7] = nrm.z;
			Gflat[b + 8] = ec.x; Gflat[b + 9] = ec.y; Gflat[b + 10] = ec.z;
		}
	}
}

//=========================================================================
// PrecomputeHexaTriangleData: Pre-compute triangle local coordinate systems
// Eliminates redundant sqrt/div operations during field computation
// Each hexahedron has 12 triangles (6 faces * 2 triangles)
//=========================================================================

void radTInteraction::PrecomputeHexaTriangleData()
{
	if(m_hexaTriDataReady || !m_hexaGeomReady) return;

	const double EPS = 1.0e-20;
	int nHex = (int)m_hexaElemIndices.size();
	if(nHex == 0) return;

	// Allocate: 12 triangles per hex, TRI_DATA_SIZE doubles per triangle
	m_hexaTriData.resize(nHex * TRIS_PER_HEX_ELEM * TRI_DATA_SIZE);

	for(int h = 0; h < nHex; h++)
	{
		const double* center = &m_hexaCenters[h * 3];

		for(int f = 0; f < 6; f++)
		{
			for(int t = 0; t < 2; t++)
			{
				int tri_idx = h * TRIS_PER_HEX_ELEM + f * 2 + t;
				double* data = &m_hexaTriData[tri_idx * TRI_DATA_SIZE];

				// Get triangle vertices from m_hexaTriVertices
				int tvIdx = ((h * 6 + f) * 2 + t) * 3 * 3;
				const double* v0 = &m_hexaTriVertices[tvIdx + 0];
				const double* v1 = &m_hexaTriVertices[tvIdx + 3];
				const double* v2 = &m_hexaTriVertices[tvIdx + 6];

				// Build local coordinate system
				double e1[3] = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]};
				double e2[3] = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]};

				// basis_c = e1 x e2 (face normal)
				double* basis_c = data + 6;
				basis_c[0] = e1[1]*e2[2] - e1[2]*e2[1];
				basis_c[1] = e1[2]*e2[0] - e1[0]*e2[2];
				basis_c[2] = e1[0]*e2[1] - e1[1]*e2[0];

				double cLen = sqrt(basis_c[0]*basis_c[0] + basis_c[1]*basis_c[1] + basis_c[2]*basis_c[2]);
				if(cLen < EPS) {
					std::memset(data, 0, TRI_DATA_SIZE * sizeof(double));
					continue;
				}
				basis_c[0] /= cLen; basis_c[1] /= cLen; basis_c[2] /= cLen;

				// basis_a = e1 normalized
				double* basis_a = data;
				basis_a[0] = e1[0]; basis_a[1] = e1[1]; basis_a[2] = e1[2];
				double aLen = sqrt(basis_a[0]*basis_a[0] + basis_a[1]*basis_a[1] + basis_a[2]*basis_a[2]);
				if(aLen < EPS) {
					std::memset(data, 0, TRI_DATA_SIZE * sizeof(double));
					continue;
				}
				basis_a[0] /= aLen; basis_a[1] /= aLen; basis_a[2] /= aLen;

				// basis_b = basis_c x basis_a
				double* basis_b = data + 3;
				basis_b[0] = basis_c[1]*basis_a[2] - basis_c[2]*basis_a[1];
				basis_b[1] = basis_c[2]*basis_a[0] - basis_c[0]*basis_a[2];
				basis_b[2] = basis_c[0]*basis_a[1] - basis_c[1]*basis_a[0];

				// origin = v0
				double* origin = data + 9;
				origin[0] = v0[0]; origin[1] = v0[1]; origin[2] = v0[2];

				// 2D coordinates (v0 = origin)
				double* XY = data + 12;  // 6 doubles: {x0,y0, x1,y1, x2,y2}
				XY[0] = 0.0; XY[1] = 0.0;  // v0
				XY[2] = e1[0]*basis_a[0] + e1[1]*basis_a[1] + e1[2]*basis_a[2];
				XY[3] = e1[0]*basis_b[0] + e1[1]*basis_b[1] + e1[2]*basis_b[2];
				XY[4] = e2[0]*basis_a[0] + e2[1]*basis_a[1] + e2[2]*basis_a[2];
				XY[5] = e2[0]*basis_b[0] + e2[1]*basis_b[1] + e2[2]*basis_b[2];

				// Edge parameters
				double* DS = data + 18;  // 3 doubles
				double* AM = data + 21;  // 3 doubles
				double* XD = data + 24;  // 3 doubles
				double* YD = data + 27;  // 3 doubles
				double EPSG = 0.0;

				for(int j = 0; j < 3; j++)
				{
					int l = (j + 1) % 3;
					double dx = XY[l*2] - XY[j*2];
					double dy = XY[l*2+1] - XY[j*2+1];
					if(fabs(dx) < EPS) dx = (dx >= 0) ? EPS : -EPS;

					DS[j] = sqrt(dx*dx + dy*dy);
					AM[j] = dy / dx;
					XD[j] = -dx / DS[j];
					YD[j] = dy / DS[j];

					if(DS[j] > EPSG) EPSG = DS[j];
				}

				// Store EPSG and sign
				data[30] = EPSG * 1.0e-12;  // EPSG

				// Get pre-computed sign from m_hexaTriSigns
				data[31] = m_hexaTriSigns[(h * 6 + f) * 2 + t];
			}
		}
	}

	m_hexaTriDataReady = true;
}

//=========================================================================
// FieldFromTrianglePrecomputed: Ultra-fast field using pre-computed data
// Uses pre-computed basis vectors and edge parameters (no sqrt/div)
// Returns field WITHOUT 4pi divisor
//=========================================================================

void radTInteraction::FieldFromTrianglePrecomputed(int hex_idx, int tri_idx, const double* obs, double sigma, double* H_out) const
{
	const double EPS = 1.0e-20;

	// Get pre-computed data
	int global_tri_idx = hex_idx * TRIS_PER_HEX_ELEM + tri_idx;
	const double* data = &m_hexaTriData[global_tri_idx * TRI_DATA_SIZE];

	const double* basis_a = data;       // [0..2]
	const double* basis_b = data + 3;   // [3..5]
	const double* basis_c = data + 6;   // [6..8]
	const double* origin = data + 9;    // [9..11]
	const double* XY = data + 12;       // [12..17] = {x0,y0, x1,y1, x2,y2}
	const double* DS = data + 18;       // [18..20]
	const double* AM = data + 21;       // [21..23]
	const double* XD = data + 24;       // [24..26]
	const double* YD = data + 27;       // [27..29]
	double EPSG = data[30];
	double sign = data[31];

	// Apply sign to sigma
	sigma *= sign;

	// Transform observation point to local coordinates
	double d[3] = {obs[0] - origin[0], obs[1] - origin[1], obs[2] - origin[2]};
	double EE1 = d[0]*basis_a[0] + d[1]*basis_a[1] + d[2]*basis_a[2];
	double EE2 = d[0]*basis_b[0] + d[1]*basis_b[1] + d[2]*basis_b[2];
	double EE3 = d[0]*basis_c[0] + d[1]*basis_c[1] + d[2]*basis_c[2];

	// Vertex-relative coordinates
	double X[3], Y[3], H[3], E[3], R[3];
	for(int j = 0; j < 3; j++)
	{
		X[j] = EE1 - XY[j*2];
		Y[j] = EE2 - XY[j*2+1];
		H[j] = Y[j] * X[j];
		E[j] = EE3*EE3 + X[j]*X[j];
		R[j] = sqrt(X[j]*X[j] + Y[j]*Y[j] + EE3*EE3);
	}

	double Z = EE3;

	// Edge contributions (log terms)
	double AL[3];
	for(int j = 0; j < 3; j++)
	{
		int jp1 = (j + 1) % 3;
		double RM = R[j] + R[jp1] - DS[j];
		double RP = R[j] + R[jp1] + DS[j];
		double RR = (RM / RP > EPS) ? (RM / RP) : EPS;
		AL[j] = log(RR);
	}

	// Tangential field components
	double HH1 = sigma * (-YD[0]*AL[0] - YD[1]*AL[1] - YD[2]*AL[2]);
	double HH2 = sigma * (-XD[0]*AL[0] - XD[1]*AL[1] - XD[2]*AL[2]);
	double HH3 = 0.0;

	// Normal component (atan terms) - only if not on surface
	if(fabs(Z) > EPSG)
	{
		double ZR[3];
		for(int j = 0; j < 3; j++)
		{
			ZR[j] = Z * R[j];
		}

		double AT[3], BT[3];
		for(int j = 0; j < 3; j++)
		{
			int jp1 = (j + 1) % 3;
			AT[j] = (AM[j]*E[j] - H[j]) / ZR[j];
			BT[j] = (AM[j]*E[jp1] - H[jp1]) / ZR[jp1];
		}

		HH3 = sigma * (-atan(AT[0]) - atan(AT[1]) - atan(AT[2])
		               +atan(BT[0]) + atan(BT[1]) + atan(BT[2]));
	}

	// Transform back to global coordinates
	H_out[0] = HH1*basis_a[0] + HH2*basis_b[0] + HH3*basis_c[0];
	H_out[1] = HH1*basis_a[1] + HH2*basis_b[1] + HH3*basis_c[1];
	H_out[2] = HH1*basis_a[2] + HH2*basis_b[2] + HH3*basis_c[2];
}

//=========================================================================
// PrecomputeWedgeGeometry: Pre-compute wedge face geometry
// Wedge: 5 faces (2 triangular + 3 quadrilateral)
// Each tri face -> 1 triangle, each quad face -> 2 triangles
// Total: up to 8 triangles per wedge (same approach as hex)
//=========================================================================

void radTInteraction::PrecomputeWedgeGeometry()
{
	if(m_wedgeGeomReady || AmOfMainElem == 0) return;

	int nWedge = 0;
	m_wedgeElemIndices.clear();
	m_globalToWedgeIdx.assign(AmOfMainElem, -1);
	for(int e = 0; e < AmOfMainElem; e++)
	{
		if(m_elemDOF[e] == 5)
		{
			m_globalToWedgeIdx[e] = nWedge;
			m_wedgeElemIndices.push_back(e);
			nWedge++;
		}
	}
	if(nWedge == 0) return;

	m_wedgeCenters.resize(nWedge * 3);
	m_wedgeFaceNormals.resize(nWedge * 5 * 3);
	m_wedgeFaceAreas.resize(nWedge * 5);
	m_wedgeFaceNumTris.resize(nWedge * 5);
	m_wedgeTriOffset.resize(nWedge * 5);
	m_wedgeTriVertices.resize(nWedge * WEDGE_MAX_TRIS * 3 * 3);
	m_wedgeTriSigns.resize(nWedge * WEDGE_MAX_TRIS);


	for(int w = 0; w < nWedge; w++)
	{
		int elemIdx = m_wedgeElemIndices[w];
		radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[elemIdx]);
		if(!poly || poly->AmOfFaces != 5) continue;

		int cIdx = w * 3;
		m_wedgeCenters[cIdx + 0] = poly->CentrPoint.x;
		m_wedgeCenters[cIdx + 1] = poly->CentrPoint.y;
		m_wedgeCenters[cIdx + 2] = poly->CentrPoint.z;

		int triCount = 0;
		for(int f = 0; f < 5; f++)
		{
			m_wedgeFaceNormals[(w*5+f)*3+0] = poly->FaceNormal[f].x;
			m_wedgeFaceNormals[(w*5+f)*3+1] = poly->FaceNormal[f].y;
			m_wedgeFaceNormals[(w*5+f)*3+2] = poly->FaceNormal[f].z;
			m_wedgeFaceAreas[w*5+f] = poly->FaceArea[f];
			const radTHandlePgnAndTrans& hpt = poly->VectHandlePgnAndTrans[f];
			radTPolygon* pgn = hpt.PgnHndl.rep;
			radTrans* tr = hpt.TransHndl.rep;
			const radTVect2dVect& verts2d = pgn->EdgePointsVector;
			int nv = (int)verts2d.size();

			m_wedgeTriOffset[w*5+f] = triCount;

			// Get face vertices and split into triangles. Stored real vertices keep warped
			// quad geometry exact; triangular faces fall through to the flattened-polygon
			// reconstruction.
			TVector3d V[4];
			TVector3d RV[4]; int rnv = 0;
			if(poly->GetRealFaceVerts(f, RV, rnv) && rnv >= 4)
			{
				nv = 4;
				V[0] = RV[0]; V[1] = RV[1]; V[2] = RV[2]; V[3] = RV[3];
			}
			else
			{
				for(int v = 0; v < nv && v < 4; v++)
					V[v] = tr->TrPoint(TVector3d(verts2d[v].x, verts2d[v].y, pgn->CoordZ));
			}

			int numTris = (nv == 3) ? 1 : 2;
			m_wedgeFaceNumTris[w*5+f] = numTris;

			// Build triangles: nv==3 -> 1 tri (V0,V1,V2), nv==4 -> 2 tris (V0,V1,V2) + (V0,V2,V3)
			TVector3d tri_v[2][3];
			tri_v[0][0] = V[0]; tri_v[0][1] = V[1]; tri_v[0][2] = V[2];
			if(numTris == 2) { tri_v[1][0] = V[0]; tri_v[1][1] = V[2]; tri_v[1][2] = V[3]; }

			for(int t = 0; t < numTris; t++)
			{
				int tvIdx = (w * WEDGE_MAX_TRIS + triCount) * 3 * 3;
				for(int vi = 0; vi < 3; vi++)
				{
					m_wedgeTriVertices[tvIdx + vi*3 + 0] = tri_v[t][vi].x;
					m_wedgeTriVertices[tvIdx + vi*3 + 1] = tri_v[t][vi].y;
					m_wedgeTriVertices[tvIdx + vi*3 + 2] = tri_v[t][vi].z;
				}

				TVector3d e1 = {tri_v[t][1].x-tri_v[t][0].x, tri_v[t][1].y-tri_v[t][0].y, tri_v[t][1].z-tri_v[t][0].z};
				TVector3d e2 = {tri_v[t][2].x-tri_v[t][0].x, tri_v[t][2].y-tri_v[t][0].y, tri_v[t][2].z-tri_v[t][0].z};
				TVector3d tn = {e1.y*e2.z-e1.z*e2.y, e1.z*e2.x-e1.x*e2.z, e1.x*e2.y-e1.y*e2.x};
				double nLen = sqrt(tn.x*tn.x + tn.y*tn.y + tn.z*tn.z);
				double sign = 1.0;
				if(nLen > 1e-20)
				{
					tn.x /= nLen; tn.y /= nLen; tn.z /= nLen;
					TVector3d tc = {(tri_v[t][0].x+tri_v[t][1].x+tri_v[t][2].x)/3.0,
					                (tri_v[t][0].y+tri_v[t][1].y+tri_v[t][2].y)/3.0,
					                (tri_v[t][0].z+tri_v[t][1].z+tri_v[t][2].z)/3.0};
					TVector3d toTri = {tc.x-poly->CentrPoint.x, tc.y-poly->CentrPoint.y, tc.z-poly->CentrPoint.z};
					sign = (tn.x*toTri.x + tn.y*toTri.y + tn.z*toTri.z >= 0) ? 1.0 : -1.0;
				}
				m_wedgeTriSigns[w * WEDGE_MAX_TRIS + triCount] = sign;
				triCount++;
			}
		}
	}


	m_wedgeGeomReady = true;
}

// Legacy dense face-coefficient block kernels were deleted.  Mesh-backed
// magnetic-material solves route through HDiv-VIM.

//-------------------------------------------------------------------------
// SetIMASymmetry: Configure IMA symmetry mode
// symmetry: IMA_X, IMA_Y, IMA_Z, etc.
// signX, signY, signZ: +1 for symmetric BC, -1 for antisymmetric BC per axis
// Returns: number of elements in IMA region
//-------------------------------------------------------------------------
int radTInteraction::SetIMASymmetry(int symmetry, int signX, int signY, int signZ)
{
	m_imaSymmetry = symmetry;
	m_imaSignX = (signX >= 0) ? 1 : -1;  // Normalize to +1 or -1
	m_imaSignY = (signY >= 0) ? 1 : -1;
	m_imaSignZ = (signZ >= 0) ? 1 : -1;
	m_imaEnabled = (symmetry != IMA_NONE);

	if(!m_imaEnabled)
	{
		m_imaNumElements = AmOfMainElem;
		m_imaToFull.clear();
		m_imaMirrorMap.clear();
		m_imaUseVirtualMirror.clear();
		return m_imaNumElements;
	}

	// Build mapping of elements in IMA region
	// Elements in IMA region have positive coordinates for each active symmetry axis
	m_imaToFull.clear();
	m_imaMirrorMap.clear();
	m_imaUseVirtualMirror.clear();

	for(int i = 0; i < AmOfMainElem; i++)
	{
		if(IsElementInIMARegion(i))
		{
			m_imaToFull.push_back(i);
		}
	}

	m_imaNumElements = (int)m_imaToFull.size();

	// Build mirror map: for each IMA element, find its mirror in full model
	// For quarter models, there may be no physical mirror - use virtual mirror
	m_imaMirrorMap.resize(m_imaNumElements);
	m_imaUseVirtualMirror.resize(m_imaNumElements, false);

	const double eps = 1e-6;  // Tolerance for matching
	int numVirtualMirrors = 0;

	for(int ima_i = 0; ima_i < m_imaNumElements; ima_i++)
	{
		int full_i = m_imaToFull[ima_i];
		TVector3d center = g3dRelaxPtrVect[full_i]->ReturnCentrPoint();

		// Compute mirrored center
		TVector3d mirror_center = center;
		if(m_imaSymmetry & IMA_X) mirror_center.x = -center.x;
		if(m_imaSymmetry & IMA_Y) mirror_center.y = -center.y;
		if(m_imaSymmetry & IMA_Z) mirror_center.z = -center.z;

		// Find element closest to mirror_center
		double min_dist_sq = 1e30;
		int best_match = -1;

		for(int i = 0; i < AmOfMainElem; i++)
		{
			TVector3d c = g3dRelaxPtrVect[i]->ReturnCentrPoint();
			double dx = c.x - mirror_center.x;
			double dy = c.y - mirror_center.y;
			double dz = c.z - mirror_center.z;
			double dist_sq = dx*dx + dy*dy + dz*dz;

			if(dist_sq < min_dist_sq)
			{
				min_dist_sq = dist_sq;
				best_match = i;
			}
		}

		// Check if we found a physical mirror element
		if(min_dist_sq > eps * eps)
		{
			// No physical mirror found - use virtual mirror
			m_imaMirrorMap[ima_i] = full_i;  // Store self for reference
			m_imaUseVirtualMirror[ima_i] = true;
			numVirtualMirrors++;
		}
		else
		{
			// Physical mirror found
			m_imaMirrorMap[ima_i] = best_match;
			m_imaUseVirtualMirror[ima_i] = false;
		}
	}

	// IMA symmetry enabled (silent in production)

	return m_imaNumElements;
}

//-------------------------------------------------------------------------
// IsElementInIMARegion: Check if element center is in positive half-space
// for all active symmetry axes
//-------------------------------------------------------------------------
bool radTInteraction::IsElementInIMARegion(int elemIdx) const
{
	if(elemIdx < 0 || elemIdx >= AmOfMainElem) return false;

	TVector3d center = g3dRelaxPtrVect[elemIdx]->ReturnCentrPoint();

	// For x-mirror: element must have x >= 0 (or x > -epsilon to handle centerline elements)
	const double eps = 1e-10;

	if(m_imaSymmetry & IMA_X)
	{
		if(center.x < -eps) return false;
	}
	if(m_imaSymmetry & IMA_Y)
	{
		if(center.y < -eps) return false;
	}
	if(m_imaSymmetry & IMA_Z)
	{
		if(center.z < -eps) return false;
	}

	return true;
}

//-------------------------------------------------------------------------
// GetMirrorElementIndex: Find element that is mirror image of given element
// For x-mirror: find element with center at (-x, y, z)
//-------------------------------------------------------------------------
int radTInteraction::GetMirrorElementIndex(int elemIdx, int symmetryAxis) const
{
	if(elemIdx < 0 || elemIdx >= AmOfMainElem) return -1;

	TVector3d center = g3dRelaxPtrVect[elemIdx]->ReturnCentrPoint();

	// Compute mirrored center
	TVector3d mirror_center = center;
	if(symmetryAxis & IMA_X) mirror_center.x = -center.x;
	if(symmetryAxis & IMA_Y) mirror_center.y = -center.y;
	if(symmetryAxis & IMA_Z) mirror_center.z = -center.z;

	// Find element closest to mirror_center
	const double eps = 1e-6;  // Tolerance for matching
	double min_dist_sq = 1e30;
	int best_match = -1;

	for(int i = 0; i < AmOfMainElem; i++)
	{
		TVector3d c = g3dRelaxPtrVect[i]->ReturnCentrPoint();
		double dx = c.x - mirror_center.x;
		double dy = c.y - mirror_center.y;
		double dz = c.z - mirror_center.z;
		double dist_sq = dx*dx + dy*dy + dz*dz;

		if(dist_sq < min_dist_sq)
		{
			min_dist_sq = dist_sq;
			best_match = i;
		}
	}

	// Verify match is close enough
	if(min_dist_sq > eps * eps)
	{
		return -1;  // No physical mirror found
	}

	return best_match;
}

// Legacy ELF face-permutation + IMA-mirror block (ApplyDOFPermutation, ApplyRowPermutation,
// Compute6x6BlockIMA, Compute6x6BlockMirrored, Compute6x6BlockMirroredTarget) deleted.

//-------------------------------------------------------------------------
// SetupInteractMatrix_IMA: Build IMA interaction matrix.
// Compact 3-component elements use Compute3x3BlockFast.  Mesh-backed
// magnetic-material solves route through HDiv-VIM.
//-------------------------------------------------------------------------
int radTInteraction::SetupInteractMatrix_IMA(bool skipDenseMatrix)
{
	if(!m_imaEnabled)
	{
		std::cerr << "[Radia] Error: IMA not enabled" << std::endl;
		return 0;
	}


	// Check all elements have valid DOF.
	bool allHex = true;
	bool allTet = true;
	bool allWedge = true;
	for(int i = 0; i < AmOfMainElem; i++)
	{
		if(m_elemDOF[i] < 3)
		{
			std::cerr << "[Radia] Error: IMA requires elements with DOF >= 3" << std::endl;
			return 0;
		}
		if(m_elemDOF[i] != 6) allHex = false;
		if(m_elemDOF[i] != 5) allWedge = false;
		if(m_elemDOF[i] != 3) allTet = false;
	}
	// Any independent face-coefficient element (DOF 4/5/6) present?
	bool hasFaceCharge = !allTet;

	// Pre-compute geometry for fast path (all element types, including mixed meshes)
	if(!m_hexaGeomReady)
		PrecomputeHexaGeometry();
	if(!m_wedgeGeomReady)
		PrecomputeWedgeGeometry();
	if(!m_tetraGeomReady)
		PrecomputeTetraGeometry();

	// Compute IMA DOF count from actual element DOFs
	int imaDOF = 0;
	for(int ima_i = 0; ima_i < m_imaNumElements; ima_i++)
	{
		int full_i = m_imaToFull[ima_i];
		imaDOF += m_elemDOF[full_i];
	}

	// Build IMA interaction matrix (imaDOF x imaDOF)

	// Allocate IMA matrix (skip for HACApK - kernel computes entries on demand)
	if(!skipDenseMatrix)
	{
		size_t matrix_size = (size_t)imaDOF * (size_t)imaDOF;
		try {
			m_flatInteractMatrix.resize(matrix_size, 0.0);
		} catch(const std::bad_alloc&) {
			std::cerr << "[Radia] Error: Memory allocation failed for IMA matrix" << std::endl;
			return 0;
		}
	}

	// Save original values
	int originalAmOfMainElem = AmOfMainElem;
	std::vector<int> originalElemDOF = m_elemDOF;
	std::vector<int> originalElemDOFOffset = m_elemDOFOffset;
	radTVectPtrg3dRelax originalG3dRelaxPtrVect = g3dRelaxPtrVect;

	// Save external field values for IMA elements BEFORE resizing
	std::vector<double> savedExternField(imaDOF, 0.0);
	if(!m_flatExternFieldArray.empty())
	{
		int ima_offset = 0;
		for(int ima_i = 0; ima_i < m_imaNumElements; ima_i++)
		{
			int full_i = m_imaToFull[ima_i];
			int full_offset = originalElemDOFOffset[full_i];
			int elem_dof = originalElemDOF[full_i];

			for(int dof = 0; dof < elem_dof; dof++)
			{
				if((size_t)(full_offset + dof) < m_flatExternFieldArray.size())
				{
					savedExternField[ima_offset + dof] = m_flatExternFieldArray[full_offset + dof];
				}
			}
			ima_offset += elem_dof;
		}
	}

	// Resize arrays
	m_flatExternFieldArray = std::move(savedExternField);
	m_flatMagnArray.resize(imaDOF, 0.0);
	m_flatFieldArray.resize(imaDOF, 0.0);

	// Update AmOfMainElem to IMA element count
	AmOfMainElem = m_imaNumElements;

	// Rebuild m_elemDOF and m_elemDOFOffset for IMA elements (variable DOF)
	m_elemDOF.resize(m_imaNumElements);
	m_elemDOFOffset.resize(m_imaNumElements + 1);
	int acc = 0;
	for(int ima_i = 0; ima_i < m_imaNumElements; ima_i++)
	{
		int full_i = m_imaToFull[ima_i];
		m_elemDOF[ima_i] = originalElemDOF[full_i];
		m_elemDOFOffset[ima_i] = acc;
		acc += m_elemDOF[ima_i];
	}
	m_elemDOFOffset[m_imaNumElements] = acc;

	m_totalDOF = imaDOF;

	// Remap g3dRelaxPtrVect to only contain IMA elements
	g3dRelaxPtrVect.resize(m_imaNumElements);
	for(int ima_i = 0; ima_i < m_imaNumElements; ima_i++)
	{
		int full_i = m_imaToFull[ima_i];
		g3dRelaxPtrVect[ima_i] = originalG3dRelaxPtrVect[full_i];
	}

	// IMA: AmOfMainElem updated, m_totalDOF set

	// For HACApK: skip dense matrix and recompute compact-element geometry for
	// the reduced IMA element set.
	if(skipDenseMatrix)
	{
		// Reset precomputed geometry so HACApK recomputes for the reduced IMA elements
		m_hexaGeomReady = false;
		m_wedgeGeomReady = false;
		m_tetraGeomReady = false;
		return 1;
	}

	// Independent face-coefficient IMA dense assembly is retired; HDiv-VIM owns
	// mesh-backed magnetic-material solves.  Compact 3-component IMA still
	// builds its 3x3 dense matrix below.
	if(hasFaceCharge)
		return 1;

	// Build mapping from IMA index to type-specific geometry index
	// Uses geometry precomputed from the full model (before system reduction)
	// NOTE: Always populate for ALL element types (not just pure meshes) to support mixed meshes
	// Tet mapping: stores full model index (Compute3x3BlockFast resolves internally)
	std::vector<int> imaToTet(m_imaNumElements, -1);
	if(m_tetraGeomReady)
	{
		for(int ima_i = 0; ima_i < m_imaNumElements; ima_i++)
		{
			if(m_elemDOF[ima_i] != 3) continue;
			imaToTet[ima_i] = m_imaToFull[ima_i];
		}
	}

	// Build IMA interaction matrix with TaskManager parallelization
	ngcore::ParallelFor(ngcore::IntRange(m_imaNumElements), [&](size_t ima_col)
	{
		int offset_col = m_elemDOFOffset[(int)ima_col];
		int dof_col = m_elemDOF[(int)ima_col];

		radTPolyhedron* poly_col = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[(int)ima_col]);
		if(!poly_col) return;

		for(int ima_row = 0; ima_row < m_imaNumElements; ima_row++)
		{
			int offset_row = m_elemDOFOffset[ima_row];
			int dof_row = m_elemDOF[ima_row];

			radTPolyhedron* poly_row = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[ima_row]);
			if(!poly_row) continue;

			double* block = &m_flatInteractMatrix[(size_t)offset_row * imaDOF + offset_col];

			// Fast path: compact 3DOF blocks; Compute3x3BlockFast handles IMA inline.
			if(dof_row == 3 && dof_col == 3)
			{
				int tet_row = imaToTet[ima_row];
				int tet_col = imaToTet[ima_col];
				if(tet_row >= 0 && tet_col >= 0)
				{
					double N_ima[9];
					Compute3x3BlockFast(tet_row, tet_col, N_ima);
					for(int i = 0; i < 3; i++)
						for(int j = 0; j < 3; j++)
							block[(size_t)i * imaDOF + j] = N_ima[i * 3 + j];
					continue;
				}
			}
		}
	});

	// IMA interaction matrix built successfully
	return 1;
}

//-------------------------------------------------------------------------
