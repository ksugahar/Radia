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
#include <cmath>
#include <algorithm>
#include <utility>

#include "rad_parallel.h"

// MSC (Magnetic Surface Charge) support for 6 DOF hexahedra
// radTPolyhedron hexahedra use 6 DOF MSC (surface charge on each face)
// MSC is always enabled (unconditional)

// Note: Dipole-dipole method for tetrahedra was tested but found numerically unstable.
// Radia production solver uses the surface charge (MSC) method.

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

	// Check if any element has variable DOF (e.g., 6 DOF MSC hexahedra)
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
		// The moment-yano H-matrix path (method 2) enumerates hexes via m_hexaElemIndices, populated by
		// PrecomputeHexaGeometry() -- normally done inside SetupInteractMatrix_VariableDOF (skipped here).
		// Call it directly (O(N) per-element geometry, idempotent) so the scalable moment solve sees the
		// hexes even with the dense matrix skipped (Phase 2 Increment 4 storage decoupling).  Harmless for
		// tet-only / EIEM2-HACApK models (nHex == 0 -> early return).
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
	// True if any relaxable element is a surface-charge (Use6DOF_MSC) polyhedron --
	// a soft-iron hexahedron (6 faces) or wedge/pyramid (5 faces).  These are solved by
	// the canonical moment-yano path in C++ unless the Python wrapper routes a mesh-backed
	// soft iron to FEEC HDiv-VIM.
	for(int i = 0; i < AmOfMainElem; i++)
	{
		radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[i]);
		if(poly && poly->Use6DOF_MSC) return true;
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
		// Check memory requirements before allocation (LU/BiCGSTAB need dense matrix)
		// Dense matrix requires N^2 * sizeof(TMatrix3df) = N^2 * 36 bytes
		size_t matrix_size = (size_t)AmOfMainElem * (size_t)AmOfMainElem;
		size_t required_bytes = matrix_size * sizeof(TMatrix3df);
		const size_t MAX_DENSE_MATRIX_BYTES = 100ULL * 1024 * 1024 * 1024;  // 100 GB limit

		if(required_bytes > MAX_DENSE_MATRIX_BYTES)
		{
			std::cerr << "[Radia] Error: Dense matrix too large for LU/BiCGSTAB solver." << std::endl;
			std::cerr << "[Radia] Elements=" << AmOfMainElem << ", required memory="
			          << (required_bytes / (1024*1024*1024)) << " GB" << std::endl;
			std::cerr << "[Radia] Use HACApK solver (method 2) for large problems." << std::endl;
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
			std::cerr << "[Radia] Use HACApK solver (method 2) for large problems." << std::endl;
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
// Variable DOF support for hybrid MSC + standard element analysis
// Reference: Yano & Sugahara, "MMM with MSC", J. Magn. Soc. Jpn., 2023
//=========================================================================
//-------------------------------------------------------------------------

void radTInteraction::ComputeDOFOffsets()
{
	// Compute DOF offsets for each element
	// This allows mixing elements with different DOF counts (3 for standard, 6 for MSC)

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
		std::cerr << "[Radia] Error: Dense matrix too large for LU/BiCGSTAB solver." << std::endl;
		std::cerr << "[Radia] DOF=" << m_totalDOF << ", required memory="
		          << (required_bytes / (1024*1024*1024)) << " GB" << std::endl;
		std::cerr << "[Radia] Use HACApK solver (method 2) for large problems (>100,000 DOF)." << std::endl;
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
		std::cerr << "[Radia] Use HACApK solver (method 2) for large problems." << std::endl;
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

	// Check if we have any MSC elements (5 DOF wedges or 6 DOF hexahedra)
	// MSC elements require Yano midpoint evaluation which is more complex
	bool hasMSCElements = false;
	for(int i = 0; i < AmOfMainElem && !hasMSCElements; i++)
	{
		if(m_elemDOF[i] >= 5) hasMSCElements = true;
	}

	// Build interaction matrix with variable-size blocks
	// For each pair (row_elem, col_elem), compute the interaction block

	// FAST PATH: Only for pure tetrahedra meshes (no MSC hexahedra) without symmetry
	if(!hasSymmetry && !hasMSCElements)
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

				// Compute the interaction block based on DOF types
				// FAST PATH: Only for 3x3 blocks (tetrahedra)
				// MSC hexahedra (6 DOF) fall through to slow path for correctness
				// (MSC requires Yano midpoint evaluation and proper transforms)
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
				// Note: MSC blocks (3x6, 6x3, 6x6) not handled in fast path
				// They will fall through and be computed in slow path below
			}
		});
		return 1;
	}

	// EIEM2 retirement (Phase 3b): the dense MSC (surface-charge) interaction blocks are no longer
	// assembled.  After the fail-loud guards in MakeAutoRelax (mixed MMM+MSC and B-input-on-MSC are
	// rejected), the ONLY path that reaches here with an MSC element (DOF>=5) is a pure surface-charge
	// (hex / wedge / pyramid) model solved by the parameter-free moment-yano formulation.  That path assembles its
	// own system (BuildMomentSystemCore, on-the-fly geometry via momentFaceGeom) and NEVER reads this
	// dense interaction matrix -- it is identical to the method-2 path, which runs with no dense matrix at
	// all (Phase 2 Increment 4).  So leave the MSC blocks zero (the matrix is already zero-initialized) and
	// return -- the EIEM2 collocation block kernels (Compute6x6/5x5/MixedBlockFast) are retired (Phase 3b).
	// Precompute the per-element geometry caches so the
	// method-2 H-matrix moment path (which enumerates hexes via m_hexaElemIndices) and any precompute-based
	// scaffold see the elements -- idempotent, O(N), mirrors the skipDenseMatrix=1 branch in Setup.
	if(hasMSCElements)
	{
		PrecomputeHexaGeometry();
		PrecomputeWedgeGeometry();
		return 1;
	}

	// EIEM2 retirement (Phase 3b): the dense MSC fast paths (pure-hex / pure-wedge via Compute6x6/5x5BlockFast)
	// and the no-symmetry MEDIUM MSC path were deleted.  The only surface-charge path now reaches the
	// hasMSCElements early-return above (moment-yano assembles its own system and never reads this dense
	// matrix).  What remains below is the SLOW PATH, now MMM-only (3x3, symmetry-aware) -- reached only by
	// an all-tetrahedron model WITH space-group symmetry (pure-tet without symmetry returned via the fast
	// path at the top of this function).

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
				// EIEM2 retirement (Phase 3b): the dof>=5 MSC branches were deleted.  In this MMM-only path
				// dof is always 3, so the 3x3 branch above always matches; the defensive else zeroes any
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
			// MSC (hex/wedge/pyramid): no per-face external field -- the moment-yano solve samples ExternFieldArray (centroid).
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
				// MSC: no per-face external field (moment-yano uses the centroid ExternFieldArray).
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
// Reference: ELF-style optimization (same as RadHACApKMMMManager::PrecomputeGeometry3DOF)
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
// Compute3x3BlockFast: Fast 3x3 interaction block for tetrahedra
// Uses pre-computed geometry (no B_comp overhead)
// Reference: ELF-style optimization (same as RadHACApKMMMManager::Compute3x3BlockFast)
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
	// Same pattern as Compute6x6BlockFast: mirror source geometry inline
	// For MMM (tet), M is a pseudovector: sign matrix S[beta] per component
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
// Reference: Yano MSC method for hexahedral elements
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

			// Get face vertices and split into 2 triangles
			const radTHandlePgnAndTrans& hpt = poly->VectHandlePgnAndTrans[f];
			radTPolygon* pgn = hpt.PgnHndl.rep;
			radTrans* tr = hpt.TransHndl.rep;

			const radTVect2dVect& verts2d = pgn->EdgePointsVector;
			if(verts2d.size() < 4) continue;

			TVector3d V[4];
			for(int v = 0; v < 4; v++)
			{
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
// BuildLoopBasis: cell-graph cycle (loop) basis of the yano-MSC operator.
//
// The loops are the field-null subspace of the retired EIEM2 surface-charge
// collocation operator N (== the HDiv ker(B) cell/internal-face cycle space).
// Runtime deflation/gauge machinery was removed; the basis remains as a
// historical geometry probe and for comparisons with the loop-free HDiv-VIM.
//
// Construction is geometry-only (no dense SVD): match internal faces by
// coincident face centers -> cell-adjacency graph -> spanning forest ->
// fundamental cycles.  Each internal face shared by hex a (DOF dA) and hex b
// (DOF dB) contributes +c to dA and -c to dB along the cycle (the field-null
// "+q / -q on the same physical face" convention).  Output Lflat is ROW-MAJOR
// (m_totalDOF rows x nLoop cols): Lflat[d * nLoop + col].  nLoop = the cell-graph
// cycle count = nInternalFaces - (nHex - nComponents).
//=========================================================================

void radTInteraction::BuildLoopBasis(std::vector<double>& Lflat, int& nLoop) const
{
	Lflat.clear(); nLoop = 0;
	int nHex = (int)m_hexaElemIndices.size();
	if(nHex == 0) return;

	// 1) face center + area + DOF index for every hex face (DOF f of hex h = offset + f)
	struct FaceRec { double cx, cy, cz, area; int hex, dof; };
	std::vector<FaceRec> faces; faces.reserve((size_t)nHex * 6);
	double scale = 0.0;
	for(int h = 0; h < nHex; h++)
	{
		int elemIdx = m_hexaElemIndices[h];
		radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[elemIdx]);
		if(!poly || poly->AmOfFaces != 6) continue;
		int off = m_elemDOFOffset[elemIdx];
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
			FaceRec fr;
			fr.cx = 0.25 * (V4[0].x + V4[1].x + V4[2].x + V4[3].x);
			fr.cy = 0.25 * (V4[0].y + V4[1].y + V4[2].y + V4[3].y);
			fr.cz = 0.25 * (V4[0].z + V4[1].z + V4[2].z + V4[3].z);
			fr.area = (area > 0 ? area : 1.0);
			fr.hex = h; fr.dof = off + f;
			faces.push_back(fr);
			scale = std::max(scale, std::fabs(fr.cx) + std::fabs(fr.cy) + std::fabs(fr.cz));
		}
	}
	double tol = 1e-9 * (scale > 0 ? scale : 1.0);

	// 2) match coincident face centers -> internal edges (conforming mesh: shared centers coincide).
	// Bucket per spatial-hash key as a VECTOR (not a single int) and match against ALL members by
	// distance, always registering the face.  The previous single-bucket map dropped any face whose
	// XOR hash COLLIDED with an unrelated earlier face (it found the wrong face, failed the distance
	// check, and never registered itself, so its true partner never matched) -- that silently
	// undercounted the internal faces and hence the cycle (loop) basis by exactly the collision count.
	struct Edge { int hexA, dofA, hexB, dofB; double areaA, areaB; };
	std::vector<Edge> edges;
	std::unordered_map<long long, std::vector<int> > buckets;   // spatial-hash key -> all face indices
	for(int i = 0; i < (int)faces.size(); i++)
	{
		long long kx = (long long)std::llround(faces[i].cx / tol);
		long long ky = (long long)std::llround(faces[i].cy / tol);
		long long kz = (long long)std::llround(faces[i].cz / tol);
		long long key = (kx * 73856093LL) ^ (ky * 19349663LL) ^ (kz * 83492791LL);
		std::vector<int>& bucket = buckets[key];
		for(size_t b = 0; b < bucket.size(); b++)
		{
			int j = bucket[b];
			double d = std::fabs(faces[i].cx - faces[j].cx) + std::fabs(faces[i].cy - faces[j].cy)
			         + std::fabs(faces[i].cz - faces[j].cz);
			if(d <= 10.0 * tol && faces[i].hex != faces[j].hex)
			{
				Edge e; e.hexA = faces[j].hex; e.dofA = faces[j].dof; e.areaA = faces[j].area;
				e.hexB = faces[i].hex; e.dofB = faces[i].dof; e.areaB = faces[i].area;
				edges.push_back(e);
			}
		}
		bucket.push_back(i);
	}
	int nInternal = (int)edges.size();
	if(nInternal == 0) return;

	// 3) spanning forest (BFS) with parent edge; non-tree edges close fundamental cycles
	std::vector<std::vector<std::pair<int,int> > > adj(nHex);   // hex -> (neighbor, edgeIdx)
	for(int e = 0; e < nInternal; e++)
	{
		adj[edges[e].hexA].push_back(std::make_pair(edges[e].hexB, e));
		adj[edges[e].hexB].push_back(std::make_pair(edges[e].hexA, e));
	}
	std::vector<int> parentNode(nHex, -1), parentEdge(nHex, -1);
	std::vector<char> visited(nHex, 0), treeEdge(nInternal, 0);
	for(int s = 0; s < nHex; s++)
	{
		if(visited[s]) continue;
		visited[s] = 1; std::vector<int> q; q.push_back(s); size_t qi = 0;
		while(qi < q.size())
		{
			int u = q[qi++];
			for(size_t a = 0; a < adj[u].size(); a++)
			{
				int v = adj[u][a].first, eidx = adj[u][a].second;
				if(!visited[v]) { visited[v] = 1; parentNode[v] = u; parentEdge[v] = eidx; treeEdge[eidx] = 1; q.push_back(v); }
			}
		}
	}
	std::vector<int> nonTree;
	for(int e = 0; e < nInternal; e++) if(!treeEdge[e]) nonTree.push_back(e);
	nLoop = (int)nonTree.size();
	if(nLoop == 0) return;

	// 4) assemble L: fundamental cycle = non-tree edge + tree paths of its two endpoints
	Lflat.assign((size_t)m_totalDOF * nLoop, 0.0);
	for(int c = 0; c < nLoop; c++)
	{
		int e = nonTree[c];
		std::map<int,int> coeff;     // edgeIdx -> integer cycle coefficient
		int x = edges[e].hexA; while(parentNode[x] != -1) { coeff[parentEdge[x]] += 1; x = parentNode[x]; }
		x = edges[e].hexB;     while(parentNode[x] != -1) { coeff[parentEdge[x]] -= 1; x = parentNode[x]; }
		coeff[e] += 1;
		for(std::map<int,int>::iterator it = coeff.begin(); it != coeff.end(); ++it)
		{
			if(it->second == 0) continue;
			const Edge& ee = edges[it->first];
			// integer coefficient = FLUX along the cycle; charge DENSITY sigma = flux / face area
			// (the field-null condition needs Sum_f area*sigma = 0 per cell -> use flux, not density;
			//  on a cube all areas are equal so 1/area is a constant, but on a distorted hex it matters).
			double q = (double)it->second;
			Lflat[(size_t)ee.dofA * nLoop + c] += q / ee.areaA;
			Lflat[(size_t)ee.dofB * nLoop + c] -= q / ee.areaB;
		}
	}
}

//=========================================================================
// BuildFaceGeom: per-DOF hex face geometry in the matrix DOF order.
// Row-major (m_totalDOF x 11): [elem_local, area, cx,cy,cz, nx,ny,nz(outward), ecx,ecy,ecz].
// Mirrors the DOF<->face mapping of BuildLoopBasis (DOF = m_elemDOFOffset[elem] + f), so the rows
// align 1:1 with GetInteractMatrix.  Lets Python form the div(B)=0 constraint (Sum_f area_f*sigma_f
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
			// quad area via two triangles (V0V1V2 + V0V2V3) -- same as BuildLoopBasis
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
// CentroidFieldGradFromFace: the reusable single-(target,source) kernel of the moment formulation -- field
// H[3] + grad gH[6] (xx,yy,zz,xy,xz,yz, RAW = before the 1/4pi factor) at target centroid ce3 from ONE unit
// surface-charge face (V4 corners, srcCenter = the face's element center, area), INCLUDING the IMA mirror
// images.  isSelf drops the original yano center charge (singularity-free); each IMA mirror always carries
// its (reflected) center charge.  Shared by BuildCentroidFieldGrad (dense Cflat) and MomentSystemEntry
// (on-demand H-matrix entry) so the two never diverge.  Caller applies 1/4pi.
//=========================================================================
void radTInteraction::CentroidFieldGradFromFace(const double ce3[3], const double V4[4][3],
        const double srcCenter[3], bool isSelf, double area, double Hout[3], double gHout[6]) const
{
	static const int NG = 8;
	static const double gl_x[NG] = {-0.9602898564975363,-0.7966664774136267,-0.5255324099163290,-0.1834346424956498,
	                                 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
	static const double gl_w[NG] = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
	                                 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
	double gp[NG], gw[NG];
	for(int i = 0; i < NG; i++) { gp[i] = 0.5*(gl_x[i]+1.0); gw[i] = 0.5*gl_w[i]; }

	Hout[0] = Hout[1] = Hout[2] = 0.0;
	for(int k = 0; k < 6; k++) gHout[k] = 0.0;

	// field+grad from a bilinear quad Vq (unit charge density, 8x8 Gauss) + optional point charge q=-area at
	// cen; scaled by sgn (IMA image sign; 1.0 original).  Computing FRESH from mirrored geometry makes the
	// rank-2 gradient transform under the reflection automatically -- sgn only carries the BC charge sign.
	auto accumQG = [&](const double Vq[4][3], bool withCenter, const double cen[3], double sgn)
	{
		for(int iu = 0; iu < NG; iu++) for(int iv = 0; iv < NG; iv++)
		{
			double u = gp[iu], v = gp[iv], wuv = gw[iu]*gw[iv];
			double a0 = (1-u)*(1-v), a1 = u*(1-v), a2 = u*v, a3 = (1-u)*v;
			double Px = a0*Vq[0][0]+a1*Vq[1][0]+a2*Vq[2][0]+a3*Vq[3][0];
			double Py = a0*Vq[0][1]+a1*Vq[1][1]+a2*Vq[2][1]+a3*Vq[3][1];
			double Pz = a0*Vq[0][2]+a1*Vq[1][2]+a2*Vq[2][2]+a3*Vq[3][2];
			double Tux = (1-v)*(Vq[1][0]-Vq[0][0])+v*(Vq[2][0]-Vq[3][0]);
			double Tuy = (1-v)*(Vq[1][1]-Vq[0][1])+v*(Vq[2][1]-Vq[3][1]);
			double Tuz = (1-v)*(Vq[1][2]-Vq[0][2])+v*(Vq[2][2]-Vq[3][2]);
			double Tvx = (1-u)*(Vq[3][0]-Vq[0][0])+u*(Vq[2][0]-Vq[1][0]);
			double Tvy = (1-u)*(Vq[3][1]-Vq[0][1])+u*(Vq[2][1]-Vq[1][1]);
			double Tvz = (1-u)*(Vq[3][2]-Vq[0][2])+u*(Vq[2][2]-Vq[1][2]);
			double jx = Tuy*Tvz-Tuz*Tvy, jy = Tuz*Tvx-Tux*Tvz, jz = Tux*Tvy-Tuy*Tvx;
			double dA = std::sqrt(jx*jx+jy*jy+jz*jz)*wuv*sgn;
			double dx = ce3[0]-Px, dy = ce3[1]-Py, dz = ce3[2]-Pz;
			double r2 = dx*dx+dy*dy+dz*dz, inv_r = 1.0/std::sqrt(r2);
			double inv_r3 = inv_r/r2, inv_r5 = inv_r3/r2;
			double c3 = inv_r3*dA, c5 = inv_r5*dA;
			Hout[0] += dx*c3; Hout[1] += dy*c3; Hout[2] += dz*c3;
			gHout[0] += c3 - 3.0*dx*dx*c5; gHout[1] += c3 - 3.0*dy*dy*c5; gHout[2] += c3 - 3.0*dz*dz*c5;
			gHout[3] += -3.0*dx*dy*c5;     gHout[4] += -3.0*dx*dz*c5;     gHout[5] += -3.0*dy*dz*c5;
		}
		if(withCenter)
		{
			double dx = ce3[0]-cen[0], dy = ce3[1]-cen[1], dz = ce3[2]-cen[2];
			double r2 = dx*dx+dy*dy+dz*dz, inv_r = 1.0/std::sqrt(r2);
			double inv_r3 = inv_r/r2, inv_r5 = inv_r3/r2, q = -area*sgn;
			Hout[0] += q*dx*inv_r3; Hout[1] += q*dy*inv_r3; Hout[2] += q*dz*inv_r3;
			gHout[0] += q*(inv_r3-3.0*dx*dx*inv_r5); gHout[1] += q*(inv_r3-3.0*dy*dy*inv_r5);
			gHout[2] += q*(inv_r3-3.0*dz*dz*inv_r5); gHout[3] += q*(-3.0*dx*dy*inv_r5);
			gHout[4] += q*(-3.0*dx*dz*inv_r5);       gHout[5] += q*(-3.0*dy*dz*inv_r5);
		}
	};

	accumQG(V4, !isSelf, srcCenter, 1.0);                  // original (center charge only for the mutual pairing)
	if(m_imaEnabled)                                       // IMA mirror images (scalar BC sign)
	{
		auto addMir = [&](int ax, double sgn)
		{
			double Vm[4][3]; double cm[3] = {srcCenter[0], srcCenter[1], srcCenter[2]};
			for(int c = 0; c < 4; c++) { Vm[c][0]=V4[c][0]; Vm[c][1]=V4[c][1]; Vm[c][2]=V4[c][2]; }
			if(ax & IMA_X) { for(int c=0;c<4;c++) Vm[c][0]=-Vm[c][0]; cm[0]=-cm[0]; }
			if(ax & IMA_Y) { for(int c=0;c<4;c++) Vm[c][1]=-Vm[c][1]; cm[1]=-cm[1]; }
			if(ax & IMA_Z) { for(int c=0;c<4;c++) Vm[c][2]=-Vm[c][2]; cm[2]=-cm[2]; }
			accumQG(Vm, true, cm, sgn);
		};
		bool hX = (m_imaSymmetry & IMA_X) != 0, hY = (m_imaSymmetry & IMA_Y) != 0, hZ = (m_imaSymmetry & IMA_Z) != 0;
		if(hX) addMir(IMA_X, (double)m_imaSignX);
		if(hY) addMir(IMA_Y, (double)m_imaSignY);
		if(hZ) addMir(IMA_Z, (double)m_imaSignZ);
		if(hX && hY) addMir(IMA_XY, (double)m_imaSignX*m_imaSignY);
		if(hX && hZ) addMir(IMA_XZ, (double)m_imaSignX*m_imaSignZ);
		if(hY && hZ) addMir(IMA_YZ, (double)m_imaSignY*m_imaSignZ);
		if(hX && hY && hZ) addMir(IMA_XYZ, (double)m_imaSignX*m_imaSignY*m_imaSignZ);
	}
}

//=========================================================================
// momentFaceGeom: per-face geometry for the moment formulation, generalized to triangular (3-vert) and
// quad (4-vert) faces.  Fills V4col (the corners, with a TRIANGLE returned as a degenerate quad V[3]=V[2]
// so the bilinear-quad kernel CentroidFieldGradFromFace integrates the triangle), the face center fc, the
// OUTWARD unit normal nf, the polygon area (fan over nv-2 triangles), and d = fc - element centroid.  For a
// hex quad (nv=4) this reproduces the previous inline geometry BIT-FOR-BIT.  File-static: used only by the
// three moment kernel functions below (BuildCentroidFieldGrad / BuildMomentSystemCore / MomentSystemEntry).
//=========================================================================
static int momentFaceGeom(radTPolyhedron* poly, int f, const TVector3d& ce,
                          double V4col[4][3], double fc[3], double nf[3], double& area, double d[3])
{
	const radTHandlePgnAndTrans& hpt = poly->VectHandlePgnAndTrans[f];
	radTPolygon* pgn = hpt.PgnHndl.rep; radTrans* tr = hpt.TransHndl.rep;
	const radTVect2dVect& v2d = pgn->EdgePointsVector;
	int nv = (int)v2d.size(); if(nv > 4) nv = 4; if(nv < 3) nv = 3;
	TVector3d V[4];
	for(int v = 0; v < nv; v++) V[v] = tr->TrPoint(TVector3d(v2d[v].x, v2d[v].y, pgn->CoordZ));
	if(nv == 3) V[3] = V[2];                       // degenerate quad -> bilinear kernel integrates the triangle
	double cx = 0, cy = 0, cz = 0;
	for(int v = 0; v < nv; v++) { cx += V[v].x; cy += V[v].y; cz += V[v].z; }
	fc[0] = cx/nv; fc[1] = cy/nv; fc[2] = cz/nv;
	double sx = 0, sy = 0, sz = 0; area = 0.0;
	for(int t = 0; t < nv-2; t++)                  // fan triangulation (quad: 2 triangles, tri: 1)
	{
		const TVector3d& P0 = V[0]; const TVector3d& P1 = V[t+1]; const TVector3d& P2 = V[t+2];
		double ux = P1.x-P0.x, uy = P1.y-P0.y, uz = P1.z-P0.z;
		double wx = P2.x-P0.x, wy = P2.y-P0.y, wz = P2.z-P0.z;
		double rx = uy*wz-uz*wy, ry = uz*wx-ux*wz, rz = ux*wy-uy*wx;
		area += 0.5*std::sqrt(rx*rx+ry*ry+rz*rz); sx += rx; sy += ry; sz += rz;
	}
	double nlen = std::sqrt(sx*sx+sy*sy+sz*sz); if(nlen < 1e-300) nlen = 1.0;
	double ox = fc[0]-ce.x, oy = fc[1]-ce.y, oz = fc[2]-ce.z;
	double sgn = (sx*ox+sy*oy+sz*oz >= 0.0) ? 1.0 : -1.0;
	nf[0] = sgn*sx/nlen; nf[1] = sgn*sy/nlen; nf[2] = sgn*sz/nlen;
	d[0] = fc[0]-ce.x; d[1] = fc[1]-ce.y; d[2] = fc[2]-ce.z;
	for(int v = 0; v < 4; v++) { V4col[v][0] = V[v].x; V4col[v][1] = V[v].y; V4col[v][2] = V[v].z; }
	return nv;
}

//=========================================================================
// momentResidualEigenmodes: orthonormal basis of the RESIDUAL subspace -- the (nF-4) charge patterns with
// ZERO monopole and ZERO dipole MOMENT (null space of the 4 functionals {Ae, Ae*d_x, Ae*d_y, Ae*d_z} in
// R^nF).  These are the element's natural QUADRUPOLE eigenmodes; using them as the moment-system quadrupole
// test directions (in place of the hand-picked dx^2-dy^2 / axial forms) is geometry-adaptive and ALWAYS
// full-rank -- it never hits the degenerate near-null mode a fixed quadratic can on a distorted element
// (the wedge dx^2-dy^2 -> M~1e9 blow-up that motivated the axial hand-fix; verified equivalent to the
// hand-pick on symmetric elements -> cube demag N=1/3 preserved).  Deterministic Gram-Schmidt (nF <= 6,
// no LAPACK needed); identical Ae/d in BuildMomentSystemCore and MomentSystemEntry -> identical modes, so
// the dense (method 0/1) and H-matrix (method 2) moment systems stay consistent.  Returns the mode count
// (nF-4 generically); phi[q][f] is the q-th orthonormal mode (sum_f phi[q][f]^2 = 1).
//=========================================================================
static int momentResidualEigenmodes(const double Ae[], const double d[][3], int nF, double phi[][6])
{
	double basis[6][6]; int nb = 0;
	auto orthoAdd = [&](double v[6]) -> bool
	{
		for(int b = 0; b < nb; b++)
		{
			double dot = 0.0; for(int f = 0; f < nF; f++) dot += basis[b][f]*v[f];
			for(int f = 0; f < nF; f++) v[f] -= dot*basis[b][f];
		}
		double nrm = 0.0; for(int f = 0; f < nF; f++) nrm += v[f]*v[f]; nrm = std::sqrt(nrm);
		if(nrm > 1.0e-9) { for(int f = 0; f < nF; f++) basis[nb][f] = v[f]/nrm; nb++; return true; }
		return false;
	};
	// 1) orthonormalize the monopole + 3 dipole MOMENT functionals -> ortho basis of the NON-residual part
	for(int i = 0; i < 4; i++)
	{
		double v[6]; for(int f = 0; f < nF; f++) v[f] = (i == 0) ? Ae[f] : Ae[f]*d[f][i-1];
		orthoAdd(v);
	}
	// 2) extend to R^nF with the standard basis; each newly-independent residual vector is a quad eigenmode
	int nq = 0;
	for(int e = 0; e < nF && nb < nF; e++)
	{
		double v[6] = {0,0,0,0,0,0}; v[e] = 1.0;
		if(orthoAdd(v)) { for(int f = 0; f < nF; f++) phi[nq][f] = basis[nb-1][f]; nq++; }
	}
	return nq;
}

//=========================================================================
// CollectMomentElems: e-ordered list of MOMENT elements (the surface-charge polyhedra: hex with 6 DOF +
// wedge with 5 DOF).  For a pure-hex model this equals m_hexaElemIndices order, so the moment Cflat / row
// layout is unchanged (hex stays bit-identical).  THE single source of element ordering for the moment
// dense path (BuildCentroidFieldGrad + BuildMomentSystemCore) and the SolveLinearStep moment branch.
//=========================================================================
void radTInteraction::CollectMomentElems(std::vector<int>& out) const
{
	out.clear();
	for(int e = 0; e < AmOfMainElem; e++)
		if(m_elemDOF[e] == 6 || m_elemDOF[e] == 5) out.push_back(e);
}

//=========================================================================
// BuildCentroidFieldGrad: per moment element demag field H and gradient gradH at the element CENTROID, as linear
// functionals of each source DOF charge -- the kernel of the parameter-free moment formulation (the fix
// for the eval-point alpha + the finite-difference conditioning noise; validated in
// examples/vim/yano_moment_analytic_selfterm.py).
//   SELF face (same element): bare charged-face field, no center charge.  The centroid is interior, at
//     finite distance from every face, so the integrand is smooth -> exact by Gauss quadrature, and it is
//     patch-test exact (single cube -> demag N=1/3).
//   MUTUAL face: yano dipole layer = (bare face) - area*(point charge @ source element center).  The
//     per-face center charge is REQUIRED (without it a regular grid develops a charge-free dipole-free
//     quadrupole near-null mode); the source center is at finite distance from this centroid -> finite.
//     Excluding it for SELF is exactly where it would be singular -> the kernel is singularity-free.
// Output Cflat ROW-MAJOR (nMom x 9 x m_totalDOF): comp k (Hx,Hy,Hz, gxx,gyy,gzz,gxy,gxz,gyz), source DOF g
// -> Cflat[(h*9+k)*m_totalDOF + g].  Field convention H = (1/4pi) int sigma (r-r')/|r-r'|^3 dA'.
//=========================================================================

void radTInteraction::BuildCentroidFieldGrad(std::vector<double>& Cflat, int& nHexOut) const
{
	const int NK = 9;
	std::vector<int> melem; CollectMomentElems(melem);   // hex (6) + wedge/pyramid (5), e-ordered
	int nMom = (int)melem.size();
	nHexOut = nMom;                                       // historical name; = moment-element count (hex+wedge+pyramid)
	Cflat.assign((size_t)nMom * NK * m_totalDOF, 0.0);
	if(nMom == 0 || m_totalDOF == 0) return;

	// gather per-DOF source-face data (corners as a possibly-degenerate quad, polygon area, source element
	// center + index) for every moment element's faces (hex 6 / wedge 5; triangle faces via momentFaceGeom).
	struct FaceRec { double V[4][3]; double area; TVector3d srcEC; int srcElem; int valid; };
	std::vector<FaceRec> faces((size_t)m_totalDOF);
	for(size_t i = 0; i < faces.size(); i++) { faces[i].valid = 0; faces[i].srcElem = -1; }

	for(int h = 0; h < nMom; h++)
	{
		int elemIdx = melem[h];
		radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[elemIdx]);
		if(!poly) continue;
		int nF = poly->AmOfFaces; if(nF != 6 && nF != 5) continue;
		int off = m_elemDOFOffset[elemIdx];
		const TVector3d ec = poly->CentrPoint;
		for(int f = 0; f < nF; f++)
		{
			const radTHandlePgnAndTrans& hpt = poly->VectHandlePgnAndTrans[f];
			if(!hpt.PgnHndl.rep || !hpt.TransHndl.rep) continue;
			double V4col[4][3], fc[3], nf[3], d[3], area;
			momentFaceGeom(poly, f, ec, V4col, fc, nf, area, d);   // tri -> degenerate quad; hex bit-identical
			FaceRec& fr = faces[off+f];
			for(int v = 0; v < 4; v++) { fr.V[v][0]=V4col[v][0]; fr.V[v][1]=V4col[v][1]; fr.V[v][2]=V4col[v][2]; }
			fr.area = area; fr.srcEC = ec; fr.srcElem = elemIdx; fr.valid = 1;
		}
	}

	const double INV4PI = 1.0/(4.0*3.14159265358979323846);

	for(int h = 0; h < nMom; h++)
	{
		int elemIdx = melem[h];
		radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[elemIdx]);
		if(!poly) continue;
		const TVector3d ce = poly->CentrPoint;
		const double ce3[3] = {ce.x, ce.y, ce.z};
		size_t base = (size_t)h * NK * m_totalDOF;
		for(int g = 0; g < m_totalDOF; g++)
		{
			const FaceRec& fr = faces[g];
			if(!fr.valid) continue;
			const double V4[4][3] = {{fr.V[0][0],fr.V[0][1],fr.V[0][2]}, {fr.V[1][0],fr.V[1][1],fr.V[1][2]},
			                         {fr.V[2][0],fr.V[2][1],fr.V[2][2]}, {fr.V[3][0],fr.V[3][1],fr.V[3][2]}};
			const double cen[3] = {fr.srcEC.x, fr.srcEC.y, fr.srcEC.z};
			double H[3], gH[6];
			CentroidFieldGradFromFace(ce3, V4, cen, fr.srcElem == elemIdx, fr.area, H, gH);  // incl. IMA mirrors
			Cflat[base + 0*m_totalDOF + g] = H[0]*INV4PI;
			Cflat[base + 1*m_totalDOF + g] = H[1]*INV4PI;
			Cflat[base + 2*m_totalDOF + g] = H[2]*INV4PI;
			for(int k = 0; k < 6; k++) Cflat[base + (size_t)(3+k)*m_totalDOF + g] = gH[k]*INV4PI;
		}
	}
}

//=========================================================================
// BuildMomentSystemCore: the parameter-free MOMENT-yano system matrix A and RHS for PER-ELEMENT linear
// susceptibility chiPerHex[h] in a PER-ELEMENT external field HextPerHex[h*3+k] (at the element centroid) -- the
// C++ port of examples/vim/yano_moment_iter_scaling.py::build, generalized for the solve (coil/source fields
// are not uniform).  Per moment element (hex 6 DOF, wedge/pyramid 5 DOF):
//   3 dipole rows : (local dipole moment of sigma)/Ve - chi*H_k(centroid) . sigma = chi*Hext_k
//   1 monopole row: sum_f area_f sigma_f = 0                       (= div B = 0)
//   2 quad rows   : (local diagonal-quadrupole moment of sigma) - chi*(Dvec . gradH(centroid)) . sigma = 0
// Global field/grad functionals H,gradH(centroid) come from BuildCentroidFieldGrad; the local geometric
// moments (face center fc, outward normal n, area, d=fc-centroid, volume Ve) from the element geometry.  Each
// row is 2-norm normalized.  A is ROW-MAJOR (dof x dof); rhs length dof.  The column index of A is the face
// DOF, so dgesv's solution is sigma in DOF order (drop-in for the retired EIEM2 LU write-back).
// The uniform BuildMomentSystem(chi,Happ,...) wrapper below broadcasts a scalar chi + uniform Happ.
//=========================================================================
void radTInteraction::BuildMomentSystemCore(const double* chiPerHex, const double* HextPerHex,
                                            std::vector<double>& A, std::vector<double>& rhs, bool normalize) const
{
	std::vector<double> Cflat; int nMom = 0;
	BuildCentroidFieldGrad(Cflat, nMom);                 // (nMom x 9 x dof): H[3], gradH[6] per moment elem
	std::vector<int> melem; CollectMomentElems(melem);   // SAME hex(6)+wedge/pyramid(5) e-order as Cflat
	const int dof = m_totalDOF;
	A.assign((size_t)dof * dof, 0.0);
	rhs.assign(dof, 0.0);
	if(nMom == 0 || dof == 0) return;

	std::vector<double> r((size_t)dof);
	int row = 0;
	for(int h = 0; h < nMom; h++)
	{
		int elemIdx = melem[h];
		radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[elemIdx]);
		if(!poly) continue;
		int nF = poly->AmOfFaces; if(nF != 6 && nF != 5) continue;     // hex 6 / wedge 5
		const int off = m_elemDOFOffset[elemIdx];
		const TVector3d ce = poly->CentrPoint;
		const double chiH = chiPerHex[h];                 // per-element susceptibility (linear or current Picard)
		const double* Hext = &HextPerHex[(size_t)h*3];    // external field at this element centroid

		// per-face local geometry (center fc, OUTWARD unit normal nf, area Ae, d = fc - centroid); the shared
		// momentFaceGeom handles tri(3-vert) + quad(4-vert) faces and is bit-identical for the hex quad case.
		double fc[6][3], nf[6][3], Ae[6], d[6][3]; int fnv[6];
		for(int f = 0; f < nF; f++)
		{
			double V4col[4][3];
			fnv[f] = momentFaceGeom(poly, f, ce, V4col, fc[f], nf[f], Ae[f], d[f]);
		}
		// residual quadrupole eigenmodes (the geometry-adaptive replacement for the hand-picked dx^2-dy^2 /
		// axial forms): the (nF-4) zero-monopole, zero-dipole-moment charge patterns the element can carry.
		double phiQ[6][6]; int nQ = momentResidualEigenmodes(Ae, d, nF, phiQ); (void)nQ;
		double Ve = 0.0;
		for(int f = 0; f < nF; f++) Ve += Ae[f]*(fc[f][0]*nf[f][0]+fc[f][1]*nf[f][1]+fc[f][2]*nf[f][2]);
		Ve *= (1.0/3.0);

		const double* Hh = &Cflat[(size_t)h*9*dof];           // Hh[k*dof..] field (k<3), grad (k 3..8)

		auto putRow = [&](double rh)
		{
			if(normalize)
			{
			double nn = 0.0; for(int g = 0; g < dof; g++) nn += r[g]*r[g];
			nn = std::sqrt(nn);
			if(nn > 1e-300) { double inv = 1.0/nn; for(int g = 0; g < dof; g++) r[g] *= inv; rh *= inv; }
			}
			double* Arow = &A[(size_t)row*dof];
			for(int g = 0; g < dof; g++) Arow[g] = r[g];
			rhs[row] = rh; row++;
		};

		// 3 dipole rows
		for(int k = 0; k < 3; k++)
		{
			std::fill(r.begin(), r.end(), 0.0);
			for(int f = 0; f < nF; f++) r[off+f] += Ae[f]*d[f][k]/Ve;
			const double* F0k = &Hh[(size_t)k*dof];
			for(int g = 0; g < dof; g++) r[g] -= chiH*F0k[g];
			putRow(chiH*Hext[k]);
		}
		// monopole row
		std::fill(r.begin(), r.end(), 0.0);
		for(int f = 0; f < nF; f++) r[off+f] = Ae[f];
		putRow(0.0);
		// residual-eigenmode quadrupole rows: (nF-4) of them (hex 2, wedge 1).  3 dipole + 1 monopole +
		// (nF-4) quad = nF rows = nF DOF (square per element).  Bm = the qq-th residual eigenmode written as
		// a per-face form so the test term Ae*Bm = phi_qq exactly; the cm-correction + Dm field-balance below
		// are unchanged (they act on Ae*Bm).
		for(int qq = 0; qq < nF-4; qq++)
		{
			double Bm[6];
			for(int f = 0; f < nF; f++) Bm[f] = (Ae[f] > 1.0e-300) ? phiQ[qq][f]/Ae[f] : 0.0;
			std::fill(r.begin(), r.end(), 0.0);
			for(int f = 0; f < nF; f++) r[off+f] += Ae[f]*Bm[f];
			double cm[3] = {0,0,0};
			for(int f = 0; f < nF; f++) for(int k = 0; k < 3; k++) cm[k] += Ae[f]*nf[f][k]*Bm[f];
			for(int f = 0; f < nF; f++)
			{
				double cd = cm[0]*d[f][0]+cm[1]*d[f][1]+cm[2]*d[f][2];
				r[off+f] -= Ae[f]*cd/Ve;
			}
			double Dm[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
			for(int f = 0; f < nF; f++) for(int ii = 0; ii < 3; ii++) for(int jj = 0; jj < 3; jj++)
				Dm[ii][jj] += Ae[f]*d[f][jj]*nf[f][ii]*Bm[f];
			double Dvec[6] = {Dm[0][0], Dm[1][1], Dm[2][2], Dm[0][1]+Dm[1][0], Dm[0][2]+Dm[2][0], Dm[1][2]+Dm[2][1]};
			for(int m = 0; m < 6; m++)
			{
				const double* Gm = &Hh[(size_t)(3+m)*dof];
				double w = chiH*Dvec[m];
				for(int g = 0; g < dof; g++) r[g] -= w*Gm[g];
			}
			putRow(0.0);
		}
	}
}

// Uniform-field wrapper: broadcast a scalar chi + uniform applied field Happ to every moment element
// (verification path + uniform-source linear solves).  See BuildMomentSystemCore.
void radTInteraction::BuildMomentSystem(double chi, const double Happ[3],
                                        std::vector<double>& A, std::vector<double>& rhs) const
{
	std::vector<int> melem; CollectMomentElems(melem);   // hex(6)+wedge/pyramid(5), matches BuildMomentSystemCore order
	int nMom = (int)melem.size();
	if(nMom <= 0) { A.clear(); rhs.assign(m_totalDOF, 0.0); return; }
	std::vector<double> chiv((size_t)nMom, chi), Hv((size_t)nMom*3);
	for(int h = 0; h < nMom; h++) { Hv[(size_t)h*3] = Happ[0]; Hv[(size_t)h*3+1] = Happ[1]; Hv[(size_t)h*3+2] = Happ[2]; }
	BuildMomentSystemCore(chiv.data(), Hv.data(), A, rhs);
}

//=========================================================================
// MomentSystemEntry: the ON-DEMAND un-normalized moment system entry A_raw[rowGlobal][colDOF] (the HACApK
// H-matrix entry; see docs/moment_yano/ACA_MOMENT_DESIGN.md).  Reproduces BuildMomentSystemCore's row math
// for a SINGLE (row,col) WITHOUT building the full system or normalizing -- the row 2-norm is a diagonal
// scaling that leaves the direct solve invariant, so the H-LU path uses A_raw.  HEX-ONLY; assumes all
// m_hexaElemIndices are valid 6-face hexes (rowGlobal = 6*h + t).
//=========================================================================
double radTInteraction::MomentSystemEntry(int rowGlobal, int colDOF, const double* chiPerHex) const
{
	const int dof = m_totalDOF;
	if(colDOF < 0 || colDOF >= dof) return 0.0;
	int nHex = (int)m_hexaElemIndices.size();
	int h = rowGlobal / 6, t = rowGlobal % 6;
	if(h < 0 || h >= nHex) return 0.0;
	int elemIdx = m_hexaElemIndices[h];
	radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[elemIdx]);
	if(!poly || poly->AmOfFaces != 6) return 0.0;
	const int off = m_elemDOFOffset[elemIdx];
	const TVector3d ce = poly->CentrPoint;
	const double chiH = chiPerHex[h];

	// per-face local geometry of the ROW element (identical to BuildMomentSystemCore)
	double fc[6][3], nf[6][3], Ae[6], d[6][3];
	for(int f = 0; f < 6; f++)
	{
		const radTHandlePgnAndTrans& hpt = poly->VectHandlePgnAndTrans[f];
		radTPolygon* pgn = hpt.PgnHndl.rep; radTrans* tr = hpt.TransHndl.rep;
		const radTVect2dVect& v2d = pgn->EdgePointsVector;
		TVector3d V4[4];
		for(int v = 0; v < 4; v++) V4[v] = tr->TrPoint(TVector3d(v2d[v].x, v2d[v].y, pgn->CoordZ));
		double cx = 0, cy = 0, cz = 0;
		for(int v = 0; v < 4; v++) { cx += V4[v].x; cy += V4[v].y; cz += V4[v].z; }
		fc[f][0] = cx*0.25; fc[f][1] = cy*0.25; fc[f][2] = cz*0.25;
		double sx = 0, sy = 0, sz = 0, area = 0;
		for(int tt = 0; tt < 2; tt++)
		{
			const TVector3d& P0 = V4[0]; const TVector3d& P1 = V4[tt+1]; const TVector3d& P2 = V4[tt+2];
			double ux = P1.x-P0.x, uy = P1.y-P0.y, uz = P1.z-P0.z;
			double wx = P2.x-P0.x, wy = P2.y-P0.y, wz = P2.z-P0.z;
			double rx = uy*wz-uz*wy, ry = uz*wx-ux*wz, rz = ux*wy-uy*wx;
			area += 0.5*std::sqrt(rx*rx+ry*ry+rz*rz); sx += rx; sy += ry; sz += rz;
		}
		double nlen = std::sqrt(sx*sx+sy*sy+sz*sz); if(nlen < 1e-300) nlen = 1.0;
		double ox = fc[f][0]-ce.x, oy = fc[f][1]-ce.y, oz = fc[f][2]-ce.z;
		double sgn = (sx*ox+sy*oy+sz*oz >= 0.0) ? 1.0 : -1.0;
		nf[f][0] = sgn*sx/nlen; nf[f][1] = sgn*sy/nlen; nf[f][2] = sgn*sz/nlen;
		Ae[f] = area;
		d[f][0] = fc[f][0]-ce.x; d[f][1] = fc[f][1]-ce.y; d[f][2] = fc[f][2]-ce.z;
	}
	double Ve = 0.0;
	for(int f = 0; f < 6; f++) Ve += Ae[f]*(fc[f][0]*nf[f][0]+fc[f][1]*nf[f][1]+fc[f][2]*nf[f][2]);
	Ve *= (1.0/3.0);

	// on-demand field/grad C[h, :, colDOF] (INV4PI-scaled): find colDOF's element + face, build its geometry,
	// then evaluate the shared kernel at THIS element's centroid (incl. IMA mirrors).
	double Hc[3] = {0,0,0}, gHc[6] = {0,0,0,0,0,0};
	{
		// O(1) DOF->(hex,face): for pure hex m_elemDOFOffset[m_hexaElemIndices[h]] == 6h, so colDOF/6 is
		// the hex position and colDOF%6 the local face (no O(nHex) search).  Guard the pure-hex premise.
		int h_col = colDOF / 6, colF = colDOF % 6, colElem = -1;
		if(h_col >= 0 && h_col < nHex)
		{
			int ei = m_hexaElemIndices[h_col];
			if(m_elemDOFOffset[ei] == 6 * h_col) colElem = ei;   // premise holds -> direct map
		}
		radTPolyhedron* pc = (colElem >= 0) ? dynamic_cast<radTPolyhedron*>(g3dRelaxPtrVect[colElem]) : nullptr;
		if(pc && pc->AmOfFaces == 6)
		{
			const radTHandlePgnAndTrans& hpt = pc->VectHandlePgnAndTrans[colF];
			radTPolygon* pgn = hpt.PgnHndl.rep; radTrans* tr = hpt.TransHndl.rep;
			if(pgn && tr)
			{
				const radTVect2dVect& v2d = pgn->EdgePointsVector;
				double V4c[4][3]; double areaC = 0.0;
				for(int v = 0; v < 4; v++)
				{
					TVector3d P = tr->TrPoint(TVector3d(v2d[v].x, v2d[v].y, pgn->CoordZ));
					V4c[v][0] = P.x; V4c[v][1] = P.y; V4c[v][2] = P.z;
				}
				for(int tt = 0; tt < 2; tt++)
				{
					double ux=V4c[tt+1][0]-V4c[0][0], uy=V4c[tt+1][1]-V4c[0][1], uz=V4c[tt+1][2]-V4c[0][2];
					double wx=V4c[tt+2][0]-V4c[0][0], wy=V4c[tt+2][1]-V4c[0][1], wz=V4c[tt+2][2]-V4c[0][2];
					double rx=uy*wz-uz*wy, ry=uz*wx-ux*wz, rz=ux*wy-uy*wx;
					areaC += 0.5*std::sqrt(rx*rx+ry*ry+rz*rz);
				}
				const TVector3d cce = pc->CentrPoint;
				const double cenC[3] = {cce.x, cce.y, cce.z};
				const double ce3[3] = {ce.x, ce.y, ce.z};
				CentroidFieldGradFromFace(ce3, V4c, cenC, colElem == elemIdx, areaC, Hc, gHc);
				const double INV4PI = 1.0/(4.0*3.14159265358979323846);
				for(int k = 0; k < 3; k++) Hc[k] *= INV4PI;
				for(int k = 0; k < 6; k++) gHc[k] *= INV4PI;
			}
		}
	}

	const int lf = colDOF - off;
	const bool localFace = (lf >= 0 && lf < 6);
	double val = 0.0;
	if(t < 3)                                   // dipole-t
	{
		if(localFace) val += Ae[lf]*d[lf][t]/Ve;
		val -= chiH*Hc[t];
	}
	else if(t == 3)                             // monopole
	{
		if(localFace) val = Ae[lf];
	}
	else                                        // diagonal-quadrupole qq = t-4
	{
		int qq = t - 4;
		// SAME residual quadrupole eigenmode as BuildMomentSystemCore (identical Ae/d -> identical modes),
		// so the H-matrix (method 2) and dense (method 0/1) moment systems stay consistent.  Bm = phi_qq/Ae.
		double phiQ[6][6]; momentResidualEigenmodes(Ae, d, 6, phiQ);
		double Bm[6];
		for(int f = 0; f < 6; f++) Bm[f] = (Ae[f] > 1.0e-300) ? phiQ[qq][f]/Ae[f] : 0.0;
		if(localFace) val += Ae[lf]*Bm[lf];
		double cm[3] = {0,0,0};
		for(int f = 0; f < 6; f++) for(int k = 0; k < 3; k++) cm[k] += Ae[f]*nf[f][k]*Bm[f];
		if(localFace)
		{
			double cd = cm[0]*d[lf][0]+cm[1]*d[lf][1]+cm[2]*d[lf][2];
			val -= Ae[lf]*cd/Ve;
		}
		double Dm[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
		for(int f = 0; f < 6; f++) for(int ii = 0; ii < 3; ii++) for(int jj = 0; jj < 3; jj++)
			Dm[ii][jj] += Ae[f]*d[f][jj]*nf[f][ii]*Bm[f];
		double Dvec[6] = {Dm[0][0], Dm[1][1], Dm[2][2], Dm[0][1]+Dm[1][0], Dm[0][2]+Dm[2][0], Dm[1][2]+Dm[2][1]};
		for(int m = 0; m < 6; m++) val -= chiH*Dvec[m]*gHc[m];
	}
	return val;
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

			TVector3d V[4];
			for(int v = 0; v < nv && v < 4; v++)
				V[v] = tr->TrPoint(TVector3d(verts2d[v].x, verts2d[v].y, pgn->CoordZ));

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

// EIEM2 surface-charge kernels retired (Phase 3b, live/dead step C): Compute5x5BlockFast (wedge),
// Compute6x6BlockFast (hex), and ComputeMixedBlockFast (cross-DOF) are deleted -- the moment-yano
// formulation (BuildMomentSystemCore) is the sole surface-charge demag, and the method-2 HACApK
// path is MMM-only (3x3) for the remaining tetrahedron solves.

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
// Compute6x6BlockIMA, Compute6x6BlockMirrored, Compute6x6BlockMirroredTarget) deleted (Phase 3b).
// It was #if 0 dead code (kernel-based IMA replaced it); the moment-yano path adds IMA mirror
// images via CentroidFieldGradFromFace.

//-------------------------------------------------------------------------
// SetupInteractMatrix_IMA: Build IMA interaction matrix.
// The tet MMM path uses Compute3x3BlockFast. Surface-charge moment-yano handles IMA in
// BuildMomentSystemCore / CentroidFieldGradFromFace and skips this dense MSC matrix.
//-------------------------------------------------------------------------
int radTInteraction::SetupInteractMatrix_IMA(bool skipDenseMatrix)
{
	if(!m_imaEnabled)
	{
		std::cerr << "[Radia] Error: IMA not enabled" << std::endl;
		return 0;
	}


	// Check all elements have valid DOF (3=tet MMM, 5=wedge MSC, 6=hex MSC)
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
	// EIEM2 retirement (Phase 3b): any MSC surface-charge element (DOF 5/6) present?  (DOF is in {3,5,6}
	// here; the loop above errors on < 3.)  Used below to skip the dense MSC IMA build.
	bool hasMSC = !allTet;

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

	// For HACApK: skip dense matrix. Pure tet MMM computes entries on demand through
	// Compute3x3BlockFast; surface-charge moment-yano uses RadHACApKMomentSystem instead.
	// Reset geometry so it gets recomputed for the reduced IMA element set.
	if(skipDenseMatrix)
	{
		// Reset precomputed geometry so HACApK recomputes for the reduced IMA elements
		m_hexaGeomReady = false;
		m_wedgeGeomReady = false;
		m_tetraGeomReady = false;
		return 1;
	}

	// EIEM2 retirement (Phase 3b): skip the dense MSC IMA matrix.  moment+IMA assembles its own system
	// (BuildMomentSystemCore + CentroidFieldGradFromFace, which adds the IMA mirror images) and never reads
	// this matrix; mixed MMM+MSC and B-input-on-MSC are rejected fail-loud in MakeAutoRelax.  Pure-tet (MMM)
	// IMA still builds its 3x3 dense matrix below.  Geometry precomputed above is left intact (moment uses
	// on-the-fly momentFaceGeom; the precompute caches stay valid for the reduced IMA element set).
	if(hasMSC)
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

			// Fast path: 3DOF-3DOF tet (MMM) - Compute3x3BlockFast handles IMA inline
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
