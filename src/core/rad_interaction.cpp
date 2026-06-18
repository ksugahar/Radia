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
	// True if any relaxable element is a collocation surface-charge (Use6DOF_MSC) polyhedron --
	// a soft-iron hexahedron (6 faces) or wedge (5 faces).  The yano-type MSC demag for these is
	// removed; the soft iron must go through the FEEC HDiv-VIM (radia.vim.soft_iron_from_mesh).
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

	// Check if all elements are hexahedra (for ultra-fast pure-hexa path)
	bool allHex = hasMSCElements;
	for(int i = 0; i < AmOfMainElem && allHex; i++)
	{
		if(m_elemDOF[i] != 6) allHex = false;
	}

	// ULTRA-FAST PATH: Pure hexahedra without symmetry
	if(!hasSymmetry && allHex)
	{
		PrecomputeHexaGeometry();
		if(m_hexaGeomReady)
		{
			int nHex = (int)m_hexaElemIndices.size();

			ngcore::ParallelFor(ngcore::IntRange(nHex), [&](size_t hex_col)
			{
				int col = m_hexaElemIndices[(int)hex_col];
				int offset_col = m_elemDOFOffset[col];

				for(int hex_row = 0; hex_row < nHex; hex_row++)
				{
					int row = m_hexaElemIndices[hex_row];
					int offset_row = m_elemDOFOffset[row];

					double K_block[36];
					Compute6x6BlockFast(hex_row, (int)hex_col, K_block);

					// Copy to row-major flat matrix (both row-major, direct copy pattern)
					// K_block is [target][source] row-major, matrix is [target][source] row-major
					// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
					double* block = &m_flatInteractMatrix[(size_t)offset_row * m_totalDOF + offset_col];
					for(int i = 0; i < 6; i++)
					{
						for(int j = 0; j < 6; j++)
						{
							// Both row-major: K[i][j] at i*6+j -> A[target_i][source_j] at i*stride+j
							block[(size_t)i * m_totalDOF + j] = K_block[i * 6 + j];
						}
					}
				}
			});

			return 1;
		}
		// Fall through to MEDIUM PATH if geometry precompute failed
	}

	// Check if all MSC elements are wedges (for pure-wedge fast path)
	bool allWedge = hasMSCElements;
	for(int i = 0; i < AmOfMainElem && allWedge; i++)
	{
		if(m_elemDOF[i] != 5) allWedge = false;
	}

	// ULTRA-FAST PATH: Pure wedges without symmetry
	if(!hasSymmetry && allWedge)
	{
		PrecomputeWedgeGeometry();
		if(m_wedgeGeomReady)
		{
			int nWedge = (int)m_wedgeElemIndices.size();

			ngcore::ParallelFor(ngcore::IntRange(nWedge), [&](size_t w_col)
			{
				int col = m_wedgeElemIndices[(int)w_col];
				int offset_col = m_elemDOFOffset[col];

				for(int w_row = 0; w_row < nWedge; w_row++)
				{
					int row = m_wedgeElemIndices[w_row];
					int offset_row = m_elemDOFOffset[row];

					double K_block[25];
					Compute5x5BlockFast(w_row, (int)w_col, K_block);

					double* block = &m_flatInteractMatrix[(size_t)offset_row * m_totalDOF + offset_col];
					for(int i = 0; i < 5; i++)
						for(int j = 0; j < 5; j++)
							block[(size_t)i * m_totalDOF + j] = K_block[i * 5 + j];
				}
			});

			return 1;
		}
	}

	// MEDIUM PATH: MSC hexahedra without symmetry - uses OpenMP
	// This handles mixed hex/tetra meshes with 6x6 blocks
	if(!hasSymmetry && hasMSCElements)
	{
		// Use unified 1/(4*pi) constant for all MSC interactions
		// (ELF-compatible sign convention: K_ij / (4*pi))

		ngcore::ParallelFor(ngcore::IntRange(AmOfMainElem), [&](size_t col)
		{
			radTg3dRelax* elem_col = g3dRelaxPtrVect[(int)col];
			int dof_col = m_elemDOF[(int)col];
			int offset_col = m_elemDOFOffset[(int)col];

			// Check if source is MSC element (wedge: 5 DOF, hexahedron: 6 DOF)
			radTPolyhedron* poly_col = nullptr;
			if(dof_col >= 5)
			{
				poly_col = dynamic_cast<radTPolyhedron*>(elem_col);
			}

			for(int row = 0; row < AmOfMainElem; row++)
			{
				radTg3dRelax* elem_row = g3dRelaxPtrVect[row];
				int dof_row = m_elemDOF[row];
				int offset_row = m_elemDOFOffset[row];

				// ROW-MAJOR: A(row, col) at index [row * m_totalDOF + col]
				// A[target][source] format: ELF-compatible, BiCGSTAB/HACApK-optimal
				// CRITICAL: Use size_t cast to avoid int32 overflow for DOF > 46340
				double* block = &m_flatInteractMatrix[(size_t)offset_row * m_totalDOF + offset_col];

				// Check if target is MSC element (wedge: 5 DOF, hexahedron: 6 DOF)
				radTPolyhedron* poly_row = nullptr;
				if(dof_row >= 5)
				{
					poly_row = dynamic_cast<radTPolyhedron*>(elem_row);
				}

				if(dof_row == 3 && dof_col == 3)
				{
					// 3x3 N-matrix block: H-field at row center from col magnetization
					// PreRelax mode computes ALL 3 unit responses in ONE call:
					//   Field.B = dH/dMx, Field.H = dH/dMy, Field.A = dH/dMz
					TVector3d ObsPoiVect = elem_row->ReturnCentrPoint();

					radTField Field(FieldKeyInteract, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
					Field.AmOfIntrctElemWithSym = AmOfElemWithSym;
					elem_col->B_comp(&Field);

					// Store N-matrix (row-major): N[row_comp][col_comp]
					// Field.B = response to Mx, Field.H = response to My, Field.A = response to Mz
					block[(size_t)0 * m_totalDOF + 0] = Field.B.x;  // dHx/dMx
					block[(size_t)0 * m_totalDOF + 1] = Field.H.x;  // dHx/dMy
					block[(size_t)0 * m_totalDOF + 2] = Field.A.x;  // dHx/dMz
					block[(size_t)1 * m_totalDOF + 0] = Field.B.y;  // dHy/dMx
					block[(size_t)1 * m_totalDOF + 1] = Field.H.y;  // dHy/dMy
					block[(size_t)1 * m_totalDOF + 2] = Field.A.y;  // dHy/dMz
					block[(size_t)2 * m_totalDOF + 0] = Field.B.z;  // dHz/dMx
					block[(size_t)2 * m_totalDOF + 1] = Field.H.z;  // dHz/dMy
					block[(size_t)2 * m_totalDOF + 2] = Field.A.z;  // dHz/dMz
				}
				else if(dof_row == 6 && dof_col == 6 && poly_row && poly_col)
				{
					// 6x6 block: MSC hexahedron to MSC hexahedron
					for(int face_i = 0; face_i < 6; face_i++)
					{
						// Yano evaluation point: midpoint between face center and element center
						TVector3d EvalPt;
						EvalPt.x = 0.5 * (poly_row->FaceCenter[face_i].x + poly_row->CentrPoint.x);
						EvalPt.y = 0.5 * (poly_row->FaceCenter[face_i].y + poly_row->CentrPoint.y);
						EvalPt.z = 0.5 * (poly_row->FaceCenter[face_i].z + poly_row->CentrPoint.z);

						for(int face_j = 0; face_j < 6; face_j++)
						{
							// Field from unit sigma on face j
							TVector3d H_face = poly_col->FieldFromQuadFace(EvalPt, face_j, 1.0);

							// Point charge contribution: m = -sigma * area
							double unit_point_charge = -1.0 * poly_col->FaceArea[face_j];
							TVector3d H_point = poly_col->FieldFromPointCharge(EvalPt, unit_point_charge);

							TVector3d H_total;
							H_total.x = H_face.x + H_point.x;
							H_total.y = H_face.y + H_point.y;
							H_total.z = H_face.z + H_point.z;

							// K_ij = normal_i dot H_total
							double K_ij = H_total.x * poly_row->FaceNormal[face_i].x +
							              H_total.y * poly_row->FaceNormal[face_i].y +
							              H_total.z * poly_row->FaceNormal[face_i].z;

							// ROW-MAJOR: A[target_face_i][source_face_j] at face_i*stride+face_j
							// Store K/(4pi) - the solver will negate when building system matrix
							block[(size_t)face_i * m_totalDOF + face_j] = K_ij * RadConst::INV_FOUR_PI;
						}
					}
				}
				else if(dof_row == 3 && dof_col == 6 && poly_col)
				{
					// 3x6 block: tetrahedron from MSC hexahedron
					TVector3d ObsPoiVect = elem_row->ReturnCentrPoint();

					for(int face_j = 0; face_j < 6; face_j++)
					{
						TVector3d H_face = poly_col->FieldFromQuadFace(ObsPoiVect, face_j, 1.0);
						double unit_point_charge = -1.0 * poly_col->FaceArea[face_j];
						TVector3d H_point = poly_col->FieldFromPointCharge(ObsPoiVect, unit_point_charge);

						TVector3d H_total;
						H_total.x = H_face.x + H_point.x;
						H_total.y = H_face.y + H_point.y;
						H_total.z = H_face.z + H_point.z;

						// ROW-MAJOR: A[target_i][source_face_j] at target_i*stride+face_j
						block[(size_t)0 * m_totalDOF + face_j] = H_total.x * RadConst::INV_FOUR_PI;
						block[(size_t)1 * m_totalDOF + face_j] = H_total.y * RadConst::INV_FOUR_PI;
						block[(size_t)2 * m_totalDOF + face_j] = H_total.z * RadConst::INV_FOUR_PI;
					}
				}
				else if(dof_row == 6 && dof_col == 3 && poly_row)
				{
					// 6x3 block: MSC hexahedron from 3DOF polyhedron (tetra/wedge)
					// K(face_i, Mj) = normal_i · H_j where H_j is H-field from unit M in direction j
					// PreRelax mode computes ALL 3 unit responses in ONE call:
					//   Field.B = dH/dMx, Field.H = dH/dMy, Field.A = dH/dMz

					for(int face_i = 0; face_i < dof_row; face_i++)
					{
						// Yano evaluation point
						TVector3d EvalPt;
						EvalPt.x = 0.5 * (poly_row->FaceCenter[face_i].x + poly_row->CentrPoint.x);
						EvalPt.y = 0.5 * (poly_row->FaceCenter[face_i].y + poly_row->CentrPoint.y);
						EvalPt.z = 0.5 * (poly_row->FaceCenter[face_i].z + poly_row->CentrPoint.z);

						TVector3d& n = poly_row->FaceNormal[face_i];

						radTField Field(FieldKeyInteract, CompCriterium, EvalPt, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
						Field.AmOfIntrctElemWithSym = AmOfElemWithSym;
						elem_col->B_comp(&Field);

						// K(face_i, m_dir) = normal · H_response
						// Field.B = response to Mx, Field.H = response to My, Field.A = response to Mz
						double K_Mx = n.x * Field.B.x + n.y * Field.B.y + n.z * Field.B.z;
						double K_My = n.x * Field.H.x + n.y * Field.H.y + n.z * Field.H.z;
						double K_Mz = n.x * Field.A.x + n.y * Field.A.y + n.z * Field.A.z;

						// ROW-MAJOR: A[target_face_i][source_m_dir]
						// No INV_FOUR_PI - already included in field computation
						block[(size_t)face_i * m_totalDOF + 0] = K_Mx;
						block[(size_t)face_i * m_totalDOF + 1] = K_My;
						block[(size_t)face_i * m_totalDOF + 2] = K_Mz;
					}
				}
				else if(dof_row >= 5 && dof_col >= 5 && poly_row && poly_col)
				{
					// NxM block: MSC-to-MSC (covers 5x5, 5x6, 6x5 wedge/hex combinations)
					// Uses FieldFromFace for generalized tri/quad face handling
					int nFacesRow = dof_row;  // 5 for wedge, 6 for hex
					int nFacesCol = dof_col;

					for(int face_i = 0; face_i < nFacesRow; face_i++)
					{
						// Yano evaluation point: midpoint between face center and element center
						TVector3d EvalPt;
						EvalPt.x = 0.5 * (poly_row->FaceCenter[face_i].x + poly_row->CentrPoint.x);
						EvalPt.y = 0.5 * (poly_row->FaceCenter[face_i].y + poly_row->CentrPoint.y);
						EvalPt.z = 0.5 * (poly_row->FaceCenter[face_i].z + poly_row->CentrPoint.z);

						for(int face_j = 0; face_j < nFacesCol; face_j++)
						{
							// Field from unit sigma on face j (handles both tri and quad)
							TVector3d H_face = poly_col->FieldFromFace(EvalPt, face_j, 1.0);

							// Point charge contribution: m = -sigma * area
							double unit_point_charge = -1.0 * poly_col->FaceArea[face_j];
							TVector3d H_point = poly_col->FieldFromPointCharge(EvalPt, unit_point_charge);

							TVector3d H_total;
							H_total.x = H_face.x + H_point.x;
							H_total.y = H_face.y + H_point.y;
							H_total.z = H_face.z + H_point.z;

							// K_ij = normal_i dot H_total
							double K_ij = H_total.x * poly_row->FaceNormal[face_i].x +
							              H_total.y * poly_row->FaceNormal[face_i].y +
							              H_total.z * poly_row->FaceNormal[face_i].z;

							block[(size_t)face_i * m_totalDOF + face_j] = K_ij * RadConst::INV_FOUR_PI;
						}
					}
				}
				else if(dof_row == 3 && dof_col == 5 && poly_col)
				{
					// 3x5 block: tetrahedron from MSC wedge
					TVector3d ObsPoiVect = elem_row->ReturnCentrPoint();

					for(int face_j = 0; face_j < 5; face_j++)
					{
						TVector3d H_face = poly_col->FieldFromFace(ObsPoiVect, face_j, 1.0);
						double unit_point_charge = -1.0 * poly_col->FaceArea[face_j];
						TVector3d H_point = poly_col->FieldFromPointCharge(ObsPoiVect, unit_point_charge);

						TVector3d H_total;
						H_total.x = H_face.x + H_point.x;
						H_total.y = H_face.y + H_point.y;
						H_total.z = H_face.z + H_point.z;

						block[(size_t)0 * m_totalDOF + face_j] = H_total.x * RadConst::INV_FOUR_PI;
						block[(size_t)1 * m_totalDOF + face_j] = H_total.y * RadConst::INV_FOUR_PI;
						block[(size_t)2 * m_totalDOF + face_j] = H_total.z * RadConst::INV_FOUR_PI;
					}
				}
				else if(dof_row == 5 && dof_col == 3 && poly_row)
				{
					// 5x3 block: MSC wedge from 3DOF polyhedron (tetra)
					radTPolyhedron* poly_col_3dof = dynamic_cast<radTPolyhedron*>(elem_col);
					TVector3d orig_magn(0., 0., 0.);
					if(poly_col_3dof) {
						orig_magn = poly_col_3dof->Magn;
					}

					for(int face_i = 0; face_i < 5; face_i++)
					{
						TVector3d EvalPt;
						EvalPt.x = 0.5 * (poly_row->FaceCenter[face_i].x + poly_row->CentrPoint.x);
						EvalPt.y = 0.5 * (poly_row->FaceCenter[face_i].y + poly_row->CentrPoint.y);
						EvalPt.z = 0.5 * (poly_row->FaceCenter[face_i].z + poly_row->CentrPoint.z);

						TVector3d& n = poly_row->FaceNormal[face_i];

						for(int m_dir = 0; m_dir < 3; m_dir++)
						{
							TVector3d unit_M(0., 0., 0.);
							if(m_dir == 0) unit_M.x = 1.0;
							else if(m_dir == 1) unit_M.y = 1.0;
							else unit_M.z = 1.0;

							if(poly_col_3dof) {
								poly_col_3dof->Magn = unit_M;
							}

							radTField Field(FieldKeyInteract, CompCriterium, EvalPt, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
							Field.AmOfIntrctElemWithSym = AmOfElemWithSym;
							elem_col->B_comp(&Field);

							double K_val = n.x * Field.H.x + n.y * Field.H.y + n.z * Field.H.z;
							block[(size_t)face_i * m_totalDOF + m_dir] = K_val;
						}
					}

					if(poly_col_3dof) {
						poly_col_3dof->Magn = orig_magn;
					}
				}
				else
				{
					// Zero out unknown blocks
					for(int i = 0; i < dof_row; i++)
					{
						for(int j = 0; j < dof_col; j++)
						{
							// ROW-MAJOR: A[i][j] at i*stride+j
							block[(size_t)i * m_totalDOF + j] = 0.0;
						}
					}
				}
			}
		});
		return 1;
	}

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
			else if(dof_row == 3 && dof_col >= 5)
			{
				// 3xN block: Field at tetrahedron (3 DOF) center from MSC element (5 or 6 DOF)
				// K(Mk, face_j) = H_field_k at tetra center due to unit sigma on face j
				// Matrix stores: H_field(k) / (4*pi) (following ELF convention)

				radTPolyhedron* poly_col = dynamic_cast<radTPolyhedron*>(elem_col);
				if(poly_col && poly_col->Use6DOF_MSC)
				{
					TVector3d InitObsPoiVect = MainTransPtrArray[row]->TrPoint(elem_row->ReturnCentrPoint());

					for(int face_j = 0; face_j < dof_col; face_j++)
					{
						TVector3d H_total(0., 0., 0.);

						for(unsigned tr = 0; tr < TransPtrVect.size(); tr++)
						{
							// Transform obs point from world frame to column element's local frame
							TVector3d ObsPoiVect = TransPtrVect[tr]->TrPoint_inv(InitObsPoiVect);

							// Field from source at observation point (both in column's local frame)
							// FieldFromFace handles both tri and quad faces
							TVector3d H_face = poly_col->FieldFromFace(ObsPoiVect, face_j, 1.0);
							double unit_point_charge = -1.0 * poly_col->FaceArea[face_j];
							TVector3d H_point = poly_col->FieldFromPointCharge(ObsPoiVect, unit_point_charge);

							// Sum field contributions (in column's local frame)
							TVector3d H_local;
							H_local.x = H_face.x + H_point.x;
							H_local.y = H_face.y + H_point.y;
							H_local.z = H_face.z + H_point.z;

							// Transform field from column's local frame back to world frame
							TVector3d H_world = TransPtrVect[tr]->TrAxialVect(H_local);

							H_total.x += H_world.x;
							H_total.y += H_world.y;
							H_total.z += H_world.z;
						}

						// Transform field from world frame to row element's local frame
						TVector3d H_final = MainTransPtrArray[row]->TrAxialVect_inv(H_total);

						// Store in block (ROW-MAJOR): A[i][j] at [i * stride + j]
						block[(size_t)0 * m_totalDOF + face_j] = H_final.x * RadConst::INV_FOUR_PI;
						block[(size_t)1 * m_totalDOF + face_j] = H_final.y * RadConst::INV_FOUR_PI;
						block[(size_t)2 * m_totalDOF + face_j] = H_final.z * RadConst::INV_FOUR_PI;
					}
				}
			}
			else if(dof_row >= 5 && dof_col == 3)
			{
				// Nx3 block: Field at MSC element (5 or 6 DOF) eval points from tetrahedron (3 DOF)
				// K(face_i, Mj) = normal_i dot N_mat(:, j)

				radTPolyhedron* poly_row = dynamic_cast<radTPolyhedron*>(elem_row);
				if(poly_row && poly_row->Use6DOF_MSC)
				{
					for(int face_i = 0; face_i < dof_row; face_i++)
					{
						// Yano evaluation point: midpoint between face center and element center
						TVector3d EvalPt;
						EvalPt.x = 0.5 * (poly_row->FaceCenter[face_i].x + poly_row->CentrPoint.x);
						EvalPt.y = 0.5 * (poly_row->FaceCenter[face_i].y + poly_row->CentrPoint.y);
						EvalPt.z = 0.5 * (poly_row->FaceCenter[face_i].z + poly_row->CentrPoint.z);

						TVector3d InitObsPoiVect = MainTransPtrArray[row]->TrPoint(EvalPt);

						// N_mat stores the 3x3 demagnetization tensor: H = -N*M/(4*pi)
						// SubMatrix.Str0 = column 0 (response to Mx)
						// SubMatrix.Str1 = column 1 (response to My)
						// SubMatrix.Str2 = column 2 (response to Mz)
						TMatrix3d SubMatrix(TVector3d(0., 0., 0.), TVector3d(0., 0., 0.), TVector3d(0., 0., 0.));
						TMatrix3d BufSubMatrix;

						for(unsigned tr = 0; tr < TransPtrVect.size(); tr++)
						{
							TVector3d ObsPoiVect = TransPtrVect[tr]->TrPoint_inv(InitObsPoiVect);

							radTField Field(FieldKeyInteract, CompCriterium, ObsPoiVect, TVector3d(0., 0., 0.),
							                TVector3d(0., 0., 0.), TVector3d(0., 0., 0.), TVector3d(0., 0., 0.), 0.);
							Field.AmOfIntrctElemWithSym = AmOfElemWithSym;

							elem_col->B_comp(&Field);

							BufSubMatrix.Str0 = Field.B;
							BufSubMatrix.Str1 = Field.H;
							BufSubMatrix.Str2 = Field.A;

							TransPtrVect[tr]->TrMatrix(BufSubMatrix);
							SubMatrix += BufSubMatrix;
						}

						MainTransPtrArray[row]->TrMatrix_inv(SubMatrix);

						// Get face normal (outward)
						TVector3d& n = poly_row->FaceNormal[face_i];

						// K(face_i, Mj) = normal · N_mat(:, j) / (4*pi)
						// SubMatrix stores the 3x3 demagnetization response
						// Str0 = [dHx/dMx, dHy/dMx, dHz/dMx] (column for Mx)
						// Str1 = [dHx/dMy, dHy/dMy, dHz/dMy] (column for My)
						// Str2 = [dHx/dMz, dHy/dMz, dHz/dMz] (column for Mz)
						double K_face_Mx = n.x * SubMatrix.Str0.x + n.y * SubMatrix.Str0.y + n.z * SubMatrix.Str0.z;
						double K_face_My = n.x * SubMatrix.Str1.x + n.y * SubMatrix.Str1.y + n.z * SubMatrix.Str1.z;
						double K_face_Mz = n.x * SubMatrix.Str2.x + n.y * SubMatrix.Str2.y + n.z * SubMatrix.Str2.z;

						// Store in block (ROW-MAJOR): A[face_i][col] at [face_i * stride + col]
						// row is face_i (0-5), col is component (0,1,2 = Mx,My,Mz)
						// Sign convention: +K/(4*pi) for hex-tetra (following ELF)
						// CRITICAL: Use size_t cast for indexing with m_totalDOF
						block[(size_t)face_i * m_totalDOF + 0] = K_face_Mx * RadConst::INV_FOUR_PI;  // (face_i, Mx)
						block[(size_t)face_i * m_totalDOF + 1] = K_face_My * RadConst::INV_FOUR_PI;  // (face_i, My)
						block[(size_t)face_i * m_totalDOF + 2] = K_face_Mz * RadConst::INV_FOUR_PI;  // (face_i, Mz)
					}
				}
			}
			else if(dof_row >= 5 && dof_col >= 5)
			{
				// NxM block: MSC-to-MSC (covers 5x5, 5x6, 6x5, 6x6 combinations)
				// K(face_i, face_j) = normal_i dot H_field(eval_pt_i, src_face_j)

				radTPolyhedron* poly_row = dynamic_cast<radTPolyhedron*>(elem_row);
				radTPolyhedron* poly_col = dynamic_cast<radTPolyhedron*>(elem_col);

				if(poly_row && poly_row->Use6DOF_MSC && poly_col && poly_col->Use6DOF_MSC)
				{
					for(int face_i = 0; face_i < dof_row; face_i++)
					{
						// Yano evaluation point: midpoint between face center and element center
						// FaceCenter and CentrPoint are already in GLOBAL frame
						TVector3d EvalPt;
						EvalPt.x = 0.5 * (poly_row->FaceCenter[face_i].x + poly_row->CentrPoint.x);
						EvalPt.y = 0.5 * (poly_row->FaceCenter[face_i].y + poly_row->CentrPoint.y);
						EvalPt.z = 0.5 * (poly_row->FaceCenter[face_i].z + poly_row->CentrPoint.z);

						// EvalPt is already in GLOBAL frame - no transform needed
						TVector3d InitObsPoiVect = EvalPt;

						for(int face_j = 0; face_j < dof_col; face_j++)
						{
							TVector3d H_total(0., 0., 0.);

							for(unsigned tr = 0; tr < TransPtrVect.size(); tr++)
							{
								// Transform obs point from GLOBAL frame to column element's local frame
								TVector3d ObsPoiVect = TransPtrVect[tr]->TrPoint_inv(InitObsPoiVect);

								// FieldFromFace handles both tri and quad faces
								TVector3d H_face = poly_col->FieldFromFace(ObsPoiVect, face_j, 1.0);
								double unit_point_charge = -1.0 * poly_col->FaceArea[face_j];
								TVector3d H_point = poly_col->FieldFromPointCharge(ObsPoiVect, unit_point_charge);

								// Sum field contributions (in column's local frame)
								TVector3d H_local;
								H_local.x = H_face.x + H_point.x;
								H_local.y = H_face.y + H_point.y;
								H_local.z = H_face.z + H_point.z;

								// Transform field from column's local frame back to GLOBAL frame
								TVector3d H_world = TransPtrVect[tr]->TrAxialVect(H_local);

								H_total.x += H_world.x;
								H_total.y += H_world.y;
								H_total.z += H_world.z;
							}

							// K_ij = normal_i dot H_total (both in GLOBAL frame)
							// FaceNormal is already in GLOBAL frame (from SetupFaceGeometry)
							double K_ij = H_total.x * poly_row->FaceNormal[face_i].x +
							              H_total.y * poly_row->FaceNormal[face_i].y +
							              H_total.z * poly_row->FaceNormal[face_i].z;

							// Store K_ij / (4*pi) (ROW-MAJOR): A[face_i][face_j] at [face_i * stride + face_j]
							double K_val = K_ij * RadConst::INV_FOUR_PI;
							block[(size_t)face_i * m_totalDOF + face_j] = K_val;
						}
					}
				}
			}
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
			else if(dof >= 5)
			{
				// MSC element (5=wedge, 6=hex): compute H_ext dot n at each face
				// External field is evaluated at element positions
				radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(elem);
				if(poly && poly->Use6DOF_MSC)
				{
					for(int face_i = 0; face_i < dof; face_i++)
					{
						// Eval point for face i (midpoint between face center and element center)
						// This is in element's local coordinates
						TVector3d EvalPt;
						EvalPt.x = 0.5 * (poly->FaceCenter[face_i].x + poly->CentrPoint.x);
						EvalPt.y = 0.5 * (poly->FaceCenter[face_i].y + poly->CentrPoint.y);
						EvalPt.z = 0.5 * (poly->FaceCenter[face_i].z + poly->CentrPoint.z);

						// Transform EvalPt from element's local coords to world coords
						TVector3d WorldEvalPt = MainTransPtrArray[StrNo]->TrPoint(EvalPt);

						// Compute H_ext at eval point from all external sources
						// Use same transform handling as 3DOF elements
						TVector3d H_world(0., 0., 0.);
						for(int ExtElNo = 0; ExtElNo < AmOfExtElem; ExtElNo++)
						{
							FillInTransPtrVectForElem(ExtElNo, 'E');
							radTg3d* ExtElPtr = g3dExternPtrVect[ExtElNo];

							TVector3d BufVect(0., 0., 0.);
							for(unsigned t = 0; t < TransPtrVect.size(); t++)
							{
								// Transform obs point to external element's local coords
								TVector3d ObsPoiVect = TransPtrVect[t]->TrPoint_inv(WorldEvalPt);
								radTField Field(FieldKeyExtern, CompCriterium, ObsPoiVect, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
								ExtElPtr->B_comp(&Field);
								// Transform field back to world coords using pseudo-vector transformation
								// H field is a pseudo-vector: H' = det(T) * T * H
								BufVect += TransPtrVect[t]->TrAxialVect(Field.H);
							}
							H_world += BufVect;
							EmptyTransPtrVect();
						}

						// FaceNormal is stored in GLOBAL coordinates (from SetupFaceGeometry)
						// H_world is already in global coordinates
						// Compute dot product directly in global coordinates (ELF-compatible)
						// Note: Do NOT transform H to local coords - FaceNormal is already global!
						double H_dot_n = H_world.x * poly->FaceNormal[face_i].x +
						                 H_world.y * poly->FaceNormal[face_i].y +
						                 H_world.z * poly->FaceNormal[face_i].z;

						m_flatExternFieldArray[offset + face_i] += H_dot_n;
					}
				}
			}
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
				else if(dof == 6)
				{
					// MSC hexahedron: compute H_ext dot n at each face
					radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(elem);
					if(poly && poly->Use6DOF_MSC)
					{
						for(int face_i = 0; face_i < 6; face_i++)
						{
							// Eval point for face i (midpoint between face center and element center)
							// This is in element's local coordinates
							TVector3d EvalPt;
							EvalPt.x = 0.5 * (poly->FaceCenter[face_i].x + poly->CentrPoint.x);
							EvalPt.y = 0.5 * (poly->FaceCenter[face_i].y + poly->CentrPoint.y);
							EvalPt.z = 0.5 * (poly->FaceCenter[face_i].z + poly->CentrPoint.z);

							// Transform EvalPt to world coords (same as 3DOF code on line 1520)
							TVector3d WorldEvalPt = MainTransPtrArray[StrNo]->TrPoint(EvalPt);

							// Compute H_ext at world eval point
							// B_genComp handles internal transforms recursively
							radTField Field(FieldKeyExtern, CompCriterium, WorldEvalPt, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
							(static_cast<radTg3d*>(MoreExtSourceHandle.rep))->B_genComp(&Field);

							// FaceNormal is stored in GLOBAL coordinates (from SetupFaceGeometry)
							// Field.H from B_genComp is already in global coordinates
							// Compute dot product directly in global coordinates (ELF-compatible)
							// Note: Do NOT transform H to local coords - FaceNormal is already global!
							double H_dot_n = Field.H.x * poly->FaceNormal[face_i].x +
							                 Field.H.y * poly->FaceNormal[face_i].y +
							                 Field.H.z * poly->FaceNormal[face_i].z;

							m_flatExternFieldArray[offset + face_i] = H_dot_n;
						}
					}
				}
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
// Reference: ELF-style optimization (same as RadHACApKMSCManager::PrecomputeGeometry3DOF)
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
// Reference: ELF-style optimization (same as RadHACApKMSCManager::Compute3x3BlockFast)
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
	m_hexaEvalPoints.resize(nHex * 6 * 3);  // 6 faces, xyz
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

			// Store Yano evaluation point: midpoint(face_center, element_center)
			int epIdx = (h * 6 + f) * 3;
			m_hexaEvalPoints[epIdx + 0] = 0.5 * (poly->FaceCenter[f].x + poly->CentrPoint.x);
			m_hexaEvalPoints[epIdx + 1] = 0.5 * (poly->FaceCenter[f].y + poly->CentrPoint.y);
			m_hexaEvalPoints[epIdx + 2] = 0.5 * (poly->FaceCenter[f].z + poly->CentrPoint.z);

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
	m_wedgeEvalPoints.resize(nWedge * 5 * 3);
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
			m_wedgeEvalPoints[(w*5+f)*3+0] = 0.5*(poly->FaceCenter[f].x + poly->CentrPoint.x);
			m_wedgeEvalPoints[(w*5+f)*3+1] = 0.5*(poly->FaceCenter[f].y + poly->CentrPoint.y);
			m_wedgeEvalPoints[(w*5+f)*3+2] = 0.5*(poly->FaceCenter[f].z + poly->CentrPoint.z);

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

//=========================================================================
// Compute5x5BlockFast: Fast 5x5 interaction block for wedges (MSC)
// Same pattern as Compute6x6BlockFast: Yano eval points,
// face-triangle decomposition, IMA inline with scalar sign.
//=========================================================================

void radTInteraction::Compute5x5BlockFast(int wedge_i, int wedge_j, double* K_mat) const
{
	std::memset(K_mat, 0, 25 * sizeof(double));
	if(!m_wedgeGeomReady) return;

	int nWedge = (int)m_wedgeElemIndices.size();
	if(wedge_i < 0 || wedge_i >= nWedge || wedge_j < 0 || wedge_j >= nWedge) return;

	const double src_center[3] = {m_wedgeCenters[wedge_j*3+0],
	                               m_wedgeCenters[wedge_j*3+1],
	                               m_wedgeCenters[wedge_j*3+2]};

	for(int fi = 0; fi < 5; fi++)
	{
		int epIdx = (wedge_i*5+fi)*3;
		const double obs[3] = {m_wedgeEvalPoints[epIdx+0], m_wedgeEvalPoints[epIdx+1], m_wedgeEvalPoints[epIdx+2]};
		int fnIdx_i = (wedge_i*5+fi)*3;
		const double n_i[3] = {m_wedgeFaceNormals[fnIdx_i+0], m_wedgeFaceNormals[fnIdx_i+1], m_wedgeFaceNormals[fnIdx_i+2]};

		for(int fj = 0; fj < 5; fj++)
		{
			double H_total[3] = {0,0,0};

			// Original source: triangles of face fj
			int triOff = m_wedgeTriOffset[wedge_j*5+fj];
			int numTris = m_wedgeFaceNumTris[wedge_j*5+fj];
			for(int t = 0; t < numTris; t++)
			{
				int tvIdx = (wedge_j*WEDGE_MAX_TRIS + triOff + t) * 3 * 3;
				const double* V0 = &m_wedgeTriVertices[tvIdx+0];
				const double* V1 = &m_wedgeTriVertices[tvIdx+3];
				const double* V2 = &m_wedgeTriVertices[tvIdx+6];
				double sign_tri = m_wedgeTriSigns[wedge_j*WEDGE_MAX_TRIS + triOff + t];
				double H_tri[3];
				FieldFromChargedTriangleLocal(obs, V0, V1, V2, sign_tri, H_tri);
				H_total[0] += H_tri[0]; H_total[1] += H_tri[1]; H_total[2] += H_tri[2];
			}

			// Point charge cancellation
			double area_j = m_wedgeFaceAreas[wedge_j*5+fj];
			{
				double r[3] = {obs[0]-src_center[0], obs[1]-src_center[1], obs[2]-src_center[2]};
				double d2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
				if(d2 > 1e-30) {
					double id3 = 1.0 / (sqrt(d2)*d2);
					double c = -area_j * id3;
					H_total[0] += c*r[0]; H_total[1] += c*r[1]; H_total[2] += c*r[2];
				}
			}

			// IMA: mirrored source contributions (scalar sign, same as hex)
			if(m_imaEnabled)
			{
				auto addMirrorWedge = [&](int mirrorAxis, int sign) {
					double imaSign = (double)sign;
					double mir_center[3] = {src_center[0], src_center[1], src_center[2]};
					if(mirrorAxis & IMA_X) mir_center[0] = -mir_center[0];
					if(mirrorAxis & IMA_Y) mir_center[1] = -mir_center[1];
					if(mirrorAxis & IMA_Z) mir_center[2] = -mir_center[2];

					int numMir = 0;
					if(mirrorAxis & IMA_X) numMir++;
					if(mirrorAxis & IMA_Y) numMir++;
					if(mirrorAxis & IMA_Z) numMir++;
					bool flipW = (numMir % 2 == 1);

					for(int t = 0; t < numTris; t++)
					{
						int tvIdx = (wedge_j*WEDGE_MAX_TRIS + triOff + t) * 3 * 3;
						double V0[3] = {m_wedgeTriVertices[tvIdx+0], m_wedgeTriVertices[tvIdx+1], m_wedgeTriVertices[tvIdx+2]};
						double V1[3] = {m_wedgeTriVertices[tvIdx+3], m_wedgeTriVertices[tvIdx+4], m_wedgeTriVertices[tvIdx+5]};
						double V2[3] = {m_wedgeTriVertices[tvIdx+6], m_wedgeTriVertices[tvIdx+7], m_wedgeTriVertices[tvIdx+8]};
						if(mirrorAxis & IMA_X) { V0[0]=-V0[0]; V1[0]=-V1[0]; V2[0]=-V2[0]; }
						if(mirrorAxis & IMA_Y) { V0[1]=-V0[1]; V1[1]=-V1[1]; V2[1]=-V2[1]; }
						if(mirrorAxis & IMA_Z) { V0[2]=-V0[2]; V1[2]=-V1[2]; V2[2]=-V2[2]; }
						if(flipW) { for(int k=0;k<3;k++) std::swap(V1[k],V2[k]); }

						double st = m_wedgeTriSigns[wedge_j*WEDGE_MAX_TRIS + triOff + t];
						double H_tri[3];
						FieldFromChargedTriangleLocal(obs, V0, V1, V2, st, H_tri);
						H_total[0] += imaSign*H_tri[0]; H_total[1] += imaSign*H_tri[1]; H_total[2] += imaSign*H_tri[2];
					}

					// Mirror point charge
					double r[3] = {obs[0]-mir_center[0], obs[1]-mir_center[1], obs[2]-mir_center[2]};
					double d2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
					if(d2 > 1e-30) {
						double id3 = 1.0/(sqrt(d2)*d2);
						double c = -area_j * id3 * imaSign;
						H_total[0] += c*r[0]; H_total[1] += c*r[1]; H_total[2] += c*r[2];
					}
				};

				bool hasX = (m_imaSymmetry & IMA_X) != 0;
				bool hasY = (m_imaSymmetry & IMA_Y) != 0;
				bool hasZ = (m_imaSymmetry & IMA_Z) != 0;
				if(hasX) addMirrorWedge(IMA_X, m_imaSignX);
				if(hasY) addMirrorWedge(IMA_Y, m_imaSignY);
				if(hasZ) addMirrorWedge(IMA_Z, m_imaSignZ);
				if(hasX && hasY) addMirrorWedge(IMA_XY, m_imaSignX*m_imaSignY);
				if(hasX && hasZ) addMirrorWedge(IMA_XZ, m_imaSignX*m_imaSignZ);
				if(hasY && hasZ) addMirrorWedge(IMA_YZ, m_imaSignY*m_imaSignZ);
				if(hasX && hasY && hasZ) addMirrorWedge(IMA_XYZ, m_imaSignX*m_imaSignY*m_imaSignZ);
			}

			double K_ij = (n_i[0]*H_total[0] + n_i[1]*H_total[1] + n_i[2]*H_total[2]) * RadConst::INV_FOUR_PI;
			K_mat[fi*5+fj] = K_ij;
		}
	}
}

//=========================================================================
// Compute6x6BlockFast: Fast 6x6 interaction block for hexahedra
// Uses pre-computed geometry (avoiding FieldFromQuadFace overhead)
// Reference: Yano MSC method
//
// When IMA is enabled (m_imaEnabled=true), this function computes:
//   K[i,j] = field at target i from original source j + field from mirrored source j
// This is the kernel-based IMA approach (user suggestion 2026-01-31):
// - Mirroring is done directly in the kernel, not via virtual elements/DOFs
// - No permutation arrays needed - coordinates are mirrored directly
// - HACApK compatible - kernel returns correct IMA value directly
//=========================================================================

void radTInteraction::Compute6x6BlockFast(int hex_i, int hex_j, double* K_mat) const
{
	std::memset(K_mat, 0, 36 * sizeof(double));

	if(!m_hexaGeomReady) return;

	int nHex = (int)m_hexaElemIndices.size();
	if(hex_i < 0 || hex_i >= nHex || hex_j < 0 || hex_j >= nHex) return;

	// Source element center (for point charge)
	const double src_center_orig[3] = {m_hexaCenters[hex_j * 3 + 0],
	                                   m_hexaCenters[hex_j * 3 + 1],
	                                   m_hexaCenters[hex_j * 3 + 2]};

	// Check if pre-computed triangle data is available
	const bool usePrecomputed = m_hexaTriDataReady;

	// For each target face i
	for(int face_i = 0; face_i < 6; face_i++)
	{
		// Yano evaluation point for target face
		int epIdx = (hex_i * 6 + face_i) * 3;
		const double obs[3] = {m_hexaEvalPoints[epIdx + 0],
		                       m_hexaEvalPoints[epIdx + 1],
		                       m_hexaEvalPoints[epIdx + 2]};

		// Target face normal
		int fnIdx_i = (hex_i * 6 + face_i) * 3;
		const double n_i[3] = {m_hexaFaceNormals[fnIdx_i + 0],
		                       m_hexaFaceNormals[fnIdx_i + 1],
		                       m_hexaFaceNormals[fnIdx_i + 2]};

		// For each source face j
		for(int face_j = 0; face_j < 6; face_j++)
		{
			// Field from unit sigma on source face j
			double H_total[3] = {0.0, 0.0, 0.0};

			// =========== Original source contribution ===========
			// Sum contributions from 2 triangles of face j
			for(int t = 0; t < 2; t++)
			{
				double H_tri[3];

				if(usePrecomputed)
				{
					// Use pre-computed triangle data (fast path)
					int tri_idx = face_j * 2 + t;
					FieldFromTrianglePrecomputed(hex_j, tri_idx, obs, 1.0, H_tri);
				}
				else
				{
					// Fallback: compute on the fly
					int tvIdx = ((hex_j * 6 + face_j) * 2 + t) * 3 * 3;
					const double* V0 = &m_hexaTriVertices[tvIdx + 0];
					const double* V1 = &m_hexaTriVertices[tvIdx + 3];
					const double* V2 = &m_hexaTriVertices[tvIdx + 6];
					double sign = m_hexaTriSigns[(hex_j * 6 + face_j) * 2 + t];
					FieldFromChargedTriangleLocal(obs, V0, V1, V2, sign, H_tri);
				}

				H_total[0] += H_tri[0];
				H_total[1] += H_tri[1];
				H_total[2] += H_tri[2];
			}

			// Point charge contribution: m = -sigma * area
			// H_point = -area * (r - p) / |r - p|^3
			double area_j = m_hexaFaceAreas[hex_j * 6 + face_j];
			{
				double r[3] = {obs[0] - src_center_orig[0],
				               obs[1] - src_center_orig[1],
				               obs[2] - src_center_orig[2]};
				double dist_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

				if(dist_sq > 1e-30)
				{
					double dist = sqrt(dist_sq);
					double inv_dist3 = 1.0 / (dist * dist_sq);
					double coef = -area_j * inv_dist3;
					H_total[0] += coef * r[0];
					H_total[1] += coef * r[1];
					H_total[2] += coef * r[2];
				}
			}

			// =========== IMA: Mirrored source contributions ===========
			// Kernel-based IMA: compute field from mirrored source coordinates directly
			// This avoids virtual elements/DOFs and permutation arrays
			if(m_imaEnabled)
			{
				// Helper lambda: add field contribution from mirrored source
				auto addMirroredSourceContribution = [&](int mirrorAxis, int sign) {
					// Mirror sign: +1 for symmetric BC, -1 for antisymmetric BC
					double imaSign = (double)sign;

					// Mirrored source center
					double src_mirror[3] = {src_center_orig[0], src_center_orig[1], src_center_orig[2]};
					if(mirrorAxis & IMA_X) src_mirror[0] = -src_mirror[0];
					if(mirrorAxis & IMA_Y) src_mirror[1] = -src_mirror[1];
					if(mirrorAxis & IMA_Z) src_mirror[2] = -src_mirror[2];

					// Count number of mirror axes for winding correction
					int numMirrors = 0;
					if(mirrorAxis & IMA_X) numMirrors++;
					if(mirrorAxis & IMA_Y) numMirrors++;
					if(mirrorAxis & IMA_Z) numMirrors++;
					bool flipWinding = (numMirrors % 2 == 1);

					// Field from mirrored triangles
					for(int t = 0; t < 2; t++)
					{
						int tvIdx = ((hex_j * 6 + face_j) * 2 + t) * 3 * 3;
						double V0[3] = {m_hexaTriVertices[tvIdx + 0],
						                m_hexaTriVertices[tvIdx + 1],
						                m_hexaTriVertices[tvIdx + 2]};
						double V1[3] = {m_hexaTriVertices[tvIdx + 3],
						                m_hexaTriVertices[tvIdx + 4],
						                m_hexaTriVertices[tvIdx + 5]};
						double V2[3] = {m_hexaTriVertices[tvIdx + 6],
						                m_hexaTriVertices[tvIdx + 7],
						                m_hexaTriVertices[tvIdx + 8]};

						// Mirror vertices
						if(mirrorAxis & IMA_X) { V0[0] = -V0[0]; V1[0] = -V1[0]; V2[0] = -V2[0]; }
						if(mirrorAxis & IMA_Y) { V0[1] = -V0[1]; V1[1] = -V1[1]; V2[1] = -V2[1]; }
						if(mirrorAxis & IMA_Z) { V0[2] = -V0[2]; V1[2] = -V1[2]; V2[2] = -V2[2]; }

						// Swap V1/V2 to restore winding if odd number of mirrors
						if(flipWinding) {
							std::swap(V1[0], V2[0]);
							std::swap(V1[1], V2[1]);
							std::swap(V1[2], V2[2]);
						}

						double sign_tri = m_hexaTriSigns[(hex_j * 6 + face_j) * 2 + t];
						double H_tri[3];
						FieldFromChargedTriangleLocal(obs, V0, V1, V2, sign_tri, H_tri);

						H_total[0] += imaSign * H_tri[0];
						H_total[1] += imaSign * H_tri[1];
						H_total[2] += imaSign * H_tri[2];
					}

					// Point charge from mirrored center
					double r[3] = {obs[0] - src_mirror[0],
					               obs[1] - src_mirror[1],
					               obs[2] - src_mirror[2]};
					double dist_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

					if(dist_sq > 1e-30)
					{
						double dist = sqrt(dist_sq);
						double inv_dist3 = 1.0 / (dist * dist_sq);
						double coef = -area_j * inv_dist3 * imaSign;
						H_total[0] += coef * r[0];
						H_total[1] += coef * r[1];
						H_total[2] += coef * r[2];
					}
				};

				// Add contributions based on active symmetry axes
				bool hasX = (m_imaSymmetry & IMA_X) != 0;
				bool hasY = (m_imaSymmetry & IMA_Y) != 0;
				bool hasZ = (m_imaSymmetry & IMA_Z) != 0;

				// Single axis mirrors
				if(hasX) addMirroredSourceContribution(IMA_X, m_imaSignX);
				if(hasY) addMirroredSourceContribution(IMA_Y, m_imaSignY);
				if(hasZ) addMirroredSourceContribution(IMA_Z, m_imaSignZ);

				// Dual axis mirrors (for quarter models)
				if(hasX && hasY) addMirroredSourceContribution(IMA_XY, m_imaSignX * m_imaSignY);
				if(hasX && hasZ) addMirroredSourceContribution(IMA_XZ, m_imaSignX * m_imaSignZ);
				if(hasY && hasZ) addMirroredSourceContribution(IMA_YZ, m_imaSignY * m_imaSignZ);

				// Triple axis mirror (for eighth models)
				if(hasX && hasY && hasZ) addMirroredSourceContribution(IMA_XYZ, m_imaSignX * m_imaSignY * m_imaSignZ);
			}

			// K_ij = n_i dot H_total / (4*pi)
			double K_ij = (n_i[0]*H_total[0] + n_i[1]*H_total[1] + n_i[2]*H_total[2])
			              * RadConst::INV_FOUR_PI;

			// Store in ROW-MAJOR format: K[i][j] at index i*6 + j
			// The solver will negate when building system matrix
			K_mat[face_i * 6 + face_j] = K_ij;
		}
	}
}

//-------------------------------------------------------------------------
// IMA (Image) Symmetry Implementation
// Reference: IMA approach - matrix construction with image summation
//-------------------------------------------------------------------------

//=========================================================================
// ComputeMixedBlockFast: Unified cross-DOF interaction block computation
// Handles ALL element type pairs: 3x3, 3x5, 3x6, 5x3, 5x5, 5x6, 6x3, 6x5, 6x6
//
// Row element (target): observation points from precomputed geometry
//   3DOF (tet): obs = element center, result = H components directly
//   5DOF (wedge): obs = Yano midpoint per face, result = n_i dot H
//   6DOF (hex): obs = Yano midpoint per face, result = n_i dot H
//
// Col element (source): field from precomputed triangles + point charge
//   3DOF (tet): for each unit M_beta, sigma = n_f dot e_beta
//   5DOF (wedge): sigma on source face (1-2 triangles per face)
//   6DOF (hex): sigma on source face (2 triangles per face)
//
// IMA: MSC source -> scalar sign, MMM source -> component sign matrix S[beta]
//=========================================================================

void radTInteraction::ComputeMixedBlockFast(
	int elem_row, int dof_row, int elem_col, int dof_col,
	double* block_out) const
{
	std::memset(block_out, 0, dof_row * dof_col * sizeof(double));

	// Helper: convert global element index to type-specific index
	// O(1) reverse lookup via m_globalToHexIdx / m_globalToWedgeIdx / m_globalToTetraIdx
	// (built in PrecomputeHexaGeometry / PrecomputeWedgeGeometry / PrecomputeTetraGeometry).
	auto globalToHexIdx = [&](int globalIdx) -> int {
		if(globalIdx < 0 || globalIdx >= (int)m_globalToHexIdx.size()) return -1;
		return m_globalToHexIdx[globalIdx];
	};
	auto globalToWedgeIdx = [&](int globalIdx) -> int {
		if(globalIdx < 0 || globalIdx >= (int)m_globalToWedgeIdx.size()) return -1;
		return m_globalToWedgeIdx[globalIdx];
	};
	auto globalToTetIdx = [&](int globalIdx) -> int {
		if(globalIdx < 0 || globalIdx >= (int)m_globalToTetraIdx.size()) return -1;
		return m_globalToTetraIdx[globalIdx];
	};

	// Get type-specific indices
	int hex_row = -1, wedge_row = -1, hex_col = -1, wedge_col = -1;
	int tet_row = -1, tet_col = -1;
	if(dof_row == 6) hex_row = globalToHexIdx(elem_row);
	else if(dof_row == 5) wedge_row = globalToWedgeIdx(elem_row);
	else if(dof_row == 3) tet_row = globalToTetIdx(elem_row);
	if(dof_col == 6) hex_col = globalToHexIdx(elem_col);
	else if(dof_col == 5) wedge_col = globalToWedgeIdx(elem_col);
	else if(dof_col == 3) tet_col = globalToTetIdx(elem_col);

	// Helper: compute H field at obs from MSC source face (scalar sigma = 1)
	// Returns H without 4pi factor. Includes IMA mirror contributions.
	auto fieldFromMSCSourceFace = [&](const double* obs, int src_elem_type, int src_elem_idx, int src_face,
	                                   double* H_out) {
		H_out[0] = H_out[1] = H_out[2] = 0.0;
		double src_center[3];
		int triOff, numTris;
		const double* triVerts;
		const double* triSigns;
		double faceArea;
		int maxTris;

		if(src_elem_type == 6 && m_hexaGeomReady) {
			src_center[0] = m_hexaCenters[src_elem_idx*3+0];
			src_center[1] = m_hexaCenters[src_elem_idx*3+1];
			src_center[2] = m_hexaCenters[src_elem_idx*3+2];
			triOff = (src_elem_idx * 6 + src_face) * 2;
			numTris = 2;
			triVerts = &m_hexaTriVertices[triOff * 3 * 3];
			triSigns = &m_hexaTriSigns[triOff];
			faceArea = m_hexaFaceAreas[src_elem_idx*6 + src_face];
			maxTris = 2;
		} else if(src_elem_type == 5 && m_wedgeGeomReady) {
			src_center[0] = m_wedgeCenters[src_elem_idx*3+0];
			src_center[1] = m_wedgeCenters[src_elem_idx*3+1];
			src_center[2] = m_wedgeCenters[src_elem_idx*3+2];
			int off = m_wedgeTriOffset[src_elem_idx*5 + src_face];
			numTris = m_wedgeFaceNumTris[src_elem_idx*5 + src_face];
			triVerts = &m_wedgeTriVertices[(src_elem_idx*WEDGE_MAX_TRIS + off) * 3 * 3];
			triSigns = &m_wedgeTriSigns[src_elem_idx*WEDGE_MAX_TRIS + off];
			faceArea = m_wedgeFaceAreas[src_elem_idx*5 + src_face];
			maxTris = numTris;
		} else return;

		// Direct contribution
		for(int t = 0; t < numTris; t++) {
			const double* V0 = &triVerts[t*9+0];
			const double* V1 = &triVerts[t*9+3];
			const double* V2 = &triVerts[t*9+6];
			double H_tri[3];
			FieldFromChargedTriangleLocal(obs, V0, V1, V2, triSigns[t], H_tri);
			H_out[0] += H_tri[0]; H_out[1] += H_tri[1]; H_out[2] += H_tri[2];
		}
		// Point charge
		double r[3] = {obs[0]-src_center[0], obs[1]-src_center[1], obs[2]-src_center[2]};
		double d2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
		if(d2 > 1e-30) { double id3=1.0/(sqrt(d2)*d2); double c=-faceArea*id3;
			H_out[0]+=c*r[0]; H_out[1]+=c*r[1]; H_out[2]+=c*r[2]; }

		// IMA mirrors
		if(m_imaEnabled) {
			auto addMir = [&](int mirrorAxis, int sign) {
				double imaSign = (double)sign;
				double mc[3] = {src_center[0], src_center[1], src_center[2]};
				if(mirrorAxis & IMA_X) mc[0]=-mc[0];
				if(mirrorAxis & IMA_Y) mc[1]=-mc[1];
				if(mirrorAxis & IMA_Z) mc[2]=-mc[2];
				int nm=0; if(mirrorAxis&IMA_X) nm++; if(mirrorAxis&IMA_Y) nm++; if(mirrorAxis&IMA_Z) nm++;
				bool fw = (nm%2==1);
				for(int t=0;t<numTris;t++) {
					double V0[3]={triVerts[t*9+0],triVerts[t*9+1],triVerts[t*9+2]};
					double V1[3]={triVerts[t*9+3],triVerts[t*9+4],triVerts[t*9+5]};
					double V2[3]={triVerts[t*9+6],triVerts[t*9+7],triVerts[t*9+8]};
					if(mirrorAxis&IMA_X){V0[0]=-V0[0];V1[0]=-V1[0];V2[0]=-V2[0];}
					if(mirrorAxis&IMA_Y){V0[1]=-V0[1];V1[1]=-V1[1];V2[1]=-V2[1];}
					if(mirrorAxis&IMA_Z){V0[2]=-V0[2];V1[2]=-V1[2];V2[2]=-V2[2];}
					if(fw) for(int k=0;k<3;k++) std::swap(V1[k],V2[k]);
					double H_tri[3]; FieldFromChargedTriangleLocal(obs, V0, V1, V2, triSigns[t], H_tri);
					H_out[0]+=imaSign*H_tri[0]; H_out[1]+=imaSign*H_tri[1]; H_out[2]+=imaSign*H_tri[2];
				}
				double r[3]={obs[0]-mc[0],obs[1]-mc[1],obs[2]-mc[2]};
				double d2=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
				if(d2>1e-30){double id3=1.0/(sqrt(d2)*d2);double c=-faceArea*id3*imaSign;
					H_out[0]+=c*r[0];H_out[1]+=c*r[1];H_out[2]+=c*r[2];}
			};
			bool hX=(m_imaSymmetry&IMA_X)!=0, hY=(m_imaSymmetry&IMA_Y)!=0, hZ=(m_imaSymmetry&IMA_Z)!=0;
			if(hX) addMir(IMA_X,m_imaSignX); if(hY) addMir(IMA_Y,m_imaSignY); if(hZ) addMir(IMA_Z,m_imaSignZ);
			if(hX&&hY) addMir(IMA_XY,m_imaSignX*m_imaSignY);
			if(hX&&hZ) addMir(IMA_XZ,m_imaSignX*m_imaSignZ);
			if(hY&&hZ) addMir(IMA_YZ,m_imaSignY*m_imaSignZ);
			if(hX&&hY&&hZ) addMir(IMA_XYZ,m_imaSignX*m_imaSignY*m_imaSignZ);
		}
	};

	// --- Case 1: MSC row x MSC col (5x5, 5x6, 6x5, 6x6) ---
	if(dof_row >= 5 && dof_col >= 5)
	{
		for(int fi = 0; fi < dof_row; fi++) {
			double obs[3], n_i[3];
			if(dof_row == 6 && m_hexaGeomReady && hex_row >= 0) {
				int ep=(hex_row*6+fi)*3; obs[0]=m_hexaEvalPoints[ep]; obs[1]=m_hexaEvalPoints[ep+1]; obs[2]=m_hexaEvalPoints[ep+2];
				int fn=(hex_row*6+fi)*3; n_i[0]=m_hexaFaceNormals[fn]; n_i[1]=m_hexaFaceNormals[fn+1]; n_i[2]=m_hexaFaceNormals[fn+2];
			} else if(dof_row == 5 && m_wedgeGeomReady && wedge_row >= 0) {
				int ep=(wedge_row*5+fi)*3; obs[0]=m_wedgeEvalPoints[ep]; obs[1]=m_wedgeEvalPoints[ep+1]; obs[2]=m_wedgeEvalPoints[ep+2];
				int fn=(wedge_row*5+fi)*3; n_i[0]=m_wedgeFaceNormals[fn]; n_i[1]=m_wedgeFaceNormals[fn+1]; n_i[2]=m_wedgeFaceNormals[fn+2];
			} else continue;

			int col_idx = (dof_col == 6) ? hex_col : wedge_col;
			for(int fj = 0; fj < dof_col; fj++) {
				double H[3];
				fieldFromMSCSourceFace(obs, dof_col, col_idx, fj, H);
				block_out[fi * dof_col + fj] = (n_i[0]*H[0]+n_i[1]*H[1]+n_i[2]*H[2]) * RadConst::INV_FOUR_PI;
			}
		}
		return;
	}

	// --- Case 2: MMM row (3DOF) x MSC col (5/6DOF) ---
	if(dof_row == 3 && dof_col >= 5)
	{
		if(!m_tetraGeomReady || tet_row < 0) return;
		const double* obs = &m_tetraCenters[tet_row * 3];
		int col_idx = (dof_col == 6) ? hex_col : wedge_col;
		for(int fj = 0; fj < dof_col; fj++) {
			double H[3];
			fieldFromMSCSourceFace(obs, dof_col, col_idx, fj, H);
			// N[alpha][fj] = H_alpha (directly, no n_i projection)
			for(int a = 0; a < 3; a++)
				block_out[a * dof_col + fj] = H[a] * RadConst::INV_FOUR_PI;
		}
		return;
	}

	// --- Case 3: MSC row (5/6DOF) x MMM col (3DOF) ---
	if(dof_row >= 5 && dof_col == 3)
	{
		// For each target face i, for each unit M_beta:
		// Compute sigma=n_f dot e_beta on each source tet face, field at obs, dot n_i
		if(!m_tetraGeomReady || tet_col < 0) return;
		const double* col_center = &m_tetraCenters[tet_col * 3];

		for(int fi = 0; fi < dof_row; fi++) {
			double obs[3], n_i[3];
			if(dof_row == 6 && m_hexaGeomReady && hex_row >= 0) {
				int ep=(hex_row*6+fi)*3; obs[0]=m_hexaEvalPoints[ep]; obs[1]=m_hexaEvalPoints[ep+1]; obs[2]=m_hexaEvalPoints[ep+2];
				int fn=(hex_row*6+fi)*3; n_i[0]=m_hexaFaceNormals[fn]; n_i[1]=m_hexaFaceNormals[fn+1]; n_i[2]=m_hexaFaceNormals[fn+2];
			} else if(dof_row == 5 && m_wedgeGeomReady && wedge_row >= 0) {
				int ep=(wedge_row*5+fi)*3; obs[0]=m_wedgeEvalPoints[ep]; obs[1]=m_wedgeEvalPoints[ep+1]; obs[2]=m_wedgeEvalPoints[ep+2];
				int fn=(wedge_row*5+fi)*3; n_i[0]=m_wedgeFaceNormals[fn]; n_i[1]=m_wedgeFaceNormals[fn+1]; n_i[2]=m_wedgeFaceNormals[fn+2];
			} else continue;

			// Compute H at obs for each unit M_beta from tet source
			for(int beta = 0; beta < 3; beta++) {
				double H_total[3] = {0,0,0};
				double total_charge = 0;
				for(int f = 0; f < 4; f++) {
					int fvIdx = (tet_col*4+f)*3*3;
					const double* V0=&m_tetraFaceVertices[fvIdx]; const double* V1=&m_tetraFaceVertices[fvIdx+3]; const double* V2=&m_tetraFaceVertices[fvIdx+6];
					double sigma = m_tetraFaceNormals[(tet_col*4+f)*3+beta];
					total_charge += sigma * m_tetraFaceAreas[tet_col*4+f];
					if(fabs(sigma) > 1e-20) {
						double H_f[3]; FieldFromChargedTriangleLocal(obs, V0, V1, V2, sigma, H_f);
						H_total[0]+=H_f[0]; H_total[1]+=H_f[1]; H_total[2]+=H_f[2];
					}
				}
				// Point charge
				double r[3]={obs[0]-col_center[0],obs[1]-col_center[1],obs[2]-col_center[2]};
				double d2=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
				if(d2>1e-30){double id3=1.0/(sqrt(d2)*d2);
					H_total[0]+=-total_charge*r[0]*id3; H_total[1]+=-total_charge*r[1]*id3; H_total[2]+=-total_charge*r[2]*id3;}

				// IMA for MMM source: sign matrix S[beta]
				if(m_imaEnabled) {
					auto addMirTet = [&](int mirrorAxis, int combinedSign) {
						double S_beta = (double)combinedSign;
						if((mirrorAxis & IMA_X) && beta==0) S_beta = -S_beta;
						if((mirrorAxis & IMA_Y) && beta==1) S_beta = -S_beta;
						if((mirrorAxis & IMA_Z) && beta==2) S_beta = -S_beta;

						double mc[3]={col_center[0],col_center[1],col_center[2]};
						if(mirrorAxis&IMA_X)mc[0]=-mc[0]; if(mirrorAxis&IMA_Y)mc[1]=-mc[1]; if(mirrorAxis&IMA_Z)mc[2]=-mc[2];
						int nm=0; if(mirrorAxis&IMA_X)nm++; if(mirrorAxis&IMA_Y)nm++; if(mirrorAxis&IMA_Z)nm++;
						bool fw=(nm%2==1);
						double mirH[3]={0,0,0}; double mirCharge=0;
						for(int f=0;f<4;f++){
							int fvIdx=(tet_col*4+f)*3*3;
							double V0[3],V1[3],V2[3];
							for(int k=0;k<3;k++){V0[k]=m_tetraFaceVertices[fvIdx+k];V1[k]=m_tetraFaceVertices[fvIdx+3+k];V2[k]=m_tetraFaceVertices[fvIdx+6+k];}
							if(mirrorAxis&IMA_X){V0[0]=-V0[0];V1[0]=-V1[0];V2[0]=-V2[0];}
							if(mirrorAxis&IMA_Y){V0[1]=-V0[1];V1[1]=-V1[1];V2[1]=-V2[1];}
							if(mirrorAxis&IMA_Z){V0[2]=-V0[2];V1[2]=-V1[2];V2[2]=-V2[2];}
							if(fw) for(int k=0;k<3;k++) std::swap(V1[k],V2[k]);
							double e1[3]={V1[0]-V0[0],V1[1]-V0[1],V1[2]-V0[2]};
							double e2[3]={V2[0]-V0[0],V2[1]-V0[1],V2[2]-V0[2]};
							double n_f[3]={e1[1]*e2[2]-e1[2]*e2[1],e1[2]*e2[0]-e1[0]*e2[2],e1[0]*e2[1]-e1[1]*e2[0]};
							double nL=sqrt(n_f[0]*n_f[0]+n_f[1]*n_f[1]+n_f[2]*n_f[2]);
							double area=0.5*nL; if(nL>1e-20){n_f[0]/=nL;n_f[1]/=nL;n_f[2]/=nL;}
							double sigma=n_f[beta]; mirCharge+=sigma*area;
							if(fabs(sigma)>1e-20){double H_f[3];FieldFromChargedTriangleLocal(obs,V0,V1,V2,sigma,H_f);
								mirH[0]+=H_f[0];mirH[1]+=H_f[1];mirH[2]+=H_f[2];}
						}
						double r[3]={obs[0]-mc[0],obs[1]-mc[1],obs[2]-mc[2]};
						double d2=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
						if(d2>1e-30){double id3=1.0/(sqrt(d2)*d2);
							mirH[0]+=-mirCharge*r[0]*id3;mirH[1]+=-mirCharge*r[1]*id3;mirH[2]+=-mirCharge*r[2]*id3;}
						for(int a=0;a<3;a++) H_total[a]+=S_beta*mirH[a];
					};
					bool hX=(m_imaSymmetry&IMA_X)!=0,hY=(m_imaSymmetry&IMA_Y)!=0,hZ=(m_imaSymmetry&IMA_Z)!=0;
					if(hX) addMirTet(IMA_X,m_imaSignX); if(hY) addMirTet(IMA_Y,m_imaSignY); if(hZ) addMirTet(IMA_Z,m_imaSignZ);
					if(hX&&hY) addMirTet(IMA_XY,m_imaSignX*m_imaSignY);
					if(hX&&hZ) addMirTet(IMA_XZ,m_imaSignX*m_imaSignZ);
					if(hY&&hZ) addMirTet(IMA_YZ,m_imaSignY*m_imaSignZ);
					if(hX&&hY&&hZ) addMirTet(IMA_XYZ,m_imaSignX*m_imaSignY*m_imaSignZ);
				}

				block_out[fi * 3 + beta] = (n_i[0]*H_total[0]+n_i[1]*H_total[1]+n_i[2]*H_total[2]) * RadConst::INV_FOUR_PI;
			}
		}
		return;
	}

	// --- Case 4: MMM row x MMM col (3x3) - delegate to Compute3x3BlockFast ---
	if(dof_row == 3 && dof_col == 3) {
		Compute3x3BlockFast(elem_row, elem_col, block_out);
		return;
	}
}

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

// Legacy IMA functions removed (2026-03-31):
// - ApplyDOFPermutation, ApplyRowPermutation: ELF-style face permutation (replaced by kernel-based IMA)
// - Compute6x6BlockIMA: used ELF permutations, never called
// - Compute6x6BlockMirrored, Compute6x6BlockMirroredTarget: helper for Compute6x6BlockIMA

// REMOVED_BLOCK_START (was ApplyDOFPermutation through Compute6x6BlockMirroredTarget)
#if 0  // Dead code preserved for reference only
void radTInteraction::ApplyDOFPermutation(const double* input, const int* perm, double* result) const
{
	// result = input @ P
	// P[j][k] = 1 if perm[j] == k, else 0
	// (input @ P)[i][k] = sum_j input[i][j] * P[j][k] = input[i][perm^-1[k]]
	//
	// Alternative interpretation: column j of input goes to column perm[j] of result
	// result[:, perm[j]] = input[:, j]

	// Initialize result to zero
	for(int k = 0; k < 36; k++) result[k] = 0.0;

	// For each column j of input, copy to column perm[j] of result
	for(int j = 0; j < 6; j++)
	{
		int k = perm[j];  // column j maps to column k
		for(int i = 0; i < 6; i++)
		{
			// input[i][j] -> result[i][k]
			result[i * 6 + k] = input[i * 6 + j];
		}
	}
}

//-------------------------------------------------------------------------
// ApplyRowPermutation: Apply row permutation Q to 6x6 block
// result[i, :] = input[perm[i], :]
// For ELF IMA: result = Q @ K_BA
// Q[i] = row index in K_BA that becomes row i in result
//-------------------------------------------------------------------------
void radTInteraction::ApplyRowPermutation(const double* input, const int* perm, double* result) const
{
	// result = Q @ input
	// Q[i][k] = 1 if perm[i] == k, else 0
	// (Q @ input)[i][j] = sum_k Q[i][k] * input[k][j] = input[perm[i]][j]

	// For each row i of result, copy from row perm[i] of input
	for(int i = 0; i < 6; i++)
	{
		int k = perm[i];  // row i comes from row k
		for(int j = 0; j < 6; j++)
		{
			// input[k][j] -> result[i][j]
			result[i * 6 + j] = input[k * 6 + j];
		}
	}
}

//-------------------------------------------------------------------------
// Compute6x6BlockIMA: Compute IMA block with image summation (ELF-compatible)
// ELF formula: K_IMA[i,j] = K[i,j] + sign * Q @ K[mirror(i), j]
// where K[mirror(i), j] = K_BA (field at mirrored target from original source)
// For single axis: K_IMA[i,j] = K[i,j] + sign * Q @ K[mirror(i), j]
// For dual axis (e.g., XZ):
//   K_IMA[i,j] = K[i,j] + sign_x * K[i, mx(j)] @ Px
//                       + sign_z * K[i, mz(j)] @ Pz
//                       + sign_x * sign_z * K[i, mxz(j)] @ Pxz
// For quarter models: if no physical mirror exists, compute virtual mirror
//-------------------------------------------------------------------------
void radTInteraction::Compute6x6BlockIMA(int ima_i, int ima_j, double* K_ima) const
{
	if(!m_hexaGeomReady)
	{
		std::cerr << "[Radia] Error: Hexahedron geometry not precomputed for IMA" << std::endl;
		for(int k = 0; k < 36; k++) K_ima[k] = 0.0;
		return;
	}

	int full_i = m_imaToFull[ima_i];
	int full_j = m_imaToFull[ima_j];

	// Find hex indices in precomputed arrays via O(1) reverse lookup
	int hex_i = -1, hex_j = -1;
	if(full_i >= 0 && full_i < (int)m_globalToHexIdx.size()) hex_i = m_globalToHexIdx[full_i];
	if(full_j >= 0 && full_j < (int)m_globalToHexIdx.size()) hex_j = m_globalToHexIdx[full_j];

	if(hex_i < 0 || hex_j < 0)
	{
		std::cerr << "[Radia] Error: Element not found in hex arrays" << std::endl;
		for(int k = 0; k < 36; k++) K_ima[k] = 0.0;
		return;
	}

	// Start with direct interaction: K[full_i][full_j] = K_AA
	Compute6x6BlockFast(hex_i, hex_j, K_ima);

	// ELF-compatible IMA formula:
	// K_IMA[i,i] = K[i,i] + sign * Q @ K[mirror(i), i]
	//            = K_AA + sign * Q @ K_BA
	//
	// where K_BA = K[mirror(i), i] is computed by Compute6x6BlockMirroredTarget
	// Q is the row permutation to reorder faces after mirroring
	//
	// For x-mirror: Q = [0, 3, 2, 1, 4, 5] (swap x+ (face 1) and x- (face 3))

	// Helper lambda to add mirror contribution using ELF formula
	auto addMirrorContributionELF = [&](int mirrorAxis, int sign, const int* rowPerm, const int* colPerm) {
		double K_BA[36];        // K[mirror(i), j] from Compute6x6BlockMirroredTarget
		double K_BA_perm[36];   // After row permutation: Q @ K_BA

		// Compute K_BA: mirrored TARGET element i
		Compute6x6BlockMirroredTarget(hex_i, hex_j, mirrorAxis, K_BA);

		// Apply row permutation Q: K_BA_perm[i, :] = K_BA[rowPerm[i], :]
		ApplyRowPermutation(K_BA, rowPerm, K_BA_perm);

		// Add contribution: K_ima += sign * Q @ K_BA
		for(int k = 0; k < 36; k++)
		{
			K_ima[k] += sign * K_BA_perm[k];
		}
	};

	// Add contributions based on active symmetry axes
	bool hasX = (m_imaSymmetry & IMA_X) != 0;
	bool hasY = (m_imaSymmetry & IMA_Y) != 0;
	bool hasZ = (m_imaSymmetry & IMA_Z) != 0;

	// Single axis contributions using ELF permutations
	if(hasX) addMirrorContributionELF(IMA_X, m_imaSignX, IMA_ROW_PERM_X, IMA_COL_PERM_X);
	if(hasY) addMirrorContributionELF(IMA_Y, m_imaSignY, IMA_ROW_PERM_Y, IMA_COL_PERM_Y);
	if(hasZ) addMirrorContributionELF(IMA_Z, m_imaSignZ, IMA_ROW_PERM_Z, IMA_COL_PERM_Z);

	// Dual axis contributions (combined permutations)
	if(hasX && hasY)
	{
		// XY mirror: apply column perms, then row perms
		double K_AB[36], K_col1[36], K_col2[36], K_row1[36], K_perm[36];
		Compute6x6BlockMirrored(hex_i, hex_j, IMA_XY, K_AB);
		ApplyDOFPermutation(K_AB, IMA_COL_PERM_X, K_col1);
		ApplyDOFPermutation(K_col1, IMA_COL_PERM_Y, K_col2);
		ApplyRowPermutation(K_col2, IMA_ROW_PERM_X, K_row1);
		ApplyRowPermutation(K_row1, IMA_ROW_PERM_Y, K_perm);
		int sign = m_imaSignX * m_imaSignY;
		for(int k = 0; k < 36; k++) K_ima[k] += sign * K_perm[k];
	}
	if(hasX && hasZ)
	{
		// XZ mirror: apply column perms, then row perms
		double K_AB[36], K_col1[36], K_col2[36], K_row1[36], K_perm[36];
		Compute6x6BlockMirrored(hex_i, hex_j, IMA_XZ, K_AB);
		ApplyDOFPermutation(K_AB, IMA_COL_PERM_X, K_col1);
		ApplyDOFPermutation(K_col1, IMA_COL_PERM_Z, K_col2);
		ApplyRowPermutation(K_col2, IMA_ROW_PERM_X, K_row1);
		ApplyRowPermutation(K_row1, IMA_ROW_PERM_Z, K_perm);
		int sign = m_imaSignX * m_imaSignZ;
		for(int k = 0; k < 36; k++) K_ima[k] += sign * K_perm[k];
	}
	if(hasY && hasZ)
	{
		// YZ mirror: apply column perms, then row perms
		double K_AB[36], K_col1[36], K_col2[36], K_row1[36], K_perm[36];
		Compute6x6BlockMirrored(hex_i, hex_j, IMA_YZ, K_AB);
		ApplyDOFPermutation(K_AB, IMA_COL_PERM_Y, K_col1);
		ApplyDOFPermutation(K_col1, IMA_COL_PERM_Z, K_col2);
		ApplyRowPermutation(K_col2, IMA_ROW_PERM_Y, K_row1);
		ApplyRowPermutation(K_row1, IMA_ROW_PERM_Z, K_perm);
		int sign = m_imaSignY * m_imaSignZ;
		for(int k = 0; k < 36; k++) K_ima[k] += sign * K_perm[k];
	}

	// Triple axis contribution (eighth model)
	if(hasX && hasY && hasZ)
	{
		double K_AB[36], K_col1[36], K_col2[36], K_col3[36];
		double K_row1[36], K_row2[36], K_perm[36];
		Compute6x6BlockMirrored(hex_i, hex_j, IMA_XYZ, K_AB);
		ApplyDOFPermutation(K_AB, IMA_COL_PERM_X, K_col1);
		ApplyDOFPermutation(K_col1, IMA_COL_PERM_Y, K_col2);
		ApplyDOFPermutation(K_col2, IMA_COL_PERM_Z, K_col3);
		ApplyRowPermutation(K_col3, IMA_ROW_PERM_X, K_row1);
		ApplyRowPermutation(K_row1, IMA_ROW_PERM_Y, K_row2);
		ApplyRowPermutation(K_row2, IMA_ROW_PERM_Z, K_perm);
		int sign = m_imaSignX * m_imaSignY * m_imaSignZ;
		for(int k = 0; k < 36; k++) K_ima[k] += sign * K_perm[k];
	}
}

//-------------------------------------------------------------------------
// Compute6x6BlockMirrored: Compute interaction with virtually mirrored element j
// For quarter model support: element j's geometry is mirrored on-the-fly
// mirrorAxis: IMA_X, IMA_Y, IMA_Z, or combinations
//-------------------------------------------------------------------------
void radTInteraction::Compute6x6BlockMirrored(int hex_i, int hex_j, int mirrorAxis, double* K_mat) const
{
	std::memset(K_mat, 0, 36 * sizeof(double));

	if(!m_hexaGeomReady) return;

	int nHex = (int)m_hexaElemIndices.size();
	if(hex_i < 0 || hex_i >= nHex || hex_j < 0 || hex_j >= nHex) return;

	// Mirror source element center
	double src_center[3] = {m_hexaCenters[hex_j * 3 + 0],
	                        m_hexaCenters[hex_j * 3 + 1],
	                        m_hexaCenters[hex_j * 3 + 2]};
	if(mirrorAxis & IMA_X) src_center[0] = -src_center[0];
	if(mirrorAxis & IMA_Y) src_center[1] = -src_center[1];
	if(mirrorAxis & IMA_Z) src_center[2] = -src_center[2];

	// For each target face i (unchanged)
	for(int face_i = 0; face_i < 6; face_i++)
	{
		// Yano evaluation point for target face
		int epIdx = (hex_i * 6 + face_i) * 3;
		const double obs[3] = {m_hexaEvalPoints[epIdx + 0],
		                       m_hexaEvalPoints[epIdx + 1],
		                       m_hexaEvalPoints[epIdx + 2]};

		// Target face normal
		int fnIdx_i = (hex_i * 6 + face_i) * 3;
		const double n_i[3] = {m_hexaFaceNormals[fnIdx_i + 0],
		                       m_hexaFaceNormals[fnIdx_i + 1],
		                       m_hexaFaceNormals[fnIdx_i + 2]};

		// For each source face j (mirrored)
		for(int face_j = 0; face_j < 6; face_j++)
		{
			// Field from unit sigma on mirrored source face j
			double H_total[3] = {0.0, 0.0, 0.0};

			// Sum contributions from 2 triangles of mirrored face j
			for(int t = 0; t < 2; t++)
			{
				// Get original triangle vertices
				int tvIdx = ((hex_j * 6 + face_j) * 2 + t) * 3 * 3;
				double V0[3] = {m_hexaTriVertices[tvIdx + 0],
				                m_hexaTriVertices[tvIdx + 1],
				                m_hexaTriVertices[tvIdx + 2]};
				double V1[3] = {m_hexaTriVertices[tvIdx + 3],
				                m_hexaTriVertices[tvIdx + 4],
				                m_hexaTriVertices[tvIdx + 5]};
				double V2[3] = {m_hexaTriVertices[tvIdx + 6],
				                m_hexaTriVertices[tvIdx + 7],
				                m_hexaTriVertices[tvIdx + 8]};

				// Mirror vertices
				if(mirrorAxis & IMA_X)
				{
					V0[0] = -V0[0]; V1[0] = -V1[0]; V2[0] = -V2[0];
				}
				if(mirrorAxis & IMA_Y)
				{
					V0[1] = -V0[1]; V1[1] = -V1[1]; V2[1] = -V2[1];
				}
				if(mirrorAxis & IMA_Z)
				{
					V0[2] = -V0[2]; V1[2] = -V1[2]; V2[2] = -V2[2];
				}

				// Original sign (for orientation)
				double sign = m_hexaTriSigns[(hex_j * 6 + face_j) * 2 + t];

				// Mirror reverses winding -> flip sign
				// Count number of mirror axes (odd = flip, even = no flip)
				int numMirrors = 0;
				if(mirrorAxis & IMA_X) numMirrors++;
				if(mirrorAxis & IMA_Y) numMirrors++;
				if(mirrorAxis & IMA_Z) numMirrors++;
				if(numMirrors % 2 == 1)
				{
					// Odd number of mirrors: swap V1 and V2 to restore winding
					std::swap(V1[0], V2[0]);
					std::swap(V1[1], V2[1]);
					std::swap(V1[2], V2[2]);
				}

				double H_tri[3];
				FieldFromChargedTriangleLocal(obs, V0, V1, V2, sign, H_tri);

				H_total[0] += H_tri[0];
				H_total[1] += H_tri[1];
				H_total[2] += H_tri[2];
			}

			// Point charge contribution from mirrored element center
			// (Same approximation as Compute6x6BlockFast: use element center for all faces)
			double area_j = m_hexaFaceAreas[hex_j * 6 + face_j];
			double r[3] = {obs[0] - src_center[0],
			               obs[1] - src_center[1],
			               obs[2] - src_center[2]};
			double dist_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

			if(dist_sq > 1e-30)
			{
				double dist = sqrt(dist_sq);
				double inv_dist3 = 1.0 / (dist * dist_sq);
				double coef = -area_j * inv_dist3;
				H_total[0] += coef * r[0];
				H_total[1] += coef * r[1];
				H_total[2] += coef * r[2];
			}

			// K_ij = n_i dot H_total / (4*pi)
			double K_ij = (n_i[0]*H_total[0] + n_i[1]*H_total[1] + n_i[2]*H_total[2])
			              * RadConst::INV_FOUR_PI;

			// Store in ROW-MAJOR format: K[i][j] at index i*6 + j
			// The solver will negate when building system matrix
			K_mat[face_i * 6 + face_j] = K_ij;
		}
	}
}

//-------------------------------------------------------------------------
// Compute6x6BlockMirroredTarget: Compute K[mirror(i), j] for ELF IMA formula
// This mirrors the TARGET element i, keeping SOURCE element j unchanged.
// K_BA[face_i, face_j] = n'_i dot H_j / (4*pi)
// where n'_i is the mirrored target face normal,
// and H_j is computed at the mirrored target evaluation point
// mirrorAxis: IMA_X, IMA_Y, IMA_Z, or combinations
//-------------------------------------------------------------------------
void radTInteraction::Compute6x6BlockMirroredTarget(int hex_i, int hex_j, int mirrorAxis, double* K_mat) const
{
	std::memset(K_mat, 0, 36 * sizeof(double));

	if(!m_hexaGeomReady) return;

	int nHex = (int)m_hexaElemIndices.size();
	if(hex_i < 0 || hex_i >= nHex || hex_j < 0 || hex_j >= nHex) return;

	// Source element center (unchanged)
	double src_center[3] = {m_hexaCenters[hex_j * 3 + 0],
	                        m_hexaCenters[hex_j * 3 + 1],
	                        m_hexaCenters[hex_j * 3 + 2]};

	// For each target face i (MIRRORED evaluation point and normal)
	for(int face_i = 0; face_i < 6; face_i++)
	{
		// Get original target evaluation point and mirror it
		int epIdx = (hex_i * 6 + face_i) * 3;
		double obs[3] = {m_hexaEvalPoints[epIdx + 0],
		                 m_hexaEvalPoints[epIdx + 1],
		                 m_hexaEvalPoints[epIdx + 2]};

		// Mirror the evaluation point
		if(mirrorAxis & IMA_X) obs[0] = -obs[0];
		if(mirrorAxis & IMA_Y) obs[1] = -obs[1];
		if(mirrorAxis & IMA_Z) obs[2] = -obs[2];

		// Get original target face normal and mirror it
		int fnIdx_i = (hex_i * 6 + face_i) * 3;
		double n_i[3] = {m_hexaFaceNormals[fnIdx_i + 0],
		                 m_hexaFaceNormals[fnIdx_i + 1],
		                 m_hexaFaceNormals[fnIdx_i + 2]};

		// Mirror the normal (flip component in mirror axis)
		if(mirrorAxis & IMA_X) n_i[0] = -n_i[0];
		if(mirrorAxis & IMA_Y) n_i[1] = -n_i[1];
		if(mirrorAxis & IMA_Z) n_i[2] = -n_i[2];

		// For each source face j (unchanged geometry)
		for(int face_j = 0; face_j < 6; face_j++)
		{
			// Field from unit sigma on ORIGINAL source face j
			double H_total[3] = {0.0, 0.0, 0.0};

			// Sum contributions from 2 triangles of original source face j
			for(int t = 0; t < 2; t++)
			{
				// Get original triangle vertices (NOT mirrored)
				int tvIdx = ((hex_j * 6 + face_j) * 2 + t) * 3 * 3;
				double V0[3] = {m_hexaTriVertices[tvIdx + 0],
				                m_hexaTriVertices[tvIdx + 1],
				                m_hexaTriVertices[tvIdx + 2]};
				double V1[3] = {m_hexaTriVertices[tvIdx + 3],
				                m_hexaTriVertices[tvIdx + 4],
				                m_hexaTriVertices[tvIdx + 5]};
				double V2[3] = {m_hexaTriVertices[tvIdx + 6],
				                m_hexaTriVertices[tvIdx + 7],
				                m_hexaTriVertices[tvIdx + 8]};

				double sign = m_hexaTriSigns[(hex_j * 6 + face_j) * 2 + t];

				double H_tri[3];
				FieldFromChargedTriangleLocal(obs, V0, V1, V2, sign, H_tri);

				H_total[0] += H_tri[0];
				H_total[1] += H_tri[1];
				H_total[2] += H_tri[2];
			}

			// Point charge contribution from original source element center
			double area_j = m_hexaFaceAreas[hex_j * 6 + face_j];
			double r[3] = {obs[0] - src_center[0],
			               obs[1] - src_center[1],
			               obs[2] - src_center[2]};
			double dist_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

			if(dist_sq > 1e-30)
			{
				double dist = sqrt(dist_sq);
				double inv_dist3 = 1.0 / (dist * dist_sq);
				double coef = -area_j * inv_dist3;
				H_total[0] += coef * r[0];
				H_total[1] += coef * r[1];
				H_total[2] += coef * r[2];
			}

			// K_ij = n'_i dot H_total / (4*pi)
			// where n'_i is the mirrored normal
			double K_ij = (n_i[0]*H_total[0] + n_i[1]*H_total[1] + n_i[2]*H_total[2])
			              * RadConst::INV_FOUR_PI;

			// Store in ROW-MAJOR format: K[i][j] at index i*6 + j
			// The solver will negate when building system matrix
			K_mat[face_i * 6 + face_j] = K_ij;
		}
	}
}
#endif  // Dead code end
// REMOVED_BLOCK_END

//-------------------------------------------------------------------------
// SetupInteractMatrix_IMA: Build IMA interaction matrix
// Both hex and tet fast paths now delegate to kernel functions that
// handle IMA mirror contributions internally (Compute6x6BlockFast, Compute3x3BlockFast).
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

	// For HACApK: skip dense matrix, let kernel functions compute on demand.
	// All element types (pure and mixed) are handled by kernel functions:
	// - Same-type: Compute6x6BlockFast, Compute5x5BlockFast, Compute3x3BlockFast
	// - Cross-DOF: ComputeMixedBlockFast (with type-specific indexing)
	// Reset geometry so it gets recomputed for the reduced IMA element set.
	if(skipDenseMatrix)
	{
		// Reset precomputed geometry so HACApK recomputes for the reduced IMA elements
		m_hexaGeomReady = false;
		m_wedgeGeomReady = false;
		m_tetraGeomReady = false;
		return 1;
	}

	// Build mapping from IMA index to type-specific geometry index
	// Uses geometry precomputed from the full model (before system reduction)
	// NOTE: Always populate for ALL element types (not just pure meshes) to support mixed meshes
	std::vector<int> imaToHex(m_imaNumElements, -1);
	if(m_hexaGeomReady)
	{
		for(int ima_i = 0; ima_i < m_imaNumElements; ima_i++)
		{
			if(m_elemDOF[ima_i] != 6) continue;
			int full_i = m_imaToFull[ima_i];
			// O(1) reverse lookup via m_globalToHexIdx
			if(full_i >= 0 && full_i < (int)m_globalToHexIdx.size())
				imaToHex[ima_i] = m_globalToHexIdx[full_i];
		}
	}

	std::vector<int> imaToWedge(m_imaNumElements, -1);
	if(m_wedgeGeomReady)
	{
		for(int ima_i = 0; ima_i < m_imaNumElements; ima_i++)
		{
			if(m_elemDOF[ima_i] != 5) continue;
			int full_i = m_imaToFull[ima_i];
			// O(1) reverse lookup via m_globalToWedgeIdx
			if(full_i >= 0 && full_i < (int)m_globalToWedgeIdx.size())
				imaToWedge[ima_i] = m_globalToWedgeIdx[full_i];
		}
	}

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

			// Fast path: 6DOF-6DOF hex with precomputed geometry
			if(dof_row == 6 && dof_col == 6)
			{
				int hex_row = imaToHex[ima_row];
				int hex_col = imaToHex[ima_col];
				if(hex_row >= 0 && hex_col >= 0)
				{
					double K_ima[36];
					Compute6x6BlockFast(hex_row, hex_col, K_ima);
					for(int i = 0; i < 6; i++)
						for(int j = 0; j < 6; j++)
							block[(size_t)i * imaDOF + j] = K_ima[i * 6 + j];
					continue;
				}
			}

			// Fast path: 5DOF-5DOF wedge - Compute5x5BlockFast handles IMA inline
			if(dof_row == 5 && dof_col == 5)
			{
				int w_row = imaToWedge[ima_row];
				int w_col = imaToWedge[(int)ima_col];
				if(w_row >= 0 && w_col >= 0)
				{
					double K_ima[25];
					Compute5x5BlockFast(w_row, w_col, K_ima);
					for(int i = 0; i < 5; i++)
						for(int j = 0; j < 5; j++)
							block[(size_t)i * imaDOF + j] = K_ima[i * 5 + j];
					continue;
				}
			}

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

			// Unified kernel path: handles ALL DOF combinations (same-type and cross-DOF)
			// including IMA mirror contributions, using precomputed geometry from full model.
			// ComputeMixedBlockFast takes full model element indices and resolves type-specific
			// geometry indices internally.
			{
				int full_row = m_imaToFull[ima_row];
				int full_col = m_imaToFull[(int)ima_col];
				double mixed_block[36];  // max 6x6
				ComputeMixedBlockFast(full_row, dof_row, full_col, dof_col, mixed_block);
				for(int i = 0; i < dof_row; i++)
					for(int j = 0; j < dof_col; j++)
						block[(size_t)i * imaDOF + j] = mixed_block[i * dof_col + j];
			}
		}
	});

	// IMA interaction matrix built successfully
	return 1;
}

//-------------------------------------------------------------------------
