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
#include "rad_subdivided_rectangle.h"
#include "rad_polyhedron.h"  // For IsTetrahedron() check in N_self fix

#ifdef _OPENMP
#include <omp.h>
#endif

// MSC (Magnetic Surface Charge) support for 6 DOF hexahedra
// radTPolyhedron hexahedra use 6 DOF MSC (surface charge on each face)
#define RADIA_MSC_SUPPORT
#ifdef RADIA_MSC_SUPPORT
// Note: radTPolyhedron is already included above (line 19)
// radTExtrPolygonMSC is deprecated - use radTPolyhedron with 6 faces instead
#endif

// Note: Dipole-dipole method for tetrahedra was tested but found numerically unstable.
// Radia production solver uses the surface charge (MSC) method.

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
		AllocateMemory(AuxOldMagnArrayIsNeeded); //In case of MPI-parallelization, this has to be executed by master only

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
		//--New
		radTSubdividedRecMag* SubdividedRecMagPtr = Cast.SubdividedRecMagCast(GroupPtr);
		if(SubdividedRecMagPtr != 0)
		{
			radTg3dRelax* g3dRelaxFromSbdRecMagPtr = (radTg3dRelax*)SubdividedRecMagPtr;

			radTRecMag* SubElRecMagPtr = Cast.RecMagCast((radTg3dRelax*)((*(SubdividedRecMagPtr->GroupMapOfHandlers.begin())).second.rep));

			if((g3dRelaxFromSbdRecMagPtr->MaterHandle.rep != 0) && (SubElRecMagPtr != 0))
			{
				int SubIntervStart = AmOfMainElem;
				if(SubdividedRecMagPtr->FldCmpMeth==1)
				{
					for(int ix=0; ix<int(SubdividedRecMagPtr->kx); ix++)
						for(int iy=0; iy<int(SubdividedRecMagPtr->ky); iy++)
							for(int iz=0; iz<int(SubdividedRecMagPtr->kz); iz++)
							{
								g3dRelaxPtrVect.push_back(g3dRelaxFromSbdRecMagPtr);
								AmOfMainElem++;

								radTlphgPtr* TotalListOfElemTransPtrPtr = new radTlphgPtr(*CurrListOfTransPtrPtr);
								PushFrontNativeElemTransList(g3dRelaxFromSbdRecMagPtr, TotalListOfElemTransPtrPtr);
								IntVectOfPtrToListsOfTransPtr.push_back(TotalListOfElemTransPtrPtr);
							}
				}
				int SubIntervFin = SubIntervStart + (int)(SubdividedRecMagPtr->GroupMapOfHandlers.size()) - 1;

				if(RelaxSubIntervConstrVect.empty())
				{
					radTRelaxSubInterval RlxSbIntrv(SubIntervStart, SubIntervFin, TRelaxSubIntervalID::RelaxTogether);
					RelaxSubIntervConstrVect.push_back(RlxSbIntrv);
				}
				else
				{
					radTRelaxSubInterval& LastEnteredSubIntrv = RelaxSubIntervConstrVect.back();
					if((SubIntervStart != LastEnteredSubIntrv.StartNo) && (SubIntervFin != LastEnteredSubIntrv.FinNo))
					{
						radTRelaxSubInterval RlxSbIntrv(SubIntervStart, SubIntervFin, TRelaxSubIntervalID::RelaxTogether);
						RelaxSubIntervConstrVect.push_back(RlxSbIntrv);
					}
				}
			}
		}
		if((SubdividedRecMagPtr == 0) || ((SubdividedRecMagPtr != 0) && (SubdividedRecMagPtr->FldCmpMeth != 1)))
		{
		//--EndNew
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
		//--New
		}
		//--EndNew
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

void radTInteraction::AllocateMemory(char AuxOldMagnArrayIsNeeded)
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

		#pragma omp parallel for if(AmOfMainElem > 20)
		for(int ColNo=0; ColNo<AmOfMainElem; ColNo++)
		{
			radTg3dRelax* g3dRelaxPtrColNo = g3dRelaxPtrVect[ColNo];

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

				InteractMatrix[StrNo][ColNo] = SubMatrix;
			}
		}

		//DEBUG
		//long long nTotMatrElem = ((long long)AmOfMainElem)*((long long)AmOfMainElem);
		//std::cout << "rank=" << m_rankMPI << ": iCntBcomp= " << iCntBcomp << "; nTotMatrElem=" << nTotMatrElem; //DEBUG
		//std::cout.flush();
		//END DEBUG

		// SELF-INTERACTION NOTE for tetrahedral elements:
		// For tetrahedral elements with correct coordinate transforms, B_comp() should
		// compute correct self-demagnetization (~-1/3). If there are issues, they may
		// be in the coordinate transformation chain, not the diagonal values.
		// The following code has been disabled pending further investigation.
		/*
		const float N_self = -1.0f / 3.0f;
		for(int diagNo = 0; diagNo < AmOfMainElem; diagNo++)
		{
			radTg3dRelax* g3dRelaxPtr = g3dRelaxPtrVect[diagNo];
			radTPolyhedron* polyPtr = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
			if(polyPtr != nullptr && polyPtr->IsTetrahedron())
			{
				TMatrix3df& diag = InteractMatrix[diagNo][diagNo];
				diag.Str0.x = N_self;
				diag.Str1.y = N_self;
				diag.Str2.z = N_self;
			}
		}
		*/

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
// DEPRECATED: Dipole-Dipole Interaction Matrix
//
// This method was tested but found numerically unstable.
// Kept for historical reference. Results varied wildly with mesh size (117K-549K A/m).
// Radia production code uses surface charge (MSC) method.
//
// The dipole-dipole approximation:
//   Diagonal (self-demagnetization): N_ii = 1/3 * I  (isotropic, sphere approx)
//   Off-diagonal: N_ij = (V_j / 4*pi) * (3*r*r^T/r^5 - I/r^3)
//
// Problems:
// 1. N_self=1/3 is only exact for spheres, not tetrahedra
// 2. Far-field approximation breaks down for adjacent elements
//=========================================================================
//-------------------------------------------------------------------------
#if 0  // DISABLED - numerically unstable, kept for reference
int radTInteraction::SetupInteractMatrix_DipoleDipole()
{
	const double PI = 3.14159265358979323846;
	const double ONE_OVER_4PI = 1.0 / (4.0 * PI);
	const double ONE_THIRD = 1.0 / 3.0;

	TVector3d ZeroVect(0., 0., 0.);

	if(m_nProcMPI < 2)
	{
		for(int ColNo = 0; ColNo < AmOfMainElem; ColNo++)
		{
			radTg3dRelax* elem_col = g3dRelaxPtrVect[ColNo];
			TVector3d center_col = MainTransPtrArray[ColNo]->TrPoint(elem_col->ReturnCentrPoint());
			double vol_col = elem_col->Volume();

			for(int StrNo = 0; StrNo < AmOfMainElem; StrNo++)
			{
				radTg3dRelax* elem_row = g3dRelaxPtrVect[StrNo];
				TVector3d center_row = MainTransPtrArray[StrNo]->TrPoint(elem_row->ReturnCentrPoint());

				TMatrix3d SubMatrix(ZeroVect, ZeroVect, ZeroVect);

				if(ColNo == StrNo)
				{
					// Diagonal: self-demagnetization N_self = 1/3 * I
					// This is the exact value for a uniformly magnetized sphere
					// and a good approximation for compact elements
					SubMatrix.Str0.x = -ONE_THIRD;
					SubMatrix.Str1.y = -ONE_THIRD;
					SubMatrix.Str2.z = -ONE_THIRD;
				}
				else
				{
					// Off-diagonal: dipole-dipole interaction
					// r = center_row - center_col (displacement from source to target)
					TVector3d r;
					r.x = center_row.x - center_col.x;
					r.y = center_row.y - center_col.y;
					r.z = center_row.z - center_col.z;

					double dist2 = r.x*r.x + r.y*r.y + r.z*r.z;
					double dist = sqrt(dist2);
					double dist3 = dist2 * dist;
					double dist5 = dist3 * dist2;

					// Coefficient: -vol_col / (4*pi)
					// NEGATIVE sign because Radia stores demagnetization tensor (H_demag = N*M)
					// and the dipole field H = +N_dipole*M enhances magnetization,
					// while demagnetization opposes it.
					double coef = -vol_col * ONE_OVER_4PI;

					// Demagnetization tensor: N_ij = -coef * (3*r*r^T/r^5 - I/r^3)
					// (negative of dipole field tensor)

					double coef_r3 = coef / dist3;
					double coef_r5_3 = 3.0 * coef / dist5;

					// Row 0: dH_demag/dMx
					SubMatrix.Str0.x = coef_r5_3 * r.x * r.x - coef_r3;
					SubMatrix.Str0.y = coef_r5_3 * r.x * r.y;
					SubMatrix.Str0.z = coef_r5_3 * r.x * r.z;

					// Row 1: dH_demag/dMy
					SubMatrix.Str1.x = coef_r5_3 * r.y * r.x;
					SubMatrix.Str1.y = coef_r5_3 * r.y * r.y - coef_r3;
					SubMatrix.Str1.z = coef_r5_3 * r.y * r.z;

					// Row 2: dH_demag/dMz
					SubMatrix.Str2.x = coef_r5_3 * r.z * r.x;
					SubMatrix.Str2.y = coef_r5_3 * r.z * r.y;
					SubMatrix.Str2.z = coef_r5_3 * r.z * r.z - coef_r3;
				}

				// Store in interaction matrix
				// Note: MainTransPtrArray[StrNo]->TrMatrix_inv would transform the matrix
				// but for dipole-dipole, we work directly in global coordinates
				InteractMatrix[StrNo][ColNo] = SubMatrix;
			}
		}

		// Update formal interaction member pointers (same as SetupInteractMatrix)
		for(int ClNo = 0; ClNo < AmOfMainElem; ClNo++)
		{
			radTg3dRelax* g3dRelaxPtrClNo = g3dRelaxPtrVect[ClNo];
			g3dRelaxPtrVect[ClNo] = g3dRelaxPtrClNo->FormalIntrctMemberPtr();
		}
	}

	return 1;
}
#endif  // DISABLED dipole-dipole method

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
	// Matrix is COLUMN-MAJOR: A(i,j) at index [j * m_totalDOF + i]
	if(m_flatInteractMatrix.empty()) return nullptr;

	int offset_row = m_elemDOFOffset[row_elem];
	int offset_col = m_elemDOFOffset[col_elem];

	// Block starts at row offset_row, column offset_col in the totalDOF x totalDOF matrix
	// Column-major: element at (row, col) is at index [col * m_totalDOF + row]
	return &m_flatInteractMatrix[offset_col * m_totalDOF + offset_row];
}

const double* radTInteraction::GetInteractBlock(int row_elem, int col_elem) const
{
	if(m_flatInteractMatrix.empty()) return nullptr;

	int offset_row = m_elemDOFOffset[row_elem];
	int offset_col = m_elemDOFOffset[col_elem];

	// Column-major: element at (row, col) is at index [col * m_totalDOF + row]
	return &m_flatInteractMatrix[offset_col * m_totalDOF + offset_row];
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

	// Allocate flattened interaction matrix
	m_flatInteractMatrix.resize(m_totalDOF * m_totalDOF, 0.0);

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

	// Build interaction matrix with variable-size blocks
	// For each pair (row_elem, col_elem), compute the interaction block
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
			// COLUMN-MAJOR: A(row, col) at index [col * m_totalDOF + row]
			double* block = &m_flatInteractMatrix[offset_col * m_totalDOF + offset_row];

			// Compute the interaction block based on DOF types
			// For now, only support 3x3 blocks (standard elements)
			// TODO: Add support for 3x6, 6x3, 6x6 blocks for MSC elements

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

				// Copy 3x3 matrix to flattened block (COLUMN-MAJOR)
				// A(i,j) at [j * stride + i] where stride = m_totalDOF
				block[0 * m_totalDOF + 0] = SubMatrix.Str0.x;  // (0,0)
				block[0 * m_totalDOF + 1] = SubMatrix.Str1.x;  // (1,0)
				block[0 * m_totalDOF + 2] = SubMatrix.Str2.x;  // (2,0)
				block[1 * m_totalDOF + 0] = SubMatrix.Str0.y;  // (0,1)
				block[1 * m_totalDOF + 1] = SubMatrix.Str1.y;  // (1,1)
				block[1 * m_totalDOF + 2] = SubMatrix.Str2.y;  // (2,1)
				block[2 * m_totalDOF + 0] = SubMatrix.Str0.z;  // (0,2)
				block[2 * m_totalDOF + 1] = SubMatrix.Str1.z;  // (1,2)
				block[2 * m_totalDOF + 2] = SubMatrix.Str2.z;  // (2,2)
			}
#ifdef RADIA_MSC_SUPPORT
			else if(dof_row == 3 && dof_col == 6)
			{
				// 3x6 block: Field at standard element (3 DOF) center from MSC hexahedron (6 DOF)
				// For each sigma_j on face j of MSC element, compute field at standard element center
				// N[row][col]_ij = component i of field at row's center due to unit sigma on col's face j

				radTPolyhedron* poly_col = dynamic_cast<radTPolyhedron*>(elem_col);
				if(poly_col && poly_col->Use6DOF_MSC)
				{
					TVector3d InitObsPoiVect = MainTransPtrArray[row]->TrPoint(elem_row->ReturnCentrPoint());

					for(int face_j = 0; face_j < 6; face_j++)
					{
						TVector3d H_total(0., 0., 0.);

						for(unsigned tr = 0; tr < TransPtrVect.size(); tr++)
						{
							TVector3d ObsPoiVect = TransPtrVect[tr]->TrPoint_inv(InitObsPoiVect);

							// Field from unit sigma on face j (quad face + point charge)
							TVector3d H_face = poly_col->FieldFromQuadFace(ObsPoiVect, face_j, 1.0);

							// Point charge contribution (m = -sigma * area)
							double unit_point_charge = -1.0 * poly_col->FaceArea[face_j];
							TVector3d H_point = poly_col->FieldFromPointCharge(ObsPoiVect, unit_point_charge);

							TVector3d H_local;
							H_local.x = H_face.x + H_point.x;
							H_local.y = H_face.y + H_point.y;
							H_local.z = H_face.z + H_point.z;

							// Transform back
							H_total.x += TransPtrVect[tr]->TrVectField(H_local).x;
							H_total.y += TransPtrVect[tr]->TrVectField(H_local).y;
							H_total.z += TransPtrVect[tr]->TrVectField(H_local).z;
						}

						// Transform by row's main transform (inverse)
						TVector3d H_final = MainTransPtrArray[row]->TrVectField_inv(H_total);

						// Store in block (COLUMN-MAJOR): A(i,j) at [j * stride + i]
						// Here: row is component (0,1,2=x,y,z), col is face_j
						block[face_j * m_totalDOF + 0] = H_final.x;  // (0, face_j)
						block[face_j * m_totalDOF + 1] = H_final.y;  // (1, face_j)
						block[face_j * m_totalDOF + 2] = H_final.z;  // (2, face_j)
					}
				}
			}
			else if(dof_row == 6 && dof_col == 3)
			{
				// 6x3 block: Field at MSC hexahedron (6 DOF) eval points from standard element (3 DOF)
				// For each eval point i on MSC element, compute field component from standard element
				// N[row][col]_ij = H dot n_i at eval point i due to unit M_j

				radTPolyhedron* poly_row = dynamic_cast<radTPolyhedron*>(elem_row);
				if(poly_row && poly_row->Use6DOF_MSC)
				{
					for(int face_i = 0; face_i < 6; face_i++)
					{
						// Eval point for face i (midpoint between face center and element center)
						TVector3d EvalPt;
						EvalPt.x = 0.5 * (poly_row->FaceCenter[face_i].x + poly_row->CentrPoint.x);
						EvalPt.y = 0.5 * (poly_row->FaceCenter[face_i].y + poly_row->CentrPoint.y);
						EvalPt.z = 0.5 * (poly_row->FaceCenter[face_i].z + poly_row->CentrPoint.z);

						TVector3d InitObsPoiVect = MainTransPtrArray[row]->TrPoint(EvalPt);

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

						// Get face normal
						TVector3d& n = poly_row->FaceNormal[face_i];

						// H dot n for each magnetization component
						// N[face_i][Mx] = (dHx/dMx*nx + dHy/dMx*ny + dHz/dMx*nz)
						double H_dot_n_Mx = SubMatrix.Str1.x * n.x + SubMatrix.Str1.y * n.y + SubMatrix.Str1.z * n.z;
						// Rows of SubMatrix correspond to: 0=Hx, 1=Hy, 2=Hz for unit magnetization
						// But actually SubMatrix stores dB/dM, not dH/dM directly
						// For now, use H = B/mu0 - M approximation for linear materials

						// Store in block (COLUMN-MAJOR): A(i,j) at [j * stride + i]
						// row is face_i, col is component (0,1,2=Mx,My,Mz)
						block[0 * m_totalDOF + face_i] = SubMatrix.Str1.x * n.x;  // (face_i, 0)
						block[1 * m_totalDOF + face_i] = SubMatrix.Str1.y * n.y;  // (face_i, 1)
						block[2 * m_totalDOF + face_i] = SubMatrix.Str1.z * n.z;  // (face_i, 2)
					}
				}
			}
			else if(dof_row == 6 && dof_col == 6)
			{
				// 6x6 block: Field at MSC hexahedron (6 DOF) eval points from MSC hexahedron (6 DOF)
				// K(face_i, face_j) = normal_i dot H_field(eval_pt_i, src_face_j)
				// Field functions return values WITHOUT 4pi divisor (ELF_MAGIC convention)
				// Matrix stores: -K_ij / (4*pi)

				static const double PI_MSC = 3.14159265358979323846;
				static const double INV_4PI_MSC = 1.0 / (4.0 * PI_MSC);

				radTPolyhedron* poly_row = dynamic_cast<radTPolyhedron*>(elem_row);
				radTPolyhedron* poly_col = dynamic_cast<radTPolyhedron*>(elem_col);

				if(poly_row && poly_row->Use6DOF_MSC && poly_col && poly_col->Use6DOF_MSC)
				{
					for(int face_i = 0; face_i < 6; face_i++)
					{
						// Yano-Sugahara evaluation point: midpoint between face center and element center
						TVector3d EvalPt;
						EvalPt.x = 0.5 * (poly_row->FaceCenter[face_i].x + poly_row->CentrPoint.x);
						EvalPt.y = 0.5 * (poly_row->FaceCenter[face_i].y + poly_row->CentrPoint.y);
						EvalPt.z = 0.5 * (poly_row->FaceCenter[face_i].z + poly_row->CentrPoint.z);

						TVector3d InitObsPoiVect = MainTransPtrArray[row]->TrPoint(EvalPt);

						for(int face_j = 0; face_j < 6; face_j++)
						{
							double K_ij = 0.0;

							for(unsigned tr = 0; tr < TransPtrVect.size(); tr++)
							{
								TVector3d ObsPoiVect = TransPtrVect[tr]->TrPoint_inv(InitObsPoiVect);

								// Field from unit sigma on face j (quad face + point charge)
								// Field functions return values WITHOUT 4pi divisor
								TVector3d H_face = poly_col->FieldFromQuadFace(ObsPoiVect, face_j, 1.0);

								// Point charge contribution: m = -sigma * area (Yano-Sugahara MSC method)
								double unit_point_charge = -1.0 * poly_col->FaceArea[face_j];
								TVector3d H_point = poly_col->FieldFromPointCharge(ObsPoiVect, unit_point_charge);

								TVector3d H_local;
								H_local.x = H_face.x + H_point.x;
								H_local.y = H_face.y + H_point.y;
								H_local.z = H_face.z + H_point.z;

								// Transform back
								TVector3d H_global = TransPtrVect[tr]->TrVectField(H_local);
								K_ij += H_global.x * poly_row->FaceNormal[face_i].x +
								        H_global.y * poly_row->FaceNormal[face_i].y +
								        H_global.z * poly_row->FaceNormal[face_i].z;
							}

							// Store -K_ij / (4*pi) (COLUMN-MAJOR): A(i,j) at [j * stride + i]
							block[face_j * m_totalDOF + face_i] = -K_ij * INV_4PI_MSC;
						}
					}
				}
			}
#endif // RADIA_MSC_SUPPORT
			else
			{
				// Unknown DOF combination - zero out the block (COLUMN-MAJOR)
				for(int i = 0; i < dof_row; i++)
				{
					for(int j = 0; j < dof_col; j++)
					{
						block[j * m_totalDOF + i] = 0.0;
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
				BufVect += TransPtrVect[i]->TrVectField(Field.H);
			}
			ExternFieldArray[StrNo] += MainTransPtrArray[StrNo]->TrVectField_inv(BufVect);
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
			else if(dof == 6)
			{
				// MSC hexahedron: compute H_ext dot n at each face
				radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(elem);
				if(poly && poly->Use6DOF_MSC)
				{
					for(int face_i = 0; face_i < 6; face_i++)
					{
						// Eval point for face i (midpoint between face center and element center)
						TVector3d EvalPt;
						EvalPt.x = 0.5 * (poly->FaceCenter[face_i].x + poly->CentrPoint.x);
						EvalPt.y = 0.5 * (poly->FaceCenter[face_i].y + poly->CentrPoint.y);
						EvalPt.z = 0.5 * (poly->FaceCenter[face_i].z + poly->CentrPoint.z);

						// Compute H_ext at eval point from all external sources
						TVector3d H_total(0., 0., 0.);
						for(int ExtElNo = 0; ExtElNo < AmOfExtElem; ExtElNo++)
						{
							radTg3d* ExtElPtr = g3dExternPtrVect[ExtElNo];
							radTField Field(FieldKeyExtern, CompCriterium, EvalPt, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
							ExtElPtr->B_comp(&Field);
							H_total.x += Field.H.x;
							H_total.y += Field.H.y;
							H_total.z += Field.H.z;
						}

						// H_ext dot n_i
						double H_dot_n = H_total.x * poly->FaceNormal[face_i].x +
						                 H_total.y * poly->FaceNormal[face_i].y +
						                 H_total.z * poly->FaceNormal[face_i].z;

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

			ExternFieldArray[StrNo] += MainTransPtrArray[StrNo]->TrVectField_inv(Field.H);
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
							TVector3d EvalPt;
							EvalPt.x = 0.5 * (poly->FaceCenter[face_i].x + poly->CentrPoint.x);
							EvalPt.y = 0.5 * (poly->FaceCenter[face_i].y + poly->CentrPoint.y);
							EvalPt.z = 0.5 * (poly->FaceCenter[face_i].z + poly->CentrPoint.z);

							// Compute H_ext at eval point
							radTField Field(FieldKeyExtern, CompCriterium, EvalPt, ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);
							(static_cast<radTg3d*>(MoreExtSourceHandle.rep))->B_genComp(&Field);

							// H_ext dot n_i
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

		ExternFieldArray[StrNo] += MainTransPtrArray[StrNo]->TrVectField_inv(Field.H);
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

//-------------------------------------------------------------------------

void radTInteraction::DumpBinVectOfPtrToListsOfTransPtr(CAuxBinStrVect& oStr, radVectPtr_lphgPtr& VectOfPtrToListsOfTransPtr, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers)
{
	int sizeVectOfPtrToListsOfTransPtr = (int)VectOfPtrToListsOfTransPtr.size();
	oStr << sizeVectOfPtrToListsOfTransPtr;
	for(int i=0; i<sizeVectOfPtrToListsOfTransPtr; i++)
	{
		radTlphgPtr* curListOfElemTransPtrPtr = VectOfPtrToListsOfTransPtr[i];
		int size_curListOfElemTransPtr = 0;
		if(curListOfElemTransPtrPtr != 0) size_curListOfElemTransPtr = (int)curListOfElemTransPtrPtr->size();
		
		oStr << size_curListOfElemTransPtr;
		if(size_curListOfElemTransPtr > 0)
		{
			for(radTlphgPtr::iterator TrIter = curListOfElemTransPtrPtr->begin();	TrIter != curListOfElemTransPtrPtr->end(); ++TrIter)
			{
				radTPair_int_hg *p_m_hg = *TrIter;
				//int mult = 0;

				if(p_m_hg != 0) 
				{
					int mult = p_m_hg->m;
					radThg &hg = p_m_hg->Handler_g;

					int existKey = 0;
					for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
					{
						if(mit->second == hg) { existKey = mit->first; break;}
					}
					oStr << mult;
					oStr << existKey;
				}
				else oStr << (int)0;
			}
		}
	}
}

//-------------------------------------------------------------------------

void radTInteraction::DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut, map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers, int& gUniqueMapKey, int elemKey)
{
	//radThg SourceHandle;
	int existKeySource = 0;
	if(SourceHandle.rep != 0)
	{
		//oStr << (char)1;
		//int existKey = 0;
		//const radThg &cur_hg = iter->second;
		for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
		{
			if(mit->second == SourceHandle) { existKeySource = mit->first; break;}
		}
		if(existKeySource == 0)
		{
			existKeySource = gUniqueMapKey; 
			gMapOfHandlers[gUniqueMapKey++] = SourceHandle;
		}
		int indExist = CAuxParse::FindElemInd(existKeySource, vElemKeysOut);
		if(indExist < 0) SourceHandle.rep->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, existKeySource);
	}
	//else oStr << (char)0;

	//radThg MoreExtSourceHandle;
	int existKeyMoreExtSource = 0;
	if(MoreExtSourceHandle.rep != 0)
	{
		//oStr << (char)1;
		//int existKey = 0;
		//const radThg &cur_hg = iter->second;
		for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
		{
			if(mit->second == MoreExtSourceHandle) { existKeyMoreExtSource = mit->first; break;}
		}
		if(existKeyMoreExtSource == 0)
		{
			existKeyMoreExtSource = gUniqueMapKey; 
			gMapOfHandlers[gUniqueMapKey++] = MoreExtSourceHandle;
		}
		int indExist = CAuxParse::FindElemInd(existKeyMoreExtSource, vElemKeysOut);
		if(indExist < 0) MoreExtSourceHandle.rep->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, existKeyMoreExtSource);
	}
	//else oStr << (char)0;

	//radTVectPtrg3dRelax g3dRelaxPtrVect;
	vector<int> vInd_g3dRelax;
	int size_g3dRelaxPtrVect = (int)g3dRelaxPtrVect.size();
	//oStr << size_g3dRelaxPtrVect;
	for(int i=0; i<size_g3dRelaxPtrVect; i++)
	{
		radTg3dRelax *p_g3dRelax = g3dRelaxPtrVect[i];
		if(p_g3dRelax != 0)
		{
			radTg *p_g = (radTg*)p_g3dRelax;
			//try to find element in the global map by pointer
			int oldKey = 0;
			for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
			{
				if(mit->second.rep == p_g) { oldKey = mit->first; break;}
			}
			if(oldKey == 0)
			{
				oldKey = gUniqueMapKey;
				radThg hg(p_g3dRelax);
				gMapOfHandlers[gUniqueMapKey++] = hg;
			}
			int indExist = CAuxParse::FindElemInd(oldKey, vElemKeysOut);
			if(indExist < 0) p_g3dRelax->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, oldKey);

			vInd_g3dRelax.push_back(oldKey);
		}
	}

	//radTVectPtr_g3d g3dExternPtrVect;
	vector<int> vInd_g3dExternPtrVect;
	int size_g3dExternPtrVect = (int)g3dExternPtrVect.size();
	for(int i=0; i<size_g3dExternPtrVect; i++)
	{
		radTg3d *p_g3d = g3dExternPtrVect[i];
		if(p_g3d != 0)
		{
			radTg *p_g = (radTg*)p_g3d;

			//try to find this element in the global map by pointer
			int oldKey = 0;
			for(radTmhg::iterator mit = gMapOfHandlers.begin(); mit != gMapOfHandlers.end(); ++mit)
			{
				if(mit->second.rep == p_g) { oldKey = mit->first; break;}
			}
			if(oldKey == 0)
			{
				oldKey = gUniqueMapKey;
				radThg hg(p_g3d);
				gMapOfHandlers[gUniqueMapKey++] = hg;
			}
			int indExist = CAuxParse::FindElemInd(oldKey, vElemKeysOut);
			if(indExist < 0) p_g3d->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, oldKey);

			vInd_g3dExternPtrVect.push_back(oldKey);
		}
	}

	//radTVectPtrTrans TransPtrVect; //not required?
	vector<int> vIndTransPtrVect;
	int size_TransPtrVect = (int)TransPtrVect.size();
	for(int i=0; i<size_TransPtrVect; i++)
	{
		radTrans *pTrans = TransPtrVect[i];
		if(pTrans != 0)
		{
			if(Cast.IdentTransCast(pTrans))
			{
				vIndTransPtrVect.push_back(-1); //indicator of IdentTrans
			}
			else
			{
				radTrans *pTransCopy = new radTrans(*pTrans);

				radThg hg(pTransCopy);
				int oldKey = gUniqueMapKey;
				gMapOfHandlers[gUniqueMapKey++] = hg;
				
				pTransCopy->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, oldKey);
				vIndTransPtrVect.push_back(oldKey);
			}
		}
		else vIndTransPtrVect.push_back(0);
	}

	//radTrans** MainTransPtrArray; //required
	vector<int> vIndMainTrans;
	if(mKeepTransData && (MainTransPtrArray != 0))
	{
		for(int i=0; i<AmOfMainElem; i++)
		{
			radTrans *pTrans = MainTransPtrArray[i];
			if(pTrans != 0)
			{
				if(Cast.IdentTransCast(pTrans))
				{
					vIndTransPtrVect.push_back(-1); //indicator of IdentTrans
				}
				else
				{
					radTrans *pTransCopy = new radTrans(*pTrans);

					radThg hg(pTransCopy);
					int oldKey = gUniqueMapKey;
					gMapOfHandlers[gUniqueMapKey++] = hg;
				
					pTransCopy->DumpBin(oStr, vElemKeysOut, gMapOfHandlers, gUniqueMapKey, oldKey);
					vIndMainTrans.push_back(oldKey);
				}
			}
			else vIndMainTrans.push_back(0);
		}
	}

	vElemKeysOut.push_back(elemKey);
	oStr << elemKey;

	//Next 5 bytes define/encode element type:
	oStr << (char)Type_g();
	oStr << (char)0;
	oStr << (char)0;
	oStr << (char)0;
	oStr << (char)0;

	//int AmOfMainElem;
	oStr << AmOfMainElem;

	//int AmOfExtElem;
	oStr << AmOfExtElem;

	//radThg SourceHandle;
	oStr << existKeySource;

	//radThg MoreExtSourceHandle;
	oStr << existKeyMoreExtSource;

	//radTVectPtrg3dRelax g3dRelaxPtrVect;
	int size_vInd_g3dRelax = (int)vInd_g3dRelax.size();
	oStr << size_vInd_g3dRelax;
	for(int i=0; i<size_vInd_g3dRelax; i++) oStr << vInd_g3dRelax[i];

	//radTVectPtr_g3d g3dExternPtrVect;
	int size_vInd_g3dExternPtrVect = (int)vInd_g3dExternPtrVect.size();
	oStr << size_vInd_g3dExternPtrVect;
	for(int i=0; i<size_vInd_g3dExternPtrVect; i++) oStr << vInd_g3dExternPtrVect[i];

	//radTVectPtrTrans TransPtrVect; //not required?
	int size_vIndTransPtrVect = (int)vIndTransPtrVect.size();
	oStr << size_vIndTransPtrVect;
	for(int i=0; i<size_vIndTransPtrVect; i++) oStr << vIndTransPtrVect[i];

	//radTCompCriterium CompCriterium;
	//short BasedOnPrecLevel; // Actually this is used nowhere at the moment
	oStr << CompCriterium.BasedOnPrecLevel;
	//double AbsPrecB;
	oStr << CompCriterium.AbsPrecB;
	//double AbsPrecA;
	oStr << CompCriterium.AbsPrecA;
	//double AbsPrecB_int;
	oStr << CompCriterium.AbsPrecB_int;
	//double AbsPrecForce;
	oStr << CompCriterium.AbsPrecForce;
	//double AbsPrecTorque;
	oStr << CompCriterium.AbsPrecTorque;
	//double AbsPrecEnergy;
	oStr << CompCriterium.AbsPrecTorque;
	//double AbsPrecTrjCoord;
	oStr << CompCriterium.AbsPrecTrjCoord;
	//double AbsPrecTrjAngle;
	oStr << CompCriterium.AbsPrecTrjAngle;
	//double MltplThresh[4]; // Threshold ratios for 4 diff. orders of multipole approx. at field computation
	oStr << CompCriterium.MltplThresh[0] << CompCriterium.MltplThresh[1] << CompCriterium.MltplThresh[2] << CompCriterium.MltplThresh[3];
	//double WorstRelPrec;
	oStr << CompCriterium.WorstRelPrec;
	//char BasedOnWorstRelPrec; // Used at energy - force computation
	oStr << CompCriterium.BasedOnWorstRelPrec;

	//radTRelaxStatusParam RelaxStatusParam;
	//double MisfitM, MaxModM, MaxModH;
	oStr << RelaxStatusParam.MisfitM;
	oStr << RelaxStatusParam.MaxModM;
	oStr << RelaxStatusParam.MaxModH;

	//short RelaxationStarted;
	oStr << RelaxationStarted;

	//TMatrix3df** InteractMatrix; //OC250504
	////TMatrix3d** InteractMatrix; //OC250504
	if(InteractMatrix != nullptr)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++)
		{
			TMatrix3df *pLineInteractMatrix = InteractMatrix[i];
			if(pLineInteractMatrix != nullptr)
			{
				oStr << (char)1;
				for(int j=0; j<AmOfMainElem; j++)
				{
					oStr << pLineInteractMatrix[j];
				}
			}
			else oStr << (char)0;
		}
	}
	else oStr << (char)0;

	//TVector3d* ExternFieldArray;
	if(ExternFieldArray != nullptr)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++) oStr << ExternFieldArray[i];
	}
	else oStr << (char)0;

	//TVector3d* NewMagnArray;
	if(NewMagnArray != nullptr)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++) oStr << NewMagnArray[i];
	}
	else oStr << (char)0;

	//TVector3d* NewFieldArray;
	if(NewFieldArray != nullptr)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++) oStr << NewFieldArray[i];
	}
	else oStr << (char)0;

	//TVector3d* AuxOldMagnArray;
	if(AuxOldMagnArray != nullptr)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++) oStr << AuxOldMagnArray[i];
	}
	else oStr << (char)0;

	//TVector3d* AuxOldFieldArray;
	if(AuxOldFieldArray != nullptr)
	{
		oStr << (char)1;
		for(int i=0; i<AmOfMainElem; i++) oStr << AuxOldFieldArray[i];
	}
	else oStr << (char)0;

	//radTVectRelaxSubInterval RelaxSubIntervConstrVect; // New
	int sizeRelaxSubIntervConstrVect = (int)RelaxSubIntervConstrVect.size();
	oStr << sizeRelaxSubIntervConstrVect;	
	if(sizeRelaxSubIntervConstrVect > 0)
	{
		for(int i=0; i<sizeRelaxSubIntervConstrVect; i++)
		{
			radTRelaxSubInterval &relaxSubInterval = RelaxSubIntervConstrVect[i];
			oStr << relaxSubInterval.StartNo;
			oStr << relaxSubInterval.FinNo;
			oStr << (int)(relaxSubInterval.SubIntervalID);
		}

		//radTRelaxSubInterval* RelaxSubIntervArray; // New 
		if(RelaxSubIntervArray != nullptr)
		{
			int MaxSubIntervArraySize = 2*sizeRelaxSubIntervConstrVect + 1;
			oStr << (int)MaxSubIntervArraySize;
			radTRelaxSubInterval *t_RelaxSubIntervArray = RelaxSubIntervArray;
			for(int i=0; i<MaxSubIntervArraySize; i++)
			{
				oStr << (t_RelaxSubIntervArray->StartNo);
				oStr << (t_RelaxSubIntervArray->FinNo);
				oStr << (int)(t_RelaxSubIntervArray->SubIntervalID);
				t_RelaxSubIntervArray++;
			}
		}
		else oStr << (int)0;
	}

	//radVectPtr_lphgPtr IntVectOfPtrToListsOfTransPtr; //required
	DumpBinVectOfPtrToListsOfTransPtr(oStr, IntVectOfPtrToListsOfTransPtr, gMapOfHandlers);

	//radVectPtr_lphgPtr ExtVectOfPtrToListsOfTransPtr; //required
	DumpBinVectOfPtrToListsOfTransPtr(oStr, ExtVectOfPtrToListsOfTransPtr, gMapOfHandlers);

	//radIdentTrans* IdentTransPtr; //required, but doesn't need to be saved
	//radTCast Cast; //no members?
	//radTSend Send; //no members?

	//short FillInMainTransOnly;
	oStr << FillInMainTransOnly;

	//char mKeepTransData;
	oStr << mKeepTransData;

	//radTrans** MainTransPtrArray; //required
	int size_vIndMainTrans = (int)vIndMainTrans.size();
	oStr << size_vIndMainTrans;
	for(int k=0; k<size_vIndMainTrans; k++) oStr << vIndMainTrans[k];
	
	//int AmOfRelaxSubInterv;
	oStr << AmOfRelaxSubInterv;

	//short SomethingIsWrong;
	oStr << SomethingIsWrong;

	//short MemAllocTotAtOnce;
	oStr << MemAllocTotAtOnce;
}

//-------------------------------------------------------------------------

//void radTInteraction::DumpBinParseSourceHandle(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers, bool do_g3dCast, bool do_g3dRelaxCast, radThg& out_hg)
int radTInteraction::DumpBinParseSourceHandle(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers, bool do_g3dCast, bool do_g3dRelaxCast, radThg& out_hg)
{//move to g3d?
	int oldKey = 0;
	inStr >> oldKey;
	if(oldKey > 0)
	{
		map<int, int>::const_iterator itKey = mKeysOldNew.find(oldKey);
		if(itKey != mKeysOldNew.end())
		{
			int newKey = itKey->second;
			if(newKey > 0)
			{
				radTmhg::const_iterator iter = gMapOfHandlers.find(newKey);
				if(iter != gMapOfHandlers.end())
				{
					radThg hg = (*iter).second;
					if(hg.rep != 0)
					{
						if(do_g3dCast || do_g3dRelaxCast)
						{
							radTg3d *g3dPtr = radTCast::g3dCast(hg.rep);
							if(g3dPtr != 0)
							{
								if(do_g3dRelaxCast)
								{
									if(radTCast::g3dRelaxCast(g3dPtr) != 0) out_hg = hg;
								}
								else out_hg = hg;
							}
						}
						else out_hg = hg;
					}
				}
			}
		}
	}
	return oldKey;
}

//-------------------------------------------------------------------------

void radTInteraction::DumpBinParseVectOfPtrToListsOfTransPtr(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers, radVectPtr_lphgPtr& VectOfPtrToListsOfTransPtr)
{
	int sizeVectOfPtrToListsOfTransPtr = 0;
	inStr >> sizeVectOfPtrToListsOfTransPtr;

	for(int i=0; i<sizeVectOfPtrToListsOfTransPtr; i++)
	{
		int size_curListOfElemTransPtr = 0;
		inStr >> size_curListOfElemTransPtr;

		if(size_curListOfElemTransPtr > 0)
		{
			radTlphgPtr *pCurListOfElemTransPtr = new radTlphgPtr();
			for(int j=0; j<size_curListOfElemTransPtr; j++)
			{
				int mult = 0;
				inStr >> mult;
				if(mult > 0)
				{
					radThg hg;
					DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, false, false, hg);
					pCurListOfElemTransPtr->push_back(new radTPair_int_hg(mult, hg));
				}
			}
			VectOfPtrToListsOfTransPtr.push_back(pCurListOfElemTransPtr);
		}
	}
}

//-------------------------------------------------------------------------

radTInteraction::radTInteraction(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers)
{
	//radIdentTrans* IdentTransPtr; //required
	IdentTransPtr = new radIdentTrans();

	//int AmOfMainElem;
	inStr >> AmOfMainElem;

	//int AmOfExtElem;
	inStr >> AmOfExtElem;

	//radThg SourceHandle;
	DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, true, false, SourceHandle);

	//radThg MoreExtSourceHandle;
	DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, true, false, MoreExtSourceHandle);

	//radTVectPtrg3dRelax g3dRelaxPtrVect;
	int size_g3dRelaxPtrVect = 0;
	inStr >> size_g3dRelaxPtrVect;
	if(g3dRelaxPtrVect.size() > 0) g3dRelaxPtrVect.erase(g3dRelaxPtrVect.begin(), g3dRelaxPtrVect.end()); //?
	for(int i=0; i<size_g3dRelaxPtrVect; i++)
	{
		radThg hg;
		DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, true, true, hg);
		if(hg.rep != 0) g3dRelaxPtrVect.push_back((radTg3dRelax*)((radTg3d*)hg.rep));
	}

	//radTVectPtr_g3d g3dExternPtrVect;
	int size_g3dExternPtrVect = 0;
	inStr >> size_g3dExternPtrVect;
	if(g3dExternPtrVect.size() > 0) g3dExternPtrVect.erase(g3dExternPtrVect.begin(), g3dExternPtrVect.end()); //?
	for(int i=0; i<size_g3dExternPtrVect; i++)
	{
		radThg hg;
		DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, true, false, hg);
		if(hg.rep != 0) g3dExternPtrVect.push_back((radTg3d*)hg.rep);
	}

	//radTVectPtrTrans TransPtrVect; //not required?
	int sizeTransPtrVect = 0;
	inStr >> sizeTransPtrVect;
	if(TransPtrVect.size() > 0) TransPtrVect.erase(TransPtrVect.begin(), TransPtrVect.end()); //?
	for(int i=0; i<sizeTransPtrVect; i++)
	{
		radThg hg;
		int oldKey = DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, false, false, hg);
		if(oldKey < 0) TransPtrVect.push_back(IdentTransPtr);
		else if(hg.rep != 0) TransPtrVect.push_back(new radTrans(*((radTrans*)hg.rep))); //will be deleted at distraction
	}

	//radTCompCriterium CompCriterium;
	//short BasedOnPrecLevel; // Actually this is used nowhere at the moment
	inStr >> CompCriterium.BasedOnPrecLevel;
	//double AbsPrecB;
	inStr >> CompCriterium.AbsPrecB;
	//double AbsPrecA;
	inStr >> CompCriterium.AbsPrecA;
	//double AbsPrecB_int;
	inStr >> CompCriterium.AbsPrecB_int;
	//double AbsPrecForce;
	inStr >> CompCriterium.AbsPrecForce;
	//double AbsPrecTorque;
	inStr >> CompCriterium.AbsPrecTorque;
	//double AbsPrecEnergy;
	inStr >> CompCriterium.AbsPrecTorque;
	//double AbsPrecTrjCoord;
	inStr >> CompCriterium.AbsPrecTrjCoord;
	//double AbsPrecTrjAngle;
	inStr >> CompCriterium.AbsPrecTrjAngle;
	//double MltplThresh[4]; // Threshold ratios for 4 diff. orders of multipole approx. at field computation
	inStr >> CompCriterium.MltplThresh[0];
	inStr >> CompCriterium.MltplThresh[1];
	inStr >> CompCriterium.MltplThresh[2];
	inStr >> CompCriterium.MltplThresh[3];
	//double WorstRelPrec;
	inStr >> CompCriterium.WorstRelPrec;
	//char BasedOnWorstRelPrec; // Used at energy - force computation
	inStr >> CompCriterium.BasedOnWorstRelPrec;

	//radTRelaxStatusParam RelaxStatusParam;
	//double MisfitM, MaxModM, MaxModH;
	inStr >> RelaxStatusParam.MisfitM;
	inStr >> RelaxStatusParam.MaxModM;
	inStr >> RelaxStatusParam.MaxModH;

	//short RelaxationStarted;
	inStr >> RelaxationStarted;

	//TMatrix3df** InteractMatrix;
	char matrixExists = 0;
	inStr >> matrixExists;
	if(matrixExists && (AmOfMainElem > 0))
	{
		vInteractMatrixPtrs.resize(AmOfMainElem, nullptr);
		vInteractMatrix.resize(AmOfMainElem);
		InteractMatrix = vInteractMatrixPtrs.data();

		for(int i=0; i<AmOfMainElem; i++)
		{
			char matrixRowExists = 0;
			inStr >> matrixRowExists;
			if(matrixRowExists)
			{
				vInteractMatrix[i].resize(AmOfMainElem);
				InteractMatrix[i] = vInteractMatrix[i].data();
				vInteractMatrixPtrs[i] = InteractMatrix[i];
				TMatrix3df *tLine = InteractMatrix[i];
				for(int j=0; j<AmOfMainElem; j++)
				{
					inStr >> *(tLine++);
				}
			}
		}
	}

	//TVector3d* ExternFieldArray;
	char externFieldArrayExists = 0;
	ExternFieldArray = 0;
	inStr >> externFieldArrayExists;
	if(externFieldArrayExists && (AmOfMainElem > 0))
	{
		vExternFieldArray.resize(AmOfMainElem);
		ExternFieldArray = vExternFieldArray.data();
		for(int i=0; i<AmOfMainElem; i++) inStr >> ExternFieldArray[i];
	}

	//TVector3d* NewMagnArray;
	char newMagnArrayExists = 0;
	NewMagnArray = 0;
	inStr >> newMagnArrayExists;
	if(newMagnArrayExists && (AmOfMainElem > 0))
	{
		vNewMagnArray.resize(AmOfMainElem);
		NewMagnArray = vNewMagnArray.data();
		for(int i=0; i<AmOfMainElem; i++) inStr >> NewMagnArray[i];
	}

	//TVector3d* NewFieldArray;
	char newFieldArrayExists = 0;
	NewFieldArray = 0;
	inStr >> newFieldArrayExists;
	if(newFieldArrayExists && (AmOfMainElem > 0))
	{
		vNewFieldArray.resize(AmOfMainElem);
		NewFieldArray = vNewFieldArray.data();
		for(int i=0; i<AmOfMainElem; i++) inStr >> NewFieldArray[i];
	}

	//TVector3d* AuxOldMagnArray;
	char auxOldMagnArrayExists = 0;
	AuxOldMagnArray = 0;
	inStr >> auxOldMagnArrayExists;
	if(auxOldMagnArrayExists && (AmOfMainElem > 0))
	{
		vAuxOldMagnArray.resize(AmOfMainElem);
		AuxOldMagnArray = vAuxOldMagnArray.data();
		for(int i=0; i<AmOfMainElem; i++) inStr >> AuxOldMagnArray[i];
	}

	//TVector3d* AuxOldFieldArray;
	char auxOldFieldArrayExists = 0;
	AuxOldFieldArray = 0;
	inStr >> auxOldFieldArrayExists;
	if(auxOldFieldArrayExists && (AmOfMainElem > 0))
	{
		vAuxOldFieldArray.resize(AmOfMainElem);
		AuxOldFieldArray = vAuxOldFieldArray.data();
		for(int i=0; i<AmOfMainElem; i++) inStr >> AuxOldFieldArray[i];
	}

	//radTVectRelaxSubInterval RelaxSubIntervConstrVect; // New
	int sizeRelaxSubIntervConstrVect = 0;
	RelaxSubIntervArray = 0;
	inStr >> sizeRelaxSubIntervConstrVect;
	if(sizeRelaxSubIntervConstrVect > 0)
	{
		for(int i=0; i<sizeRelaxSubIntervConstrVect; i++)
		{
			radTRelaxSubInterval relaxSubInterval;
			inStr >> relaxSubInterval.StartNo;
			inStr >> relaxSubInterval.FinNo;
			int subIntervalID = 0;
			inStr >> subIntervalID;
			relaxSubInterval.SubIntervalID = (TRelaxSubIntervalID)subIntervalID;

			RelaxSubIntervConstrVect.push_back(relaxSubInterval);
		}

		//radTRelaxSubInterval* RelaxSubIntervArray; // New 
		int MaxSubIntervArraySize = 0;
		inStr >> MaxSubIntervArraySize;
		if(MaxSubIntervArraySize > 0)
		{
			vRelaxSubIntervArray.resize(MaxSubIntervArraySize);
			RelaxSubIntervArray = vRelaxSubIntervArray.data();
			radTRelaxSubInterval *t_RelaxSubIntervArray = RelaxSubIntervArray;
			for(int i=0; i<MaxSubIntervArraySize; i++)
			{
				inStr >> (t_RelaxSubIntervArray->StartNo);
				inStr >> (t_RelaxSubIntervArray->FinNo);
				int subIntervalID = 0;
				inStr >> subIntervalID;
				t_RelaxSubIntervArray->SubIntervalID = (TRelaxSubIntervalID)subIntervalID;
				t_RelaxSubIntervArray++;
			}
		}
	}

	//radVectPtr_lphgPtr IntVectOfPtrToListsOfTransPtr; //required
	DumpBinParseVectOfPtrToListsOfTransPtr(inStr, mKeysOldNew, gMapOfHandlers, IntVectOfPtrToListsOfTransPtr);

	//radVectPtr_lphgPtr ExtVectOfPtrToListsOfTransPtr; //required
	DumpBinParseVectOfPtrToListsOfTransPtr(inStr, mKeysOldNew, gMapOfHandlers, ExtVectOfPtrToListsOfTransPtr);

	//radTCast Cast; //no members?
	//radTSend Send; //no members?

	//short FillInMainTransOnly;
	inStr >> FillInMainTransOnly;

	//char mKeepTransData;
	inStr >> mKeepTransData;

	//radTrans** MainTransPtrArray; //required
	MainTransPtrArray= 0;
	int size_vIndMainTrans = 0;
	inStr >> size_vIndMainTrans;
	if(size_vIndMainTrans > 0)
	{
		vMainTransPtrArray.resize(AmOfMainElem);
		MainTransPtrArray = vMainTransPtrArray.data();

		for(int i=0; i<AmOfMainElem; i++)
		{
			radThg hg;
			int oldKey = DumpBinParseSourceHandle(inStr, mKeysOldNew, gMapOfHandlers, false, false, hg);
			if(oldKey < 0) MainTransPtrArray[i] = IdentTransPtr;
			else if(hg.rep != 0) MainTransPtrArray[i] = new radTrans(*((radTrans*)hg.rep)); //will be deleted at distraction
		}
	}

	//int AmOfRelaxSubInterv;
	inStr >> AmOfRelaxSubInterv;

	//short SomethingIsWrong;
	inStr >> SomethingIsWrong;

	//short MemAllocTotAtOnce;
	inStr >> MemAllocTotAtOnce;
}

//-------------------------------------------------------------------------
