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
//-------------------------------------------------------------------------

void radTSimpleRelaxation::DefineNewMagnetizations()
{
	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	// Dense matrix-vector multiplication
	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d H_atElemStrNo(0.,0.,0.);
		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
			H_atElemStrNo += (IntrctPtr->InteractMatrix[StrNo][ColNo])*(IntrctPtr->NewMagnArray[ColNo]);

		IntrctPtr->NewFieldArray[StrNo] = H_atElemStrNo + IntrctPtr->ExternFieldArray[StrNo];
	}

	double One_mi_RelaxParam = 1.- RelaxParam;
	radTg3dRelax* g3dRelaxPtr = nullptr;
	for(int StNo=0; StNo<LocAmOfMainElem; StNo++)
	{
		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StNo];

		IntrctPtr->NewMagnArray[StNo] = One_mi_RelaxParam*g3dRelaxPtr->Magn;
		g3dRelaxPtr->Magn = ((radTMaterial*)(g3dRelaxPtr->MaterHandle.rep))->M(IntrctPtr->NewFieldArray[StNo]);
		IntrctPtr->NewMagnArray[StNo] += RelaxParam*g3dRelaxPtr->Magn;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRelaxationMethNo_2::DefineNewMagnetizations()
{
	TMatrix3d InstantKsiTensor;
	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;
	int AmOfMainElem_mi_One = LocAmOfMainElem - 1;

	TVector3d* OldField = IntrctPtr->NewMagnArray;
	TVector3d* NewField = IntrctPtr->NewFieldArray;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504

	// Dense matrix multiplication
	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d H_atElemStrNo(0.,0.,0.);
		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
			H_atElemStrNo += (IntrcMat[StrNo][ColNo])*((IntrctPtr->g3dRelaxPtrVect[ColNo])->Magn);

		OldField[StrNo] = NewField[StrNo];
		NewField[StrNo] = H_atElemStrNo + IntrctPtr->ExternFieldArray[StrNo];
	}

	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.), MagnFromMaterRel, InstantMr; // The later is not actually used here
	TMatrix3d E(E_Str0, E_Str1, E_Str2), mi_Eta, E_pl_Eta, InvE_pl_Eta;
	radTg3dRelax* g3dRelaxPtr = nullptr;
	double One_mi_RelaxParam = 1.- RelaxParam;
	for(int StNo=0; StNo<LocAmOfMainElem; StNo++)
	{
		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StNo];

		((radTMaterial*)(g3dRelaxPtr->MaterHandle.rep))->
			DefineInstantKsiTensor(OldField[StNo], InstantKsiTensor, InstantMr);
		mi_Eta = InstantKsiTensor*IntrcMat[StNo][StNo];

		E_pl_Eta = E - mi_Eta;
		Matrix3d_inv(E_pl_Eta, InvE_pl_Eta);

		MagnFromMaterRel = ((radTMaterial*)(g3dRelaxPtr->MaterHandle.rep))->M(NewField[StNo]);
		g3dRelaxPtr->Magn = RelaxParam*(InvE_pl_Eta*(MagnFromMaterRel - (mi_Eta*g3dRelaxPtr->Magn)))
							+ One_mi_RelaxParam*g3dRelaxPtr->Magn;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRelaxationMethNo_3::DefineNewMagnetizations()
{
	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.);
	TMatrix3d E(E_Str0, E_Str1, E_Str2), BufMatr, InvBufMatr;
	TMatrix3d MultByInstKsi;
	TVector3d MultByInstMr, Mnew_mi_MoldVect;
	double BufMisfitM=0.;

	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = nullptr;
	radTMaterial* MaterPtr = nullptr;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	int AmOfMainElem_mi_One = LocAmOfMainElem - 1;
	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);
		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
			if(ColNo!=StrNo) QuasiExtFieldAtElemStrNo += IntrcMat[StrNo][ColNo] * MagnAr[ColNo];
		QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		MaterPtr->MultMatrByInstKsiAndMr(NewFieldAr[StrNo], IntrcMat[StrNo][StrNo], MultByInstKsi, MultByInstMr);

		BufMatr = E - MultByInstKsi;
		Matrix3d_inv(BufMatr, InvBufMatr);
		NewFieldAr[StrNo] = InvBufMatr * (MultByInstMr + QuasiExtFieldAtElemStrNo);

		MagnAr[StrNo] = MaterPtr->M(NewFieldAr[StrNo]);

		Mnew_mi_MoldVect = MagnAr[StrNo] - g3dRelaxPtr->Magn;
		BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y 
					+ Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;

		g3dRelaxPtr->Magn = MagnAr[StrNo];
	}
	InstMisfitM = sqrt(BufMisfitM/LocAmOfMainElem);
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_3::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM(); // Consider removing
	}

	int IterCount = 0;
	while(InstMisfitM > PrecOnMagnetiz)
	{
		if(++IterCount > MaxIterNumber) break;
		DefineNewMagnetizations();

		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this
	}

	IntrctPtr->RelaxStatusParam.MisfitM = -1.;
	ComputeRelaxStatusParam(IntrctPtr->NewMagnArray, nullptr, IntrctPtr->NewFieldArray);
	IntrctPtr->RelaxStatusParam.MisfitM = InstMisfitM;

	return IterCount-1;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTRelaxationMethNo_a5::radTRelaxationMethNo_a5(radTInteraction* InInteractionPtr) : radTIterativeRelaxMeth(InInteractionPtr) 
{ 
	IntrctPtr = InInteractionPtr; InstMisfitM = 1.E+23;

	std::vector<radTRelaxSubInterval> vTmpSubIntervArray(IntrctPtr->AmOfRelaxSubInterv);
	radTRelaxSubInterval* TmpSubIntervArray = vTmpSubIntervArray.data();

	int RelaxTogetherCount = 0;
	int MaxSize = 0;

	for(int i=0; i<IntrctPtr->AmOfRelaxSubInterv; i++)
	{
		radTRelaxSubInterval& LocSubInterv = IntrctPtr->RelaxSubIntervArray[i];
		if(LocSubInterv.SubIntervalID == TRelaxSubIntervalID::RelaxTogether)
		{
			TmpSubIntervArray[RelaxTogetherCount] = LocSubInterv;
			RelaxTogetherCount++;

			int CurSize = LocSubInterv.FinNo - LocSubInterv.StartNo + 1;
			if(CurSize>MaxSize) MaxSize = CurSize;
		}
	}
	AmOfRelaxTogether = RelaxTogetherCount;

	SizeOfAuxs = 3*MaxSize;

	MathMethPtr = new radTMathLinAlgEq(SizeOfAuxs);

	if(AmOfRelaxTogether != 0)
	{
		vAuxMatr1Storage.resize(SizeOfAuxs);
		vAuxMatr2Storage.resize(SizeOfAuxs);
		vAuxMatr1.resize(SizeOfAuxs);
		vAuxMatr2.resize(SizeOfAuxs);

		vAuxArray.resize(SizeOfAuxs);
		AuxArray = vAuxArray.data();

		for(int m=0; m<SizeOfAuxs; m++)
		{
			vAuxMatr1Storage[m].resize(SizeOfAuxs);
			vAuxMatr2Storage[m].resize(SizeOfAuxs);
			vAuxMatr1[m] = vAuxMatr1Storage[m].data();
			vAuxMatr2[m] = vAuxMatr2Storage[m].data();
		}
		AuxMatr1 = vAuxMatr1.data();
		AuxMatr2 = vAuxMatr2.data();
	}
	// Automatic cleanup via RAII 
}

//-------------------------------------------------------------------------

radTRelaxationMethNo_a5::~radTRelaxationMethNo_a5()
{
	// Automatic cleanup via RAII (std::vector)
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_a5::DefineNewMagnetizations()
{

	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.), ZeroVect(0.,0.,0.);
	TMatrix3d E(E_Str0, E_Str1, E_Str2), BufMatr, InvBufMatr, ZeroMatr(ZeroVect, ZeroVect, ZeroVect);
	TMatrix3d MultByInstKsi;
	TVector3d MultByInstMr, Mnew_mi_MoldVect;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = nullptr;
	radTMaterial* MaterPtr = nullptr;

	double BufMisfitM=0.;

	int StrNo = 0;
	int RelaxTogetherCount = -1;

	for(int IntrvNo=0; IntrvNo<IntrctPtr->AmOfRelaxSubInterv; IntrvNo++)
	{
		radTRelaxSubInterval& CurrentSubInterv = IntrctPtr->RelaxSubIntervArray[IntrvNo];

		if(CurrentSubInterv.SubIntervalID == TRelaxSubIntervalID::RelaxTogether)
		{
			RelaxTogetherCount++;
			for(StrNo = CurrentSubInterv.StartNo; StrNo <= CurrentSubInterv.FinNo; StrNo++)
			{
				TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);
				int ColNo=0;
				for(ColNo = 0; ColNo < CurrentSubInterv.StartNo; ColNo++)
					QuasiExtFieldAtElemStrNo += IntrcMat[StrNo][ColNo] * MagnAr[ColNo];
				for(ColNo = CurrentSubInterv.FinNo+1; ColNo < LocAmOfMainElem; ColNo++)
					QuasiExtFieldAtElemStrNo += IntrcMat[StrNo][ColNo] * MagnAr[ColNo];
				QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

				int AuxMatrStrNo = 3*(StrNo - CurrentSubInterv.StartNo);
				int AuxMatrStrNo_p1=AuxMatrStrNo+1, AuxMatrStrNo_p2=AuxMatrStrNo+2;

				double* AuxMatr1StrNoPtr = AuxMatr1[AuxMatrStrNo];
				double* AuxMatr1StrNo_p1Ptr = AuxMatr1[AuxMatrStrNo_p1];
				double* AuxMatr1StrNo_p2Ptr = AuxMatr1[AuxMatrStrNo_p2];

				g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
				MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

				TVector3d ContribFromMr(0.,0.,0.);
				for(ColNo = CurrentSubInterv.StartNo; ColNo <= CurrentSubInterv.FinNo; ColNo++)
				{
					MaterPtr->MultMatrByInstKsiAndMr(NewFieldAr[ColNo], IntrcMat[StrNo][ColNo], MultByInstKsi, MultByInstMr);

					if(ColNo == StrNo) BufMatr = E - MultByInstKsi;
					else BufMatr = ZeroMatr - MultByInstKsi;
					
					ContribFromMr += MultByInstMr;

					int AuxMatrColNo = 3*(ColNo - CurrentSubInterv.StartNo);
					int AuxMatrColNo_p1=AuxMatrColNo+1, AuxMatrColNo_p2=AuxMatrColNo+2;

					TVector3d& BufMatrStr0 = BufMatr.Str0;
					AuxMatr1StrNoPtr[AuxMatrColNo] = BufMatrStr0.x; AuxMatr1StrNoPtr[AuxMatrColNo_p1] = BufMatrStr0.y; AuxMatr1StrNoPtr[AuxMatrColNo_p2] = BufMatrStr0.z;
					TVector3d& BufMatrStr1 = BufMatr.Str1;
					AuxMatr1StrNo_p1Ptr[AuxMatrColNo] = BufMatrStr1.x; AuxMatr1StrNo_p1Ptr[AuxMatrColNo_p1] = BufMatrStr1.y; AuxMatr1StrNo_p1Ptr[AuxMatrColNo_p2] = BufMatrStr1.z;
					TVector3d& BufMatrStr2 = BufMatr.Str2;
					AuxMatr1StrNo_p2Ptr[AuxMatrColNo] = BufMatrStr2.x; AuxMatr1StrNo_p2Ptr[AuxMatrColNo_p1] = BufMatrStr2.y; AuxMatr1StrNo_p2Ptr[AuxMatrColNo_p2] = BufMatrStr2.z;
				}

				QuasiExtFieldAtElemStrNo += ContribFromMr;
				AuxArray[AuxMatrStrNo] = QuasiExtFieldAtElemStrNo.x;
				AuxArray[AuxMatrStrNo_p1] = QuasiExtFieldAtElemStrNo.y;
				AuxArray[AuxMatrStrNo_p2] = QuasiExtFieldAtElemStrNo.z;
			}
			int SizeMatr = 3*(CurrentSubInterv.FinNo - CurrentSubInterv.StartNo + 1);
			MathMethPtr->InverseMatrix(AuxMatr1, SizeMatr, AuxMatr2);

			int TriplCount = -1;
			int NewFieldArIndx = CurrentSubInterv.StartNo;

			for(int i=0; i<SizeMatr; i++)
			{
				TriplCount++;

				double Sum=0.;
				double* AuxMatr2_i = AuxMatr2[i];
				for(int j=0; j<SizeMatr; j++) Sum += AuxMatr2_i[j] * AuxArray[j];

				if(TriplCount==0) NewFieldAr[NewFieldArIndx].x = Sum;
				else if(TriplCount==1) NewFieldAr[NewFieldArIndx].y = Sum;
				else
				{
					NewFieldAr[NewFieldArIndx].z = Sum;

					NewFieldArIndx++;
					TriplCount = -1;
				}
			}
			for(StrNo = CurrentSubInterv.StartNo; StrNo <= CurrentSubInterv.FinNo; StrNo++)
			{
				g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
				MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

				MagnAr[StrNo] = MaterPtr->M(NewFieldAr[StrNo]);
				Mnew_mi_MoldVect = MagnAr[StrNo] - g3dRelaxPtr->Magn;
				BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y 
							+ Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;
				g3dRelaxPtr->Magn = MagnAr[StrNo];
			}
		}
		
		if(CurrentSubInterv.SubIntervalID == TRelaxSubIntervalID::RelaxApart)
		{
			for(StrNo = CurrentSubInterv.StartNo; StrNo <= CurrentSubInterv.FinNo; StrNo++)
			{
				TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);
				for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
					if(ColNo!=StrNo) QuasiExtFieldAtElemStrNo += IntrcMat[StrNo][ColNo] * MagnAr[ColNo];
				QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

				g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
				MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
				MaterPtr->MultMatrByInstKsiAndMr(NewFieldAr[StrNo], IntrcMat[StrNo][StrNo], MultByInstKsi, MultByInstMr);

				BufMatr = E - MultByInstKsi;
				Matrix3d_inv(BufMatr, InvBufMatr);
				NewFieldAr[StrNo] = InvBufMatr * (MultByInstMr + QuasiExtFieldAtElemStrNo);

				MagnAr[StrNo] = MaterPtr->M(NewFieldAr[StrNo]);

				Mnew_mi_MoldVect = MagnAr[StrNo] - g3dRelaxPtr->Magn;
				BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y 
							+ Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;
				g3dRelaxPtr->Magn = MagnAr[StrNo];
			}
		}
	}
	InstMisfitM = sqrt(BufMisfitM/LocAmOfMainElem);
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_a5::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{ // Absolutely the same as for MethNo_3

	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();  // Consider removing
	}

	int IterCount = 0;
	while(InstMisfitM > PrecOnMagnetiz)
	{
		if(++IterCount > MaxIterNumber) break;
		DefineNewMagnetizations();

		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this
	}

	IntrctPtr->RelaxStatusParam.MisfitM = -1.;
	ComputeRelaxStatusParam(IntrctPtr->NewMagnArray, nullptr, IntrctPtr->NewFieldArray);
	IntrctPtr->RelaxStatusParam.MisfitM = InstMisfitM;

	return IterCount-1;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// Good operational version; failed to converge in case of fine mesh of Soleil booster quads
// restored 150505, because another version appears to fail to converge even in the Radia dipole magnet example
void radTRelaxationMethNo_4::DefineNewMagnetizations()
{
	//TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.);
	//TMatrix3d E(E_Str0, E_Str1, E_Str2), BufMatr, InvBufMatr;
	//TMatrix3d MultByInstKsi;
	//TVector3d MultByInstMr, Mnew_mi_MoldVect;
	TVector3d Mnew_mi_MoldVect;
	double BufMisfitM=0.;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d** IntrcMat = IntrctPtr->InteractMatrix; //OC250504

	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = nullptr;
	radTMaterial* MaterPtr = nullptr;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	double NormFact = 1./double(LocAmOfMainElem);
	double BestPrecMagnE2 =  DesiredPrecOnMagnetizE2*NormFact;
	double LocPrecMagnE2 = (InstMisfitMe2 > 1.E+20)? BestPrecMagnE2 : 0.25*NormFact*InstMisfitMe2;
	//double LocPrecMagnE2 = (InstMisfitMe2 > 1.E+20)? BestPrecMagnE2 : 0.1*NormFact*InstMisfitMe2;

	if(LocPrecMagnE2 < BestPrecMagnE2) LocPrecMagnE2 = BestPrecMagnE2;
	radTRelaxAuxData *tRelaxAuxData = mpRelaxAuxData; //OC06112003
	const int MaxConseqBadPasses = 1; //OC06112003

	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);
		TMatrix3df* MatrArrayPtr = IntrcMat[StrNo]; //OC250504
		//TMatrix3d* MatrArrayPtr = IntrcMat[StrNo]; //OC250504

		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
		{
			if(ColNo!=StrNo) QuasiExtFieldAtElemStrNo += MatrArrayPtr[ColNo] * MagnAr[ColNo];
		}
		QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		TVector3d& InstantH = NewFieldAr[StrNo];
		TVector3d& InstantM = MagnAr[StrNo];

		TVector3d PrevH = InstantH;

		MaterPtr->FindNewH(InstantH, MatrArrayPtr[StrNo], QuasiExtFieldAtElemStrNo, LocPrecMagnE2);
		//MaterPtr->FindNewH(InstantH, MatrArrayPtr[StrNo], QuasiExtFieldAtElemStrNo, LocPrecMagnE2, g3dRelaxPtr);
		//MaterPtr->FindNewH(InstantH, MatrArrayPtr[StrNo], QuasiExtFieldAtElemStrNo, LocPrecMagnE2, g3dRelaxPtr, gpRelaxAuxData + StrNo); //OC140103

		//InstantH = (tRelaxAuxData->mRelaxPar)*InstantH + (1. - tRelaxAuxData->mRelaxPar)*PrevH; //OC131103
		//OC: commented out 150304

		//InstantM = MaterPtr->M(InstantH);
		
		TVector3d PureNewM = MaterPtr->M(InstantH);
		InstantM = PureNewM;

		Mnew_mi_MoldVect = PureNewM - g3dRelaxPtr->Magn;
		double NewDifMe2 = Mnew_mi_MoldVect.AmpE2(); //OC06112003
		BufMisfitM += NewDifMe2;

		//tRelaxAuxData->Update(NewDifMe2, MaxConseqBadPasses, mRelaxParModFact, mRelaxParMin); //OC06112003
		//OC: commented out 150304

		//InstantM = mRelaxPar*PureNewM + (1. - mRelaxPar)*(g3dRelaxPtr->Magn); //OC041103
		//InstantM = (tRelaxAuxData->mRelaxPar)*PureNewM + (1. - tRelaxAuxData->mRelaxPar)*(g3dRelaxPtr->Magn); //OC041103
		//InstantM = PureNewM;

		//tRelaxAuxData++; //OC041103
		//OC: commented out 150304

		//Mnew_mi_MoldVect = InstantM - g3dRelaxPtr->Magn;
		//BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y + Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;

		//double CurRelaxPar = (gpRelaxAuxData + StrNo)->RelaxPar; //OC140103
		//BufMisfitM += (Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y + Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z)/(CurRelaxPar*CurRelaxPar); //OC140103

		g3dRelaxPtr->Magn = InstantM; 
	}
	double NewInstMisfitMe2 = BufMisfitM/LocAmOfMainElem;

	//if(NewInstMisfitMe2 > mMisfitE2RatToStartModifRelaxPar*InstMisfitMe2) //OC041103
	//{
	//	if(mRelaxPar > mRelaxParMin) mRelaxPar *= mRelaxParModFact;
	//	int Aha = 1;
	//}
	InstMisfitMe2 = NewInstMisfitMe2;
}

//-------------------------------------------------------------------------
// Strange version; failed to converge in case of Radia dipole magnet example
// commented out 150505
/**
void radTRelaxationMethNo_4::DefineNewMagnetizations()
{
	TVector3d Mnew_mi_MoldVect;
	double BufMisfitM=0.;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;

	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = nullptr;
	radTMaterial* MaterPtr = nullptr;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	double NormFact = 1./double(LocAmOfMainElem);
	double BestPrecMagnE2 =  DesiredPrecOnMagnetizE2*NormFact;
	double LocPrecMagnE2 = (InstMisfitMe2 > 1.E+20)? BestPrecMagnE2 : 0.25*NormFact*InstMisfitMe2;

	if(LocPrecMagnE2 < BestPrecMagnE2) LocPrecMagnE2 = BestPrecMagnE2;
	//radTRelaxAuxData *tRelaxAuxData = mpRelaxAuxData; //OC06112003
	//const int MaxConseqBadPasses = 1; //OC06112003
	const int MaxInitialBadPasses = 10; //OC06112003

	//if((!mKeepPrevOldValues) && ((mNumConvergPasses > 0) || ((mNumConvergPasses == 0) && (mNumDivergPasses == 0)) || ((mIterCount < MaxInitialBadPasses)))) IntrctPtr->StoreAuxOldArrays(); //OC300504
	if(!mKeepPrevOldValues) IntrctPtr->StoreAuxOldArrays(); //OC300504

	double *tElemVolumes = mElemVolumeArray; //OC010604
	const double ConstForM = -1./(400.*(3.14159265358979));
	double SumEnergy = 0;

	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TMatrix3df* MatrArrayPtr = IntrcMat[StrNo]; //OC250504
		TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);

		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
		{
			if(ColNo!=StrNo) QuasiExtFieldAtElemStrNo += MatrArrayPtr[ColNo] * MagnAr[ColNo];
		}
		QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		TVector3d& InstantH = NewFieldAr[StrNo];
		TVector3d& InstantM = MagnAr[StrNo];
		TVector3d PrevH = InstantH;

		if(StrNo == 42)
		{
			int aha = 1;
		}

		MaterPtr->FindNewH(InstantH, MatrArrayPtr[StrNo], QuasiExtFieldAtElemStrNo, LocPrecMagnE2);

		TVector3d PureNewM = MaterPtr->M(InstantH);
		Mnew_mi_MoldVect = PureNewM - g3dRelaxPtr->Magn;
		double NewDifMe2 = Mnew_mi_MoldVect.AmpE2(); //OC06112003
		BufMisfitM += NewDifMe2;

		//InstantH = mRelaxPar*InstantH + (1. - mRelaxPar)*PrevH; //OC260504
		//InstantM = MaterPtr->M(InstantH);

		//InstantM = mRelaxPar*PureNewM + (1. - mRelaxPar)*(g3dRelaxPtr->Magn); //OC041103
		//InstantM = 0.5*(MaterPtr->M(InstantH) + (mRelaxPar*PureNewM + (1. - mRelaxPar)*(g3dRelaxPtr->Magn))); //OC280504

		//Mnew_mi_MoldVect = InstantM - g3dRelaxPtr->Magn;
		//BufMisfitM += Mnew_mi_MoldVect.x*Mnew_mi_MoldVect.x + Mnew_mi_MoldVect.y*Mnew_mi_MoldVect.y + Mnew_mi_MoldVect.z*Mnew_mi_MoldVect.z;

		InstantM = PureNewM;
		g3dRelaxPtr->Magn = InstantM; 

		SumEnergy += (InstantM*(InstantH + InstantM))*(*(tElemVolumes++)); //OC010604
	}
	SumEnergy *= ConstForM; //OC010604
	
	//InstMisfitMe2 = BufMisfitM*NormFact;
	double NewInstMisfitMe2 = BufMisfitM*NormFact;

	//if((NewInstMisfitMe2 > InstMisfitMe2) && (mIterCount > MaxInitialBadPasses)) //OC053004
	//{
	//       //if(mNumConvergPasses > 0) IntrctPtr->StoreAuxOldArrays();
	//       if(mRelaxPar > mRelaxParMin) 
	//	{
	//		mRelaxPar *= mRelaxParModFact; //OC300504
	//           mKeepPrevOldValues = true;
	//	}
	//	else mKeepPrevOldValues = false;

	//       mNumConvergPasses = 0;
	//       mNumDivergPasses++;

	////	//if(mNumDivergPasses > mNumDivergPassesMax)
	////	//{
	////       //    if(mRelaxPar > mRelaxParMin) mRelaxPar *= mRelaxParModFact;
	////	//           int Aha = 1;
	////	//		////test
	////	//		//char ErrorMesTitle[] = "Radia Debug";
	////	//		//char ErrorStr[100];
	////	//		//int j = sprintf(ErrorStr, "mRelaxPar: %g ", mRelaxPar);
	////	//		//j += sprintf(ErrorStr + j, "          NewInstMisfitMe2: %g   ", NewInstMisfitMe2);
	////	//		//UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	////	//		//int MesBoxInf = MessageBox(nullptr, ErrorStr, ErrorMesTitle, DlgStyle); 
	////	//		////end test
	////	//}
	//}
	//else
	//{
	//       mNumDivergPasses = 0;
	//       mNumConvergPasses++;
	//	mKeepPrevOldValues = false;
	//	mSysEnergyMin = SumEnergy;

	//	//mRelaxPar = 1.; //?????

	////	//      if(mNumConvergPasses > mNumConvergPassesMax)
	////	//{
	////	//          if(mRelaxPar < 1.) 
	////	//	{
	////	//		mRelaxPar /= mRelaxParModFact;
	////	//		if(mRelaxPar > 1.) mRelaxPar = 1.;
	////	//		mNumConvergPasses = 0;
	////	//	}
	////	//          //int Aha = 1;
	////	//}
	////	//      InstMisfitMe2 = NewInstMisfitMe2;
	//}

	if((NewInstMisfitMe2 > InstMisfitMe2) && (mIterCount > MaxInitialBadPasses))
	{
		mNumConvergPasses = 0;
		mNumDivergPasses++;
		if(mNumDivergPasses >= MaxInitialBadPasses) mBadConverg = true;
	}
	else
	{
		mNumDivergPasses = 0;
		mNumConvergPasses++;
	}

	if(mBadConverg)
	{
		if(SumEnergy < mSysEnergyMin)
		{
			StoreOptimValuesFromOldArrays();
			mRelaxPar = 1.; //?????
			mKeepPrevOldValues = false;
			mSysEnergyMin = SumEnergy;
		}
		else
		{
			if(SumEnergy < mSysEnergy)
			{
				mKeepPrevOldValues = false;
			}
			else
			{
				//if(mRelaxPar > mRelaxParMin)
				//{
				mRelaxPar *= mRelaxParModFact;
				mKeepPrevOldValues = true;
				RestoreOptimValuesToOldArrays();
				//}
				//else
				//{
				//	mKeepPrevOldValues = false;
				//}
			}
		}
		if(mRelaxPar < 1.) CorrectMagnAndFieldArraysWithRelaxPar();
	}

	mSysEnergy = SumEnergy;
	InstMisfitMe2 = NewInstMisfitMe2;
}
**/
//-------------------------------------------------------------------------

int radTRelaxationMethNo_4::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	DesiredPrecOnMagnetizE2 = PrecOnMagnetiz * PrecOnMagnetiz;

	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	//SetupAuxArrays(); //OC140103
	//SetupElemVolumeArray(); //OC150505 //OC010604
	//SetupOptimValuesArrays(); //OC150505 //OC020604

	mMethNo = 1;
	mIterCount = 0;
	mSysEnergy = 0;
	mSysEnergyMin = 1e+23;
	mKeepPrevOldValues = false;
	mBadConverg = false;

	double MinInstMisfitMe2 = 1.e+30;
	while(InstMisfitMe2 > DesiredPrecOnMagnetizE2)
	{
		if(++mIterCount > MaxIterNumber) break;
		DefineNewMagnetizations();

		//if(MinInstMisfitMe2 > InstMisfitMe2) 
		//{
		//	MinInstMisfitMe2 = InstMisfitMe2;
		//	IntrctPtr->StoreAuxOldArrays();
		//}

		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this

		//test
		//mRelaxPar = 1./pow((double)IterCount + 1., 0.35);
		//mRelaxPar *= pow(((double)mIterCount)/((double)mIterCount + 1.), 0.4);
		//end test
	}

	//if(mIterCount > MaxIterNumber)
	//{
	//	IntrctPtr->RestoreAuxOldArrays();
	//	InstMisfitMe2 = MinInstMisfitMe2;
	//}

	IntrctPtr->RelaxStatusParam.MisfitM = -1.;
	ComputeRelaxStatusParam(IntrctPtr->NewMagnArray, nullptr, IntrctPtr->NewFieldArray);
	
	IntrctPtr->RelaxStatusParam.MisfitM = sqrt(InstMisfitMe2);

	//DeleteAuxArrays(); //OC140103
	//DeleteElemVolumeArray(); //OC150505 //OC010604
	//DeleteOptimValuesArrays(); //OC150505 //OC020604

	return mIterCount-1;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_4::CorrectMagnAndFieldArraysWithRelaxPar() //OC300504
{
	if((mRelaxPar < 0.) || (mRelaxPar > 1.)) return;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	TVector3d *tMagnAr = IntrctPtr->NewMagnArray;
	TVector3d *tFieldAr = IntrctPtr->NewFieldArray;
	TVector3d *tOldMagnAr = IntrctPtr->AuxOldMagnArray;
	TVector3d *tOldFieldAr = IntrctPtr->AuxOldFieldArray;

	radTg3dRelax* g3dRelaxPtr = nullptr;
	radTMaterial* MaterPtr = nullptr;

	for(int k=0; k<LocAmOfMainElem; k++)
	{
		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[k];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		*tFieldAr = mRelaxPar*(*tFieldAr) + (1. - mRelaxPar)*(*tOldFieldAr);
		*tMagnAr = MaterPtr->M(*tFieldAr);
		//*tMagnAr = mRelaxPar*(*tMagnAr) + (1. - mRelaxPar)*(*tOldMagnAr);

		g3dRelaxPtr->Magn = *tMagnAr;
		tMagnAr++; tFieldAr++; tOldMagnAr++; tOldFieldAr++;
	}
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_4::DefineNewMagnetizationsTest()
{
	double MaxNumDivergPassesToSwitchMeth = 5;

	TVector3d Mnew_mi_MoldVect;
	double BufMisfitM=0.;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	//TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = nullptr;
	radTMaterial* MaterPtr = nullptr;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;

	double NormFact = 1./double(LocAmOfMainElem);
	//double BestPrecMagnE2 =  DesiredPrecOnMagnetizE2*NormFact;
	//double LocPrecMagnE2 = (InstMisfitMe2 > 1.E+20)? BestPrecMagnE2 : 0.25*NormFact*InstMisfitMe2;
	//if(LocPrecMagnE2 < BestPrecMagnE2) LocPrecMagnE2 = BestPrecMagnE2;

	if((mMethNo == 1) && (mNumDivergPasses >= MaxNumDivergPassesToSwitchMeth)) mMethNo = 2;

	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d& InstantH = NewFieldAr[StrNo];
		TVector3d PrevH = InstantH;

		InstantH = FindNewFieldTwoSteps(StrNo, mMethNo);

		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		TVector3d& InstantM = MagnAr[StrNo];

		InstantM = MaterPtr->M(InstantH);
		Mnew_mi_MoldVect = InstantM - g3dRelaxPtr->Magn;
		double NewDifMe2 = Mnew_mi_MoldVect.AmpE2();
		BufMisfitM += NewDifMe2;

		g3dRelaxPtr->Magn = InstantM; 
	}
	double NewInstMisfitMe2 = BufMisfitM*NormFact;

	if(NewInstMisfitMe2 > InstMisfitMe2)
	{
		mNumConvergPasses = 0;
		mNumDivergPasses++;
		//if(mNumDivergPasses > mNumDivergPassesMax)
		//{
		//          if(mRelaxPar > mRelaxParMin) mRelaxPar *= mRelaxParModFact;
		//          int Aha = 1;
		//	////test
		//	//char ErrorMesTitle[] = "Radia Debug";
		//	//char ErrorStr[100];
		//	//int j = sprintf(ErrorStr, "mRelaxPar: %g ", mRelaxPar);
		//	//j += sprintf(ErrorStr + j, "          NewInstMisfitMe2: %g   ", NewInstMisfitMe2);
		//	//UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
		//	//int MesBoxInf = MessageBox(nullptr, ErrorStr, ErrorMesTitle, DlgStyle); 
		//	////end test
		//}
	}
	else
	{
		mNumDivergPasses = 0;
		mNumConvergPasses++;
		//if(mNumConvergPasses > mNumConvergPassesMax)
		//{
		//          if(mRelaxPar < 1.) 
		//	{
		//		mRelaxPar /= mRelaxParModFact;
		//		if(mRelaxPar > 1.) mRelaxPar = 1.;
		//	}
		//          int Aha = 1;
		//}
	}

	InstMisfitMe2 = BufMisfitM*NormFact;
}

//-------------------------------------------------------------------------

TVector3d radTRelaxationMethNo_4::FindNewFieldTwoSteps(int i, int MethNo)
{
	const double DesiredPrecOnMagn = 1.e-08;
	//const double DesiredPrecOnMagnE2 = DesiredPrecOnMagn*DesiredPrecOnMagn;
	int MaxLoopPasses = 10;
	double AbsMisfitM = 5;

	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.);
	TMatrix3d E(E_Str0, E_Str1, E_Str2), BufMatr, InvBufMatr, LocBufMatr, InvLocBufMatr, ExtBufMatr, InvExtBufMatr;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;

	radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
	radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

	int N = IntrctPtr->AmOfMainElem;
	TVector3d& Hi = NewFieldAr[i];
	TVector3d& Mi = MagnAr[i];
	TMatrix3d Qii = IntrcMat[i][i];

	TMatrix3d SumMatr1;
	TVector3d SumVect1 = ExternFieldAr[i];
	for(int j=0; j<N; j++)
	{
		if(j == i) continue;

		TMatrix3d Qij = IntrcMat[i][j];

		if(MethNo == 1)
		{
			SumVect1 += Qij*MagnAr[j];
		}
		else if(MethNo == 2)
		{
			TMatrix3d Qjj = IntrcMat[j][j];
			TMatrix3d Qji = IntrcMat[j][i];

			TVector3d &Hj = NewFieldAr[j], Mr_j;
			TMatrix3d Ksi_j;
			radTg3dRelax* Loc_g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[j];
			radTMaterial* Loc_MaterPtr = (radTMaterial*)(Loc_g3dRelaxPtr->MaterHandle.rep);
			Loc_MaterPtr->DefineInstantKsiTensor(Hj, Ksi_j, Mr_j);

			//LocBufMatr = E - Qjj*Ksi_j;
			//Matrix3d_inv(LocBufMatr, InvLocBufMatr);
			//TMatrix3d Ksi_j_InvLocBufMatr = Ksi_j*InvLocBufMatr;
			//TMatrix3d Qij_Ksi_j_InvLocBufMatr = Qij*Ksi_j_InvLocBufMatr;
			//SumMatr1 += Qij_Ksi_j_InvLocBufMatr*Qji;

			SumMatr1 += (Qij*Ksi_j)*Qji;

			TVector3d LocSumVect(0,0,0);
			for(int k=0; k<N; k++)
			{
				if((k == j) || (k == i)) continue;
				LocSumVect += (IntrcMat[j][k])*MagnAr[k];
			}
			TVector3d Qjj_Mr_j_p_Qji_Mr_i_p_He_j_p_LocSumVect = (Qjj*Mr_j) + ExternFieldAr[j] + LocSumVect;

			//TVector3d Qjj_Mr_j_p_Qji_Mr_i_p_He_j_p_LocSumVect = (Qjj*Mr_j) + (Qji*Mr_i) + ExternFieldAr[j] + LocSumVect;
			//SumVect1 += Qij*((Ksi_j_InvLocBufMatr*Qjj_Mr_j_p_Qji_Mr_i_p_LocSumVect) + Mr_j);

			SumVect1 += Qij*((Ksi_j*Qjj_Mr_j_p_Qji_Mr_i_p_He_j_p_LocSumVect) + Mr_j);
		}
	}

	for(int p=0; p<MaxLoopPasses; p++)
	{
		TMatrix3d Ksi_i;
		TVector3d Mr_i;
		MaterPtr->DefineInstantKsiTensor(Hi, Ksi_i, Mr_i);
		TVector3d PrevHi = Hi;

		TVector3d AuxVect1 = SumVect1 + (Qii*Mr_i);

		if(MethNo == 1)
		{
			ExtBufMatr = E - (Qii*Ksi_i);
		}
		else if(MethNo == 2)
		{
			AuxVect1 += (SumMatr1*Mr_i);
			ExtBufMatr = E - ((Qii + SumMatr1)*Ksi_i);
		}

		Matrix3d_inv(ExtBufMatr, InvExtBufMatr);
		Hi = InvExtBufMatr*AuxVect1;

		Mi = MaterPtr->M(Hi);
		TVector3d MiLin = (Ksi_i*Hi) + Mr_i;
		TVector3d Mi_mi_MiLin = Mi - MiLin;
		double AbsNewMisfitM = sqrt(Mi_mi_MiLin.AmpE2());

		if(AbsNewMisfitM <= DesiredPrecOnMagn) break;

		double Alpha = AbsMisfitM/(AbsMisfitM + AbsNewMisfitM);
		//AbsInstantH = Alpha*NewAbsInstantH + (1 - Alpha)*AbsInstantH;

		AbsMisfitM = AbsNewMisfitM;

		Hi = Alpha*Hi + (1 - Alpha)*PrevHi;
	}

	return Hi;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_4::LpTau(int i, double& q)
{/**
	double a = i;
	int m = 1 + int(log(a)/0.693147);

	double s = 0.;
	for(int k = 1; k <= m; k++)
	{
		int ns = 0;
		for(int l = k; l <= m; l++)
		{
			ns += int(2*D(a/pow(2,l)))*int(2*D(1./pow(2,l+1-k)));
		}
		s += D(0.5*ns)/pow(2,k-1);
	}
	q = s;
**/
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int radTRelaxationMethNo_8::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsMotNeeded)
{
	if(IntrctPtr == 0) return 0;
	mDesiredPrecOnMagnetizE2 = PrecOnMagnetiz * PrecOnMagnetiz;

	if(!MagnResetIsMotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	double MinInstMisfitMe2 = 1.e+30;
	int ItCnt=0;
	for(ItCnt=0; ItCnt<MaxIterNumber; ItCnt++)
	{
		if(ItCnt > MaxIterNumber) break;
		DefineNewMagnetizations();

		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this
	}

/**
	//SetupAuxArrays(); //OC140103
	SetupElemVolumeArray(); //OC010604
	SetupOptimValuesArrays(); //OC020604

	mMethNo = 1;
	mSysEnergy = 0;
	mSysEnergyMin = 1e+23;
	mKeepPrevOldValues = false;
	mBadConverg = false;

	while(InstMisfitMe2 > DesiredPrecOnMagnetizE2)
	{
		if(++mIterCount > MaxIterNumber) break;
		DefineNewMagnetizations();

		//if(MinInstMisfitMe2 > InstMisfitMe2) 
		//{
		//	MinInstMisfitMe2 = InstMisfitMe2;
		//	IntrctPtr->StoreAuxOldArrays();
		//}

		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this

		//test
		//mRelaxPar = 1./pow((double)IterCount + 1., 0.35);
		//mRelaxPar *= pow(((double)mIterCount)/((double)mIterCount + 1.), 0.4);
		//end test
	}

	//if(mIterCount > MaxIterNumber)
	//{
	//	IntrctPtr->RestoreAuxOldArrays();
	//	InstMisfitMe2 = MinInstMisfitMe2;
	//}

	IntrctPtr->RelaxStatusParam.MisfitM = -1.;
	ComputeRelaxStatusParam(IntrctPtr->NewMagnArray, nullptr, IntrctPtr->NewFieldArray);
	
	IntrctPtr->RelaxStatusParam.MisfitM = sqrt(InstMisfitMe2);

	//DeleteAuxArrays(); //OC140103
	DeleteElemVolumeArray(); //OC010604
	DeleteOptimValuesArrays(); //OC020604

**/

	return ItCnt;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_8::DefineNewMagnetizations()
{
	TVector3d E_Str0(1.,0.,0.), E_Str1(0.,1.,0.), E_Str2(0.,0.,1.), MatrElemByInstMr, Mnew_mi_MoldVect;
	TMatrix3d E(E_Str0, E_Str1, E_Str2), BufMatr, InvBufMatr, MatrElemByInstKsi;

	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;
	radTg3dRelax* g3dRelaxPtr = nullptr;
	radTMaterial* MaterPtr = nullptr;

	int LocAmOfMainElem = IntrctPtr->AmOfMainElem;
	double NormFact = 1./double(LocAmOfMainElem);
	double BufMisfitM=0.;

	for(int StrNo=0; StrNo<LocAmOfMainElem; StrNo++)
	{
		TVector3d QuasiExtFieldAtElemStrNo(0.,0.,0.);
		TMatrix3df* MatrArrayPtr = IntrcMat[StrNo];

		for(int ColNo=0; ColNo<LocAmOfMainElem; ColNo++)
		{
			if(ColNo!=StrNo) QuasiExtFieldAtElemStrNo += MatrArrayPtr[ColNo] * MagnAr[ColNo];
		}
		QuasiExtFieldAtElemStrNo += ExternFieldAr[StrNo];

		g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
		MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		TVector3d& InstantH = NewFieldAr[StrNo];
		TVector3d& InstantM = MagnAr[StrNo];

		MaterPtr->MultMatrByInstKsiAndMr(InstantH, MatrArrayPtr[StrNo], MatrElemByInstKsi, MatrElemByInstMr);

		BufMatr = E - MatrElemByInstKsi;
		Matrix3d_inv(BufMatr, InvBufMatr);
		InstantH = InvBufMatr*(QuasiExtFieldAtElemStrNo + MatrElemByInstMr);

		InstantM = MaterPtr->M(InstantH);

	//BufMatr = E - Matr*InstantKsiTensor;
	//Matrix3d_inv(BufMatr, InvBufMatr);
	//H = InvBufMatr*(H_Ext + Matr*RemMagn);

		Mnew_mi_MoldVect = InstantM - g3dRelaxPtr->Magn;
		double NewDifMe2 = Mnew_mi_MoldVect.AmpE2();
		BufMisfitM += NewDifMe2;

		g3dRelaxPtr->Magn = InstantM; 
	}
	mInstMisfitMe2 = BufMisfitM*NormFact;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRelaxationMethNo_6::SetupInteractionMatrices(const radThg& hg, const radTCompCriterium& CompCrit)
{
	mAmOfParts = 0;

	radTg3d* g3dPtr = Cast.g3dCast(hg.rep); 
	if(g3dPtr==0) { Send.ErrorMessage("Radia::Error003"); throw 0;}
	radTGroup* GroupPtr = Cast.GroupCast(g3dPtr); 
	if(GroupPtr==0) { Send.ErrorMessage("Radia::Error091"); throw 0;}
	mhGroup = hg;

	radThg hEmpty;
	mAmOfParts = (int)(GroupPtr->GroupMapOfHandlers.size());
	vIntrctPtr.resize(mAmOfParts);
	IntrctPtr = vIntrctPtr.data();
	radTInteraction *tIntrct = IntrctPtr;

	int LocMapCount = 0;
	for(radTmhg::const_iterator iter = GroupPtr->GroupMapOfHandlers.begin(); iter != GroupPtr->GroupMapOfHandlers.end(); ++iter)
	{
		tIntrct->Setup((*iter).second, hEmpty, CompCrit, 0, 1, 1);
		if(tIntrct->SomethingIsWrong)
		{
			// RAII: automatic cleanup via vIntrctPtr
			Send.ErrorMessage("Radia::Error116"); throw 0;
		}

		radThg hgGroup(GroupPtr->CreateGroupIncludingAllMembersExceptIt(iter));
		mMapOfPartHandlers[LocMapCount++] = hgGroup;
		tIntrct++;
	}
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_6::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, double* RelaxStatusParamArray)
{
	int TotAmOfRelaxObj = 0;
	radTInteraction *tIntrct = IntrctPtr;
	for(radTmhg::const_iterator iter = mMapOfPartHandlers.begin(); iter != mMapOfPartHandlers.end(); ++iter)
	{
		tIntrct->ResetM();
		tIntrct->ResetAuxParam();
		tIntrct->AddMoreExternField((*iter).second);
		TotAmOfRelaxObj += tIntrct->OutAmOfRelaxObjs();
		tIntrct++;
	}

	radTGroup* GroupPtr = Cast.GroupCast(Cast.g3dCast(mhGroup.rep)); 
	int ActualIterNum = 0, ActualOuterIterNum = 0;
	double DesiredPrecOnMagnetizE2 = PrecOnMagnetiz*PrecOnMagnetiz;
	double TotMisfitMe2 = 1.E+23;

	for(int i=0; i<MaxIterNumber; i++)
	{
		tIntrct = IntrctPtr;
		int PartCount = 0;
		double BufTotMisfitMe2 = 0;

		for(radTmhg::const_iterator it = GroupPtr->GroupMapOfHandlers.begin(); it != GroupPtr->GroupMapOfHandlers.end(); ++it)
		{
			if(tIntrct->OutAmOfRelaxObjs() > 0)
			{
				//tIntrct->StoreAuxOldMagnArray();
				tIntrct->StoreAuxOldArrays();
				radTRelaxationMethNo_4 RelaxMethNo_4(tIntrct);
				ActualIterNum = RelaxMethNo_4.AutoRelax(PrecOnMagnetiz, MaxIterNumber, 1);
				if(ActualIterNum >= MaxIterNumber) { Send.WarningMessage("Radia::Warning015");}

				tIntrct->SubstractOldMagn();
				UpdateExternFiledInAllIntrctExceptOne(tIntrct, (*it).second);
				tIntrct->AddOldMagn();
				BufTotMisfitMe2 += tIntrct->CalcQuadNewOldMagnDif();
			}
			tIntrct++;
		}
		TotMisfitMe2 = BufTotMisfitMe2/TotAmOfRelaxObj;
		if(TotMisfitMe2 <= DesiredPrecOnMagnetizE2)
		{
			ActualOuterIterNum = i; break;
		}
	}
	if(ActualOuterIterNum == 0) ActualOuterIterNum = MaxIterNumber;

	if(RelaxStatusParamArray != 0)
	{
		RelaxStatusParamArray[0] = sqrt(TotMisfitMe2);

		double MaxModM = 0, MaxModH = 0;
		tIntrct = IntrctPtr;
		for(int j=0; j<mAmOfParts; j++)
		{
			double LocMaxModM = 0, LocMaxModH = 0;
			(tIntrct++)->FindMaxModMandH(LocMaxModM, LocMaxModH);
			if(MaxModM < LocMaxModM) MaxModM = LocMaxModM;
			if(MaxModH < LocMaxModH) MaxModH = LocMaxModH;
		}
		RelaxStatusParamArray[1] = MaxModM;
		RelaxStatusParamArray[2] = MaxModH;
	}
	return ActualOuterIterNum;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_6::UpdateExternFiledInAllIntrctExceptOne(radTInteraction* pIntrctToSkip, const radThg& hgSrc)
{
	if((mAmOfParts <= 0) || (IntrctPtr == 0)) return;

	radTInteraction *tIntrct = IntrctPtr;
	for(int j=0; j<mAmOfParts; j++)
	{
		if(tIntrct == pIntrctToSkip) 
		{
			tIntrct++; continue;
		}
		(tIntrct++)->AddMoreExternField(hgSrc);
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRelaxationMethNo_7::SetupMainInteractionData(const radThg& hg, const radTCompCriterium& CompCrit)
{
	radTg3d* g3dPtr = radTCast::g3dCast(hg.rep); 
	if(g3dPtr==0) { radTSend::ErrorMessage("Radia::Error003"); throw 0;}

	mArrAuxQuasiExtField = nullptr;

	DeleteInterMatrData();
	radThg hEmpty;
	IntrctPtr = new radTInteraction(hg, hEmpty, CompCrit, 1, 1, 1);
	if(IntrctPtr->SomethingIsWrong)
	{
		delete[] IntrctPtr; IntrctPtr = nullptr;
	}
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_7::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, double* RelaxStatusParamArray)
{
	if(IntrctPtr == nullptr) return 0;
	int AmOfRelaxElem = IntrctPtr->OutAmOfRelaxObjs();
	if(AmOfRelaxElem <= 0) return 0;

	int AmOfSubMatr = 0;
	std::vector<int> vTotArrSubMatrNos(AmOfRelaxElem);
	int *TotArrSubMatrNos = vTotArrSubMatrNos.data();
	std::vector<int> vSubMatrLengths;
	int *SubMatrLengths=nullptr;
	int AmOfInitIter = FillInSubMatrixArrays(PrecOnMagnetiz, TotArrSubMatrNos, SubMatrLengths, AmOfSubMatr, vSubMatrLengths);
	if((TotArrSubMatrNos == nullptr) || (SubMatrLengths == nullptr) || (AmOfSubMatr <= 0)) return 0;

	vArrAuxQuasiExtField.resize(AmOfRelaxElem);
	mArrAuxQuasiExtField = vArrAuxQuasiExtField.data();

	CalcQuasiExtFieldForAll(TotArrSubMatrNos, SubMatrLengths, AmOfSubMatr);

	double DesiredPrecOnMagnetizE2 = PrecOnMagnetiz*PrecOnMagnetiz, InstMisfitMe2 = 1.e+30;
	int IterCount=0;
	try
	{
		double InvNumRelaxElem = 1./((double)AmOfRelaxElem);
		for(IterCount=0; IterCount<MaxIterNumber; IterCount++)
		{
			//IntrctPtr->StoreAuxOldMagnArray();
			IntrctPtr->StoreAuxOldArrays();
			int *tSubMatrLengths = SubMatrLengths;
			int *tTotArrSubMatrNos = TotArrSubMatrNos;
			int OffsetCurrentSubMatr = 0;

			for(int i=0; i<AmOfSubMatr; i++)
			{
				int SizeCurrentSubMatr = *(tSubMatrLengths++);
				int CurIterNum = RelaxCurrentSubMatrix(tTotArrSubMatrNos, SizeCurrentSubMatr, DesiredPrecOnMagnetizE2, MaxIterNumber);
				if(CurIterNum >= MaxIterNumber) radTSend::WarningMessage("Radia::Warning015");
				UpdateQuasiExtFieldFromCurrentSubMatrix(TotArrSubMatrNos, OffsetCurrentSubMatr, SizeCurrentSubMatr);

				OffsetCurrentSubMatr += SizeCurrentSubMatr;
				tTotArrSubMatrNos += SizeCurrentSubMatr;
			}
			InstMisfitMe2 = (IntrctPtr->CalcQuadNewOldMagnDif())*InvNumRelaxElem;
			if(InstMisfitMe2 <= DesiredPrecOnMagnetizE2) break;
			if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this
		}

		if(RelaxStatusParamArray != 0)
		{
			RelaxStatusParamArray[0] = sqrt(InstMisfitMe2);
			IntrctPtr->FindMaxModMandH(RelaxStatusParamArray[1], RelaxStatusParamArray[2]);
		}
	}
	catch(int ErrNo)
	{
		// RAII: automatic cleanup via vTotArrSubMatrNos and vSubMatrLengths
		throw ErrNo;
	}

	// RAII: automatic cleanup via vTotArrSubMatrNos and vSubMatrLengths
	return IterCount;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_7::CalcQuasiExtFieldForAll(int* TotArrSubMatrNos, int* SubMatrLengths, int AmOfSubMatr)
{
	if((TotArrSubMatrNos == nullptr) || (SubMatrLengths == nullptr) || (AmOfSubMatr == 0) || (mArrAuxQuasiExtField == nullptr)) return;

	int *tTotArrSubMatrNos = TotArrSubMatrNos;
	int *tSubMatrLengths = SubMatrLengths;

	int StartOffset = 0;
	for(int i=0; i<AmOfSubMatr; i++)
	{
		int CurSubMatrSize = *(tSubMatrLengths++);
		for(int j=0; j<CurSubMatrSize; j++)
		{
			int CurInd = *(tTotArrSubMatrNos++);
			CalcQuasiExtFieldForOneElem(CurInd, TotArrSubMatrNos, StartOffset, CurSubMatrSize);
		}
		StartOffset += CurSubMatrSize;
	}
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_7::CalcQuasiExtFieldForOneElem(int CurInd, int* TotArrSubMatrNos, int StartOffsetSkip, int SkipLength)
{
	int AmOfRelaxElem = IntrctPtr->OutAmOfRelaxObjs();
	if(StartOffsetSkip > AmOfRelaxElem) StartOffsetSkip = AmOfRelaxElem;

	TVector3d &QuasiExtField = *(mArrAuxQuasiExtField + CurInd);
	QuasiExtField.Zero();
	TMatrix3df *MatrArrayPtr = (IntrctPtr->InteractMatrix)[CurInd]; //OC250504
	//TMatrix3d *MatrArrayPtr = (IntrctPtr->InteractMatrix)[CurInd]; //OC250504
	TVector3d *MagnAr = IntrctPtr->NewMagnArray;

	int *tTotArrSubMatrNos = TotArrSubMatrNos;
	for(int i=0; i<StartOffsetSkip; i++)
	{
		int CurElemInd = *(tTotArrSubMatrNos++);
		QuasiExtField += MatrArrayPtr[CurElemInd] * MagnAr[CurElemInd];
	}
	int OffsetStart = StartOffsetSkip + SkipLength;
	if(OffsetStart > AmOfRelaxElem) OffsetStart = AmOfRelaxElem;
	tTotArrSubMatrNos = TotArrSubMatrNos + OffsetStart;
	for(int j=OffsetStart; j<AmOfRelaxElem; j++)
	{
		int CurElemInd = *(tTotArrSubMatrNos++);
		QuasiExtField += MatrArrayPtr[CurElemInd] * MagnAr[CurElemInd];
	}
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_7::AddQuasiExtFieldFromOneElem(int SrcElemInd, int* TotArrSubMatrNos, int StartOffsetSkip, int SkipLength)
{
	int AmOfRelaxElem = IntrctPtr->OutAmOfRelaxObjs();
	if(StartOffsetSkip > AmOfRelaxElem) StartOffsetSkip = AmOfRelaxElem;

	TMatrix3df **InteractMatr = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d **InteractMatr = IntrctPtr->InteractMatrix; //OC250504
	TVector3d SrcMagn = *(IntrctPtr->NewMagnArray + SrcElemInd);
	int *tTotArrSubMatrNos = TotArrSubMatrNos;

	for(int i=0; i<StartOffsetSkip; i++)
	{
		int CurElemInd = *(tTotArrSubMatrNos++);
		mArrAuxQuasiExtField[CurElemInd] += (*(InteractMatr + CurElemInd))[SrcElemInd] * SrcMagn;
	}
	int OffsetStart = StartOffsetSkip + SkipLength;
	if(OffsetStart > AmOfRelaxElem) OffsetStart = AmOfRelaxElem;
	tTotArrSubMatrNos = TotArrSubMatrNos + OffsetStart;
	for(int j=OffsetStart; j<AmOfRelaxElem; j++)
	{
		int CurElemInd = *(tTotArrSubMatrNos++);
		mArrAuxQuasiExtField[CurElemInd] += (*(InteractMatr + CurElemInd))[SrcElemInd] * SrcMagn;
	}
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_7::RelaxCurrentSubMatrix(int* pTotArrSubMatrNos, int SubMatrSize, double PrecOnMagnetizE2, int MaxIterNumForSubMatr)
{
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* ExternFieldAr = IntrctPtr->ExternFieldArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;

	radTg3dRelax* g3dRelaxPtr = nullptr;
	radTMaterial* MaterPtr = nullptr;
	TVector3d Mnew_mi_MoldVect, QuasiExtFieldAtElemStrNo;

	double NormFact = 1./double(SubMatrSize);
	double BestPrecMagnE2 =  PrecOnMagnetizE2*NormFact;
	double InstMisfitMe2 = 1.E+23;
	int IterCount;
	for(IterCount = 0; IterCount < MaxIterNumForSubMatr; IterCount++)
	{
		double BufMisfitM=0.;
		double LocPrecMagnE2 = (InstMisfitMe2 > 1.E+20)? BestPrecMagnE2 : 0.25*NormFact*InstMisfitMe2;
		if(LocPrecMagnE2 < BestPrecMagnE2) LocPrecMagnE2 = BestPrecMagnE2;

		int *tTotArrSubMatrNos = pTotArrSubMatrNos;
		for(int RelStrNo=0; RelStrNo<SubMatrSize; RelStrNo++)
		{
			int StrNo = *(tTotArrSubMatrNos++);
			QuasiExtFieldAtElemStrNo.Zero();
			TMatrix3df* MatrArrayPtr = IntrcMat[StrNo]; //OC250504
			//TMatrix3d* MatrArrayPtr = IntrcMat[StrNo]; //OC250504

			int *tColTotArrSubMatrNos = pTotArrSubMatrNos;
			for(int RelColNo=0; RelColNo<SubMatrSize; RelColNo++)
			{
				int ColNo = *(tColTotArrSubMatrNos++);
				if(RelColNo != RelStrNo) QuasiExtFieldAtElemStrNo += (MatrArrayPtr[ColNo] * MagnAr[ColNo]);
			}
			QuasiExtFieldAtElemStrNo += (ExternFieldAr[StrNo] + mArrAuxQuasiExtField[StrNo]);

			g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[StrNo];
			MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

			TVector3d *pInstantH = NewFieldAr + StrNo;
			MaterPtr->FindNewH(*pInstantH, MatrArrayPtr[StrNo], QuasiExtFieldAtElemStrNo, LocPrecMagnE2);

			TVector3d *pInstantM = MagnAr + StrNo;
			*pInstantM = MaterPtr->M(*pInstantH);

			Mnew_mi_MoldVect = *pInstantM - g3dRelaxPtr->Magn;
			BufMisfitM += Mnew_mi_MoldVect.AmpE2();
			g3dRelaxPtr->Magn = *pInstantM; 
		}
		InstMisfitMe2 = BufMisfitM/SubMatrSize;
		if(InstMisfitMe2 <= PrecOnMagnetizE2) break;
		if(radYield.Check()==0) return 0; // To allow multitasking on Mac: consider better places for this
	}
	return IterCount;
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_7::FillInSubMatrixArrays(double PrecOnMagnetiz, int*& TotArrSubMatrNos, int*& SubMatrLengths, int& AmOfSubMatr, std::vector<int>& vSubMatrLengths)
{
	if(IntrctPtr == nullptr) return 0;
	int AmOfRelaxElem = IntrctPtr->OutAmOfRelaxObjs();
	if(AmOfRelaxElem <= 0) return 0;

	int AmOfInitIter = 10; //to tune
	double ApproxPortionElemInSubMatr = 0.1; //0.03 //to tune

	radTRelaxationMethNo_4 RelaxMethNo_4(IntrctPtr);
	int ActualIterNum = RelaxMethNo_4.AutoRelax(PrecOnMagnetiz, AmOfInitIter);
	if(ActualIterNum < AmOfInitIter) return ActualIterNum;

	std::vector<radTlAuxIndNorm> vArrAuxIndNorm(AmOfRelaxElem);
	radTlAuxIndNorm *ArrAuxIndNorm = vArrAuxIndNorm.data();
	radTlAuxIndNorm *tAuxIndNorm = ArrAuxIndNorm;

	//radTg3dRelax *g3dRelaxPtr = nullptr;
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	//TMatrix3d** IntrcMat = IntrctPtr->InteractMatrix; //OC250504
	TMatrix3df *MatrArrayPtr = nullptr; //OC250504
	//TMatrix3d *MatrArrayPtr = nullptr; //OC250504

	TVector3d ContribH;
	for(int StrNo=0; StrNo<AmOfRelaxElem; StrNo++)
	{
		MatrArrayPtr = IntrcMat[StrNo];
		TVector3d *tMagnAr = IntrctPtr->NewMagnArray;

		for(int ColNo=0; ColNo<AmOfRelaxElem; ColNo++)
		{
			//g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[ColNo];
			if(ColNo != StrNo)
			{
				ContribH = (*MatrArrayPtr) * (*tMagnAr);
				radTAuxIndNorm AuxIndNorm(ColNo, ContribH.AmpE2());
				tAuxIndNorm->push_back(AuxIndNorm);
			}
			tMagnAr++;
			MatrArrayPtr++;
		}
		tAuxIndNorm->sort(radTAuxIndNorm::greater);

		//DEBUG
		//for(radTlAuxIndNorm::iterator it = tAuxIndNorm->begin(); it != tAuxIndNorm->end(); ++it)
		//{
		//	radTAuxIndNorm test = *it;
		//	int aha = 1;
		//}
		//END DBUG

		tAuxIndNorm++;
	}

	int ApproxAmOfElemInSubMatr = (int)(AmOfRelaxElem*ApproxPortionElemInSubMatr);
	if(ApproxAmOfElemInSubMatr < 1) ApproxAmOfElemInSubMatr = 1;

	AmOfSubMatr = (int)(AmOfRelaxElem/ApproxAmOfElemInSubMatr);
	if(AmOfSubMatr <= 0) AmOfSubMatr = 1;

	std::vector<radTvInt> vArrVectSubMatrNos(AmOfSubMatr);
	radTvInt* ArrVectSubMatrNos = vArrVectSubMatrNos.data();

	radTlAuxIndNorm::iterator it1, it2;
	tAuxIndNorm = ArrAuxIndNorm;
	int SubMatrixCount = 0;

	int SubMatrCount = 0;
	radTvInt *tArrVectSubMatrNos = ArrVectSubMatrNos;

	for(int i=0; i<AmOfRelaxElem; i++)
	{
		AddElemToCurrentSubMatrixIfNecessary(i, ArrVectSubMatrNos, SubMatrCount + 1);

		int LocElemCount = 0;
		for(it1 = tAuxIndNorm->begin(); it1 != tAuxIndNorm->end(); ++it1)
		{
			if(++LocElemCount > ApproxAmOfElemInSubMatr) break;
			int IndToCheck = (*it1).mInd;

			radTlAuxIndNorm *pAuxIndNormToCheck = ArrAuxIndNorm + IndToCheck;
			int LocElemCount2 = 0;
			for(it2 = pAuxIndNormToCheck->begin(); it2 != pAuxIndNormToCheck->end(); ++it2)
			{
				if(++LocElemCount2 > ApproxAmOfElemInSubMatr) break;
				if((*it2).mInd == i)
				{
					AddElemToCurrentSubMatrixIfNecessary(IndToCheck, ArrVectSubMatrNos, SubMatrCount + 1);
					break;
				}
			}
		}
		tAuxIndNorm++;

		int AmOfElemsInCurSubMatr = (int)(tArrVectSubMatrNos->size());
		if(AmOfElemsInCurSubMatr >= ApproxAmOfElemInSubMatr)
		{
			if(SubMatrCount < AmOfSubMatr)
			{
				tArrVectSubMatrNos++;
				SubMatrCount++;
			}
		}
	}
	AmOfSubMatr = SubMatrCount;
	if(!tArrVectSubMatrNos->empty()) AmOfSubMatr++;

	tAuxIndNorm = ArrAuxIndNorm;
	for(int j=0; j<AmOfRelaxElem; j++)
	{
		if(CheckIfElemIsPresentInAnySubMatr(j, ArrVectSubMatrNos, AmOfSubMatr)) continue;
		AddElemToAppropriateSubMatrix(j, tAuxIndNorm, ApproxAmOfElemInSubMatr, ArrVectSubMatrNos, AmOfSubMatr);
		tAuxIndNorm++;
	}

	vSubMatrLengths.resize(AmOfSubMatr);
	SubMatrLengths = vSubMatrLengths.data();
	CopyVectSubMatrDataToArrays(ArrVectSubMatrNos, AmOfSubMatr, TotArrSubMatrNos, SubMatrLengths);

	// RAII: automatic cleanup via vArrAuxIndNorm and vArrVectSubMatrNos
	// Lists are automatically cleaned up when vectors go out of scope
	return ActualIterNum;
}

//-------------------------------------------------------------------------

void radTRelaxationMethNo_7::FindSubMatricesToWhichElemCanBeAdded(radTlAuxIndNorm* pAuxIndNorm, int ApproxAmOfElemInSubMatr, radTvInt* ArrVectSubMatrNos, int AmOfSubMatr, radTvInt& VectIndPossibleSubMatr)
{
	if((ArrVectSubMatrNos == 0) || (AmOfSubMatr == 0) || (ApproxAmOfElemInSubMatr <= 0)) return;
	if(pAuxIndNorm == 0) return;

	radTlAuxIndNorm::iterator it;
	int ElemCount = 0;
	for(it = pAuxIndNorm->begin(); it != pAuxIndNorm->end(); ++it)
	{
		if(++ElemCount > ApproxAmOfElemInSubMatr) break;
		int IndToCheck = (*it).mInd;

		radTvInt *tArrVectSubMatrNos = ArrVectSubMatrNos;
		for(int i=0; i<AmOfSubMatr; i++)
		{
			if(CheckIfElemIsPresentInThisSubMatr(IndToCheck, tArrVectSubMatrNos))
			{
				VectIndPossibleSubMatr.push_back(i);
				break;
			}
			tArrVectSubMatrNos++;
		}
	}
}

//-------------------------------------------------------------------------

int radTRelaxationMethNo_7::FindSubMatrWithSmallestNumOfElem(radTvInt& VectIndPossibleSubMatr, radTvInt* ArrVectSubMatrNos, int AmOfSubMatr)
{
	if((ArrVectSubMatrNos == 0) || (AmOfSubMatr == 0)) return -1;

	int AmOfSubMatrToWhichElemCanBeAdded = (int)(VectIndPossibleSubMatr.size());
	if(AmOfSubMatrToWhichElemCanBeAdded == 0) return FindSubMatrWithSmallestNumOfElem(ArrVectSubMatrNos, AmOfSubMatr);

	int IndMinSize = -1, MinSize = IntrctPtr->OutAmOfRelaxObjs();
	radTvInt::iterator it;
	for(it = VectIndPossibleSubMatr.begin(); it != VectIndPossibleSubMatr.end(); ++it)
	{
		int IndSubMatr = (*it);
		if(IndSubMatr >= AmOfSubMatr) continue;
		radTvInt *tSubMatr = ArrVectSubMatrNos + IndSubMatr;
		int CurSize = (int)(tSubMatr->size());
		if(MinSize > CurSize) { MinSize = CurSize; IndMinSize = IndSubMatr;}
	}
	return IndMinSize;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//=========================================================================
// Method 9: Direct LU solver for tetrahedral elements
//=========================================================================

int radTRelaxationMethNo_9::SolveLU(std::vector<std::vector<double>>& A, std::vector<double>& b, int n)
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

int radTRelaxationMethNo_9::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
{
	if(IntrctPtr == nullptr) return 0;

	// Reset magnetization if needed
	if(!MagnResetIsNotNeeded)
	{
		IntrctPtr->ResetM();
		IntrctPtr->ResetAuxParam();
	}

	int AmOfMainElem = IntrctPtr->AmOfMainElem;
	
	// Debug: Check pointers and sizes
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

	// Build system matrix using ELF-compatible form:
	// (-N_physical + 1/chi) * M = H_ext
	// where N_physical is the physical demagnetization tensor (H = N_physical * M)
	//
	// NOW (after fix): Radia stores N_stored = N_physical directly
	// (both hexahedra and polygons use the same convention)
	// N_physical is NEGATIVE for demagnetization (e.g., -1/3 for cube, -0.29 for tetra)
	// Therefore: A = -N_stored + 1/chi = -N_physical + 1/chi (ELF form!)
	//
	// We use: A = -N_stored + diag(1/chi)
	//         RHS = H_ext
	
	std::vector<std::vector<double>> SystemMatrix(ndof, std::vector<double>(ndof, 0.0));
	std::vector<double> RHS(ndof, 0.0);

	// Build matrix: A = -N_stored + 1/chi (N_stored is now N_physical, which is negative)
	for(int i = 0; i < AmOfMainElem; i++)
	{
		radTg3dRelax* g3dRelaxPtr_i = IntrctPtr->g3dRelaxPtrVect[i];
		if(g3dRelaxPtr_i == nullptr) return 0;
		radTMaterial* MaterPtr_i = (radTMaterial*)(g3dRelaxPtr_i->MaterHandle.rep);
		if(MaterPtr_i == nullptr) return 0;

		// Get susceptibility tensor for this element
		TVector3d InstH(0., 0., 0.);
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

			// Add N_stored contributions (N_stored = N_physical, which is NEGATIVE for demagnetization)
			for(int j = 0; j < AmOfMainElem; j++)
			{
				TMatrix3df& Nij = IntrcMat[i][j];

				for(int comp_j = 0; comp_j < 3; comp_j++)
				{
					int col = 3 * j + comp_j;

					// Get N_ij component (stored value = N_physical from PreRelax)
					double Nij_val = 0.0;
					if(comp_i == 0)      Nij_val = (comp_j == 0) ? Nij.Str0.x : ((comp_j == 1) ? Nij.Str0.y : Nij.Str0.z);
					else if(comp_i == 1) Nij_val = (comp_j == 0) ? Nij.Str1.x : ((comp_j == 1) ? Nij.Str1.y : Nij.Str1.z);
					else                 Nij_val = (comp_j == 0) ? Nij.Str2.x : ((comp_j == 1) ? Nij.Str2.y : Nij.Str2.z);

					// A = -N_stored + 1/chi (negate N_stored since it's N_physical = negative)
					SystemMatrix[row][col] -= Nij_val;
				}
			}
			
			// Add 1/chi to diagonal
			SystemMatrix[row][row] += inv_chi;

			// RHS = H_ext + Mr/chi (ELF form: (-N + 1/chi)*M = H_ext + Mr/chi)
			double Hext_comp = 0.0;
			if(comp_i == 0)      Hext_comp = ExternFieldAr[i].x;
			else if(comp_i == 1) Hext_comp = ExternFieldAr[i].y;
			else                 Hext_comp = ExternFieldAr[i].z;

			double Mr_comp = 0.0;
			if(comp_i == 0)      Mr_comp = MrVect.x;
			else if(comp_i == 1) Mr_comp = MrVect.y;
			else                 Mr_comp = MrVect.z;

			// ELF form: RHS = H_ext + Mr/chi
			double Mr_over_chi = (chi_val > 1.0e-10) ? (Mr_comp / chi_val) : 0.0;
			RHS[row] = Hext_comp + Mr_over_chi;
		}
	}

	// Solve the linear system using LU decomposition
	int ierr = SolveLU(SystemMatrix, RHS, ndof);

	if(ierr != 0)
	{
		// Solver failed - singular matrix
		return 0;
	}

	// Extract solution (M values) and update arrays
	for(int i = 0; i < AmOfMainElem; i++)
	{
		MagnAr[i].x = RHS[3 * i + 0];
		MagnAr[i].y = RHS[3 * i + 1];
		MagnAr[i].z = RHS[3 * i + 2];

		// Compute H from M using material relation
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
		radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		// Update the object's magnetization
		g3dRelaxPtr->Magn = MagnAr[i];

		// H = M / chi (for linear material)
		TVector3d InstH(0., 0., 0.);
		TMatrix3d KsiTensor;
		TVector3d MrVect;
		MaterPtr->DefineInstantKsiTensor(InstH, KsiTensor, MrVect);

		// H = (M - Mr) / chi
		NewFieldAr[i].x = (KsiTensor.Str0.x > 1.e-10) ? (MagnAr[i].x - MrVect.x) / KsiTensor.Str0.x : 0.0;
		NewFieldAr[i].y = (KsiTensor.Str1.y > 1.e-10) ? (MagnAr[i].y - MrVect.y) / KsiTensor.Str1.y : 0.0;
		NewFieldAr[i].z = (KsiTensor.Str2.z > 1.e-10) ? (MagnAr[i].z - MrVect.z) / KsiTensor.Str2.z : 0.0;
	}

	// Update relaxation status
	// Set MisfitM to -1 to skip old magnetization comparison (no OldMagnArray available for direct solver)
	IntrctPtr->RelaxStatusParam.MisfitM = -1.0;
	ComputeRelaxStatusParam(MagnAr, nullptr, NewFieldAr);
	IntrctPtr->RelaxStatusParam.MisfitM = 0.0;  // Reset to 0 after (direct solver has no misfit)

	return 1;  // Single "iteration" for direct solve
}

//=========================================================================
// Method 10: BiCGSTAB with H-matrix acceleration
//=========================================================================

double radTRelaxationMethNo_10::Dot(const std::vector<double>& a, const std::vector<double>& b, int n)
{
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum) if(n > 100)
	for(int i = 0; i < n; i++)
	{
		sum += a[i] * b[i];
	}
	return sum;
}

double radTRelaxationMethNo_10::Norm2(const std::vector<double>& a, int n)
{
	return std::sqrt(Dot(a, a, n));
}

void radTRelaxationMethNo_10::Axpy(double alpha, const std::vector<double>& x, std::vector<double>& y, int n)
{
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		y[i] += alpha * x[i];
	}
}

void radTRelaxationMethNo_10::Copy(const std::vector<double>& src, std::vector<double>& dst, int n)
{
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		dst[i] = src[i];
	}
}

void radTRelaxationMethNo_10::Scale(double alpha, std::vector<double>& x, int n)
{
	#pragma omp parallel for if(n > 100)
	for(int i = 0; i < n; i++)
	{
		x[i] *= alpha;
	}
}

void radTRelaxationMethNo_10::GetDiagonalElements(std::vector<double>& diag, int n_elem)
{
	// Extract diagonal elements from interaction matrix for Jacobi preconditioner
	// Diagonal block [i][i] is a 3x3 matrix, we extract the diagonal of that
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;

	for(int i = 0; i < n_elem; i++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
		radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		TVector3d InstH(0., 0., 0.);
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

void radTRelaxationMethNo_10::DenseMatVec(const std::vector<double>& x, std::vector<double>& y, int ndof)
{
	// Computes y = A * x where A = -N + 1/chi
	// Uses H-matrix if available, otherwise dense matrix
	int n_elem = ndof / 3;
	TMatrix3df** IntrcMat = IntrctPtr->InteractMatrix;

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

			TVector3d InstH(0., 0., 0.);
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

int radTRelaxationMethNo_10::SolveBiCGSTAB(int ndof, double tol, int max_iter, double& residual)
{
	// BiCGSTAB with Jacobi preconditioner
	// Reference: van der Vorst, SIAM J. Sci. Stat. Comput. 13 (1992)

	int n_elem = ndof / 3;

	// Allocate work vectors
	std::vector<double> r(ndof), r0(ndof), p(ndof), v(ndof), s(ndof), t(ndof);
	std::vector<double> p_hat(ndof), s_hat(ndof), diag_inv(ndof);

	// Get RHS vector (H_external)
	std::vector<double> rhs(ndof);
	if(IntrctPtr->ExternFieldArray == nullptr) return 0;
	for(int i = 0; i < n_elem; i++)
	{
		rhs[3*i + 0] = IntrctPtr->ExternFieldArray[i].x;
		rhs[3*i + 1] = IntrctPtr->ExternFieldArray[i].y;
		rhs[3*i + 2] = IntrctPtr->ExternFieldArray[i].z;
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

int radTRelaxationMethNo_10::AutoRelax(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)
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

	// Solve using BiCGSTAB
	double residual = 0.0;
	int n_iter = SolveBiCGSTAB(ndof, PrecOnMagnetiz, MaxIterNumber, residual);

	// Update object magnetizations and compute H-field
	TVector3d* MagnAr = IntrctPtr->NewMagnArray;
	TVector3d* NewFieldAr = IntrctPtr->NewFieldArray;

	for(int i = 0; i < AmOfMainElem; i++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[i];
		radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

		// Update the object's magnetization
		g3dRelaxPtr->Magn = MagnAr[i];

		// H = (M - Mr) / chi
		TVector3d InstH(0., 0., 0.);
		TMatrix3d KsiTensor;
		TVector3d MrVect;
		MaterPtr->DefineInstantKsiTensor(InstH, KsiTensor, MrVect);

		NewFieldAr[i].x = (KsiTensor.Str0.x > 1.e-10) ? (MagnAr[i].x - MrVect.x) / KsiTensor.Str0.x : 0.0;
		NewFieldAr[i].y = (KsiTensor.Str1.y > 1.e-10) ? (MagnAr[i].y - MrVect.y) / KsiTensor.Str1.y : 0.0;
		NewFieldAr[i].z = (KsiTensor.Str2.z > 1.e-10) ? (MagnAr[i].z - MrVect.z) / KsiTensor.Str2.z : 0.0;
	}

	// Update relaxation status
	IntrctPtr->RelaxStatusParam.MisfitM = residual;
	ComputeRelaxStatusParam(MagnAr, nullptr, NewFieldAr);

	return n_iter;
}
