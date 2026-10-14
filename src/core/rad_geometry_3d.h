/*-------------------------------------------------------------------------
*
* File name:      radg3d.h
*
* Project:        RADIA
*
* Description:    Base class for 3D objects - magnetic field sources;
*                 auxiliary classes/structures for field computation
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

//-------------------------------------------------------------------------
//	Definition of class radTg3d - a class of objects capable
//	to generate magnetic field.
//	radTg3d is a parent for radTg3dRelax, radTArcCur, ...
//-------------------------------------------------------------------------

#ifndef __RADG3D_H
#define __RADG3D_H

#include "rad_material_def.h"
#include "rad_geometry_base.h"
#include "gmvect.h"
#include "rad_math_methods.h"

#include <list>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

using radTlphg = list<radTPair_int_hg>;

//-------------------------------------------------------------------------

class radTField;

//-------------------------------------------------------------------------

// radTAuxCompDataG3D / radTSubdivOptions / radTCylindricSubdivSpec REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

class radTApplication;
class radTField;
struct radTFieldKey;
class radTrans;
class radTGroup;

//-------------------------------------------------------------------------

class radTg3d : public radTg {
public:
	radTlphg g3dListOfTransform; // Don't make it private!!!
	int IsGroupMember;
	char ConsiderOnlyWithTrans;
	TVector3d CentrPoint; //moved from derived classes OC061008

	// HandleAuxCompData / MessageChar REMOVED (Phase C, 2026-04-16, energy/subdivide gone)
	//double gCurrentScaleCoef; //required for current-carrying objects?

	radTg3d() 
	{ 
		IsGroupMember = 0; ConsiderOnlyWithTrans = 0; //gCurrentScaleCoef = 1;
	}
	~radTg3d() {}

	int Type_g() { return 1;}
	virtual int Type_g3d() { return 0;}

	inline void AddTransform(int, const radThg&);
	inline void AddTransform_OtherSide(int Multiplicity, const radThg& hg);

	inline void FindResTransfWithMultOne(radTrans&, short&);
	inline void FindInnerTransfWithMultOne(radTrans&, short&);
	inline void EraseOuterTransform();
	inline void EraseInnerTransform();
	inline void EraseAllTransformations();

	void FlattenSpaceTransforms(radTvhg&);
	int TotalMultiplicity();

	virtual void B_genComp(radTField*);
	virtual void B_comp(radTField*) {}
	void NestedFor_B(radTField*, const radTlphg::iterator&);
	inline void B_comp_Or_NestedFor(radTField*, const radTlphg::iterator&);

	virtual void B_intComp(radTField*) {}
	void B_intCompFinNum(radTField*);

	// EnergyForceTorqueComp / NestedFor_Energy / SimpleEnergyComp / MarkFurtherSubdNeed* / SetupFurtherSubdInd
	// CreateAuxCompData / CheckIfMoreEnrFrcTrqCompNeededAndUpdate / Next+ProceedNextStepEnergyForceTorqueComp
	// REMOVED (Phase C, 2026-04-16)

	void NormStressTensor(radTField*);

	void IntOverShape_Or_NestedFor(radTField* FieldPtr, const radTlphg::iterator& Iter)
	{
		if(Iter == g3dListOfTransform.end()) IntOverShape(FieldPtr);
		else NestedFor_IntOverShape(FieldPtr, Iter);
	}
	void NestedFor_IntOverShape(radTField*, const radTlphg::iterator&);
	virtual void IntOverShape(radTField*) {}

	// SubdivideItself* / CutItself / CheckAxesExchangeForSubdInLabFrame / TransferSubdivisionStructToLocalFrame
	// FindLowestAndUppestVertices REMOVED (Phase C, 2026-04-16)
	virtual int ConvertToPolyhedron(radThg&, radTApplication*, char) { return 0;} // 0 is essential

	int FinishDuplication(radTg3d* g3dPtr, radThg& hg)
	{
		radTSend Send;
		if(g3dPtr == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
		g3dPtr->IsGroupMember = 0;
		radThg hgLoc(g3dPtr); hg = hgLoc; return 1;
	}

	virtual int CreateFromSym(radThg&, radTApplication*, char);
	int NestedFor_CreateFromSym(radTGroup*, radTApplication*, char, radTrans*, const radTlphg::iterator&);
	int CreateAndAddToGroupOrNestedFor(radTGroup*, radTApplication*, char, radTrans*, const radTlphg::iterator&);

	virtual int SetMaterial(radThg&, radTApplication*) { return 0;}
	virtual void SetM(TVector3d&) {}
	virtual int ScaleCurrent(double) { return 0;} //implemented in current-carying objects

	virtual int NumberOfDegOfFreedom() { return 0;}
	virtual int ItemIsNotFullyInternalAfterCut() { return 1;}

	virtual double Volume() { return 0.;}
	double VolumeWithSym();

	virtual void VerticesInLocFrame(radTVectorOfVector3d&, bool) {}

	virtual void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans*, radTg3d*, radTApplication*) {}
	inline void GetTrfAndCenPointInLabFrame(radTrans* pBaseTrans, radTrans& bufTrans, radTrans*& pResTrans, TVector3d& cenPointInLabFr);

	inline double TransAtans(double, double, double&);
	inline double Argument(double x, double y); 

	constexpr double Abs(double x) { return (x<0.)? -x : x;}
	constexpr double Max(double x1, double x2) { return (x1<x2)? x2 : x1;}
	constexpr double Sign(double x) { return (x<0.)? -1. : 1.;}
	constexpr double Step(double x) { return (x>0.)? 1. : 0.;}

	// FindEllipticCoordOfPoint / EstimateLengthAlongEllipse REMOVED (Phase C, 2026-04-16)
	inline char AngleIsBetween(double, double, double);
	inline double AngularDifference(double, double);

	short CheckIfPosEven(int ii)
	{
		double x = 0.5*(double(ii) + 1.E-08); return ((ii > 0) && ((x - int(x)) < 0.1));
	}

	// Dump / DumpBin / DumpBinParse / DumpTransApplied REMOVED (Phase B2b/B2c, 2026-04-15)

	static int Round(double dVal) 
	{
		int iBuf = int(dVal);
		return ((dVal - iBuf) < 0.499)? iBuf : (iBuf + 1);
	}

	static bool CheckIfThreePointsAreOnOneLine(const TVector3d& P1, const TVector3d& P2, const TVector3d& P3, double absTol)
	{//to test!
		TVector3d v1 = P2 - P1, v2 = P3 - P1;
		double v1e2 = v1.x*v1.x + v1.y*v1.y + v1.z*v1.z;
		double v2e2 = v2.x*v2.x + v2.y*v2.y + v2.z*v2.z;
		double absTolE2 = absTol*absTol;
		if((v1e2 < absTolE2) || (v2e2 < absTolE2)) return true;

		double maxAbsVe2 = (v1e2 > v2e2)? v1e2 : v2e2;
		TVector3d vTestVectProd = v1^v2;
		if(vTestVectProd.Abs() < absTol*sqrt(maxAbsVe2)) return true;
		else return false;
	}

	// SetMessageChar REMOVED (Phase C, 2026-04-16, MessageChar gone)
};

//-------------------------------------------------------------------------

struct radTFieldKey { 
	short B_, H_, A_, M_, J_, Phi_, PreRelax_, Ib_, Ih_, FinInt_, Force_, ForceEnr_, Torque_, Energy_, Q_;
	radTFieldKey(short InB_ =0, short InH_ =0, short InA_ =0, short InM_ =0, short InJ_ =0, short InPhi_ =0, short InPreRelax_ =0, short InIb_ =0, short InIh_ =0, short InFinInt_ =0, short InForce_ =0, short InForceEnr_ =0, short InTorque_ =0, short InEnergy_ =0, short InQ_ =0)
	{ 
		B_=InB_; H_=InH_; A_=InA_; M_=InM_; J_=InJ_; Phi_=InPhi_; PreRelax_=InPreRelax_; Ib_=InIb_; Ih_=InIh_; FinInt_=InFinInt_; Force_=InForce_; ForceEnr_=InForceEnr_; Torque_=InTorque_; Energy_=InEnergy_; Q_=InQ_;
		if(Q_)
		{
			B_= H_= PreRelax_= 1;
		}
	}

	short NumValues() //OC02012020
	{//Number of scalar values describing the field
		short nVal = 0;
		if(B_) nVal += 3;
		if(H_) nVal += 3;
		if(A_) nVal += 3;
		if(M_) nVal += 3;
		if(J_) nVal += 3;
		if(Phi_) nVal++;
		if(Ib_) nVal += 3;
		if(Ih_) nVal += 3;
		if(Force_) nVal += 3;
		if(Torque_) nVal += 3;
		if(Energy_) nVal++; //?
		return nVal;
	}
};

//-------------------------------------------------------------------------

struct radTCompCriterium {

	short BasedOnPrecLevel; // Actually this is used nowhere at the moment
	double AbsPrecB;
	double AbsPrecA;
	double AbsPrecB_int;
	double AbsPrecForce;
	double AbsPrecTorque;
	double AbsPrecEnergy;
	double AbsPrecTrjCoord;
	double AbsPrecTrjAngle;
	double MltplThresh[4]; // Threshold ratios for 4 diff. orders of multipole approx. at field computation

	double WorstRelPrec;

	char BasedOnWorstRelPrec; // Used at energy - force computation

	radTCompCriterium() 
	{
		//Default values for all the Project:
		AbsPrecB = 0.0001; // Tesla
		AbsPrecA = 0.001;  // Tesla * mm
		AbsPrecB_int = 0.001;  // Tesla * mm
		AbsPrecForce = 1.;  // Newton
		AbsPrecTorque = 10.;  // Newton * mm
		AbsPrecEnergy = 10.;
		AbsPrecTrjCoord = AbsPrecTrjAngle = -1.;

		WorstRelPrec = 0.1;  // Used for Force computation through energy
		BasedOnWorstRelPrec = 0;

		MltplThresh[0] = 0.; // No Computatation
		MltplThresh[1] = 0.; // Dipole only
		MltplThresh[2] = 0.; // Dipole + Quadr.
		MltplThresh[3] = 0.; // Dipole + Quadr. + Next order
	}
};

//-------------------------------------------------------------------------

struct radTStructForShapeInt {
	radThg HandleOfSource, HandleOfShape;
	int IntegrandLength; // Number of elements in TVector3d* to be integrated over a Shape
	TVector3d Normal;
	TVector3d* VectArray;
	char* VectTypeArray; // 'a' - axial, 'r' - regular
	void (radTg3d::*IntegrandFunPtr)(radTField*);
	short IntOverLine_, IntOverSurf_, IntOverVol_;
	double* AbsPrecArray;

	radTStructForShapeInt(const radTStructForShapeInt&);
	radTStructForShapeInt() {}
	~radTStructForShapeInt() {}
};

//-------------------------------------------------------------------------

// radTStructForEnergyForceTorqueComp / radTHandleStructForEnergyForceTorqueComp REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

class radTField {
public:
	TVector3d P, B, H, A, M, J, Ib, Ih, NextP, Force, Torque;
	double Phi, Energy;
	int AmOfIntrctElemWithSym; // Place this into separate Relaxation dedicated structure, if more Relax specific data appear
	short PointIsInsideFrame; // This is only used with Polyhedrons (or not used at all ?)

	radTFieldKey FieldKey;
	radTCompCriterium CompCriterium;

	// HandleEnergyForceTorqueCompData REMOVED (Phase C, 2026-04-16)
	radTStructForShapeInt* ShapeIntDataPtr;

	radTField(const radTFieldKey& InFieldKey, const radTCompCriterium& InCompCriterium, 
			  const TVector3d& InP, const TVector3d& InB, const TVector3d& InH, 
			  const TVector3d& InA, const TVector3d& InM, const TVector3d& InJ, double InPhi =0.)
	{ 
		FieldKey = InFieldKey;
		P = InP; B = InB; H = InH; A = InA; M = InM; J = InJ; Phi = InPhi;
		CompCriterium = InCompCriterium;
		PointIsInsideFrame = 0;
	}
	radTField(const radTFieldKey& InFieldKey, const radTCompCriterium& InCompCriterium, 
			  const TVector3d& InP, const TVector3d& InNextP, const TVector3d& InIb, const TVector3d& InIh)
	{// This is for Field Integral
		FieldKey = InFieldKey;
		P = InP; NextP = InNextP; Ib = InIb; Ih = InIh;
		CompCriterium = InCompCriterium;
		PointIsInsideFrame = 0;
	}
	radTField(const radTFieldKey& InFieldKey, const radTCompCriterium& InCompCriterium, 
			  radTStructForShapeInt* InShapeIntDataPtr =nullptr)
	{// This is used in NestedFor_IntOverShape
		FieldKey = InFieldKey; CompCriterium = InCompCriterium; ShapeIntDataPtr = InShapeIntDataPtr;
		TVector3d ZeroVect(0.,0.,0.);
		B = H = A = M = J = Ib = Ih = Force = Torque = ZeroVect; // Add here more members should they appear
		Phi = Energy = 0.;
		PointIsInsideFrame = 0;
	}
	// radTField(const radTFieldKey&, const radTCompCriterium&, radTHandleStructForEnergyForceTorqueComp&) REMOVED (Phase C, 2026-04-16)
	radTField(const radTFieldKey& InFieldKey, const TVector3d& InP, const TVector3d& InB,
			  const TVector3d& InH, const TVector3d& InA, const TVector3d& InM, const TVector3d& InJ, double InPhi =0.)
	{ 
		FieldKey = InFieldKey;
		P = InP; B = InB; H = InH; A = InA; M = InM; J = InJ; Phi = InPhi;
		CompCriterium.BasedOnPrecLevel = 0;
		PointIsInsideFrame = 0;
	}
	radTField(const radTFieldKey& InFieldKey, const TVector3d& InP, const TVector3d& InVect, double InPhi =0.)
	{ 
		P = InP;
		FieldKey = InFieldKey;
		CompCriterium.BasedOnPrecLevel = 0;
		if(FieldKey.B_) B = InVect;
		else if(FieldKey.H_) H = InVect;
		else if(FieldKey.A_) A = InVect;
		else if(FieldKey.M_) M = InVect;
		else if(FieldKey.J_) J = InVect;
		else if(FieldKey.Phi_) Phi = InPhi;
		PointIsInsideFrame = 0;
	}
	radTField(const TVector3d& InP, const TVector3d& InB) 
	{ 
		P = InP; B = InB;
		CompCriterium.BasedOnPrecLevel = 0;
		PointIsInsideFrame = 0;
	}
	radTField() { PointIsInsideFrame = 0;}

	radTField& operator +=(const radTField& AnotherField)
	{
		if(FieldKey.B_) B+=AnotherField.B;
		if(FieldKey.H_) H+=AnotherField.H;
		if(FieldKey.A_) A+=AnotherField.A;
		if(FieldKey.M_) M+=AnotherField.M;
		if(FieldKey.J_) J+=AnotherField.J;
		if(FieldKey.Phi_) Phi+=AnotherField.Phi;
		if(FieldKey.Ib_) Ib+=AnotherField.Ib;
		if(FieldKey.Ih_) Ih+=AnotherField.Ih;
		
		if(FieldKey.Q_) //matrix for H calculation //OC191005
		{
			//B+=AnotherField.B; //assumed to be done already, because Q_ enforces B_ and H_
			//H+=AnotherField.H;
			if(!FieldKey.A_) A+=AnotherField.A;
		}
		
		return *this;
	}

	void OutVals(const radTFieldKey& key, double*& t_ar) //OC02012020
	{
		if(key.B_) { *(t_ar++) = B.x; *(t_ar++) = B.y; *(t_ar++) = B.z;}
		if(key.H_) { *(t_ar++) = H.x; *(t_ar++) = H.y; *(t_ar++) = H.z;}
		if(key.A_) { *(t_ar++) = A.x; *(t_ar++) = A.y; *(t_ar++) = A.z;}
		if(key.M_) { *(t_ar++) = M.x; *(t_ar++) = M.y; *(t_ar++) = M.z;}
		if(key.J_) { *(t_ar++) = J.x; *(t_ar++) = J.y; *(t_ar++) = J.z;}
		if(key.Phi_) { *(t_ar++) = Phi;}
		if(key.Ib_) { *(t_ar++) = Ib.x; *(t_ar++) = Ib.y; *(t_ar++) = Ib.z;}
		if(key.Ih_) { *(t_ar++) = Ih.x; *(t_ar++) = Ih.y; *(t_ar++) = Ih.z;}
		if(key.Force_) { *(t_ar++) = Force.x; *(t_ar++) = Force.y; *(t_ar++) = Force.z;}
		if(key.Torque_) { *(t_ar++) = Torque.x; *(t_ar++) = Torque.y; *(t_ar++) = Torque.z;}
		if(key.Energy_) { *(t_ar++) = Energy;}
	}
	void InVals(const radTFieldKey& key, double*& t_ar) //OC02012020
	{
		if(key.B_) { B.x = *(t_ar++); B.y = *(t_ar++); B.z = *(t_ar++);}
		if(key.H_) { H.x = *(t_ar++); H.y = *(t_ar++); H.z = *(t_ar++);}
		if(key.A_) { A.x = *(t_ar++); A.y = *(t_ar++); A.z = *(t_ar++);}
		if(key.M_) { M.x = *(t_ar++); M.y = *(t_ar++); M.z = *(t_ar++);}
		if(key.J_) { J.x = *(t_ar++); J.y = *(t_ar++); J.z = *(t_ar++);}
		if(key.Phi_) { Phi = *(t_ar++);}
		if(key.Ib_) { Ib.x = *(t_ar++); Ib.y = *(t_ar++); Ib.z = *(t_ar++);}
		if(key.Ih_) { Ih.x = *(t_ar++); Ih.y = *(t_ar++); Ih.z = *(t_ar++);}
		if(key.Force_) { Force.x = *(t_ar++); Force.y = *(t_ar++); Force.z = *(t_ar++);}
		if(key.Torque_) { Torque.x = *(t_ar++); Torque.y = *(t_ar++); Torque.z = *(t_ar++);}
		if(key.Energy_) { Energy = *(t_ar++);}
	}


	inline friend radTField operator +(const radTField&, const radTField&);
	inline friend radTField operator -(const radTField&, const radTField&);

};

//-------------------------------------------------------------------------

inline radTField operator +(const radTField& F1, const radTField& F2)
{
	radTField resF(F1);
	if(F1.FieldKey.B_ && F2.FieldKey.B_) resF.B += F2.B;
	if(F1.FieldKey.H_ && F2.FieldKey.H_) resF.H += F2.H;
	if(F1.FieldKey.A_ && F2.FieldKey.A_) resF.A += F2.A;
	if(F1.FieldKey.M_ && F2.FieldKey.M_) resF.M += F2.M;
	if(F1.FieldKey.J_ && F2.FieldKey.J_) resF.J += F2.J;
	if(F1.FieldKey.Phi_ && F2.FieldKey.Phi_) resF.Phi += F2.Phi;
	if(F1.FieldKey.Ib_ && F2.FieldKey.Ib_) resF.Ib += F2.Ib;
	if(F1.FieldKey.Ih_ && F2.FieldKey.Ih_) resF.Ih += F2.Ih;
	
	if(F1.FieldKey.Q_ && F2.FieldKey.Q_) //matrix for H calculation //OC191005
	{
		//resF.B += F2.B; //assumed to be done already, because Q_ enforces B_ and H_
		//resF.H += F2.H;
		if(!resF.FieldKey.A_) resF.A += F2.A;
	}
	return resF;
}

//-------------------------------------------------------------------------

inline radTField operator -(const radTField& F1, const radTField& F2)
{
	radTField resF(F1);
	if(F1.FieldKey.B_ && F2.FieldKey.B_) resF.B -= F2.B;
	if(F1.FieldKey.H_ && F2.FieldKey.H_) resF.H -= F2.H;
	if(F1.FieldKey.A_ && F2.FieldKey.A_) resF.A -= F2.A;
	if(F1.FieldKey.M_ && F2.FieldKey.M_) resF.M -= F2.M;
	if(F1.FieldKey.J_ && F2.FieldKey.J_) resF.J -= F2.J;
	if(F1.FieldKey.Phi_ && F2.FieldKey.Phi_) resF.Phi -= F2.Phi;
	if(F1.FieldKey.Ib_ && F2.FieldKey.Ib_) resF.Ib -= F2.Ib;
	if(F1.FieldKey.Ih_ && F2.FieldKey.Ih_) resF.Ih -= F2.Ih;
	
	if(F1.FieldKey.Q_ && F2.FieldKey.Q_) //matrix for H calculation //OC191005
	{
		//resF.B -= F2.B; //assumed to be done already, because Q_ enforces B_ and H_
		//resF.H -= F2.H;
		if(!resF.FieldKey.A_) resF.A -= F2.A;
	}
	return resF;
}

//-------------------------------------------------------------------------

class radTg3dRelax : public radTg3d {
protected:
	TMatrix3d* pM_LinCoef;

public:
	radThg MaterHandle;

	//TVector3d CentrPoint; //moved to radTg3d
	TVector3d Magn;

	float AuxFloat1, AuxFloat2, AuxFloat3;

	radTg3dRelax(const TVector3d& InCPoiVect, const TVector3d& InMagnVect, const radThg& InMaterHandle)
	{
		CentrPoint=InCPoiVect; Magn=InMagnVect; MaterHandle=InMaterHandle;
		pM_LinCoef=0;
	}
	radTg3dRelax(const TVector3d& InMagnVect, const radThg& InMaterHandle)
	{
		Magn=InMagnVect; MaterHandle=InMaterHandle;
		pM_LinCoef=0;
	}
	radTg3dRelax(const TVector3d& InMagnVect)
	{
		Magn=InMagnVect;
		pM_LinCoef=0;
	}
	radTg3dRelax(const TVector3d& InMagnVect, TMatrix3d& InM_LinCoef)
	{
		Magn=InMagnVect; 
		if(!InM_LinCoef.isZero()) pM_LinCoef = new TMatrix3d(InM_LinCoef);
		else pM_LinCoef = 0;
	}
	//radTg3dRelax(const TVector3d& InMagnVect, TMatrix3d& InM_LinCoef, const radThg& InMaterHandle)
	radTg3dRelax(const TVector3d* pInMagnVect, TMatrix3d* pInM_LinCoef, const radThg& InMaterHandle)
	{
		if(pInMagnVect != 0) Magn = *pInMagnVect;

		pM_LinCoef = 0;
		if(pInM_LinCoef != 0)
		{
			if(!pInM_LinCoef->isZero()) pM_LinCoef = new TMatrix3d(*pInM_LinCoef);
		}
		MaterHandle = InMaterHandle;
	}
	radTg3dRelax() 
	{
		MaterHandle.rep = 0;
		pM_LinCoef=0;
	}
	radTg3dRelax(radTg3dRelax& aRelax) 
	{
		*this = aRelax;
		if(aRelax.pM_LinCoef != 0) pM_LinCoef = new TMatrix3d(*(aRelax.pM_LinCoef));
	}

	~radTg3dRelax()
	{//check if this is called
		if(pM_LinCoef != 0) delete pM_LinCoef;
	}

	int Type_g3d() { return 1;}
	virtual int Type_g3dRelax() { return 0;}

	// SimpleEnergyComp REMOVED (Phase C, 2026-04-16, energy-based API gone)

	int SetMaterial(radThg& InMatHandle, radTApplication*)
	{
		MaterHandle = InMatHandle;
		radTMaterial* MaterPtr = (radTMaterial*)(MaterHandle.rep);
		return MaterPtr->FinishSetup(Magn);
	}
	void SetM(TVector3d& M) //virtual
	{
		Magn = M;
	}	

	void Push_backCenterPointAndField(radTFieldKey*, radTVectPairOfVect3d*, radTrans*, radTg3d*, radTApplication*);

	virtual TVector3d& ReturnCentrPoint() { return CentrPoint;}
	virtual radTg3dRelax* FormalIntrctMemberPtr() { return this;}

	void CheckCenPtPositionWithRespectToPlane(TVector3d* CuttingPlane, char& CenPtPositionChar)
	{
		TVector3d& PointOnPlane = *CuttingPlane;
		TVector3d& PlaneNormal = CuttingPlane[1];
		TVector3d V = CentrPoint - PointOnPlane;
		double ScalProd = PlaneNormal*V;
		if(ScalProd >= 0.) CenPtPositionChar = 'H';
		else CenPtPositionChar = 'L';
	}

	// Dump / DumpBin / DumpBinParse / DumpMaterApplied REMOVED (Phase B2b/B2c, 2026-04-15)
};

//-------------------------------------------------------------------------

inline void radTg3d::AddTransform(int Multiplicity, const radThg& hg)
{
	radTPair_int_hg aPair(Multiplicity, hg);
	g3dListOfTransform.push_front(aPair);
}

//-------------------------------------------------------------------------

inline void radTg3d::AddTransform_OtherSide(int Multiplicity, const radThg& hg)
{
	radTPair_int_hg aPair(Multiplicity, hg);
	g3dListOfTransform.push_back(aPair);
}

//-------------------------------------------------------------------------

inline void radTg3d::EraseOuterTransform()
{
	g3dListOfTransform.pop_front();
}

//-------------------------------------------------------------------------

inline void radTg3d::EraseInnerTransform()
{
	g3dListOfTransform.pop_back();
}

//-------------------------------------------------------------------------

inline void radTg3d::EraseAllTransformations()
{
	g3dListOfTransform.erase(g3dListOfTransform.begin(), g3dListOfTransform.end());
}

//-------------------------------------------------------------------------
// B_genComp moved to rad_geometry_3d.cpp (now virtual, cannot be inline)
//-------------------------------------------------------------------------

inline void radTg3d::B_comp_Or_NestedFor(radTField* FieldPtr, const radTlphg::iterator& Iter)
{
	if(Iter == g3dListOfTransform.end())
	{
		if(FieldPtr->FieldKey.Ib_ || FieldPtr->FieldKey.Ih_) B_intComp(FieldPtr);
		else B_comp(FieldPtr);
	}
	else NestedFor_B(FieldPtr, Iter);
}

//-------------------------------------------------------------------------

inline double radTg3d::TransAtans(double x, double y, double& PiMult) 
{// To optimally compute sums of atans in derived classes
	double Buf = 1.-x*y;

	if(Buf == 0.) Buf = 1.e-50; //OC040504

	PiMult = (((Buf > 0)? 0.:1.) * ((x < 0)? -1.:1.));
	return (x+y)/Buf;
}

//-------------------------------------------------------------------------
/** Calculates principal value of argument of a complex number (-Pi < Phi <= Pi)  
	@param [out] x real part
	@param [out] y imaginary part
 	@return	calculated argument value
 	@see		... */
inline double radTg3d::Argument(double x, double y)
{

	const double Pi = 3.1415926535897932;
	if(x == 0)
	{
		if(y < 0) return -0.5*Pi;
		else if(y == 0) return 0;
		else return 0.5*Pi;
	}
	if(y == 0)
	{
		if(x >= 0) return 0.;
		else return Pi;
	}
	if(y < 0)
	{
		if(x < 0) return -Pi + atan(y/x);
		else return atan(y/x); // x > 0
	}
	else // y > 0
	{
		if(x < 0) return Pi + atan(y/x);
		else return atan(y/x); // x > 0
	}
}

//-------------------------------------------------------------------------

inline double radTg3d::AngularDifference(double Phi1, double Phi2)
{// This assumes Phi1 and Phi2 between 0 and TwoPI
	if(Phi1 > Phi2) return Phi1 - Phi2;
	else return Phi1 - (Phi2 - 6.28318530717959);
}

//-------------------------------------------------------------------------

inline char radTg3d::AngleIsBetween(double Phi, double PhiSt, double PhiFi)
{// This assumes Phi, PhiSt, PhiFi between 0 and TwoPI
	if(PhiSt > PhiFi)
	{
		if((Phi > PhiSt) && (Phi < 6.28318530717959)) return 1;
		if(Phi < PhiFi) return 1;
		return 0;
	}
	else
	{
		if((Phi > PhiSt) && (Phi < PhiFi)) return 1;
		return 0;
	}
}

//-------------------------------------------------------------------------

#include "rad_transform_def.h" //to allow making subsequent inline

//-------------------------------------------------------------------------

inline void radTg3d::GetTrfAndCenPointInLabFrame(radTrans* pInBaseTrans, radTrans& bufTrans, radTrans*& pOutResTrans, TVector3d& outCenPointInLabFr)
{//assumes that pInBaseTrans and pResTrans were allocated by calling function(s)
 //in any case, it doesn't (re-)allocate pOutResTrans
	pOutResTrans = 0;
	radTrans* pTrans = (g3dListOfTransform.empty())? 0 : (radTrans*)((*(g3dListOfTransform.begin())).Handler_g.rep);
	if(pTrans != 0)
	{
		if(pInBaseTrans != 0) 
		{
			TrProduct(pInBaseTrans, pTrans, bufTrans);
			pOutResTrans = &bufTrans;
		}
		else //OC04082010
		{
			pOutResTrans = pTrans; //?
		}
	}
	else
	{
		if(pInBaseTrans != 0) pOutResTrans = pInBaseTrans;
	}
	outCenPointInLabFr = CentrPoint;
	if(pOutResTrans != 0) outCenPointInLabFr = pOutResTrans->TrPoint(outCenPointInLabFr);
}

//-------------------------------------------------------------------------

// radTAuxCompDataG3D::StoreDataFromField / PutStoredDataToField REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

#endif
