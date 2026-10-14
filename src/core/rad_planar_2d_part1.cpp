/*-------------------------------------------------------------------------
*
* File name:      radplnr1.cpp
*
* Project:        RADIA
*
* Description:    Auxiliary 2D objects
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#include "rad_planar_2d.h"
#include "rad_application.h"

//-------------------------------------------------------------------------

extern radTConvergRepair& radCR;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void radTRectangle::IntOverSurf(radTField* FieldPtr)
{
	int LenVal = FieldPtr->ShapeIntDataPtr->IntegrandLength;
	int LenValp2 = LenVal+2;
	TVector3d ZeroVect(0.,0.,0.);

	TVector3d* InnerIntegVal[6];
	TVector3d* OuterIntegVal[6];

	std::vector<std::vector<TVector3d>> vInnerIntegValStorage(6);
	std::vector<std::vector<TVector3d>> vOuterIntegValStorage(6);
	std::vector<short> vInnerElemCompNotFinished(LenVal);
	std::vector<short> vOuterElemCompNotFinished(LenVal);
	std::vector<double> vInnerAbsPrecAndLimitsArray(LenValp2);
	std::vector<double> vOuterAbsPrecAndLimitsArray(LenValp2);
	std::vector<TVector3d> vLocalVectArray(LenVal);

	int j;
	for(j=0; j<6; j++)
	{
		vInnerIntegValStorage[j].resize(LenVal);
		vOuterIntegValStorage[j].resize(LenVal);
		InnerIntegVal[j] = vInnerIntegValStorage[j].data();
		OuterIntegVal[j] = vOuterIntegValStorage[j].data();
	}
	short* InnerElemCompNotFinished = vInnerElemCompNotFinished.data();
	short* OuterElemCompNotFinished = vOuterElemCompNotFinished.data();
	double* InnerAbsPrecAndLimitsArray = vInnerAbsPrecAndLimitsArray.data();
	double* OuterAbsPrecAndLimitsArray = vOuterAbsPrecAndLimitsArray.data();
	TVector3d* LocalVectArray = vLocalVectArray.data();

		SurfIntDataPtr = new radTRectangleSurfIntData();
	//}
	//catch (radTException* radExceptionPtr)
	//{
	//	Send.ErrorMessage(radExceptionPtr->what()); return;
	//}
	//catch (...)
	//{
	//	Send.ErrorMessage("Radia::Error999"); return;
	//}

	SurfIntDataPtr->IntegrandLen = LenVal;
	SurfIntDataPtr->IntegrandFunPtr = FieldPtr->ShapeIntDataPtr->IntegrandFunPtr;
	SurfIntDataPtr->InnerAbsPrecAndLimitsArray = InnerAbsPrecAndLimitsArray;
	SurfIntDataPtr->InnerElemCompNotFinished = InnerElemCompNotFinished;
	SurfIntDataPtr->InnerIntegVal = InnerIntegVal;
	
	SurfIntDataPtr->Field = *FieldPtr;
	radTStructForShapeInt LocShapeIntData = *(FieldPtr->ShapeIntDataPtr);

	TVector3d* OutVectArray = FieldPtr->ShapeIntDataPtr->VectArray;

	LocShapeIntData.VectArray = LocalVectArray;
	SurfIntDataPtr->Field.ShapeIntDataPtr = &LocShapeIntData;

	double SmallPositive = 1.E-10;

	int i;
	for(i=0; i<LenVal; i++)
		OuterAbsPrecAndLimitsArray[i] = (FieldPtr->ShapeIntDataPtr->AbsPrecArray)[i];

	OuterAbsPrecAndLimitsArray[LenVal] = CentrPoint.y - 0.5*Dimensions.y + SmallPositive;
	OuterAbsPrecAndLimitsArray[LenVal+1] = CentrPoint.y + 0.5*Dimensions.y;

	SurfIntDataPtr->PointOnSurface.z = CentrPoint.z + SmallPositive;
	SurfIntDataPtr->Field.ShapeIntDataPtr->Normal = TVector3d(0.,0.,1.);
	FormalOneFoldInteg(this, &radTRectangle::FunForOuterIntAtSurfInt, LenVal, OuterAbsPrecAndLimitsArray, OuterElemCompNotFinished, OuterIntegVal);
	for(i=0; i<LenVal; i++) OutVectArray[i] += (OuterIntegVal[0])[i];

// Automatic cleanup via RAII
	delete SurfIntDataPtr;
}

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTPolygon::radTPolygon(TVector2d* InEdgePointsArray, int InAmOfEdgePoints)
{
	CoordZ =0.; 
	AmOfEdgePoints = InAmOfEdgePoints;

	CheckAndRearrangeEdgePoints(InEdgePointsArray, InAmOfEdgePoints); // This seems to be needed only here

	for(int i=0; i<AmOfEdgePoints; i++) EdgePointsVector.push_back(InEdgePointsArray[i]);

	radTSend Send;
	if(!CheckIfNotSelfIntersecting()) SomethingIsWrong = 1;
	else SomethingIsWrong = 0;
	if(SomethingIsWrong) return;

	IsConvex = CheckIfConvex();
	//if(!IsConvex) Send.WarningMessage("Radia::Warning011");  // Disabled: too strict for slightly non-planar hexahedral faces

	short InsideBlock=1;
	char PointsAreDifferent=1;
	//SimpleComputeCentrPoint(InsideBlock); //OC 220902
	SimpleComputeCentrPoint(InsideBlock, PointsAreDifferent); // Test!

	if(!InsideBlock) 
	{
		Send.WarningMessage("Radia::Warning010");
	}
	//if(!PointsAreDifferent) Send.WarningMessage("Radia::Warning016"); //OC 220902
}

//-------------------------------------------------------------------------

radTPolygon::radTPolygon(const radTVect2dVect& InEdgePointsVector)
{
	CoordZ =0.; 
	AmOfEdgePoints = (int)(InEdgePointsVector.size());
	for(int i=0; i<AmOfEdgePoints; i++) EdgePointsVector.push_back(InEdgePointsVector[i]);

	radTSend Send;
	IsConvex = CheckIfConvex();
	//if(!IsConvex) Send.WarningMessage("Radia::Warning011");  // Disabled: too strict for slightly non-planar hexahedral faces

	short InsideBlock=1;
	char PointsAreDifferent=1;
	//SimpleComputeCentrPoint(InsideBlock); //OC 220902
	SimpleComputeCentrPoint(InsideBlock, PointsAreDifferent); // Test!

	//if(!InsideBlock) Send.WarningMessage("Radia::Warning010");  // Disabled: too strict for slightly non-planar faces
	//if(!PointsAreDifferent) Send.WarningMessage("Radia::Warning016"); //OC 220902

	SomethingIsWrong = 0;
}

//-------------------------------------------------------------------------

radTPolygon::radTPolygon(double InCoordZ, TVector2d* InEdgePointsArray, int InAmOfEdgePoints, const TVector3d& InMagn)
{
	CoordZ = InCoordZ; Magn = InMagn;
	AmOfEdgePoints = InAmOfEdgePoints;

	CheckAndRearrangeEdgePoints(InEdgePointsArray, InAmOfEdgePoints);

	for(int i=0; i<AmOfEdgePoints; i++) EdgePointsVector.push_back(InEdgePointsArray[i]);

	radTSend Send;
	if(!CheckIfNotSelfIntersecting()) SomethingIsWrong = 1;
	else SomethingIsWrong = 0;
	if(SomethingIsWrong) return;

	IsConvex = CheckIfConvex();
	//if(!IsConvex) Send.WarningMessage("Radia::Warning011");  // Disabled: too strict for slightly non-planar hexahedral faces

	short InsideBlock=1;
	char PointsAreDifferent=1;
	//SimpleComputeCentrPoint(InsideBlock); //OC 220902
	SimpleComputeCentrPoint(InsideBlock, PointsAreDifferent); // Test!

	//if(!InsideBlock) Send.WarningMessage("Radia::Warning010");  // Disabled: too strict for slightly non-planar faces
	//if(!PointsAreDifferent) Send.WarningMessage("Radia::Warning016"); //OC 220902
}

//-------------------------------------------------------------------------

radTPolygon::radTPolygon(radTVect2dVect& InEdgePointsVector, double InZ, const TVector3d& InMagn)
{
	CoordZ = InZ; 
	AmOfEdgePoints = (int)(InEdgePointsVector.size());

	CheckAndRearrangeEdgePoints(InEdgePointsVector);

	for(int i=0; i<AmOfEdgePoints; i++) EdgePointsVector.push_back(InEdgePointsVector[i]);

	Magn = InMagn;

	radTSend Send;
	if(!CheckIfNotSelfIntersecting()) SomethingIsWrong = 1;
	else SomethingIsWrong = 0;
	if(SomethingIsWrong) return;

	IsConvex = CheckIfConvex();
	//if(!IsConvex) Send.WarningMessage("Radia::Warning011");  // Disabled: too strict for slightly non-planar hexahedral faces

	short InsideBlock=1;
	char PointsAreDifferent=1;
	//SimpleComputeCentrPoint(InsideBlock); //OC 220902
	SimpleComputeCentrPoint(InsideBlock, PointsAreDifferent); // Test!

	//if(!InsideBlock) Send.WarningMessage("Radia::Warning010");  // Disabled: too strict for slightly non-planar faces
	//if(!PointsAreDifferent) Send.WarningMessage("Radia::Warning016"); //OC 220902

	SomethingIsWrong = 0;
}

//-------------------------------------------------------------------------

// radTPolygon(CAuxBinStrVect&) REMOVED (Phase B2c, 2026-04-15)

//-------------------------------------------------------------------------


//-------------------------------------------------------------------------

void radTPolygon::IntrsctOfTwoLines(const TVector2d& V1, const TVector2d& R01, const TVector2d& R02, const TVector2d& R12, TVector2d& IntrsctPo, TLinesIntrsctCase& IntrsctCase)
{// This is used at subdivision
	const double t_Toler = 5.E-13;
	const double V_Toler = 5.E-13;

	TVector2d V2 = R12 - R02;

	double AbsV2x = fabs(V2.x), AbsV2y = fabs(V2.y);
	double MaxR = (AbsV2x > AbsV2y)? AbsV2x : AbsV2y;
	double D_Toler = MaxR*V_Toler;

	double V1yV2x = V1.y*V2.x;
	double V1xV2y = V1.x*V2.y;
	double D = V1xV2y - V1yV2x;

	if(fabs(D) > D_Toler)
	{
		IntrsctPo.x = -(-R01.y*V1.x*V2.x + R02.y*V1.x*V2.x + R01.x*V1yV2x - R02.x*V1xV2y)/D;
		IntrsctPo.y = -(R02.y*V1yV2x - R01.y*V1xV2y + R01.x*V1.y*V2.y - R02.x*V1.y*V2.y)/D;
		double t_Intrsct = (Abs(V2.x) > V_Toler)? (IntrsctPo.x - R02.x)/V2.x : (IntrsctPo.y - R02.y)/V2.y;
		IntrsctCase = ((t_Intrsct > t_Toler) && (t_Intrsct + t_Toler < 1.))? TLinesIntrsctCase::PointWithinBound : (((Abs(t_Intrsct) < t_Toler) || (Abs(t_Intrsct-1.) < t_Toler))? TLinesIntrsctCase::PointOnBoundEdge : TLinesIntrsctCase::PointOutsideBound);
	}
	else 
	{
		double LineCoinsToler = MaxR*MaxR*V_Toler;
		IntrsctCase = (fabs(V2.x*(R01.y-R02.y) - V2.y*(R01.x-R02.x)) < LineCoinsToler)? TLinesIntrsctCase::LineIsIntrsct : TLinesIntrsctCase::Zero; // To check !!!
		IntrsctPo = R02;
	}
}

//-------------------------------------------------------------------------

void radTPolygon::IntrsctOfTwoLines2(const TVector2d& R01, const TVector2d& R11, const TVector2d& R02, const TVector2d& R12, TVector2d& IntrsctPo, TLinesIntrsctCase& IntrsctCase)
{// This is used to determine self-intersection
	const double t_Toler = 5.E-12;

	TVector2d V1 = R11 - R01;
	TVector2d V2 = R12 - R02;

	double AbsV2x = fabs(V2.x), AbsV2y = fabs(V2.y);
	double V2Norm = (AbsV2x > AbsV2y)? AbsV2x : AbsV2y;
	double AbsV1x = fabs(V1.x), AbsV1y = fabs(V1.y);
	double V1Norm = (AbsV1x > AbsV1y)? AbsV1x : AbsV1y;
	double V_Toler = V1Norm*t_Toler;
	double D_Toler = V2Norm*V_Toler;

	double V1yV2x = V1.y*V2.x;
	double V1xV2y = V1.x*V2.y;
	double D = V1xV2y - V1yV2x;

	if(fabs(D) > D_Toler)
	{
		IntrsctPo.x = -(-R01.y*V1.x*V2.x + R02.y*V1.x*V2.x + R01.x*V1yV2x - R02.x*V1xV2y)/D;
		IntrsctPo.y = -(R02.y*V1yV2x - R01.y*V1xV2y + R01.x*V1.y*V2.y - R02.x*V1.y*V2.y)/D;

		double t_Intrsct2 = (Abs(V2.x) > Abs(V2.y))? (IntrsctPo.x - R02.x)/V2.x : (IntrsctPo.y - R02.y)/V2.y;
		IntrsctCase = ((t_Intrsct2 > t_Toler) && (t_Intrsct2 + t_Toler < 1.))? TLinesIntrsctCase::PointWithinBound : TLinesIntrsctCase::PointOutsideBound;
		if(IntrsctCase == TLinesIntrsctCase::PointWithinBound)
		{
			double t_Intrsct1 = (Abs(V1.x) > Abs(V1.y))? (IntrsctPo.x - R01.x)/V1.x : (IntrsctPo.y - R01.y)/V1.y;
			IntrsctCase = ((t_Intrsct1 > t_Toler) && (t_Intrsct1 + t_Toler < 1.))? TLinesIntrsctCase::PointWithinBound : TLinesIntrsctCase::PointOutsideBound;
		}

		if(IntrsctCase == TLinesIntrsctCase::PointWithinBound)
		{
			TVector2d IP_mi_R = IntrsctPo - R01;
			double AbsIP_mi_Rx = fabs(IP_mi_R.x), AbsIP_mi_Ry = fabs(IP_mi_R.y);
			double NormIP_mi_R = (AbsIP_mi_Rx > AbsIP_mi_Ry)? AbsIP_mi_Rx : AbsIP_mi_Ry;
			double AbsRx = fabs(R01.x), AbsRy = fabs(R01.x);
			double RNorm = ((AbsRx > AbsRy)? AbsRx : AbsRy);
			if(NormIP_mi_R < RNorm*t_Toler) { IntrsctCase = TLinesIntrsctCase::PointOnBoundEdge; return;}

			IP_mi_R = IntrsctPo - R11;
			AbsIP_mi_Rx = fabs(IP_mi_R.x); AbsIP_mi_Ry = fabs(IP_mi_R.y);
			NormIP_mi_R = (AbsIP_mi_Rx > AbsIP_mi_Ry)? AbsIP_mi_Rx : AbsIP_mi_Ry;
			AbsRx = fabs(R11.x); AbsRy = fabs(R11.x);
			RNorm = ((AbsRx > AbsRy)? AbsRx : AbsRy);
			if(NormIP_mi_R < RNorm*t_Toler) { IntrsctCase = TLinesIntrsctCase::PointOnBoundEdge; return;}

			IP_mi_R = IntrsctPo - R02;
			AbsIP_mi_Rx = fabs(IP_mi_R.x); AbsIP_mi_Ry = fabs(IP_mi_R.y);
			NormIP_mi_R = (AbsIP_mi_Rx > AbsIP_mi_Ry)? AbsIP_mi_Rx : AbsIP_mi_Ry;
			AbsRx = fabs(R02.x); AbsRy = fabs(R02.x);
			RNorm = ((AbsRx > AbsRy)? AbsRx : AbsRy);
			if(NormIP_mi_R < RNorm*t_Toler) { IntrsctCase = TLinesIntrsctCase::PointOnBoundEdge; return;}

			IP_mi_R = IntrsctPo - R12;
			AbsIP_mi_Rx = fabs(IP_mi_R.x), AbsIP_mi_Ry = fabs(IP_mi_R.y);
			NormIP_mi_R = (AbsIP_mi_Rx > AbsIP_mi_Ry)? AbsIP_mi_Rx : AbsIP_mi_Ry;
			AbsRx = fabs(R12.x); AbsRy = fabs(R12.x);
			RNorm = ((AbsRx > AbsRy)? AbsRx : AbsRy);
			if(NormIP_mi_R < RNorm*t_Toler) { IntrsctCase = TLinesIntrsctCase::PointOnBoundEdge; return;}
		}
	}
	else 
	{
		IntrsctPo = R02;

		double V3x = R01.x-R02.x, V3y = R01.y-R02.y;
		double AbsV3x = fabs(V3x), AbsV3y = fabs(V3y);
		double V3Norm = ((AbsV3x > AbsV3y)? AbsV3x : AbsV3y);

		double AbsR01x = fabs(R01.x), AbsR01y = fabs(R01.y);
		double R01Norm = ((AbsR01x > AbsR01y)? AbsR01x : AbsR01y);
		double AbsR02x = fabs(R02.x), AbsR02y = fabs(R02.y);
		double R02Norm = ((AbsR02x > AbsR02y)? AbsR02x : AbsR02y);
		double MaxNormR01R02 = (R01Norm > R02Norm)? R01Norm : R02Norm;

		if(V3Norm < MaxNormR01R02*t_Toler)
		{
			IntrsctCase = TLinesIntrsctCase::LineIsIntrsct; return;
		}

		double LineCoinsBufToler = V1Norm*t_Toler;
		double LineCoinsToler = LineCoinsBufToler*V3Norm;
		double CompareVal = fabs(V1.x*V3y - V1.y*V3x);
		IntrsctCase = (CompareVal < LineCoinsToler)? TLinesIntrsctCase::LineIsIntrsct : TLinesIntrsctCase::Zero; // To check !!!
	}
}

//-------------------------------------------------------------------------

// radTPolygon::FindStPointsForIntrsctLines REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolygon::FillInIntrsctInfoStruct REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolygon::SubdivideBySetOfParallelLines REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

// radTPolygon::SubdivideItself REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

void radTPolygon::ComputeCentrPoint(short& Out_InsideBlock)
{
	const double Max_k = 1.E+10;

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	if(AbsRandX == 0.) AbsRandX = 1.E-12;
	if(AbsRandY == 0.) AbsRandY = 1.E-12;

	int AmOfEdPoInBase = AmOfEdgePoints;
	int AmOfEdPoInBase_mi_1 = AmOfEdPoInBase - 1;

	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	TVector2d First2d = *BaseIter;
	
	double x1 = First2d.x;
	double y1 = First2d.y;
	double x2, y2;
	double x1e2 = x1*x1, x2e2;

	double SS=0., SXc=0., SYc=0.;

	int i;
	for(i=0; i<AmOfEdPoInBase; i++)
	{
		++BaseIter;
		if(i!=AmOfEdPoInBase_mi_1)
		{
			x2 = (*BaseIter).x; y2 = (*BaseIter).y;
		}
		else
		{
			x2 = First2d.x; y2 = First2d.y;
		}
		x2e2 = x2*x2;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);

		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			double k = (y2-y1)/(x2-x1), b = y1 - k*x1;
			double x1px2 = x1 + x2;

			SS += x2mx1*(2.*b + k*x1px2);
			SXc += -3.*b*x1e2 - 2.*k*x1e2*x1 + x2e2*(3.*b + 2.*k*x2);
			SYc += x2mx1*(3.*(b*b + b*k*x1px2) + k*k*(x1e2 + x1*x2 + x2e2));
		}
		x1 = x2; y1 = y2;
		x1e2 = x2e2;
	}
	double Square = 0.5*SS;
	double One_d_6Square = 1./Square/6.;
	CentrPoint.x = SXc*One_d_6Square;
	CentrPoint.y = SYc*One_d_6Square;

	int X_In_Count=0, X_Out_Count=0;
	int Y_In_Count=0, Y_Out_Count=0;
	short OutsideGateX=1, OutsideGateY=1;

	First2d = First2d - CentrPoint;
	x1 = First2d.x;
	y1 = First2d.y;

	// Artificial shift (to correctly determine if center point is inside the polygon):
	if(x1==0.) x1 = AbsRandX;
	if(y1==0.) y1 = AbsRandY;

	BaseIter = EdgePointsVector.begin();

	for(i=0; i<AmOfEdPoInBase; i++)
	{
		++BaseIter;
		if(i!=AmOfEdPoInBase_mi_1)
		{
			x2 = (*BaseIter).x - CentrPoint.x; 
			y2 = (*BaseIter).y - CentrPoint.y;
		}
		else
		{
			x2 = First2d.x; y2 = First2d.y;
		}

		// Artificial shift (to correctly determine if center point is inside the polygon):
		if(x2==0.) x2 = AbsRandX;
		if(y2==0.) y2 = AbsRandY;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);
		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			double k = (y2-y1)/(x2-x1), b = y1 - k*x1;
			if(x1*x2 <= 0.)
			{
				short LocInside = ((x2>x1)? (b<=0) : (b>=0));
				X_In_Count += LocInside? 1 : 0;
				X_Out_Count += (!LocInside)? 1 : 0;

				OutsideGateX = 0;
			}
			if(y1*y2 <= 0.)	OutsideGateY = 0;
		}
		else
		{
			if(y1*y2 <= 0.)
			{
				short LocInside = ((y2>y1)? (x1>=0) : (x1<=0));
				Y_In_Count += LocInside? 1 : 0;
				Y_Out_Count += (!LocInside)? 1 : 0;

				OutsideGateY = 0;
			}
		}
		x1 = x2; y1 = y2;
	}
	int X_In = X_In_Count - X_Out_Count;
	int Y_In = X_In_Count - X_Out_Count;
	short InsideBlock = (CheckIfPosEven(X_In) && CheckIfPosEven(Y_In) && (!OutsideGateX) && (!OutsideGateY));

	Out_InsideBlock = InsideBlock;
}

//-------------------------------------------------------------------------

//void radTPolygon::SimpleComputeCentrPoint(short& Out_InsideBlock)
void radTPolygon::SimpleComputeCentrPoint(short& Out_InsideBlock, char& Out_PointsAreDifferent) //OC 220902
{
	Out_PointsAreDifferent = 1; //OC 220902
	Out_InsideBlock = 1; //OC 080108

	const double Max_k = 1.E+10;

	double AbsRandX = radCR.AbsRandMagnitude(CentrPoint.x);
	double AbsRandY = radCR.AbsRandMagnitude(CentrPoint.y);
	if(AbsRandX == 0.) AbsRandX = 1.E-12;
	if(AbsRandY == 0.) AbsRandY = 1.E-12;

	int AmOfEdPoInBase = AmOfEdgePoints;
	int AmOfEdPoInBase_mi_1 = AmOfEdPoInBase - 1;

	double SumX = 0., SumY = 0.;
	double xFirstP=0, yFirstP=0, xPrevP=0, yPrevP=0; //OC 220902

	for(int kk=0; kk<AmOfEdgePoints; kk++)
	{
		TVector2d& CurrentP = EdgePointsVector[kk];

		//SumX += CurrentP.x;
		//SumY += CurrentP.y;
		double& xCurrentP = CurrentP.x; //OC 220902
		double& yCurrentP = CurrentP.y; //OC 220902
		SumX += xCurrentP;
		SumY += yCurrentP;

		if(kk==0) //OC 220902
		{ 
			xFirstP = xCurrentP; yFirstP = yCurrentP;
		}
		else 
		{ 
			if(((Abs(xCurrentP - xFirstP) < AbsRandX) && (Abs(yCurrentP - yFirstP) < AbsRandY)) ||
			   ((Abs(xCurrentP - xPrevP) < AbsRandX) && (Abs(yCurrentP - yPrevP) < AbsRandY)))
			   Out_PointsAreDifferent = 0;
		}
		xPrevP = xCurrentP; yPrevP = yCurrentP;
	}
	double InvAmOfEdgePoints = 1./AmOfEdgePoints;
	CentrPoint.x = SumX*InvAmOfEdgePoints;
	CentrPoint.y = SumY*InvAmOfEdgePoints;

	int X_In_Count=0, X_Out_Count=0;
	int Y_In_Count=0, Y_Out_Count=0;
	short OutsideGateX=1, OutsideGateY=1;

	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	TVector2d First2d = *BaseIter;

	First2d = First2d - CentrPoint;
	double x1 = First2d.x;
	double y1 = First2d.y;
	double x2, y2;

	// Artificial shift (to correctly determine if center point is inside the polygon):
	if(x1==0.) x1 = AbsRandX;
	if(y1==0.) y1 = AbsRandY;

	BaseIter = EdgePointsVector.begin();

	for(int i=0; i<AmOfEdPoInBase; i++)
	{
		++BaseIter;
		if(i!=AmOfEdPoInBase_mi_1)
		{
			x2 = (*BaseIter).x - CentrPoint.x; 
			y2 = (*BaseIter).y - CentrPoint.y;
		}
		else
		{
			x2 = First2d.x; y2 = First2d.y;
		}

		// Artificial shift (to correctly determine if center point is inside the polygon):
		if(x2==0.) x2 = AbsRandX;
		if(y2==0.) y2 = AbsRandY;

		double x2mx1 = x2-x1;
		double y2my1 = y2-y1;
		double abs_x2mx1 = Abs(x2mx1), abs_y2my1 = Abs(y2my1);
		if(abs_x2mx1*Max_k > abs_y2my1)
		{
			double k = (y2-y1)/(x2-x1), b = y1 - k*x1;
			if(x1*x2 <= 0.)
			{
				short LocInside = ((x2>x1)? (b<=0) : (b>=0));
				X_In_Count += LocInside? 1 : 0;
				X_Out_Count += (!LocInside)? 1 : 0;

				OutsideGateX = 0;
			}
			if(y1*y2 <= 0.)	OutsideGateY = 0;
		}
		else
		{
			if(y1*y2 <= 0.)
			{
				short LocInside = ((y2>y1)? (x1>=0) : (x1<=0));
				Y_In_Count += LocInside? 1 : 0;
				Y_Out_Count += (!LocInside)? 1 : 0;

				OutsideGateY = 0;
			}
		}
		x1 = x2; y1 = y2;
	}
	if(AmOfEdgePoints <= 3) //OC 080108
	{
		Out_InsideBlock = 1; //OC 080108
		return;
	}

	int X_In = X_In_Count - X_Out_Count;
	int Y_In = X_In_Count - X_Out_Count;
	short InsideBlock = (CheckIfPosEven(X_In) && CheckIfPosEven(Y_In) && (!OutsideGateX) && (!OutsideGateY));

	Out_InsideBlock = InsideBlock;
}

//-------------------------------------------------------------------------

int radTPolygon::CheckIfConvex()
{
	radTVect2dVect::iterator BaseIter = EdgePointsVector.begin();
	TVector2d p1 = *BaseIter;
	TVector2d p2 = *(++BaseIter);

	TVector2d v1 = p2 - p1;
	TVector2d v2;
	p1 = p2;

	short IsNonCovex = 0;

	int AmOfEdgePoints_m_1 = AmOfEdgePoints - 1;
	for(int i=1; i<=AmOfEdgePoints; i++)
	{
		if(i==AmOfEdgePoints_m_1) BaseIter = EdgePointsVector.begin();
		else ++BaseIter;

		p2 = *BaseIter;
		v2 = p2 - p1;

		double VectProd = v1.x*v2.y - v2.x*v1.y;
		if(VectProd < -1.E-13)
		{
			IsNonCovex = 1;
		}

		p1 = p2; v1 = v2;
	}
	return (IsNonCovex)? 0 : 1;
}

//-------------------------------------------------------------------------

int radTPolygon::RandomizeNonConvexEdgePoints(radTVect2dVect& LocEdgePointsVector)
{
	radTVect2dVect::iterator BaseIter = LocEdgePointsVector.begin();
	TVector2d p1 = *BaseIter;
	TVector2d p2 = *(++BaseIter);

	TVector2d v1 = p2 - p1;
	TVector2d v2;
	p1 = p2;

	short IsNonCovex = 0;

	int LocAmOfEdgePoints = (int)(LocEdgePointsVector.size());

	int LocAmOfEdgePoints_m_1 = LocAmOfEdgePoints - 1;
	for(int i=1; i<=LocAmOfEdgePoints; i++)
	{
		if(i==LocAmOfEdgePoints_m_1) BaseIter = LocEdgePointsVector.begin();
		else ++BaseIter;

		p2 = *BaseIter;
		v2 = p2 - p1;

		double VectProd = v1.x*v2.y - v2.x*v1.y;
		if(VectProd < -1.E-13)
		{
			IsNonCovex = 1;
		}

		p1 = p2; v1 = v2;
	}
	return (IsNonCovex)? 0 : 1;
}

//-------------------------------------------------------------------------

void radTPolygon::CheckAndRearrangeEdgePoints(TVector2d* InEdgePointsArray, int InAmOfEdgePoints)
{
	int CheckCount = 0;

	TVector2d p1 = InEdgePointsArray[0];
	TVector2d p2 = InEdgePointsArray[1];

	TVector2d v1 = p2 - p1;
	TVector2d v2;
	p1 = p2;

	int AmOfEdgePoints_p_1 = InAmOfEdgePoints + 1;
	for(int i=2; i<=AmOfEdgePoints_p_1; i++)
	{
		p2 = (i==InAmOfEdgePoints)? InEdgePointsArray[0] : ((i!=AmOfEdgePoints_p_1)? InEdgePointsArray[i] : InEdgePointsArray[1]);
		v2 = p2 - p1;

		CheckCount += (v1.x*v2.y - v2.x*v1.y < -1.E-13)? -1 : 1;

		p1 = p2; v1 = v2;
	}

	if(CheckCount<0)
	{
		std::vector<TVector2d> vBufArray(InAmOfEdgePoints);
		TVector2d* BufArray = vBufArray.data();
		*BufArray = *InEdgePointsArray;
		for(int i=1; i<InAmOfEdgePoints; i++) BufArray[i] = InEdgePointsArray[InAmOfEdgePoints-i];
		for(int k=0; k<InAmOfEdgePoints; k++) InEdgePointsArray[k] = BufArray[k];
		// RAII: automatic cleanup
	}
}

//-------------------------------------------------------------------------

void radTPolygon::CheckAndRearrangeEdgePoints(radTVect2dVect& InEdgePointsArray)
{
	int InAmOfEdgePoints = (int)(InEdgePointsArray.size());
	int CheckCount = 0;

	TVector2d p1 = InEdgePointsArray[0];
	TVector2d p2 = InEdgePointsArray[1];

	TVector2d v1 = p2 - p1;
	TVector2d v2;
	p1 = p2;

	int AmOfEdgePoints_p_1 = InAmOfEdgePoints + 1;
	for(int i=2; i<=AmOfEdgePoints_p_1; i++)
	{
		p2 = (i==InAmOfEdgePoints)? InEdgePointsArray[0] : ((i!=AmOfEdgePoints_p_1)? InEdgePointsArray[i] : InEdgePointsArray[1]);
		v2 = p2 - p1;

		CheckCount += (v1.x*v2.y - v2.x*v1.y < -1.E-13)? -1 : 1;
		p1 = p2; v1 = v2;
	}

	if(CheckCount<0) // Crasy bad style! Change sometime.
	{
		radTVect2dVect BufArray;
		BufArray.push_back(InEdgePointsArray[0]);

		for(int i=1; i<InAmOfEdgePoints; i++) BufArray.push_back(InEdgePointsArray[InAmOfEdgePoints-i]);
		for(int k=0; k<InAmOfEdgePoints; k++) InEdgePointsArray[k] = BufArray[k];
	}
}

//-------------------------------------------------------------------------

int radTPolygon::CheckIfNotSelfIntersecting()
{
	if(AmOfEdgePoints <= 3) return 1; //OC291003

	TVector2d r0Gen = EdgePointsVector[0];
	TVector2d IntrsctPo;
	int AmOfEdgePoints_m_1 = AmOfEdgePoints-1;
	TLinesIntrsctCase IntrsctCase;

	radTSend Send;
	for(int i=0; i<AmOfEdgePoints-2; i++)
	{
		int i_Aux = i+1;
		TVector2d r1Gen = EdgePointsVector[i_Aux++];

		TVector2d r0Tmp = EdgePointsVector[i_Aux];
		for(int k=i_Aux; k<AmOfEdgePoints; k++)
		{
			TVector2d r1Tmp = EdgePointsVector[(k!=AmOfEdgePoints_m_1)? k+1 : 0];
			
			IntrsctOfTwoLines2(r0Gen, r1Gen, r0Tmp, r1Tmp, IntrsctPo, IntrsctCase);
			if(IntrsctCase == TLinesIntrsctCase::PointWithinBound) 
			{ 
				Send.ErrorMessage("Radia::Error104"); return 0;
			}

			r0Tmp = r1Tmp;
		}
		r0Gen = r1Gen;
	}
	return 1;
}

//-------------------------------------------------------------------------
