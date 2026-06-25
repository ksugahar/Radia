/*-------------------------------------------------------------------------
*
* File name:      radapl.h
*
* Project:        RADIA
*
* Description:    Wrapping RADIA application function calls
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*
-------------------------------------------------------------------------*/

#ifndef __RADAPPL_H
#define __RADAPPL_H

// Forward declaration for Python callback support
struct _object;
using PyObject = _object;

#include "rad_serialization.h"
#include "rad_geometry_base.h"
#include "rad_type_cast.h"
#include "rad_yield.h"
#include "rad_convergence.h"
#include "rad_geom_types.h"
#include "rad_auxiliary_structures.h"

#include <sstream> // Porting

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct TVector2d;

//-------------------------------------------------------------------------

class radTApplication {

	radTmhg GlobalMapOfHandlers;
	int GlobalUniqueMapKey;
	radTCast Cast;
	radTCompCriterium CompCriterium;

	int m_nProcMPI, m_rankMPI; //OC01012020

public:

	short SendingIsRequired;

	radTSend Send;
	radTConvergRepair CnRep, CnRepAux;
	short TreatRecMagsAsExtrPolygons, TreatRecMagsAsPolyhedrons, RecognizeRecMagsInPolyhedrons, TreatExtrPgnsAsPolyhedrons;
	short MemAllocForIntrctMatrTotAtOnce;

	// Nonlinear solver method: 0=mucal1 (chi-change), 1=mucal2 (B-change/Newton)
	int NonlinearMethod;

	// BiCGSTAB inner loop tolerance (default: 1e-4, ELF-compatible)
	// Can be set via Python API: rad.SolverPar("bicg_tol", value)
	double m_bicg_tol;

	// Multipole-moment method-2 linear/nonlinear acceleration controls.
	// moment_krylov_solver: 0 = BiCGSTAB, 1 = restarted GMRES.
	int m_moment_krylov_solver;
	int m_moment_gmres_restart;
	int m_moment_anderson_depth;

	// Relaxation coefficient for nonlinear iteration (default: 0.0 = full step)
	// 0.0 = full step (no under-relaxation)
	// 0.0-1.0 = under-relaxation: chi_new = chi_new*(1-relax) + chi_old*relax
	// Can be set via Python API: rad.SolverConfig(relax_param=value)
	double m_relax;

	// Keep magnetization from previous solve (default: false = reset to zero)
	// When true, rad.Solve() continues from the previous magnetization state.
	// Enables Newton->Picard hybrid:
	//   rad.SolverConfig(newton_method=True)
	//   rad.Solve(obj, 1e-3, 10, 2)   # Newton phase
	//   rad.SolverConfig(newton_method=False, keep_magnetization=True)
	//   rad.Solve(obj, 1e-3, 100, 2)  # Picard continues from Newton state
	bool m_keep_magnetization;

	// Newton-Raphson nonlinear iteration (default: false = Picard/fixed-point)
	// When true, uses differential susceptibility chi_d = (dB/dH)/mu_0 - 1 for system matrix
	// and adds correction term to RHS for quadratic convergence
	// Can be set via Python API: rad.SolverConfig(newton_method=True/False)
	bool m_use_newton;

	// Newton line search damping parameters
	// Enables adaptive backtracking line search to improve nonlinear convergence
	// Can be configured via Python API: rad.SolverConfig(newton_damping=True, ...)
	bool m_newton_damping_enabled;   // Enable line search damping (default: true when Newton active)
	int m_newton_ls_max_iter;        // Max line search backtracks (default: 5)
	double m_newton_ls_min_omega;    // Minimum omega threshold (default: 0.01)

	// B-input Newton-Raphson for energy-based hysteresis
	// When true and all materials are radTEnergyHysteresisMaterial,
	// uses Inverse(B) + analytical Jacobian instead of chi-based Picard/Newton.
	// Can be set via Python API: rad.SolverConfig(b_input_newton=True)
	bool m_b_input_newton;

	// B-input Hantila polarization method for energy-based hysteresis
	// Constant LHS (I - alpha*N), LU factored once → O(N^2) per iteration
	// Can be set via Python API: rad.SolverConfig(b_input_hantila=True)
	bool m_b_input_hantila;
	double m_hantila_alpha;   // Polarization parameter (0 = auto-compute from initial susceptibility)
	double m_hantila_relax;   // Under-relaxation (0 = full step)

	// Solve statistics (always available)
	double m_solve_t_matrix_build;   // Interaction matrix build time [s]
	double m_solve_t_moment_fieldgrad;      // Dense moment centroid field/gradient build time [s]
	double m_solve_t_moment_system_build;   // Dense moment system assembly time [s]
	double m_solve_t_lu_decomp;      // LU decomposition time [s] (Method 0 only)
	double m_solve_t_linear_solve;   // Total linear solver time [s]
	int m_solve_linear_iterations;   // Total linear iterations (BiCGSTAB only)
	int m_solve_nonl_iterations;     // Total nonlinear iterations
	int m_solve_num_threads;         // Number of threads used during solve
	int m_solve_defl_nplaq;          // Loop-deflation cycles installed (HACApK)
	double m_solve_defl_alpha;       // Loop-deflation shift alpha actually used (auto-scaled)
	bool m_solve_stats_valid;        // Whether stats are available

	// Interaction matrix cache for avoiding rebuild on repeated Solve() calls
	// The interaction matrix N only depends on geometry, not on material properties (chi)
	// Caching avoids the expensive O(N^2) matrix construction on every Solve()
	int m_cached_interact_key;       // Key of cached interaction object (0 = no cache)
	int m_cached_obj_key;            // Geometry key that the cache is valid for
	std::string m_cached_image_spec; // Image symmetry string that the cache is valid for

#ifdef RADIA_USE_HACAPK
	// HACApK parameters for H-matrix solver
	double m_hacapk_eps;       // ACA+ compression tolerance (default: 1e-4)
	int m_hacapk_leaf_size;    // Minimum cluster size (default: 32)
	double m_hacapk_eta;       // Admissibility parameter (default: 2.0)

	// HACApK statistics from last solve
	int m_hacapk_n_lowrank;
	int m_hacapk_n_dense;
	int m_hacapk_max_rank;
	int m_hacapk_n_leaves;
	int m_hacapk_n_dof;
	double m_hacapk_compression;
	double m_hacapk_build_time;
	double m_hacapk_memory_mb;        // H-matrix memory usage [MB]
	double m_hacapk_dense_memory_mb;  // Dense matrix memory [MB]
	bool m_hacapk_stats_valid;


	// Detailed timing statistics (ELF-compatible)
	double m_timing_hmatrix_build;   // H-matrix construction time
	double m_timing_linear_solve;    // Total BiCGSTAB solve time
	int m_linear_iterations;         // Total BiCGSTAB iterations (cumulative)
#endif

	radTApplication()
	{
		Initialize();
	}
	~radTApplication() {}

	void Initialize()
	{
		GlobalUniqueMapKey = 1;
		CompCriterium.BasedOnPrecLevel = 0;
		SendingIsRequired = 1;
		TreatRecMagsAsExtrPolygons = TreatExtrPgnsAsPolyhedrons = 0;
		TreatRecMagsAsPolyhedrons = 1;  // Default ON: ObjRecMag uses 6 DOF MSC hexahedron
		RecognizeRecMagsInPolyhedrons = 0; // Disable: Keep hexahedra as polyhedra for 6 DOF MSC solver
		MemAllocForIntrctMatrTotAtOnce = 0;
		NonlinearMethod = 1;  // Default: mucal2 (B-change/Newton) for faster convergence
		m_bicg_tol = 1.0e-4;  // Default: 1e-4 (ELF-compatible)
		m_moment_krylov_solver = 0;
		m_moment_gmres_restart = 40;
		m_moment_anderson_depth = 0;
		m_relax = 0.0;        // Default: 0.0 (full step, no under-relaxation)
		m_keep_magnetization = false; // Default: reset M to zero before each Solve
		m_use_newton = false; // Default: Picard iteration (backward compatible)

		// Newton line search damping init
		m_newton_damping_enabled = true;  // Default: enabled when Newton is active
		m_newton_ls_max_iter = 5;
		m_newton_ls_min_omega = 0.01;

		// B-input Newton init
		m_b_input_newton = false;

		// B-input Hantila init
		m_b_input_hantila = false;
		m_hantila_alpha = 0.0;   // 0 = auto-compute
		m_hantila_relax = 0.0;   // 0 = full step

		// Solve statistics init
		m_solve_t_matrix_build = 0.0;
		m_solve_t_moment_fieldgrad = 0.0;
		m_solve_t_moment_system_build = 0.0;
		m_solve_t_lu_decomp = 0.0;
		m_solve_t_linear_solve = 0.0;
		m_solve_linear_iterations = 0;
		m_solve_nonl_iterations = 0;
		m_solve_num_threads = 1;
		m_solve_defl_nplaq = 0;
		m_solve_defl_alpha = 0.0;
		m_solve_stats_valid = false;

		// Interaction matrix cache init
		m_cached_interact_key = 0;
		m_cached_obj_key = 0;
		m_cached_image_spec.clear();

		m_nProcMPI = 0; m_rankMPI = -1; //OC01012020

#ifdef RADIA_USE_HACAPK
		// HACApK default parameters
		m_hacapk_eps = 1.0e-4;
		m_hacapk_leaf_size = 32;
		m_hacapk_eta = 2.0;
		m_hacapk_stats_valid = false;
		m_hacapk_n_lowrank = 0;
		m_hacapk_n_dense = 0;
		m_hacapk_max_rank = 0;
		m_hacapk_n_leaves = 0;
		m_hacapk_n_dof = 0;
		m_hacapk_compression = 1.0;
		m_hacapk_build_time = 0.0;
		m_hacapk_memory_mb = 0.0;
		m_hacapk_dense_memory_mb = 0.0;
		m_timing_hmatrix_build = 0.0;
		m_timing_linear_solve = 0.0;
		m_linear_iterations = 0;
#endif

	}

	int ValidateVector3d(double* ArrayToCheck, long LenArray, TVector3d* VectorPtr);
	int ValidateVector2d(double* ArrayToCheck, long LenArray, TVector2d* VectorPtr);
	int ValidateMatrix3d(double* arToCheck, long LenAr, TMatrix3d* MatrixPtr);

	int ValidateElemKey(long ElemKey, radThg& hg);
	int ValidateFieldChar(char* FieldChar, radTFieldKey* FieldKeyPtr, bool LocSendRequired = true);
	int ValidateFieldIntChar(char* FieldIntChar, char* FinOrInfChar, radTFieldKey* FieldKeyPtr, bool LocSendRequired = true);
	// ValidateFieldEnergyForceChar REMOVED (Phase C, 2026-04-16, energy-based API gone)
	int ValidateMagnChar(char* MagnChar);
	// ValidateForceChar / ValidateTorqueChar REMOVED (Phase C, 2026-04-16, energy-based API gone)
	int ValidateIsotropMaterDescrByPoints(TVector2d* ArrayHB, int LenArrayArrayHB);

	inline int AddElementToContainer(radThg& hg);

	int SetRecMag(double* CPoi, long lenCPoi, double* Dims, long lenDims, double* Magn, long lenMagn, double* J, long lenJ, short J_IsZero);
	int SetArcCur(double* CPoi, long lenCPoi, double* Radii, long lenRadii, double* Angles, long lenAngles, double InHeight, double InJ_azim, int NumberOfSegm, char* ManOrAuto, char* Orient);
	int SetArcMag(double* CPoi, long lenCPoi, double* Radii, long lenRadii, double* Angles, long lenAngles, double InHeight, int InNumberOfSegm, double* Magn, long lenMagn, char* Orient);
	int SetCylMag(double* CPoi, long lenCPoi, double r, double h, int NumberOfSegm, double* Magn, long lenMagn, char* Orient);
	int FindSpaceTransToOrientObjAlongMainAxis(double* CPoi, char DefOrient, char Orient);
	void TransformBackMagnOrCurDensArr(int IndTr, double* Magn, long lenMagn);
	void TransformBackPointArr(int IndTr, double* arP, long lenP);

	int SetExtrudedPolygon(double* FirstPoi, long lenFirstPoi, double, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* Magn, long lenMagn, const char* OrientStr);
	inline int CheckIfExtrudedPolygonIsRecMag(TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d);
	int SetPlanarPolygon(double CoordZ, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* Magn, long lenMagn);
	
	int SetPolyhedron1(TVector3d* ArrayOfPoints, int lenArrayOfPoints, int** ArrayOfFaces, int* ArrayOfNumOfPoInFaces, int lenArrayOfFaces, double* Magn, double* arM_LinCoef=0, double* J=0, double* arJ_LinCoef=0, const char** OptionNames=0, const char** OptionValues=0, int OptionCount=0);
	int SetPolyhedron2(TVector3d** ArrayOfFaces, int* ArrayOfNumOfPoInFaces, long lenArrayOfFaces, double* Magn, long lenMagn);
	int SetArcPolygon(double* CenP, const char* OrientStr, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* Angles, int NumberOfSegm, const char* SymOrNoSymStr, double* Magn);
	int SetSolidRevolutionPolyhedron(const TVector3d& CPoiVect, const TVector3d& AxisVect, const TVector3d& AzAxVect0, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double StartAngle, int NumberOfSegm, double* Magn);

	int SetMultGenExtrPolygon(TVector2d** LayerPolygons, int* PtsNumbersInLayerPgns, double* CoordsZ, int AmOfLayerPolygons, double* Magn, long lenMagn);
	int SetMultGenExtrPolygonCur(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double* arSubdData, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, double* arMagnCompInSteps, const char** arOptionNames=0, const char** arOptionValues=0, int numOptions=0);
	
	int SetMultGenExtrRectangle(TVector3d* RectCenPoints, TVector2d* RectDims, int AmOfLayerRect, double* Magn, long lenMagn);
	int SetMultGenExtrTriangle(double* FirstPoi, long lenFirstPoi, double Lx, TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* arSubdData, double* Magn, long lenMagn, const char* OrientStr, const char** OptionNames, const char** OptionValues, int OptionCount);
	int TriangulatePolygon(TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d, double* arSubdData, char triSubdParamBorderCode, double triAngMin, double triAreaMax, const char* sTriExtOpt, TVector2d*& arTriVertPt, int& numTriVertPt, int*& arTriVertInd, int& numTri);
	int SetUpPolyhedronsFromBaseFacePolygons(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, double* arMagnCompInSteps, char frame, radThg& hgOut);
	int SetUpPolyhedronsFromBaseFacePolygonsTri(double zc, const char* strOrient, TVector2d* arPoints2d, int lenArPoints2d, TVector2d* arTriVertPt, int numTriVertPt, int* arTriVertInd, int numTri, double*** arPtrTrfParInExtrSteps, char** arStrTrfOrderInExtrSteps, int* arNumTrfInExtrSteps, int NumSteps, double avgCur, double* arMagnCompInExtrSteps, char frame, radThg& hgOut);

	int SetUpPolyhedronsFromLayerPolygons(TVector2d** LayerPolygons, int* PtsNumbersInLayerPgns, double* CoordsZ, int AmOfLayerPolygons, TVector3d& Magn, radThg& hg);
	int SetUpPolyhedronsFromLayerRectangles(TVector3d* RectCenPoints, TVector2d* RectDims, int AmOfLayerRect, TVector3d& MagnVect, radThg& hg);
	int CheckLayerPolygonStructures(TVector2d** LayerPolygons, int* PtsNumbersInLayerPgns, double* CoordsZ, int AmOfLayerPolygons);
	int CheckLayerRectangleStructures(TVector3d* RectCenPoints, TVector2d* RectDims, int AmOfLayerRect);
	int SetUpOnePolyhedronSegment(radTPtrsToPgnAndVect2d* pPtrsToPgnAndVect2d, double* z1z2, char StageChar, double* RelAbsTol, TVector3d& Magn, radTVectVect3d* pVertexPointsVect, radTVectIntPtrAndInt* pFacesVect, radThg& hgLoc);
	int CheckIfGroupIsNeeded(radThg& In_hg);
	int SetUpTetrahedronBasedOnTwoLinSegm(radTVect2dVect* pFirstVect2dVect, double z1, radTVect2dVect* pSecondVect2dVect, double z2, double* RelAbsTol, TVector3d& Magn, radThg& In_hg);
	int FindLowestPoint(radTVect2dVect* pVectP2d, TVector2d& V, double* RelAbsTol, int& LowestPointInd, char& AmOfPo);
	int FindTwoAdjacentFaces(int OneVertPoInd, int AnotherVertPoInd, radTVectIntPtrAndInt* pFacesVect, int& OneFaceInd, int& IndOfPoOnOneFace, int& AnotherFaceInd, int& IndOfPoOnAnotherFace);
	int NextCircularNumber(int CurrentNo, int Total) { return (CurrentNo == Total-1)? 0 : CurrentNo + 1;}
	int ShiftVertexPointNumbersInFaces(radTVectIntPtrAndInt* pFacesVect, int AmOfPoToDelete, char);

	int RecMagsAsExtrPolygons(char* OnOrOff);
	int RecMagsAsPolyhedrons(char* OnOrOff);
	int RecognizeRecMags(char* OnOrOff);
	int ExtPgnsAsPolyhedrons(char* OnOrOff);

	int SetGroup(int* ArrayOfKeys, long lenArrayOfKeys);
	int AddToGroup(int GroupKey, int* ArrayOfKeys, long lenArrayOfKeys);
	int OutGroupSize(int ElemKey);
	int OutGroupSubObjectKeys(int ElemKey);

	int SetRaceTrack(double* CPoi, long lenCPoi, double* Radii, long lenRadii, double* StrPartDims, long lenStrPartDims, double InHeight, double InJ_azim, int NumberOfSegm, char* ManOrAuto, char* Orient);
	int SetFlmCur(double I, TVector3d* ArrayOfPoints, int lenArrayOfPoints);
	int SetRectangle(double* CPoi, long lenCPoi, double* Dims, long lenDims);
	int SetBackgroundFieldSource(double* B, long lenB);
	int SetCoefficientFunctionFieldSource(PyObject* callback);

	int ComputeNumberOfDegOfFreedom(int ElemKey);
	int ComputeGeometricalVolume(int ElemKey);
	void ComputeMagnOrJ_InCenter(int ElemKey, char MorJ);
	int ScaleCurrent(int ElemKey, double scaleCoef);
	int SetObjMagn(int ElemKey, double Mx, double My, double Mz);
	void OutCenFieldCompRes(radTVectPairOfVect3d*);

	// FieldCompMethForSubdividedRecMag / SetLocMgnInSbdRecMag REMOVED (Phase C, 2026-04-16, radTSubdividedRecMag gone)

	int SetTranslation(double* Transl, long lenTransl);
	int SetRotation(double* PoiOnAx, long lenPoiOnAx, double* AxVect, long lenAxVect, double Angle);
	int SetPlaneSym(double* PoiOnPlane, long lenPoiOnPlane, double* PlaneNormal, long lenPlaneNormal, int s);
	int SetFieldInversion();

	int CombineTransformations(int ThisElemKey, int AnotherElemKey, char L_or_R);
	int ApplySymmetry(int g3dElemKey, int TransElemKey, int Multiplicity);

	int SetLinearMaterial(double* KsiArray, long lenKsiArray, double* RemMagnArray, long lenRemMagnArray);
	int SetLinearIsotropicMaterial(double Ksi);
	int SetLinearAnisotropicMaterial(double* KsiArray, long lenKsiArray, double* EasyAxisArray, long lenEasyAxisArray);
	int SetPermanentMagnet(double Br, double Hc, double* MagAxisArray, long lenMagAxisArray);

	// Permanent magnet materials
	int SetMagFixed(double* MagnArray, long lenMagnArray);
	int SetMagLinear(double Br, double Hc, double* MagAxisArray, long lenMagAxisArray);
	int SetMagCurve(double* pCurveData, int numPoints, double* MagAxisArray, long lenMagAxisArray);

	int SetNonlinearIsotropMaterial(double* Ms, long lenMs, double* ks, long len_ks);
	int SetNonlinearIsotropMaterial(TVector2d* ArrayHB, int LenArrayArrayHB);
	int SetNonlinearLaminatedMaterial(TVector2d* ArrayOfPoints2d, int lenArrayOfPoints2d, double PackFactor, double* dN);

	int SetNonlinearAnisotropMaterial(double** Ksi, double** Ms, double* Hc, int lenHc, char* DependenceIsNonlinear);
	int SetNonlinearAnisotropMaterial0(double* pDataPar, int lenDataPar, double* pDataPer, int lenDataPer);

	int SetEnergyHysteresisMaterial(int K, const double* chi,
	                                const std::vector<std::vector<double>>& r_tables,
	                                const std::vector<std::vector<double>>& f_tables,
	                                double eps);

	int SetPlayHysteresisMaterial(int K, const double* eta,
	                              const std::vector<std::vector<double>>& r_tables,
	                              const std::vector<std::vector<double>>& f_tables);

	int ApplyMaterial(int g3dRelaxElemKey, int MaterElemKey);
	void ComputeMvsH(int g3dRelaxOrMaterElemKey, char* MagnChar, double* H, long lenH);
	int MatHysSaveState(int MaterElemKey, double* pState, int* pLen);
	int MatHysRestoreState(int MaterElemKey, const double* pState, int Len);
	int MatHysCommitState(int MaterElemKey);
	int MatHysGetNuRev(int MaterElemKey, double* pNuRev);
	int MatHysIrreversible(int MaterElemKey, double* pB, double* pHirr);

	int PreRelax(int ElemKey, int SrcElemKey, char skipDenseMatrix=0);
	int SetRelaxSubInterval(int InteractElemKey, int StartNo, int FinNo, int RelaxTogether);
	void ShowInteractMatrix(int InteractElemKey);
	int GetInteractMatrix(int InteractElemKey, double* pMatrix, int* pDOF);
	int HMatrixDensify(int InteractElemKey, double* pMatrix, int* pDOF);  // Densify actual HACApK ACA+ operator (validation)
	int GetLoopBasis(int InteractElemKey, double* pL, int* pNLoop, int* pDOF);  // surface-charge MSC cell-graph cycle (loop) basis
	int GetFaceGeom(int InteractElemKey, double* pG, int* pDOF);  // per-DOF hex face geometry (area/centroid/normal/elem-center)
	int GetCentroidFieldGrad(int InteractElemKey, double* pC, int* pNHex, int* pDOF);  // per moment-element centroid demag field+gradient functionals
	int BuildMomentSystem(int InteractElemKey, double chi, const double* Happ, double* pA, double* pRhs, int* pDOF);  // multipole-moment MMM system matrix + rhs (Step-1 verification of the EIEM2->moment upgrade)
	int MomentSystemDenseRaw(int InteractElemKey, double chi, double* pA, int* pDOF);  // dense UN-normalized A_raw built ENTRY-BY-ENTRY via MomentSystemEntry (ACA H-matrix entry validation, Phase 2)
	int MomentHMatrixProbe(int InteractElemKey, double chi, double eps, int leaf, double eta, double* out);  // build A_raw as a HACApK H-matrix + probe H-matvec vs dense; out[8] (Phase 2 Increment 2)
	double HLUTestOnHACApK(int InteractElemKey);  // Phase 4: H-LU smoke test on real HACApK tree (returns max rel err vs MatVec round-trip)
	int HLUDebugMaterialize(int InteractElemKey, double *A_perm_out, int *lod_out, int *nd_out);  // Phase 4 debug: materialize post-convert tree
	void ShowInteractVector(int InteractElemKey, char* FieldVectID);
	int MakeManualRelax(int InteractElemKey, int MethNo, int IterNumber, double RelaxParam);
	int MakeAutoRelax(int InteractElemKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo, const char** arOptionNames=0, const char** arOptionValues=0, int numOptions=0);
	int UpdateSourcesForRelax(int InteractElemKey);
	int SolveGen(int ObjKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo, const char* image = nullptr);
	int SolveGenNonl(int ObjKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo, int NonlMethod, const char* image = nullptr);
	int BuildMatrix(int ObjKey, const char* image = nullptr);  // Build matrix without solving

	// Helper: Parse image string and apply IMA symmetry to interaction object
	// Format: "+x", "-z", "+x-z", "+x+y-z", etc.
	// Returns true on success
	bool ApplyIMASymmetryToInteraction(radTInteraction* pIntrc, const char* imageSpec, bool skipDenseIMA = false);

#ifdef RADIA_USE_HACAPK
	// HACApK parameter setting and statistics retrieval
	void SetHACApKParams(double eps, int leaf_size, double eta);
	void GetHACApKStats(double* dOut, int* nOut);
#endif

	// Solve statistics retrieval (always available)
	void GetSolveStats(double* dOut, int* nOut);

	// Point classification and batch field computation
	void ClassifyPoints(int* classification, int* nearest_elem, int n_points,
	                    double* points, int container_handle, double near_threshold);
	void ComputeFieldBatch(double* B_out, double* H_out, int n_points,
	                       double* points, int container_handle);
	void ComputeScalarPotentialBatch(double* phi_out, int n_points,
	                                 double* points, int container_handle);
	void ComputeVectorPotentialBatch(double* A_out, int n_points,
	                                 double* points, int container_handle);

	void ComputeField(int ElemKey, char* FieldChar, double* StObsPoi, long lenStObsPoi, double* FiObsPoi, long lenFiObsPoi, int Np, char* ShowArgFlag, double StrtArg);
	void ComputeField(int ElemKey, char* FieldChar, radTVectorOfVector3d& VectorOfVector3d, radTVectInputCell& VectInputCell);
	void ComputeField(int ElemKey, char* FieldChar, double** Points, long LenPoints);

	void ComputeFieldInt(int ElemKey, char* IntID, char* FieldIntChar, double* StPoi, long lenStPoi, double* FiPoi, long lenFiPoi);
	void ComputeFieldForce(int ElemKey, int ShapeElemKey);
	// ComputeFieldEnergy / ComputeFieldForceThroughEnergy / ComputeFieldTorqueThroughEnergy REMOVED (Phase C, 2026-04-16)
	// CheckForAutoDestSubdivision REMOVED (Phase C, 2026-04-16, ComputeField*ThroughEnergy gone)

	void ComputeFocusPotent(int ElemKey, double* StPoi, long lenStPoi, double* FiPoi, long lenFiPoi, int Np);
	void ComputeFocusKickPer(int ElemKey, double* P1, double* Nlong, double per, double nper, double* N1, double r1, int np1, double r2, int np2, const char* Comment, int nharm, int ns, double d1, double d2, const char* strKickUnit, double inEnergyGeV=0, const char* strOutFormat=0);
	void ComposeFocusKickPerFormStrRep(double* pKickData1, double* pKickData2, double* pBtE2Int, double* pCoordDir1, double* pCoordDir2, int np1, int np2, double per, int nper, const char* Comment);
	void ComputeFocusKick(int ElemKey, double* P1, double* Nlong, double* ArrLongDist, int lenArrLongDist, int ns, double* Ntr1, double r1, int np1, double r2, int np2, const char* StrComment, double d1, double d2);
	void ComputeShimSignature(int ElemKey, char* FldID, double* V, double* StPoi, double* FiPoi, int Np, double* Vi);

	void OutFieldCompRes(char* FieldChar, radTField* FieldArray, long Np, radTVectInputCell& VectInputCell);
	void OutFieldCompRes(char* FieldChar, radTField* FieldArray, long Np);
	void ParseAndSendOneFieldValue(radTField* tField, char* BufChar, int AmOfItem);

	// OutFieldEnergyForceCompRes REMOVED (Phase C, 2026-04-16, energy-based API gone)

	int SetCompPrecisions(const char** ValNames, double* Values, int ValCount);
	int SetCompCriterium(double InAbsPrecB, double InAbsPrecA, double InAbsPrecB_int, double InAbsPrecFrc, double InAbsPrecTrjCoord, double InAbsPrecTrjAngle);
	int SetMltplThresh(double* InMltplThresh); // Maybe to be removed later
	int SetTolForConvergence(double AbsRandMagnitude, double RelRandMagnitude, double ZeroRandMagnitude);
	int RandomizationOnOrOff(char* OnOrOff);

	// DumpElem / DumpElemParse / GenDump REMOVED (Phase B2a, 2026-04-15) -
	// .rad save/load is no longer supported.
	int RetrieveElemKey(const radTg* gPtr);

	inline double ReturnVersionID();
	void ReturnInput(double Input, int NumTimes);
	int SetMemAllocMethForIntrctMatr(char* TotOrParts);

	// Graphics3D functions removed (ApplyDrawAttrToElem_g3d, RemoveDrawAttrFromElem_g3d,
	// GraphicsForElem_g3d, GraphicsForElem_g3d_VTK, DecodeViewingOptions,
	// DeallocateAuxPgnViewData, GraphicsForAll_g3d)

	// SubdivideElement_g3d / CutElement_g3d / SubdivideElement_g3dByParPlanes REMOVED (Phase C, 2026-04-16)
	int DuplicateElement_g3d(int ElemKey, const char** OptionNames, const char** OptionValues, int OptionCount);

	int CreateFromObj_g3dWithSym(int ElemKey);
	inline int DeleteElement(int ElemKey);
	int DeleteAllElements(int DeletionMethNo);

	void ReplaceInAllGroups(radThg& OldHandle, radThg& NewHandle);
	void ReplaceInGlobalMap(radThg& OldHandle, radThg& NewHandle);

	int ProcMPI(const char* OnOrOff, double* arData=0, long* pnData=0, long* pRankFrom=0, long* pRankTo=0);

	// UnsafeGetElemByKey REMOVED (Phase B2a, 2026-04-15) - no callers

	// Get interaction pointer by key (for CplMag solver)
	// Returns nullptr if not found or not an interaction
	radTInteraction* GetInteractionByKey(int interactKey);
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

inline double radTApplication::ReturnVersionID()
{
	//double VersionID = 4.29; // Modified May 12, 2009
	//double VersionID = 4.30; // Modified June 24, 2012
	//double VersionID = 4.31; // Modified July 07, 2013
	double VersionID = 4.32; // Modified March 04, 2017

	if(SendingIsRequired) Send.Double(VersionID);
	return VersionID;
}

//-------------------------------------------------------------------------

inline int radTApplication::AddElementToContainer(radThg& hg)
{
	//CheckMemoryAvailable();
	GlobalMapOfHandlers[GlobalUniqueMapKey++] = hg;
	return GlobalUniqueMapKey - 1;
}

//-------------------------------------------------------------------------

inline int radTApplication::DeleteElement(int ElemKey)
{
	radTmhg::iterator iter = GlobalMapOfHandlers.find(ElemKey);
	if(iter == GlobalMapOfHandlers.end())
	{
		if(SendingIsRequired) Send.ErrorMessage("Radia::Error002");
		return 0;
	}
	GlobalMapOfHandlers.erase(iter);

	// Invalidate solve cache if the deleted element is the cached geometry or interaction
	if(ElemKey == m_cached_obj_key || ElemKey == m_cached_interact_key)
	{
		m_cached_interact_key = 0;
		m_cached_obj_key = 0;
	}

	if(SendingIsRequired) Send.Int(0);
	return 1;
}

//-------------------------------------------------------------------------

// CheckForAutoDestSubdivision impl REMOVED (Phase C, 2026-04-16)

//-------------------------------------------------------------------------

inline int radTApplication::CheckIfExtrudedPolygonIsRecMag(TVector2d* ArrayOfPoints2d, long lenArrayOfPoints2d)
{
	if((ArrayOfPoints2d == 0) || (lenArrayOfPoints2d != 4)) return 0;

	TVector2d &P2d_0 = ArrayOfPoints2d[0], &P2d_1 = ArrayOfPoints2d[1], &P2d_2 = ArrayOfPoints2d[2], &P2d_3 = ArrayOfPoints2d[3];
	double x0 = P2d_0.x, y0 = P2d_0.y;
	double x1 = P2d_1.x, y1 = P2d_1.y;
	double x2 = P2d_2.x, y2 = P2d_2.y;
	double x3 = P2d_3.x, y3 = P2d_3.y;
	
	double x2_mi_x0 = x2 - x0, y2_mi_y0 = y2 - y0;
	double charactLen1 = sqrt(x2_mi_x0*x2_mi_x0 + y2_mi_y0*y2_mi_y0);
	double x3_mi_x1 = x3 - x1, y3_mi_y1 = y3 - y1;
	double charactLen2 = sqrt(x3_mi_x1*x3_mi_x1 + y3_mi_y1*y3_mi_y1);
	double charactLen = (charactLen1 > charactLen2)? charactLen1 : charactLen2;
	double absTol = 2*radCR.AbsRandMagnitude(charactLen); //*2 precaution
	if(fabs(charactLen2 - charactLen1) > absTol) return 0;

	double abs_x10 = fabs(x1 - x0), abs_y10 = fabs(y1 - y0);
	double abs_x21 = fabs(x2 - x1), abs_y21 = fabs(y2 - y1);
	double abs_x32 = fabs(x3 - x2), abs_y32 = fabs(y3 - y2);
	double abs_x03 = fabs(x0 - x3), abs_y03 = fabs(y0 - y3);

	if(((abs_x10 < absTol) && (abs_y10 > absTol) &&
	    (abs_x21 > absTol) && (abs_y21 < absTol) &&
	    (abs_x32 < absTol) && (abs_y32 > absTol) &&
	    (abs_x03 > absTol) && (abs_y03 < absTol)) ||
	   ((abs_x10 > absTol) && (abs_y10 < absTol) &&
	    (abs_x21 < absTol) && (abs_y21 > absTol) &&
		(abs_x32 > absTol) && (abs_y32 < absTol) &&
		(abs_x03 < absTol) && (abs_y03 > absTol))) return 1;

	return 0;
}

//-------------------------------------------------------------------------

#endif
