/*-------------------------------------------------------------------------
*
* File name:      rad_interaction.h
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

#ifndef __RAD_INTERACTION_H
#define __RAD_INTERACTION_H

#include "rad_type_cast.h"
#include "gmvectf.h"
//#include "rad_transform_def.h"
#include "gmtrans.h"
#include "rad_geometry_3d.h"

#include <atomic>
#include <sstream>
#include <vector>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

using radTlphgPtr = list<radTPair_int_hg*>;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct radTRelaxStatusParam {
	double MisfitM, MaxModM, MaxModH;
	radTRelaxStatusParam(double InMisfitM =0., double InMaxModM =0., double InMaxModH =0.) 
	{ 
		MisfitM=InMisfitM; MaxModM=InMaxModM; MaxModH=InMaxModH;
	}
};

//-------------------------------------------------------------------------

enum class TRelaxSubIntervalID { RelaxTogether, RelaxApart };

//-------------------------------------------------------------------------

struct radTRelaxSubInterval {
	int StartNo, FinNo;
	TRelaxSubIntervalID SubIntervalID;

	radTRelaxSubInterval(int InStartNo, int InFinNo, TRelaxSubIntervalID InSubIntervalID)
	{
		StartNo = InStartNo; FinNo = InFinNo; SubIntervalID = InSubIntervalID;
	}
	radTRelaxSubInterval() {}

	inline friend int operator <(const radTRelaxSubInterval&, const radTRelaxSubInterval&);
	inline friend int operator ==(const radTRelaxSubInterval&, const radTRelaxSubInterval&);
};

//-------------------------------------------------------------------------

inline int operator <(const radTRelaxSubInterval&, const radTRelaxSubInterval&) { return 1;}

//-------------------------------------------------------------------------

inline int operator ==(const radTRelaxSubInterval& i1, const radTRelaxSubInterval& i2) 
{ 
	return (i1.StartNo == i2.StartNo) && (i1.FinNo == i2.FinNo) && (i1.SubIntervalID == i2.SubIntervalID);
}

//-------------------------------------------------------------------------

using radTVectPtrg3dRelax = vector<radTg3dRelax*>;
using radTVectPtr_g3d = vector<radTg3d*>;
using radTVectPtrTrans = vector<radTrans*>;
using radVectPtr_lphgPtr = vector<radTlphgPtr*>;
using radTVectRelaxSubInterval = vector<radTRelaxSubInterval>;

#ifdef __MWERKS__
/*
null_template
struct iterator_traits <radTg3dRelax**> {
	 typedef ptrdiff_t difference_type;
	 typedef radTg3dRelax* value_type;
	 typedef radTg3dRelax** pointer;
	 typedef radTg3dRelax*& reference;
	 typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTg3d**> {
	 typedef ptrdiff_t difference_type;
	 typedef radTg3d* value_type;
	 typedef radTg3d** pointer;
	 typedef radTg3d*& reference;
	 typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTrans**> {
	 typedef ptrdiff_t difference_type;
	 typedef radTrans* value_type;
	 typedef radTrans** pointer;
	 typedef radTrans*& reference;
	 typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTlphgPtr**> {
	 typedef ptrdiff_t difference_type;
	 typedef radTlphgPtr* value_type;
	 typedef radTlphgPtr** pointer;
	 typedef radTlphgPtr*& reference;
	 typedef random_access_iterator_tag iterator_category;
};
null_template
struct iterator_traits <radTRelaxSubInterval*> {
	 typedef ptrdiff_t difference_type;
	 typedef radTRelaxSubInterval value_type;
	 typedef radTRelaxSubInterval* pointer;
	 typedef radTRelaxSubInterval& reference;
	 typedef random_access_iterator_tag iterator_category;
};
*/
#endif

//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Global IMA Context for Field Computation
// When IMA is active, field computation needs to include contributions
// from virtual mirror elements. This context stores the IMA settings
// and is checked during B_comp() to add mirror contributions.
//-------------------------------------------------------------------------
class RadIMAFieldContext {
public:
	// Current IMA settings for field computation
	// All static members are std::atomic for thread-safety with OpenMP
	static std::atomic<bool> s_active;     // True if IMA field computation is enabled
	static std::atomic<int> s_symmetry;    // IMA symmetry flags (IMA_X, IMA_Z, etc.)
	static std::atomic<int> s_signX;       // Sign for X mirror: +1 (symmetric) or -1 (antisymmetric)
	static std::atomic<int> s_signY;       // Sign for Y mirror
	static std::atomic<int> s_signZ;       // Sign for Z mirror

	// Set IMA context before field computation (called from main thread)
	static void Set(int symmetry, int signX, int signY, int signZ) {
		s_symmetry.store(symmetry, std::memory_order_relaxed);
		s_signX.store(signX, std::memory_order_relaxed);
		s_signY.store(signY, std::memory_order_relaxed);
		s_signZ.store(signZ, std::memory_order_relaxed);
		// Release fence: all stores above are visible before s_active becomes true
		s_active.store((symmetry != 0), std::memory_order_release);
	}

	// Clear IMA context after field computation
	static void Clear() {
		s_active.store(false, std::memory_order_release);
		s_symmetry.store(0, std::memory_order_relaxed);
		s_signX.store(1, std::memory_order_relaxed);
		s_signY.store(1, std::memory_order_relaxed);
		s_signZ.store(1, std::memory_order_relaxed);
	}

	// Check if IMA is active (called from OpenMP worker threads)
	static bool IsActive() { return s_active.load(std::memory_order_acquire); }

	// Get symmetry flags (called after IsActive() returns true)
	static int GetSymmetry() { return s_symmetry.load(std::memory_order_relaxed); }
	static int GetSignX() { return s_signX.load(std::memory_order_relaxed); }
	static int GetSignY() { return s_signY.load(std::memory_order_relaxed); }
	static int GetSignZ() { return s_signZ.load(std::memory_order_relaxed); }
};

//-------------------------------------------------------------------------
// RadMomentKernelConfig: switch for the multipole-moment surface-charge kernel.
//   ON (default): each face is fan-triangulated and integrated with the analytic closed form
//                 (FieldGradFromChargedTriangleLocal) -- H = van Oosterom-Strackee, gH = its
//                 Mathematica-verified symbolic gradient (the quadrupole field-gradient).
//                 ~64x fewer kernel evals/face; EXACT for planar faces (a small flat-triangulation
//                 modeling diff only on non-planar quads).
//   OFF        : CentroidFieldGradFromFace integrates each face with a 64pt (8x8) Gauss
//                bilinear-quad quadrature for both the field H and the gradient gH; retained
//                as an explicit cross-check path via SolverConfig(moment_analytic_kernel=false).
//-------------------------------------------------------------------------
class RadMomentKernelConfig {
public:
	static std::atomic<bool> s_analyticKernel;
	static void SetAnalytic(bool on) { s_analyticKernel.store(on, std::memory_order_release); }
	static bool UseAnalytic() { return s_analyticKernel.load(std::memory_order_acquire); }
};

// Free-function accessors (defined in rad_interaction.cpp; forward-declared in radia_pybind.cpp)
void RadSetMomentAnalyticKernel(bool on);
bool RadGetMomentAnalyticKernel();

//-------------------------------------------------------------------------

class radTInteraction : public radTg {
	friend class radTHMatrixACA;    // Allow H-matrix to access interaction data
	friend class RadHACApKMMMManager;  // Allow HACApK manager to access interaction data

	// Allow unified nonlinear iteration helpers to access interaction data
	friend bool InitializeNonlinearContext(struct NonlinearContext&, radTInteraction*, bool);
	friend bool BuildBaseMatrix(struct NonlinearContext&, radTInteraction*);
	friend void StoreOldValuesAndComputeBnorm(struct NonlinearContext&, radTInteraction*);
	friend void UpdateMagnAndComputeH(struct NonlinearContext&, radTInteraction*);
	friend void ComputeActualHFieldFromSigma(struct NonlinearContext&, radTInteraction*);
	friend double UpdateChiAndCheckConvergence(struct NonlinearContext&, radTInteraction*);
	friend double ApplyLineSearchDamping(struct NonlinearContext&, radTInteraction*, const std::vector<double>&);

	int AmOfMainElem;
	int AmOfExtElem;
	radThg SourceHandle;
	radThg MoreExtSourceHandle;
	radTCompCriterium CompCriterium;
	radTRelaxStatusParam RelaxStatusParam;
	short RelaxationStarted;

	std::vector<std::vector<TMatrix3df>> vInteractMatrix;
	std::vector<TMatrix3df*> vInteractMatrixPtrs;
	std::vector<TMatrix3df> vGenMatrStorage; // Storage for MemAllocTotAtOnce mode
	TMatrix3df** InteractMatrix; //OC250504
	//TMatrix3d** InteractMatrix; //OC250504

	std::vector<TVector3d> vExternFieldArray;
	std::vector<TVector3d> vNewMagnArray;
	std::vector<TVector3d> vNewFieldArray;
	std::vector<TVector3d> vAuxOldMagnArray;
	std::vector<TVector3d> vAuxOldFieldArray;
	TVector3d* ExternFieldArray;
	TVector3d* NewMagnArray;
	TVector3d* NewFieldArray;
	TVector3d* AuxOldMagnArray;
	TVector3d* AuxOldFieldArray;

	std::vector<radTRelaxSubInterval> vRelaxSubIntervArray;
	radTRelaxSubInterval* RelaxSubIntervArray; // New 
	radTVectPtrg3dRelax g3dRelaxPtrVect;
	radTVectPtr_g3d g3dExternPtrVect;
	radTVectPtrTrans TransPtrVect;
	radVectPtr_lphgPtr IntVectOfPtrToListsOfTransPtr;
	radVectPtr_lphgPtr ExtVectOfPtrToListsOfTransPtr;
	radTVectRelaxSubInterval RelaxSubIntervConstrVect; // New
	std::vector<radTrans*> vMainTransPtrArray;
	radTrans** MainTransPtrArray;

	radTCast Cast;
	radTSend Send;
	radIdentTrans* IdentTransPtr;
	short FillInMainTransOnly;
	char mKeepTransData;

	int m_rankMPI; //21122019 (to set from Application?)
	int m_nProcMPI;

	//-------------------------------------------------------------------------
	// Variable DOF support for hybrid collocation MMMM + standard element analysis.
	//-------------------------------------------------------------------------
	int m_totalDOF;                           // Total degrees of freedom (sum of all element DOFs)
	std::vector<int> m_elemDOF;               // DOF for each element (3 for standard, 6 for MSC hexahedra)
	std::vector<int> m_elemDOFOffset;         // Starting offset in flattened arrays for each element
	bool m_hasVariableDOF;                    // True if any element has DOF != 3
	bool m_hasPEECElements;                   // True if any PEEC conductor elements are coupled

	// Variable-size interaction matrix (COLUMN-MAJOR for LAPACK compatibility)
	// Size: m_totalDOF x m_totalDOF
	// A(i,j) is stored at index [j * m_totalDOF + i] (Fortran/LAPACK convention)
	// Used when m_hasVariableDOF is true
	std::vector<double> m_flatInteractMatrix;

	// Variable-size field arrays (used when m_hasVariableDOF is true)
	std::vector<double> m_flatExternFieldArray;  // Size: m_totalDOF
	std::vector<double> m_flatMagnArray;         // Size: m_totalDOF (M or sigma values)
	std::vector<double> m_flatFieldArray;        // Size: m_totalDOF
	mutable double m_lastMomentFieldGradTime;    // Last BuildMomentSystemCore centroid field/gradient time [s]
	mutable double m_lastMomentSystemBuildTime;  // Last BuildMomentSystemCore total assembly time [s]

	//-------------------------------------------------------------------------
	// Pre-computed tetrahedron geometry for fast 3x3 block computation
	// Reference: ELF-style optimization avoiding B_comp() overhead
	//-------------------------------------------------------------------------
	bool m_tetraGeomReady;                        // True if geometry is pre-computed
	std::vector<double> m_tetraCenters;           // n_tet * 3: element centers
	std::vector<double> m_tetraFaceVertices;      // n_tet * 4 * 3 * 3: face vertices (4 faces, 3 verts, xyz)
	std::vector<double> m_tetraFaceNormals;       // n_tet * 4 * 3: outward face normals
	std::vector<double> m_tetraFaceAreas;         // n_tet * 4: face areas
	std::vector<int> m_tetraElemIndices;          // Maps tet index to element index (like m_hexaElemIndices)
	std::vector<int> m_globalToTetraIdx;          // Maps global element index to tet index (-1 if not tet)

	//-------------------------------------------------------------------------
	// Pre-computed hexahedron geometry for fast 6x6 block computation
	// Hexahedron faces are quads split into 2 triangles each
	//-------------------------------------------------------------------------
	bool m_hexaGeomReady;                         // True if geometry is pre-computed
	std::vector<double> m_hexaCenters;            // n_hex * 3: element centers
	std::vector<double> m_hexaFaceNormals;        // n_hex * 6 * 3: outward face normals
	std::vector<double> m_hexaFaceAreas;          // n_hex * 6: face areas
	std::vector<double> m_hexaTriVertices;        // n_hex * 6 * 2 * 3 * 3: 2 triangles per face, 3 verts, xyz
	std::vector<double> m_hexaTriSigns;           // n_hex * 6 * 2: sign correction for each triangle
	std::vector<int> m_hexaElemIndices;           // Maps hex index to element index
	std::vector<int> m_globalToHexIdx;            // Maps global element index to hex index (-1 if not hex)

	//-------------------------------------------------------------------------
	// Pre-computed wedge geometry for fast 5x5 block computation
	// Wedge: 5 faces (2 tri + 3 quad), each quad split into 2 triangles
	// Max triangles per wedge: 2*1 + 3*2 = 8
	//-------------------------------------------------------------------------
	bool m_wedgeGeomReady;
	std::vector<double> m_wedgeCenters;          // n_wedge * 3
	std::vector<double> m_wedgeFaceNormals;      // n_wedge * 5 * 3
	std::vector<double> m_wedgeFaceAreas;        // n_wedge * 5
	std::vector<int> m_wedgeFaceNumTris;         // n_wedge * 5: 1 for tri face, 2 for quad face
	std::vector<double> m_wedgeTriVertices;      // n_wedge * 8 * 3 * 3: up to 8 tris, 3 verts, xyz
	std::vector<double> m_wedgeTriSigns;         // n_wedge * 8: sign correction per triangle
	std::vector<int> m_wedgeTriOffset;           // n_wedge * 5: start index into TriVertices for each face
	std::vector<int> m_wedgeElemIndices;         // Maps wedge index to element index
	std::vector<int> m_globalToWedgeIdx;         // Maps global element index to wedge index (-1 if not wedge)
	static constexpr int WEDGE_MAX_TRIS = 8;    // Max triangles per wedge element

	//-------------------------------------------------------------------------
	// Pre-computed triangle local coordinate systems (for fast field computation)
	// Eliminates redundant sqrt/div operations during FieldFromChargedTriangle
	// Layout per triangle (32 doubles):
	//   basis_a[3], basis_b[3], basis_c[3]  - local coordinate system (9)
	//   origin[3]                           - triangle origin (v0) (3)
	//   XY[3][2]                            - 2D vertex coordinates (6)
	//   DS[3], AM[3], XD[3], YD[3]          - edge parameters (12)
	//   EPSG                                - geometric epsilon (1)
	//   sign                                - normal orientation sign (1)
	//-------------------------------------------------------------------------
	static constexpr int TRI_DATA_SIZE = 32;      // doubles per triangle
	static constexpr int TRIS_PER_HEX_ELEM = 12;  // 6 faces * 2 triangles per hexahedron
	std::vector<double> m_hexaTriData;            // n_hex * 12 * 32: pre-computed triangle data
	bool m_hexaTriDataReady;                      // True if triangle data is pre-computed

	//-------------------------------------------------------------------------
	// IMA (Image) Symmetry for MSC hexahedra (private member variables)
	//-------------------------------------------------------------------------
	int m_imaSymmetry;                            // IMA symmetry flags
	int m_imaSignX;                               // IMA sign for X-axis: +1 (symmetric) or -1 (antisymmetric)
	int m_imaSignY;                               // IMA sign for Y-axis: +1 (symmetric) or -1 (antisymmetric)
	int m_imaSignZ;                               // IMA sign for Z-axis: +1 (symmetric) or -1 (antisymmetric)
	bool m_imaEnabled;                            // True if IMA mode is active
	int m_imaNumElements;                         // Number of elements in IMA region (reduced set)
	std::vector<int> m_imaToFull;                 // IMA element index -> full element index
	std::vector<int> m_imaMirrorMap;              // For each IMA element, its mirror element index in full model
	std::vector<bool> m_imaUseVirtualMirror;      // True if element needs virtual mirror (quarter model)

public:

	//-------------------------------------------------------------------------
	// IMA (Image) Symmetry for MSC hexahedra
	// Reference: ELF_MAGIC IMA approach with image summation during matrix assembly
	// Enables half-model (x-mirror), quarter-model (xy), eighth-model (xyz) analysis
	//-------------------------------------------------------------------------
	enum IMASymmetryFlags {
		IMA_NONE = 0,
		IMA_X = 1,    // Mirror in X=0 plane (YZ plane)
		IMA_Y = 2,    // Mirror in Y=0 plane (XZ plane)
		IMA_Z = 4,    // Mirror in Z=0 plane (XY plane)
		IMA_XY = 3,   // Quarter model (X and Y mirrors)
		IMA_XZ = 5,   // Quarter model (X and Z mirrors)
		IMA_YZ = 6,   // Quarter model (Y and Z mirrors)
		IMA_XYZ = 7   // Eighth model (all three mirrors)
	};

	// IMA (Image Method of Analysis) implementation notes:
	// Netgen face ordering: DOFs [0, 1, 2, 3, 4, 5] = z-, x+, y-, x-, y+, z+
	// (Matches NETGEN_FACES in radia_pybind.cpp)
	//
	// IMA mirror contributions are computed inline in the active kernels:
	//   multipole-moment MMM centroid field/grad: scalar IMA sign (MSC surface charge is scalar)
	//   Compute3x3BlockFast: component sign matrix S[beta] (MMM magnetization is pseudovector)
	// No DOF permutation matrices are needed.

	int AmOfRelaxSubInterv;

	short SomethingIsWrong;
	short MemAllocTotAtOnce;

	radTInteraction(const radThg&, const radThg&, const radTCompCriterium&, short =0, char =0, char =0, int =-1, int =0, char =0); //OC08012020 + skipDenseMatrix
	//radTInteraction(const radThg&, const radThg&, const radTCompCriterium&, short =0, char =0, char =0);
	radTInteraction();
	~radTInteraction();

	int Setup(const radThg& In_hg, const radThg& In_hgMoreExtSrc, const radTCompCriterium& InCompCriterium, short InMemAllocTotAtOnce, char AuxOldMagnArrayIsNeeded, char KeepTransData, int rankMPI=-1, int nProcMPI=0, char skipDenseMatrix=0); //OC08012020 + skipDenseMatrix
	//int Setup(const radThg& In_hg, const radThg& In_hgMoreExtSrc, const radTCompCriterium& InCompCriterium, short InMemAllocTotAtOnce, char AuxOldMagnArrayIsNeeded, char KeepTransData);

	void CountMainRelaxElems(radTg3d*, radTlphgPtr*);
	void AllocateMemory(char ExtraExternFieldArrayIsNeeded, char skipInteractMatrix = 0);
	void DeallocateMemory(); //OC27122019

	int SetupInteractMatrix(); //OC26122019
	//void SetupInteractMatrix();
	void AllocateInteractMatrix();  // Allocate InteractMatrix only (for HACApK 3DOF case)

	//-------------------------------------------------------------------------
	// Variable DOF methods for hybrid MSC + standard element analysis
	//-------------------------------------------------------------------------
	void ComputeDOFOffsets();  // Compute DOF offsets for each element
	int SetupInteractMatrix_VariableDOF();  // Build interaction matrix with variable DOF blocks
	void SetupVariableDOFArrays();  // Set up flat arrays only (without building interaction matrix)

	// DOF accessors (inline for performance)
	int GetTotalDOF() const { return m_totalDOF; }
	int GetElementDOF(int elemIdx) const { return m_elemDOF[elemIdx]; }
	int GetElementDOFOffset(int elemIdx) const { return m_elemDOFOffset[elemIdx]; }
	bool HasVariableDOF() const { return m_hasVariableDOF; }
	bool HasPEECElements() const { return m_hasPEECElements; }
	void SetHasPEECElements(bool flag) { m_hasPEECElements = flag; }

	// Element accessors
	int GetAmOfMainElem() const { return AmOfMainElem; }
	radTg3dRelax* GetElement(int idx) { return g3dRelaxPtrVect[idx]; }

	// True if any RELAXABLE element is a surface-charge (Use6DOF_MSC) polyhedron -- hex (6 faces) or
	// 5-face wedge/pyramid.  These are solved by the canonical multipole-moment MMM path. Defined in
	// rad_interaction.cpp (needs the full radTPolyhedron type). Tet (MMM) / PM elements -> false.
	bool HasSurfaceChargeElements() const;

	// Triangle precomputation accessors (for HACApK sharing)
	const double* GetHexaTriData() const { return m_hexaTriData.data(); }
	bool IsTetraGeomReady() const { return m_tetraGeomReady; }
	bool IsHexaTriDataReady() const { return m_hexaTriDataReady; }
	int GetNumHexElements() const { return (int)m_hexaElemIndices.size(); }
	const std::vector<int>& GetHexaElemIndices() const { return m_hexaElemIndices; }

	// Triangle vertex/sign accessors
	const std::vector<double>& GetHexaTriVertices() const { return m_hexaTriVertices; }
	const std::vector<double>& GetHexaTriSigns() const { return m_hexaTriSigns; }
	const std::vector<double>& GetHexaCenters() const { return m_hexaCenters; }
	const std::vector<double>& GetHexaFaceAreas() const { return m_hexaFaceAreas; }
	const std::vector<double>& GetHexaFaceNormals() const { return m_hexaFaceNormals; }

	// Flat matrix/array accessors for solvers
	double* GetFlatInteractMatrix() { return m_flatInteractMatrix.data(); }
	const double* GetFlatInteractMatrix() const { return m_flatInteractMatrix.data(); }
	double* GetFlatMagnArray() { return m_flatMagnArray.data(); }
	double* GetFlatFieldArray() { return m_flatFieldArray.data(); }
	double* GetFlatExternFieldArray() { return m_flatExternFieldArray.data(); }

	// Element access for coupled solver (PEEC-MMM)
	int GetNumElements() const { return AmOfMainElem; }
	TVector3d GetElementCenter(int elemIdx) const {
		if (elemIdx >= 0 && elemIdx < AmOfMainElem) {
			return g3dRelaxPtrVect[elemIdx]->CentrPoint;
		}
		return TVector3d(0, 0, 0);
	}
	TVector3d GetElementMagnetization(int elemIdx) const {
		if (elemIdx >= 0 && elemIdx < AmOfMainElem) {
			return g3dRelaxPtrVect[elemIdx]->Magn;
		}
		return TVector3d(0, 0, 0);
	}
	double GetElementVolume(int elemIdx) const {
		if (elemIdx >= 0 && elemIdx < AmOfMainElem) {
			return g3dRelaxPtrVect[elemIdx]->Volume();
		}
		return 0.0;
	}

	// Helper: Get pointer to interaction block (row_elem, col_elem)
	// Matrix is stored COLUMN-MAJOR for LAPACK compatibility
	double* GetInteractBlock(int row_elem, int col_elem);
	const double* GetInteractBlock(int row_elem, int col_elem) const;

	//-------------------------------------------------------------------------
	// Nonlinear iteration helpers (unified across LU, BiCGSTAB, HACApK solvers)
	// Reference: ELF mucal0/mucal1/mucal2 pattern (Modern Fortran implementation)
	//-------------------------------------------------------------------------

	// Initialize chi array for all elements using ELF mucal0 style
	// Returns: true if all materials are linear (single solve suffices)
	// polyCache: output vector of cached radTPolyhedron* for 6DOF elements
	// chiArray: output chi value for each element
	bool InitializeChiArray(std::vector<double>& chiArray,
	                        std::vector<radTPolyhedron*>& polyCache);

	// Initialize H field for first iteration (ELF mucal0 style)
	// Uses H_init_mag as initial field magnitude
	void InitializeHField(double H_init_mag = 100.0);

	// Store old B-field norms for convergence check (ELF mucal2 style)
	// Also stores old magnetization values in OldMagn
	// chiArray: current chi for each element (from InitializeChiArray or UpdateChiAndCheckConvergence)
	void StoreOldBnorm(std::vector<double>& OldBnorm,
	                   const std::vector<double>& chiArray,
	                   const std::vector<radTPolyhedron*>& polyCache);

	// Update element magnetization from FlatMagn (after linear solve)
	// Copies FlatMagn to element Magn fields (3DOF) or poly->FaceSigma (6DOF)
	// Also updates FlatField with H = M / chi for chi update loop
	void UpdateElementMagnetization(const std::vector<radTPolyhedron*>& polyCache,
	                                const std::vector<double>& chiArray);

	// Update chi from H and check B-field convergence (ELF mucal1+mucal2 combined)
	// chiArray: input/output chi values (updated with new chi)
	// OldBnorm: input old B-field norms (from StoreOldBnorm)
	// relax: under-relaxation parameter (0.0 = full step, 0.0-1.0 damping)
	// max_B_rel_change: output max |B_new - B_old| / B_sat across all elements
	// Returns: true if converged (max_B_rel_change < tolerance)
	void UpdateChiAndCheckConvergence(std::vector<double>& chiArray,
	                                   const std::vector<double>& OldBnorm,
	                                   const std::vector<radTPolyhedron*>& polyCache,
	                                   double relax,
	                                   double& max_B_rel_change);

	//-------------------------------------------------------------------------
	// Fast tetrahedron matrix build (avoiding B_comp overhead)
	//-------------------------------------------------------------------------
	void PrecomputeTetraGeometry();  // Pre-compute face vertices/normals/areas
	void Compute3x3BlockFast(int elem_i, int elem_j, double* N_mat) const;  // Fast 3x3 block

	//-------------------------------------------------------------------------
	// Fast hexahedron matrix build (avoiding FieldFromQuadFace overhead)
	//-------------------------------------------------------------------------
	void PrecomputeHexaGeometry();  // Pre-compute face triangles/normals/eval points
	void PrecomputeHexaTriangleData();  // Pre-compute triangle local coordinate systems
	// Per-DOF hex face geometry in the matrix DOF order (for div(B)=0 / RHS / moment studies in Python).
	// Gflat is ROW-MAJOR (m_totalDOF x 11): [elem_local, area, cx,cy,cz, nx,ny,nz(outward), ecx,ecy,ecz].
	// Non-hex DOFs (tet/wedge/pyramid) get elem_local=-1 and zeros.  See .cpp.
	void BuildFaceGeom(std::vector<double>& Gflat) const;
	// Per moment element demag field H and gradient gradH at the element CENTROID, as linear functionals of each
	// source DOF charge -- the analytic self-term kernel of the parameter-free moment formulation.
	// SELF face: bare charged-face field (centroid is interior -> finite, no center charge needed).
	// MUTUAL face: collocation MMMM dipole layer = bare face - area*(point charge @ source center) (finite distance).
	// Cflat is ROW-MAJOR (nMom x 9 x m_totalDOF): comp k (Hx,Hy,Hz, gxx,gyy,gzz,gxy,gxz,gyz), source DOF g
	// -> Cflat[(h*9+k)*m_totalDOF + g].  H = (1/4pi) int sigma (r-r')/|r-r'|^3 dA'.  See .cpp.
	void BuildCentroidFieldGrad(std::vector<double>& Cflat, int& nHexOut) const;
	// e-ordered list of MOMENT elements (hex 6-DOF + wedge/pyramid 5-DOF surface-charge polyhedra); for a pure-hex
	// model this equals m_hexaElemIndices order.  Single source of moment element ordering (dense path +
	// the SolveLinearStep moment branch).  Generalizes the moment formula from hex-only to all 5/6-face MSC elements.
	void CollectMomentElems(std::vector<int>& out) const;
	// Field H[3] + grad gH[6] (xx,yy,zz,xy,xz,yz) at a target centroid ce3 from ONE unit-charge face (V4
	// corners, srcCenter = the face's element center, area), INCLUDING the IMA mirror images.  isSelf =
	// (the source face's element == the target element) drops the original center charge (singularity-free).
	// The reusable single-(target,source) kernel shared by BuildCentroidFieldGrad and MomentSystemEntry.
	void CentroidFieldGradFromFace(const double ce3[3], const double V4[4][3], const double srcCenter[3],
	                               bool isSelf, double area, double Hout[3], double gHout[6]) const;
	// O(N) geometry/sample cache for the HACApK moment callback path.  This is the multipole-moment
	// counterpart of the legacy 6-DoF PrecomputeHexaGeometry + cached 6x6 block path.
	void PrecomputeMomentGeometry() const;
	// O(N) geometry/sample cache for 5/6-DOF moment elements in CollectMomentElems order.
	void PrecomputeMomentAnyGeometry() const;
	// On-demand UN-normalized moment system entry A_raw[rowGlobal][colDOF] (the H-matrix entry; see
	// docs/multipole_moment_mmm/ACA_MOMENT_DESIGN.ipynb).  rowGlobal = 6*hpos + t over the valid 6-face hexes
	// (t: 0,1,2 dipole; 3 monopole; 4,5 diagonal-quadrupole); colDOF = global face DOF.  Computed ONLY from
	// the row element's local geometry + the on-demand centroid field/grad from face colDOF -- no full-system
	// build, no row normalization (the row 2-norm is a diagonal scaling that leaves the direct solve invariant).
	double MomentSystemEntry(int rowGlobal, int colDOF, const double* chiPerHex) const;
	void MomentSystemBlock6x6(int rowHexPos, int colHexPos, const double* chiPerHex, double* block, bool kernelOnly = false) const;
	// Pure-hex hot path for method-1 matrix-free BiCGSTAB.  Computes
	// y = diag(chi) * K_geometry * x using the same moment quadrature as
	// MomentSystemBlock6x6(..., kernelOnly=true), but without materializing a
	// temporary 6x6 block for every element pair.
	void MomentKernelMatVec6x6(const double* x, const double* chiPerHex, double* y) const;
	// Padded 6x6 on-demand moment block in CollectMomentElems order.  The active block is
	// row_dof x col_dof, where row/col dof are 5 or 6 from the corresponding moment elements.
	void MomentSystemBlockAny(int rowMomPos, int colMomPos, const double* chiPerMom, double* block, bool kernelOnly = false) const;
	// multipole-moment MMM system matrix (parameter-free replacement for EIEM2): per moment element assemble
	// 3 dipole + 1 monopole + residual quadrupole rows = moment of sigma matched to chi*{H,gradH}(centroid)
	// (global field/grad from BuildCentroidFieldGrad, local moments from the element geometry).  Rows 2-norm
	// normalized.  A row-major (dof x dof), rhs length dof.  Uniform linear chi + uniform applied field Happ
	// (Step-1 verification path vs examples/vim moment prototype; hex/wedge/pyramid).
	void BuildMomentSystem(double chi, const double Happ[3], std::vector<double>& A, std::vector<double>& rhs) const;  // uniform-field wrapper
	// Per-element chi + per-element external field (at moment-element centroid, HextPerHex[h*3+k]); A column = face DOF
	// so dgesv's solution is sigma in DOF order.  The solve path uses this (coil sources are not uniform).
	void BuildMomentSystemCore(const double* chiPerHex, const double* HextPerHex, std::vector<double>& A, std::vector<double>& rhs, bool normalize = true) const;
	double LastMomentFieldGradTime() const { return m_lastMomentFieldGradTime; }
	double LastMomentSystemBuildTime() const { return m_lastMomentSystemBuildTime; }
	// EIEM2 surface-charge block kernels (Compute6x6BlockFast / Compute5x5BlockFast /
	// ComputeMixedBlockFast) retired Phase 3b -- multipole-moment MMM (BuildMomentSystemCore) is the sole
	// surface-charge demag. Surface-charge (MMMM) solves use the dense LU / matrix-free moment path -- MMMM
	// does NOT connect to HACApK; the HACApK manager is MMM-only (3x3 tet).
	void FieldFromTrianglePrecomputed(int hex_idx, int tri_idx, const double* obs, double sigma, double* H_out) const;

	//-------------------------------------------------------------------------
	// Fast wedge geometry precompute (face triangles/normals; shared by moment assembly)
	//-------------------------------------------------------------------------
	void PrecomputeWedgeGeometry();  // Pre-compute face triangles/normals/eval points

	//-------------------------------------------------------------------------
	// IMA (Image) Symmetry Methods
	// Reference: ELF_MAGIC IMA approach - matrix construction with image summation
	//-------------------------------------------------------------------------

	// Configure IMA symmetry mode
	// symmetry: IMA_X, IMA_Y, IMA_Z, IMA_XY, etc.
	// signX, signY, signZ: +1 for symmetric BC, -1 for antisymmetric BC per axis
	// Returns: number of elements in IMA region (reduced set)
	int SetIMASymmetry(int symmetry, int signX = 1, int signY = 1, int signZ = 1);

	// Check if element is in IMA region (positive half-space for all active mirrors)
	bool IsElementInIMARegion(int elemIdx) const;

	// Get mirror element index for a given element and symmetry axis
	int GetMirrorElementIndex(int elemIdx, int symmetryAxis) const;

	// Build IMA interaction matrix (with image summation)
	// N_IMA[i,j] = N[i,j] ± N[i, mirror_j] @ P
	// Sign determined by m_imaSign: +1 (symmetric BC) or -1 (antisymmetric BC)
	int SetupInteractMatrix_IMA(bool skipDenseMatrix = false);

	// Legacy IMA functions removed (2026-03-31): ApplyDOFPermutation, ApplyRowPermutation,
	// Compute6x6BlockIMA, Compute6x6BlockMirrored, Compute6x6BlockMirroredTarget
	// IMA mirror logic is now inline in the multipole-moment MMM centroid field/grad path and Compute3x3BlockFast.

	// IMA accessors
	bool IsIMAEnabled() const { return m_imaEnabled; }
	int GetIMASymmetry() const { return m_imaSymmetry; }
	int GetIMASignX() const { return m_imaSignX; }
	int GetIMASignY() const { return m_imaSignY; }
	int GetIMASignZ() const { return m_imaSignZ; }
	int GetIMANumElements() const { return m_imaNumElements; }

	void SetupExternFieldArray();
	void AddExternFieldFromMoreExtSource();
	void AddMoreExternField(const radThg& hExtraExtSrc);
	//void ZeroAuxOldMagnArray();
	//void StoreAuxOldMagnArray();
	void SubstractOldMagn();
	void AddOldMagn();
	double CalcQuadNewOldMagnDif();
	int CountRelaxElemsWithSym();
	int OutAmOfRelaxObjs() { return AmOfMainElem;}
	void FindMaxModMandH(double& MaxModM, double& MaxModH);

	inline void PushFrontNativeElemTransList(radTg3d*, radTlphgPtr*);
	inline void EmptyVectOfPtrToListsOfTrans();

	inline void FillInTransPtrVectForElem(int, char);
	inline void EmptyTransPtrVect();

	void NestedFor_Trans(radTrans*, const radTlphgPtr::const_iterator&, int, char);
	inline void AddTransOrNestedFor(radTrans*, const radTlphgPtr::const_iterator&, int, char);

	void FillInMainTransPtrArray();
	inline void DestroyMainTransPtrArray();

	void FillInRelaxSubIntervArray(); //New
	void AddRelaxSubInterval(int StartNo, int FinNo, TRelaxSubIntervalID SubIntervalID); //New

	int NotEmpty() { return (AmOfMainElem==0)? 0 : 1;}
	// Dump / DumpBin REMOVED (Phase B2b/B2c, 2026-04-15)

	int Type_g() { return 4;}

	inline void ResetM();
	inline void ResetAuxParam();
	inline void InitAuxArrays();

	void ZeroAuxOldArrays(); //OC300504
	inline void StoreAuxOldArrays(); //OC300504
	inline void RestoreAuxOldArrays(); //OC300504

	inline void OutRelaxStatusParam(double*);
	inline void ShowInteractVector(char);
	inline void ShowInteractMatrix();

	inline int SizeOfThis();

	inline void UpdateExternalField();

	inline void LongLongToFloatAr(long long, float*); //OC24122019
	inline long long FloatArToLongLong(float*); //OC27122019

	template<class T> inline long OutMagnVals(T*& arMagnVals, int nExtraVals=0); //OC02012020
	template<class T> inline void SetRelaxObjMagnVals(T* arMagnVals); //OC02012020

	friend class radTIterativeRelaxMeth;
	friend class radTRelaxationMethNo_0;   // LU direct solver
	friend class radTRelaxationMethNo_1;   // BiCGSTAB
#ifdef RADIA_USE_HACAPK
	friend class radTRelaxationMethNo_2;   // BiCGSTAB + HACApK
	friend class RadHACApKMMMManager;      // HACApK H-matrix manager (MSC kernel)
#endif
};

//-------------------------------------------------------------------------

inline void radTInteraction::PushFrontNativeElemTransList(radTg3d* g3dPtr, radTlphgPtr* ListOfPtrToTransPtr)
{
	for(radTlphg::iterator TrIter = g3dPtr->g3dListOfTransform.begin();
		TrIter != g3dPtr->g3dListOfTransform.end(); ++TrIter)
		ListOfPtrToTransPtr->push_back(&(*TrIter)); // Improve dereferentiation?
}

//-------------------------------------------------------------------------

inline void radTInteraction::EmptyVectOfPtrToListsOfTrans()
{
	for(unsigned i=1; i<IntVectOfPtrToListsOfTransPtr.size(); i++)
	//for(unsigned i=0; i<IntVectOfPtrToListsOfTransPtr.size(); i++) //OC30122019: this correction was suggested by per-gron
	{
		radTlphgPtr*& p_lphgPtr = IntVectOfPtrToListsOfTransPtr[i];
		if(p_lphgPtr != 0) delete p_lphgPtr;
		p_lphgPtr = 0;
	}
	IntVectOfPtrToListsOfTransPtr.erase(IntVectOfPtrToListsOfTransPtr.begin(), IntVectOfPtrToListsOfTransPtr.end());
	for(unsigned k=1; k<ExtVectOfPtrToListsOfTransPtr.size(); k++) 
	//for(unsigned k=0; k<ExtVectOfPtrToListsOfTransPtr.size(); k++) //OC30122019: this correction was suggested by per-gron
	{
		radTlphgPtr*& p_lphgPtr = ExtVectOfPtrToListsOfTransPtr[k];
		if(p_lphgPtr != 0) delete p_lphgPtr;
		p_lphgPtr = 0;
	}
	ExtVectOfPtrToListsOfTransPtr.erase(ExtVectOfPtrToListsOfTransPtr.begin(), ExtVectOfPtrToListsOfTransPtr.end());
}

//-------------------------------------------------------------------------

inline void radTInteraction::FillInTransPtrVectForElem(int ElemLocInd, char I_or_E)
{
	radTlphgPtr* PtrToListOfPtrToTrans = nullptr;
	if(I_or_E == 'I') PtrToListOfPtrToTrans = IntVectOfPtrToListsOfTransPtr[ElemLocInd];
	else PtrToListOfPtrToTrans = ExtVectOfPtrToListsOfTransPtr[ElemLocInd];

	if(PtrToListOfPtrToTrans->empty()) TransPtrVect.push_back(IdentTransPtr);
	else NestedFor_Trans(IdentTransPtr, PtrToListOfPtrToTrans->begin(), ElemLocInd, I_or_E);
}

//-------------------------------------------------------------------------

inline void radTInteraction::EmptyTransPtrVect()
{
	if(Cast.IdentTransCast(TransPtrVect[0])==0) delete TransPtrVect[0];
	for(unsigned i=1; i<TransPtrVect.size(); i++) delete TransPtrVect[i];
	TransPtrVect.erase(TransPtrVect.begin(), TransPtrVect.end());
}

//-------------------------------------------------------------------------

inline void radTInteraction::AddTransOrNestedFor(radTrans* BaseTransPtr, const radTlphgPtr::const_iterator& Iter, int ElemLocInd, char I_or_E)
{
	radTlphgPtr* PtrToListOfPtrToTrans = nullptr;
	if(I_or_E == 'I') PtrToListOfPtrToTrans = IntVectOfPtrToListsOfTransPtr[ElemLocInd];
	else PtrToListOfPtrToTrans = ExtVectOfPtrToListsOfTransPtr[ElemLocInd];

	if(Iter == PtrToListOfPtrToTrans->end()) 
	{
		if(Cast.IdentTransCast(BaseTransPtr) == 0) TransPtrVect.push_back(new radTrans(*BaseTransPtr));
		else TransPtrVect.push_back(BaseTransPtr);
	}
	else NestedFor_Trans(BaseTransPtr, Iter, ElemLocInd, I_or_E);
}

//-------------------------------------------------------------------------

inline void radTInteraction::DestroyMainTransPtrArray()
{
	if(MainTransPtrArray == 0) return;

	for(int i=0; i<AmOfMainElem; i++)
	{
		radTrans* MainTransPtr = MainTransPtrArray[i];
		if(MainTransPtr != 0)
		{
			if(Cast.IdentTransCast(MainTransPtr)==0)
			{
				delete (MainTransPtr);
				MainTransPtr = 0;
			}
		}
	}
	// RAII: automatic cleanup via vMainTransPtrArray
	MainTransPtrArray = 0;
}

//-------------------------------------------------------------------------

inline void radTInteraction::ResetM()
{
	for(int i=0; i<AmOfMainElem; i++)
	{
		radTg3dRelax* g3dRelaxPtr = g3dRelaxPtrVect[i];

		// Initialize magnetization with remnant magnetization (usually zero for soft materials)
		g3dRelaxPtr->Magn = ((radTMaterial*)(g3dRelaxPtrVect[i]->MaterHandle.rep))->RemMagn;
		NewMagnArray[i] = g3dRelaxPtr->Magn;

		// Initialize field with external field (from ObjBckg or ObjBckgCF)
		// This allows relaxation to iterate from external field as starting point
		NewFieldArray[i] = ExternFieldArray[i];
	}
}

//-------------------------------------------------------------------------

inline void radTInteraction::ResetAuxParam()
{
	for(int i=0; i<AmOfMainElem; i++)
	{
		radTg3dRelax* g3dRelaxPtr = g3dRelaxPtrVect[i];

		g3dRelaxPtr->AuxFloat1 = 0;
		g3dRelaxPtr->AuxFloat2 = 0;
		g3dRelaxPtr->AuxFloat3 = 0;
	}
}

//-------------------------------------------------------------------------

inline void radTInteraction::InitAuxArrays()
{
	for(int i=0; i<AmOfMainElem; i++)
	{
		NewMagnArray[i] = g3dRelaxPtrVect[i]->Magn;
		NewFieldArray[i] = TVector3d(0.,0.,0.); // Or make it TVector3df
	}
}

//-------------------------------------------------------------------------

//inline void radTInteraction::StoreOldMagnData() //OC300504
//{
//	if(AuxOldMagnArray == nullptr) return;
//
//    TVector3d *tAuxOldMagnArray = AuxOldMagnArray;
//	for(int i=0; i<AmOfMainElem; i++)
//	{
//		*(tAuxOldMagnArray++) = g3dRelaxPtrVect[i]->Magn;
//	}
//}

//-------------------------------------------------------------------------

inline void radTInteraction::StoreAuxOldArrays()
{
	if(AmOfMainElem <= 0) return;
	
	if((AuxOldMagnArray != nullptr) && (AuxOldFieldArray != nullptr))
	{
		TVector3d *tAuxOldMagn = AuxOldMagnArray;
		TVector3d *tAuxOldField = AuxOldFieldArray;
		TVector3d *tNewFieldArray = NewFieldArray;

		for(int StNo=0; StNo<AmOfMainElem; StNo++)
		{
			TVector3d &M = (g3dRelaxPtrVect[StNo])->Magn; 
			*(tAuxOldMagn++) = M;
			*(tAuxOldField++) = *(tNewFieldArray++);
		}
	}
	else
	{
		if(AuxOldMagnArray != nullptr)
		{
			TVector3d *tAuxOldMagn = AuxOldMagnArray;
			for(int StNo=0; StNo<AmOfMainElem; StNo++)
			{
				TVector3d &M = (g3dRelaxPtrVect[StNo])->Magn; 
				*(tAuxOldMagn++) = M;
			}
		}
		if(AuxOldFieldArray != nullptr)
		{
			TVector3d *tAuxOldField = AuxOldFieldArray;
			TVector3d *tNewFieldArray = NewFieldArray;
			for(int StNo=0; StNo<AmOfMainElem; StNo++)
			{
				*(tAuxOldField++) = *(tNewFieldArray++);
			}
		}
	}
}

//-------------------------------------------------------------------------

inline void radTInteraction::RestoreAuxOldArrays() //OC300504
{
	if((AuxOldMagnArray == nullptr) && (AuxOldFieldArray == nullptr)) return;

	TVector3d *tAuxOldMagnArray = AuxOldMagnArray;
	TVector3d *tNewMagnArray = NewMagnArray;
	TVector3d *tAuxOldFieldArray = AuxOldFieldArray;
	TVector3d *tNewFieldArray = NewFieldArray;

	if((AuxOldMagnArray != nullptr) && (AuxOldFieldArray != nullptr))
	{
		for(int i=0; i<AmOfMainElem; i++)
		{
			g3dRelaxPtrVect[i]->Magn = *tAuxOldMagnArray;
			*(tNewMagnArray++) = *(tAuxOldMagnArray++);
			*(tNewFieldArray++) = *(tAuxOldFieldArray++);
		}
	}
	else
	{
		if(AuxOldMagnArray != nullptr)
		{
			for(int i=0; i<AmOfMainElem; i++)
			{
				g3dRelaxPtrVect[i]->Magn = *tAuxOldMagnArray;
				*(tNewMagnArray++) = *(tAuxOldMagnArray++);
			}
		}
		if(AuxOldFieldArray != nullptr)
		{
			for(int i=0; i<AmOfMainElem; i++)
			{
				*(tNewFieldArray++) = *(tAuxOldFieldArray++);
			}
		}
	}
}

//-------------------------------------------------------------------------

inline void radTInteraction::ShowInteractVector(char Ch)
{
	TVector3d* Vect3dPtr = nullptr;
	switch(Ch) 
	{
		case 'E':
			Vect3dPtr = ExternFieldArray; break;
		case 'T':
			Vect3dPtr = NewFieldArray; break;
		case 'M':
			Vect3dPtr = NewMagnArray; break;
		default :
			Vect3dPtr = NewFieldArray; break;
	}
	Send.ArrayOfVector3d(Vect3dPtr, AmOfMainElem);
}

//-------------------------------------------------------------------------

inline void radTInteraction::ShowInteractMatrix()
{
	Send.MatrixOfMatrix3d(InteractMatrix, AmOfMainElem, AmOfMainElem);
}

//-------------------------------------------------------------------------

// radTInteraction::Dump REMOVED (Phase B2b, 2026-04-15)

//-------------------------------------------------------------------------

inline int radTInteraction::SizeOfThis()
{
	long GenSize = sizeof(*this);
	GenSize += AmOfMainElem*AmOfMainElem*sizeof(TMatrix3d);
	GenSize += AmOfMainElem*sizeof(TMatrix3d*);
	GenSize += 3*AmOfMainElem*sizeof(TVector3d);
	GenSize += AmOfMainElem*sizeof(radTg3dRelax*);
	GenSize += AmOfRelaxSubInterv*sizeof(radTRelaxSubInterval);
	return GenSize;
}

//-------------------------------------------------------------------------

inline void radTInteraction::OutRelaxStatusParam(double* RelaxStatusParamArray)
{
	RelaxStatusParamArray[0] = RelaxStatusParam.MisfitM;
	RelaxStatusParamArray[1] = RelaxStatusParam.MaxModM;
	RelaxStatusParamArray[2] = RelaxStatusParam.MaxModH;
	// Add more members of radTRelaxStatusParam here, should they appear in Future
}

//-------------------------------------------------------------------------

inline void radTInteraction::UpdateExternalField()
{
	SetupExternFieldArray(); //zeros and then sets the ExternFieldArray from g3dExternPtrVect
	AddExternFieldFromMoreExtSource(); //adds field to ExternFieldArray from MoreExtSourceHandle
}

//-------------------------------------------------------------------------

template<class T> inline long radTInteraction::OutMagnVals(T*& arMagnVals, int nExtraVals) //OC02012020
{//Allocates arMagnVals!
 //This assumes that after the Relaxation, the resulting Magnetization data is stored in radTInteraction::g3dRelaxPtrVect
 //(but can be changed to take the Magnetization data from radTInteraction::NewMagnArray)
	//arMagnVals = 0;
	if(AmOfMainElem <= 0) return 0;
	long nMagnVals = AmOfMainElem*3 + nExtraVals; //nExtraVals can be used e.g. for attaching extra data to the Magnetization array

	if(arMagnVals == 0) arMagnVals = new T[nMagnVals];
	if(arMagnVals == 0) { Send.ErrorMessage("Radia::Error900"); return 0;}
	T *t_arMagnVals = arMagnVals + nExtraVals;
	for(long i=0; i<AmOfMainElem; i++)
	{
		TVector3d &M = (g3dRelaxPtrVect[i])->Magn;
		*(t_arMagnVals++) = (T)M.x; *(t_arMagnVals++) = (T)M.y; *(t_arMagnVals++) = (T)M.z;
	}
	return nMagnVals;
}

//-------------------------------------------------------------------------

template<class T> inline void radTInteraction::SetRelaxObjMagnVals(T* arMagnVals)
{
	if((AmOfMainElem <= 0) || (arMagnVals == 0)) return;

	T *t_arMagnVals = arMagnVals;
	for(long i=0; i<AmOfMainElem; i++)
	{
		TVector3d &M = (g3dRelaxPtrVect[i])->Magn;
		 M.x = *(t_arMagnVals++);  M.y = *(t_arMagnVals++); M.z = *(t_arMagnVals++);
	}
}

//-------------------------------------------------------------------------

inline void radTInteraction::LongLongToFloatAr(long long inVal, float outAr[4]) //24122019
{//Aux. function to encode one long long number to 4 floats (consider moving to some parser class)
	long long mask = 0xFFFF;
	for(int j=0; j<4; j++)
	{
		outAr[j] = (float)((inVal & mask) >> (j*16));
		mask <<= 16;
	}
}

//-------------------------------------------------------------------------

inline long long radTInteraction::FloatArToLongLong(float inAr[4]) //OC27122019
{//Aux. function to decode one long long number from 4 floats (consider moving to some parser class)
	long long res=0, aux;
	for(int j=0; j<4; j++)
	{
		aux = (long long)inAr[j];
		aux <<= (j*16);
		res |= aux;
	}
	return res;
}

//-------------------------------------------------------------------------

#endif
