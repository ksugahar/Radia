/////////////////////////////////////////////////////////////////////////////
// Name:        radentry.h
// Purpose:     Radia DLL header file
// Authors:     O.Chubar, P.Elleaume
// Version:     DLL 4.115
// Modified:    30.05.04
/////////////////////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------
// Platform- and compiler-dependent macro definitions
//-------------------------------------------------------------------------
#if !(defined(ALPHA_NONE) || defined(ALPHA__LIB__))
/*---------------For CodeWarrior PowerMac---------------*/
#if defined __POWERPC__
#if defined ALPHA__DLL__ || defined MATLAB_MEX_FILE
#define EXP __declspec(export)
#endif
/*---------------For CodeWarrior PC and Visual C++---------------*/
#elif defined __INTEL__ || defined WIN32
#if defined ALPHA__DLL__ || defined MATLAB_MEX_FILE
#define EXP __declspec(dllexport)
#else
#define EXP __declspec(dllimport)
#endif
//#define CALL __stdcall
#define CALL __cdecl //use one of these calling convntions at preparing DLL
/*---------------For HP-UX, gcc---------------*/
#else
#endif
#endif /*ALPHA_NONE*/
// Forward declaration for Python callback support
struct _object;
typedef _object PyObject;


#ifndef EXP
#define EXP
#endif

#ifndef CALL
#define CALL
#endif

//-------------------------------------------------------------------------

#ifdef __cplusplus  
extern "C" {
#endif

/* OK = 0 */
#define OK 0

/** Returns the error message associated to error number. 
This function cannot be called from Visual Basic For Applications. 
@param er [in] error number 
@return the chain of characters representing the error string
@author P. Elleaume 
@version 1.0 
@see GetErrorSize and GetErrorText 
*/ 
EXP const char* CALL RadErrGet(int er); 

/** Returns the length of the error message string not counting "\0". 
@param siz [out] length of the error message  
@param er [in] error number 
@return integer error code (0 : No Error, >0 : Error Number, <0 : Warning Number) 
@author P. Elleaume 
@version 1.0 
*/ 
EXP int CALL RadErrGetSize(int* siz, int er);

/** Returns the text of error message associated to error number. 
@param t [out] error message string
@param er [in] error number 
@return integer error code (0 : No Error, >0 : Error Number, <0 : Warning Number)
@author P. Elleaume 
@version 1.0 
*/ 
EXP int CALL RadErrGetText(char* t, int er);

/** Returns the warning message associated to warning number. 
This function cannot be called from Visual Basic For Applications. 
@param er [in] warning number 
@return the chain of characters representing the warning string
@author P. Elleaume 
@version 1.0
@see GetWarningSize and GetWarningText 
*/ 
EXP const char * CALL RadWarGet(int er); 

/** Returns the length of the warning message not counting "\0". 
@param siz [out] length of the warning message  
@param er [in] warning number 
@return	integer error code (0 : No Error, >0 : Error Number, <0 : Warning Number)
@author P. Elleaume 
@version 1.0 
@see dllGetWarningText
*/ 
EXP int CALL RadWarGetSize(int* siz, int er);

/** Returns the text of warning message associated to error number. 
@param t [out] warning message string
@param er [in] error number 
@return integer error code (0 : No Error, >0 : Error Number, <0 : Warning Number)
@author P. Elleaume 
@version 1.0 
@see GetWarningSize
*/ 
EXP int CALL RadWarGetText(char* t, int er);

/** Creates a uniformly magnetized rectangular parallelepiped.
The parallelepiped block is defined through its center point P[3], dimensions L[3], and magnetization M[3]."
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the block center of gravity
@param L [in] array of 3 dimensions of the block
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjRecMag(int* n, double* P, double* L, double* M);

/** Creates a uniformly magnetized extruded polygon.
The extrusion axis is directed along X axis of the laboratory frame.
@param n [out] reference number of the object created
@param xc [in] the horizontal coordinate of the block center of gravity
@param lx [in] the thickness (extrusion size)
@param FlatVert [in] flat array of y and z coordinates (y1, z1, y2, z2,...) of vertex points describing the polygon in 2D
@param nv [in] number of vertex points of the 2D polygon (the length of the FlatVert array is 2*nv)
@param a [in] character identifying extrusion direction (can be 'x', 'y' or 'z')
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjThckPgn(int* n, double xc, double lx, double* FlatVert, int nv, char a, double* M);

/** Creates a uniformly magnetized polyhedron (closed volume limited by planes).
@param n [out] reference number of the object created
@param FlatVert [in] flat array of x, y and z coordinates (x1, y1, z1, x2, y2, z2,...) of the polyhedron vertex points
@param nv [in] number of vertex points of the polyhedron (the length of the FlatVert array is 3*nv)
@param FlatFaces [in] flat array of indexes of vertex points defining the polyhedron faces (f1i1, f1i2,..., f2i1, f2i2,...)
@param FacesLen [in] array of integer numbers equal to the numbers of vertex points in each face of the polyhedron; the order of the faces is the same as in the FlatFaces array (the length of the FlatFaces array is equal to the sum of elements of the FacesLen array) 
@param nf [in] number of faces of the polyhedron (or the length of the FacesLen array)
@param M [in] array of 3 cartesian components of magnetization vector inside the block
@param M_LinCoef [in] array of 9 coefficients of linearly-varying magnetization vector inside the block
@param J [in] array of 3 cartesian components of current density vector inside the block
@param J_LinCoef [in] array of 9 coefficients of linearly-varying current density vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjPolyhdr(int* n, double* FlatVert, int nv, int* FlatFaces, int* FacesLen, int nf, double* M, double* M_LinCoef, double* J, double* J_LinCoef);

/**Creates a uniformly magnetized finite-length arc of polygonal cross-section.
@param n [out] reference number of the object created
@param P [in] array of 2 cartesian coordinates defining the position of the rotation axis in the plane perpendicular to this axis
@param a [in] character defining the orientation of the rotation axis in 3D space; it can be either to 'x', 'y' or 'z'
@param FlatVert [in] flat array of radial and axial coordinates (r1, z1, r2, z2,...) of vertex points of the cross-section polygon
@param nv [in] number of vertex points of the cross-section polygon (the length of the FlatVert array is 2*nv)
@param Phi [in] array of 2 numbers - initial and final azimuth angles
@param nseg [in] number of segments in the arc
@param sym_no [in] character which can be either 's' or 'n'; depending on the value of this switch, the magnetization vectors in nseg sector polyhedrons either will be linked by rotational symmetry ('s'), or will behave as independent magnetization vectors at any subsequent relaxation
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjArcPgnMag(int* n, double* P, char a, double* FlatVert, int nv, double* Phi, int nseg, char sym_no, double* M);

/** ?? TO IMPLEMENT ??
Creates a uniformly magnetized volume obtained by rotation of a planar convex polygon over 2 Pi around pre-defined axis.
@param n [out] reference number of the object created
@param P [in] array of 2 cartesian coordinates defining the position of the rotation axis in the plane perpendicular to this axis
@param FlatVert [in] flat array of radial and axial coordinates (r1, z1, r2, z2,...) of vertex points of the cross-section polygon
@param nv [in] number of vertex points of the cross-section polygon (the length of the FlatVert array is 2*nv)
@param nseg [in] number of segments in the arc
@param sym_no [in] character which can be either 's' or 'n'; depending on the value of this switch, the magnetization vectors in nseg sector polyhedrons either will be linked by rotational symmetry ('s'), or will behave as independent magnetization vectors at any subsequent relaxation
@param a [in] character defining the orientation of the arc rotation axis; it can be either to 'x', 'y' or 'z'
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
//EXP int CALL RadObjCircPgnMag(int* n, double* P, double* FlatVert, int nv, int nseg, char sym_no, char a, double* M);

/** Attempts to create one uniformly magnetized convex polyhedron or a set of convex polyhedrons based on slices.
The slices are assumed to be convex planar polygons parallel to the XY plane.
@param n [out] reference number of the object created
@param FlatVert [in] flat array of x and y coordinates (x11, y11, x12, y12,..., x21, y21, x22, y22,...) of vertex points of the slices (planar polygons)
@param SlicesLen [in] array of integer numbers specifying numbers of vertex points in each slice; the order of the slices is the same as in the FlatVert array (the length of the FlatVert array is twice the sum of elements of the SlicesLen array) 
@param Attitudes [in] array of vertical coordinates (or attitudes) of the slices, in the ascending order (which is the same as for the SlicesLen array)
@param ns [in] number of slices (or the length of the Attitudes and SlicesLen arrays)
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the whole block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjMltExtPgn(int* n, double* FlatVert, int* SlicesLen, double* Attitudes, int ns, double* M);

/** Attempts to create one uniformly magnetized convex polyhedron or a set of convex polyhedrons based on rectangular slices.
The rectangular slices are assumed to be parallel to the XY plane.
@param n [out] reference number of the object created
@param FlatCenPts [in] flat array of x, y and z coordinates (x1, y1, z1, x2, y2, z2,...) of the slices center points
@param FlatRtgSizes [in] flat array of sizes of the slice rectangles along x and y (wx1, wy1, wx2, wy2,...)
@param ns [in] number of slices (the length of the FlatCenPts array is 3*ns, the length of the FlatRtgSizes array is 2*ns)
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the whole block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjMltExtRtg(int* n, double* FlatCenPts, double* FlatRtgSizes, int ns, double* M);

/** Disabled legacy triangulated extruded polygon API.

This function is kept for C API compatibility only. Radia no longer bundles
the Triangle triangulation library; use the Netgen/Cubit mesh workflow instead.

@param n [out] reference number of the object created
@param xc [in] the horizontal coordinate of the block center of gravity
@param lx [in] the thickness (extrusion size)
@param FlatVert [in] flat array of y and z coordinates (y1, z1, y2, z2,...) of vertex points describing the polygon in 2D
@param FlatSubd [in] flat array of subdivision parameters for base polygon segments, two numbers for each segment: the first defining number of sub-segments, the second - "gradient" of the segmentation
@param nv [in] number of vertex points of the 2D polygon (the length of the FlatVert and FlatSubd arrays is 2*nv)
@param a [in] character identifying extrusion direction (can be 'x', 'y' or 'z')
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the block
@param opt [in] pointer to options string, which can be e.g. "ki->...,TriAngMin->...,TriAreaMax->...,ExtOpt->..." or each of these tokens separately, or 0; default values are "ki->Numb" (rather than "ki->Size"), "TriAngMin->20"
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjMltExtTri(int* n, double xc, double lx, double* FlatVert, double* FlatSubd, int nv, char a, double* M, char* opt);

/** Creates a finite-length arc magnet of rectangular cross-section.
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the arc center point
@param R [in] array of 2 numbers - inner and outer radii
@param Phi [in] array of 2 numbers - initial and final azimuth angles
@param h [in] height
@param nseg [in] number of segments
@param a [in] character defining the orientation of the rotation axis of the arc; it can be either to 'x', 'y' or 'z'
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the whole block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
//EXP int CALL RadObjArcMag(int* n, double* P, double* R, double* Phi, double h, int nseg, char a, double* M);

/** Creates a cylindrical magnet.
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the cylinder center point
@param r [in] cylinder radius
@param h [in] cylinder height
@param nseg [in] number of segments
@param a [in] character defining the orientation of the cylinder axis; it can be either to 'x', 'y' or 'z'
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the whole block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjCylMag(int* n, double* P, double r, double h, int nseg, char a, double* M);

/** Creates a current carrying rectangular parallelepiped block.
The parallelepiped block is defined through its center point P[3], dimensions L[3], and current density vector J[3]."
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the block center of gravity
@param L [in] array of 3 dimensions of the block
@param J [in] array of 3 cartesian coordinates of the current density vector inside the block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjRecCur(int* n, double* P, double* L, double* J);

/** Creates a current carrying finite-length arc of rectangular cross-section.
The arc rotation axis is directed along Z.
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the arc center point
@param R [in] array of 2 numbers - inner and outer radii
@param Phi [in] array of 2 numbers -  initial and final azimuth angles 
@param h [in] height
@param nseg [in] number of segments
@param man_auto [in] character which can be either 'm' or 'a'; the magnetic field from the arc is then computed based on the number of segments nseg ('m', i.e. "manual"), or on the general absolute precision level specified by the functions RadFldCmpCrt or RadFldCmpPrc ('a', i.e. "automatic")
@param a [in] character defining the orientation of the rotation axis of the arc; it can be either 'x', 'y' or 'z'
@param j [in] azimuthal current density
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjArcCur(int* n, double* P, double* R, double* Phi, double h, int nseg, char man_auto, char a, double j);

/** Creates a current carrying racetrack coil.
The coil consists of four 90-degree bents connected by four straight parts parallel to the XY plane.
@param n [out] reference number of the object created
@param P [in] array of 3 cartesian coordinates of the racetrack center point
@param R [in] array of 2 numbers - inner and outer radii
@param L [in] array of 2 numbers - straight section lengths
@param h [in] height
@param nseg [in] number of segments
@param man_auto [in] character which can be either 'm' or 'a'; the magnetic field from the arc is then computed based on the number of segments nseg ('m', i.e. "manual"), or on the general absolute precision level specified by the functions RadFldCmpCrt or RadFldCmpPrc ('a', i.e. "automatic")
@param a [in] character defining the orientation of the rotation axis of the arc; it can be either 'x', 'y' or 'z'
@param j [in] azimuthal current density
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjRaceTrk(int* n, double* P, double* R, double* L, double h, int nseg, char man_auto, char a, double j);

/** Creates a filament polygonal line conductor with current.
The line conductor is defined by sequence of points in 3D space.
@param n [out] reference number of the object created
@param FlatPts [in] flat array of x, y and z coordinates of the points (x1, y1, z1, x2, y2, z2,...)
@param np [in] number of points (the length of the array FlatPts is 3*np)
@param i [in] current
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjFlmCur(int* n, double* FlatPts, int np, double i);

/** Attempts to create a set of current-carrying convex polyhedron objects by applying a generalized extrusion to the initial planar convex polygon.
@param n [out] reference number of the object created
@param FlatVert [in] flat array of x and y coordinates (x1, y1, x2, y2,...) of vertex points of the base planar polygon
@param nv [in] number of vertex points of the 2D polygon (the length of the FlatVert array is 2*nv)

@param SlicesLen [in] array of integer numbers equal to the numbers of vertex points in each slice; the order of the slices is the same as in the FlatVert array (the length of the FlatVert array is twice the sum of elements of the SlicesLen array)
@param Attitudes [in] array of vertical coordinates (or attitudes) of the slices, in the ascending order (which is the same as for the SlicesLen array)
@param ns [in] number of slices (or the length of the Attitudes and SlicesLen arrays)
@param M [in] array of 3 cartesian coordinates of the magnetization vector inside the whole block
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
//EXP int CALL RadObjMltExtPgnCur(int* n, double* FlatVert, int nv, 
//	double* FlatVert, int* SlicesLen, double* Attitudes, int ns, double* M,      double z, char a, double i, char* opt);

/** Scales current (density) in a 3D object by multiplying it by a constant.
@param n [out] reference number of the object with current (density) to be scaled
@param scaleCoef [in] scaling coefficient
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjScaleCur(int n, double scaleCoef);

/** Creates a source of uniform background magnetic field.
@param n [out] reference number of the object created
@param B [in] array of 3 cartesian coordinates of the magnetic field vector
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjBckg(int* n, double* B);
EXP int CALL RadObjBckgCF(int* n, PyObject* callback);

/** Creates a container of magnetic field sources.
@param n [out] reference number of the object created
@param Objs [in] array of reference numbers of the objects (n1, n2, n3,...) to put to the container
@param nobj [in] number of objects (the length of the array Elems)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjCnt(int* n, int* Objs, int nobj);

/** Adds objects to the container object cnt.
@param cnt [in] reference number of the container object
@param Objs [in] array of reference numbers of the objects (n1, n2, n3,...) to put to the container
@param nobj [in] number of objects to put to the container (the length of the array Elems)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjAddToCnt(int cnt, int* Objs, int nobj);

/** Calculates the number of objects in the container cnt.
@param n [out] calculated number of objects
@param cnt [in] reference number of the container object
@param deep [in] switch specifying whether all atomic elements of eventual member containers have to be counted (1) or not (0, default)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjCntSize(int* n, int cnt);

/** Fills-in an array with the reference numbers of objects present in the container cnt. The array should be allocated in the calling application. The necessary size of the array can be determined using the function RadObjCntSize.
@param Objs [out] array of reference numbers of the objects present in the container
@param cnt [in] reference number of the container object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjCntStuf(int* Objs, int cnt);

/** Builds interaction matrix for magnetostatic problem.
@param n [out] reference number of the interaction matrix (for GetInteractMatrix)
@param ElemKey [in] reference number of the object to build matrix for
@param image [in] image symmetry string (e.g., "+x", "-z", "+x-z") or nullptr for no symmetry
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@note This replaces PreRelax + Image + BuildImageMatrix workflow
*/
EXP int CALL RadBuildMatrix(int* n, int ElemKey, const char* image);

/** Gets the interaction matrix built by BuildMatrix or Solve.
@param pMatrix [out] flat array to receive matrix data (column-major), or nullptr to query size
@param pDOF [out] number of DOF (matrix is pDOF x pDOF)
@param InteractElemKey [in] interaction handle from BuildMatrix, or 0 for last cached matrix
@return integer error code
*/
EXP int CALL RadGetInteractMatrix(double* pMatrix, int* pDOF, int InteractElemKey);

/** Gets per-DOF hex face geometry in the matrix DOF order (for div(B)=0 / RHS / moment studies).
Two-call pattern: pass pG=nullptr to read back DOF, allocate DOF*11 doubles, then call again to fill.
Each DOF row (ROW-MAJOR, stride 11): [elem_local, area, cx,cy,cz, nx,ny,nz(outward), ecx,ecy,ecz];
non-hex DOFs get elem_local=-1 and zeros.
@param pG [out] flat geometry array (DOF x 11, row-major), or nullptr to query size
@param pDOF [out] number of DOF (G has DOF rows)
@param InteractElemKey [in] interaction handle from BuildMatrix
@return integer error code
*/
EXP int CALL RadGetFaceGeom(double* pG, int* pDOF, int InteractElemKey);

/** Per moment-element centroid demag field+gradient functionals (the parameter-free moment formulation kernel).
For each moment element, the demag field H and gradient gradH at the element CENTROID as linear functionals
of every source DOF charge (SELF = bare charged face; MUTUAL = collocation MMMM dipole layer, singularity-free).
ROW-MAJOR (nMom x 9 x DOF): comp k (Hx,Hy,Hz, gxx,gyy,gzz,gxy,gxz,gyz), source DOF g -> C[(h*9+k)*DOF+g].
@param pC [out] flat array (nMom*9*DOF, row-major), or nullptr to query nMom/DOF
@param pNHex [out] number of moment elements (historical parameter name)
@param pDOF [out] total number of DOF
@param InteractElemKey [in] interaction handle from BuildMatrix
@return integer error code
*/
EXP int CALL RadGetCentroidFieldGrad(double* pC, int* pNHex, int* pDOF, int InteractElemKey);
EXP int CALL RadBuildMomentSystem(double chi, const double* Happ, double* pA, double* pRhs, int* pDOF, int InteractElemKey);
EXP int CALL RadMomentSystemDenseRaw(double chi, double* pA, int* pDOF, int InteractElemKey);

/** Densify the actual HACApK (ACA+) system operator into a dense matrix.
Builds the MSC H-matrix for the interaction handle and applies it to unit
vectors (A = -N + diag(1/chi), original DOF ordering). Use to validate the
H-matrix against the exact dense matrix (eigenvalues / deflation).
@param pMatrix [out] flat (pDOF x pDOF) row-major, or nullptr to query size
@param pDOF [out] number of DOF
@param InteractElemKey [in] interaction handle from BuildMatrix
@return integer error code
*/
EXP int CALL RadHMatrixDensify(double* pMatrix, int* pDOF, int InteractElemKey);

EXP int CALL RadHLUDebugMaterialize(int InteractElemKey, double *A_perm_out, int *lod_out, int *nd_out);  // Phase 4 debug

// RadPreRelax REMOVED (2026-01-31) - Use RadBuildMatrix instead
// The new API is: int handle = RadBuildMatrix(obj, image);
// where image is "+x", "-z", "+x-z", etc. for IMA symmetry

// DEPRECATED: Use RadSolve with method parameter instead
// EXP int CALL RadSetRelaxSubInterval(int InteractElemKey, int StartNo, int FinNo, int RelaxTogether);

/** Duplicates the object obj. 
@param n [out] reference number of the object created
@param obj [in] reference number of the object to duplicate
@param opt [in] pointer to an option string, which can be "FreeSym->False", "FreeSym->True" or 0. This specifies whether the symmetries (transformations with multiplicity more than one) previously applied to the object obj should be simply copied at the duplication ("FreeSym->False" or 0), or a container of new independent objects should be created in place of any symmetry previously applied to the object obj. In both cases the final object created by the duplication has exactly the same geometry as the initial object obj.
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjDpl(int* n, int obj, char* opt);

/** Provides coordinates of geometrical center point(s) and magnetization(s) of the object obj.
@param M [out] flat array of resulting magnetization array. If this pointer is 0 at input, the actual output data array can be obtained by calling RadUtiDataGet function. This aray will contain coordinates of geometrical center point(s) and magnetic field components. If obj is a container, the array will include the container members' center points and their magnetic field components. 
@param arMesh [out] flat array defining the structure of resulting magnetization array, i.e. number of dimensions (arMesh[0]) and number of values per dimension. The actual output data array can be obtained by calling RadUtiDataGet function. This aray will contain coordinates of geometrical center point(s) and magnetic field components. If obj is a container, the array will include the container members' center points and their magnetic field components. 
@param obj [in] reference number of a magnetic field source object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjM(double* M, int* arMesh, int obj); //OC27092018
//EXP int CALL RadObjM(double* M, int obj);

/** Provides coordinates of geometrical center point and magnetic field at that point.
@param B [out] flat array of resulting center point coordinates and field values. If this pointer is 0 at input, the actual output data array can be obtained by calling RadUtiDataGet function. This aray will contain coordinates of geometrical center point(s) and magnetic field components. If obj is a container, the array will include the container members' center points and their magnetic field components. 
@param arMesh [out] flat array defining the structure of resulting field array, i.e. number of dimensions (arMesh[0]) and number of values per dimension. The actual output data array can be obtained by calling RadUtiDataGet function. This aray will contain coordinates of geometrical center point(s) and magnetic field components. If obj is a container, the array will include the container members' center points and their magnetic field components. 
@param obj [in] reference number of a magnetic field source object
@param type [in] character identifying type of a magnetic field characteristic to return (can be 'A' or 'B' or 'H' or 'J' or 'M')
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjCenFld(double* B, int* arMesh, int obj, char type); //OC27092018
//EXP int CALL RadObjCenFld(double* B, int obj, char type);

/** Sets magnetization of the object obj.
@param obj [in] reference number of a magnetic field source object
@param M [in] flat array of 3 components of the magnetization vector
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjSetM(int obj, double* M);

// RadObjCutMag REMOVED (2026-01-14) - Use Cubit/Netgen for mesh operations
// See CLAUDE.md "Mesh Operations Policy" for details

/** Computes geometrical volume of a 3D object.
@param v [out] volume (in mm^3)
@param obj [in] reference number of a 3D object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadObjGeoVol(double* v, int obj);

/** Gives number of degrees of freedom for the relaxation of an object.
@param num [out] number of degrees of freedom
@param obj [in] reference number of a 3D object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjDegFre(int* num, int obj);

/** Creates a parallelepiped block with center point {P[0],P[1],P[2]}, dimensions {L[0],L[1],L[2]} and color {RGB[0],RGB[1],RGB[2]}.
The block is magnetized according to {M[0],M[1],M[2]} then subdivided according to {K[0],K[1],K[2]} and added into the container grp. grp should be defined in advance by calling RadObjCnt().
@param n [out] reference number of the object created
@param P [in] three cartesian coordinates of the block center point
@param L [in] block dimensions
@param M [in] three components of the magnetization vector
@param K [in] array of subdivision parameters
@param nK [in] length of array of subdivision parameters
@param grp [in] reference number of the container object
@param mat [in] reference number of the magnetic material
@param RGB [in] array of 3 numbers from 0 to 1 specifying intensities of red, green and blue colors
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadObjFullMag(int* n, double* P, double* L, double* M, double* K, int nK, int grp, int mat, double* RGB);

// RadTrfPlSym REMOVED (2026-01-31) - Use Image symmetry instead
// Use: RadSolve(..., image="+x") or RadBuildMatrix(obj, image="+x")

/** Creates a rotation.
@param trf [out] reference number of the symmetry object created
@param P [in] array of 3 numbers representing cartesian coordinates of a point belonging to the rotation axis
@param V [in] array of 3 numbers representing cartesian coordinates of a vector parallel to the rotation axis
@param phi [in] rotation angle
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfRot(int* trf, double* P, double* V, double phi);

/** Creates a translation.
@param trf [out] reference number of the symmetry object created
@param V [in] array of 3 numbers representing cartesian coordinates of the translation vector
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfTrsl(int* trf, double* V);

/** Creates a field inversion.
@param trf [out] reference number of the symmetry object created
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfInv(int* trf);

/** Multiplies original space transformation origtrf by another transformation trf from left.
@param fintrf [out] reference number of the final space transformation
@param origtrf [in] reference number of the original space transformation to be multiplied
@param trf [in] reference number of another space transformation
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfCmbL(int* fintrf, int origtrf, int trf);

/** Multiplies original space transformation origtrf by another transformation trf from right.
@param fintrf [out] reference number of the final space transformation
@param origtrf [in] reference number of the original space transformation to be multiplied
@param trf [in] reference number of another space transformation
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadTrfCmbR(int* fintrf, int origtrf, int trf);

// RadTrfMlt REMOVED (2026-01-31) - Use Image symmetry instead
// The shared-DOF approach was fundamentally incompatible with MSC 6DOF hexahedra
// Use: RadSolve(..., image="+x-z") or RadBuildMatrix(obj, image="+x-z")

/** Orients object obj by applying transformation trf to it once.
@param objout [out] reference number of the final object with space transformation applied
@param obj [in] reference number of the original object
@param trf [in] reference number of a space transformation
*/
EXP int CALL RadTrfOrnt(int* objout, int obj, int trf);

// RadTrfZerPara REMOVED (2026-01-31) - Use Image symmetry instead
// RadTrfZerPerp REMOVED (2026-01-31) - Use Image symmetry instead
// These functions used RadTrfMlt internally, which has fundamental issues with MSC 6DOF hexahedra
// Use: RadSolve(..., image="+x-z") or RadBuildMatrix(obj, image="+x-z")

/** Applies material mat to object obj.
@param objout [out] reference number of the final object with material applied
@param obj [in] reference number of the original object
@param mat [in] reference number of the material to be applied
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatApl(int* objout, int obj, int mat);

/** Creates a linear anisotropic magnetic material.
@param mat [out] reference number of the material created
@param Ksi [in] array of 2 magnetic susceptibility values for the directions parallel and perpendicular to the easy magnetization axis
@param Mr [in] array of 3 cartesian coordinates of the remanent magnetization vector
@param nMr [in] number of components of the remanent magnetization vector. If nMr = 1, Mr[0] specifies absolute value of the remanent magnetization; the direction of the easy magnetisation axis is set up by the magnetization vector in the object to which the material is applied (the magnetization vector is specified at the object creation).
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatLin(int* mat, double* Ksi, double* Mr, int nMr);

/** Creates an isotropic linear magnetic material with single susceptibility.
@param mat [out] reference number of the material created
@param ksi [in] magnetic susceptibility (same in all directions)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatLinIso(int* mat, double ksi);

/** Creates an anisotropic linear magnetic material with easy axis.
@param mat [out] reference number of the material created
@param Ksi [in] array of 2 numbers specifying parallel and perpendicular susceptibilities [ksi_par, ksi_perp]
@param EasyAxis [in] array of 3 numbers specifying easy magnetization axis direction [ex, ey, ez]
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatLinAniso(int* mat, double* Ksi, double* EasyAxis);

/** Creates a permanent magnet material with demagnetization curve (Br/Hc model).
@param mat [out] reference number of the material created
@param Br [in] residual flux density [T]
@param Hc [in] coercivity [A/m]
@param MagAxis [in] array of 3 numbers specifying easy magnetization axis direction [mx, my, mz]
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatPM(int* mat, double Br, double Hc, double* MagAxis);

/** Creates a pre-defined magnetic material.
The material is identified by its name/formula (e.g. \"NdFeB\"). 
@param mat [out] reference number of the material created
@param id [in] null-terminated string identifying the material
@param m [in] amplitude of the remanent magnetization (for permanent-magnet type materials)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadMatStd(int* mat, char* id, double m);

/** Computes magnetization from magnetic field strength vector.
@param M [out] array of magnetization components calculated
@param nM [out] number of magnetization components calculated (length of the array M)
@param obj [in] reference number of a material or of an object with material applied
@param id [in] string specifying the magnetization components to be calculated (e.g. \"mz\")
@param H [in] magnetic field strength vector in Tesla (mu0*H)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatMvsH(double* M, int* nM, int obj, char* id, double* H);

/** Creates a nonlinear isotropic magnetic material with the magnetization magnitude equal M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3), where H is the magnitude of the magnetic field strength vector (in Tesla).
@param mat [out] reference number of the material created
@param KsiMs1 [in] array of 2 numbers specifying the parameters ms1 and ksi1 (see the formula above)
@param KsiMs2 [in] array of 2 numbers specifying the parameters ms2 and ksi2 (see the formula above)
@param KsiMs3 [in] array of 2 numbers specifying the parameters ms3 and ksi3 (see the formula above)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatSatIsoFrm(int* mat, double* KsiMs1, double* KsiMs2, double* KsiMs3);

/** Creates a nonlinear isotropic magnetic material with the M versus H curve defined by the list of pairs corresponding values of H and M (H1,M1,H2,M2,...).
@param mat [out] reference number of the material created
@param MatData [in] flat array of material data points (H1,M1,H2,M2,H3,M3,...)
@param np [in] number of material data points
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatSatIsoTab(int* mat, double* MatData, int np);

/** Creates an energy-based vector hysteresis material with table-based shape functions.
K play operators with tabulated f_k(r) = U_k'(r) shape functions.
@param mat [out] reference number of the material created
@param K [in] number of partial polarizations (play operators)
@param chi [in] array of K pinning strengths chi_k in A/m
@param r_flat [in] concatenated r grid points for all K operators
@param f_flat [in] concatenated f_k values for all K operators
@param table_sizes [in] array of K integers giving the size of each table
@param eps [in] regularization parameter for smoothed norm (default 1e-8)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
*/
EXP int CALL RadMatEnergyHysteresis(int* mat, int K, double* chi,
	double* r_flat, double* f_flat, int* table_sizes, double eps);

/** Creates a direct B-input play hysteresis material with K play operators
and tabulated f_k(r) shape functions (sign-unconstrained).
Forward: B->H is O(K) direct evaluation, Inverse: H->B uses Newton.
@param mat [out] reference number of the material created
@param K [in] number of play operators
@param eta [in] array of K play thresholds [Tesla]
@param r_flat [in] concatenated |p| grid points for all K operators
@param f_flat [in] concatenated f_k values for all K operators
@param table_sizes [in] array of K integers giving the size of each table
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
*/
EXP int CALL RadMatPlayHysteresis(int* mat, int K, double* eta,
	double* r_flat, double* f_flat, int* table_sizes);

/** Saves the internal state of an energy hysteresis material to an array.
Call with pState=NULL first to get the required array size in *pLen.
@param mat [in] reference number of the material
@param pState [out] array to receive state (or NULL to query size)
@param pLen [out] number of doubles in the state array
@return integer error code (0 : no error, <0 : error)
*/
EXP int CALL RadMatHysSaveState(int mat, double* pState, int* pLen);
EXP int CALL RadMatHysGetNuRev(int mat, double* pNuRev);
EXP int CALL RadMatHysIrreversible(int mat, double* pB, double* pHirr);

/** Restores the internal state of an energy hysteresis material from an array.
@param mat [in] reference number of the material
@param pState [in] array of state values (from RadMatHysSaveState)
@param Len [in] number of doubles in the state array
@return integer error code (0 : no error, <0 : error)
*/
EXP int CALL RadMatHysRestoreState(int mat, const double* pState, int Len);

/** Commits the current state as the reference for the next time step.
Call after Picard convergence before moving to the next NI step.
@param mat [in] reference number of the material
@return integer error code (0 : no error, <0 : error)
*/
EXP int CALL RadMatHysCommitState(int mat);

/** Creates a fixed magnetization permanent magnet material (no demagnetization).
The magnetization is constant and does not change during Solve.
@param mat [out] reference number of the material created
@param Magn [in] array of 3 numbers specifying the magnetization vector [Mx, My, Mz] in A/m
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
*/
EXP int CALL RadMatMagFixed(int* mat, double* Magn);

/** Creates a linear demagnetization permanent magnet material (Br/Hc model).
Note: Current implementation behaves as fixed magnetization (no demagnetization).
@param mat [out] reference number of the material created
@param Br [in] residual flux density in Tesla
@param Hc [in] coercivity in A/m
@param MagAxis [in] array of 3 numbers specifying the easy magnetization axis direction
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
*/
EXP int CALL RadMatMagLinear(int* mat, double Br, double Hc, double* MagAxis);

/** Creates a user-defined demagnetization curve permanent magnet material.
Note: Current implementation behaves as fixed magnetization (no demagnetization).
@param mat [out] reference number of the material created
@param CurveData [in] flat array of B-H curve data points (H1,B1,H2,B2,...) where H is in A/m and B is in Tesla
@param np [in] number of data points in the curve
@param MagAxis [in] array of 3 numbers specifying the easy magnetization axis direction
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
*/
EXP int CALL RadMatMagCurve(int* mat, double* CurveData, int np, double* MagAxis);

/** Creates laminated nonlinear anisotropic magnetic material with packing factor p and the lamination planes perpendicular to the vector N. The magnetization magnitude vs magnetic field strength for the corresponding isotropic material is defined by the formula M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3), where H is the magnitude of the magnetic field strength vector (in Tesla); ksi1, ms1, ksi2, ms2, ksi3, ms3 constants are given by parameters KsiMs1, KsiMs2, KsiMs3.
@param mat [out] reference number of the material created
@param KsiMs1 [in] array of 2 numbers specifying the parameters ms1 and ksi1 (see the formula above)
@param KsiMs2 [in] array of 2 numbers specifying the parameters ms2 and ksi2 (see the formula above)
@param KsiMs3 [in] array of 2 numbers specifying the parameters ms3 and ksi3 (see the formula above)
@param p [in] lamination stacking factor
@param N [in] array of 3 numbers specifying cartesian coordinates of a vector normal to the lamination planes; if the pointer N is 0, the lamination planes are assumed to be perpendicular to the magnetization vector in the object to which the material is applied (the magnetization vector should be specified at the object creation)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatSatLamFrm(int* mat, double* KsiMs1, double* KsiMs2, double* KsiMs3, double p, double* N);

/** Creates laminated nonlinear anisotropic magnetic material with packing factor p and the lamination planes perpendicular to the vector N. The magnetization magnitude vs magnetic field strength for the corresponding isotropic material is defined by pairs of values H, M in Tesla.
@param mat [out] reference number of the material created
@param MatData [in] flat array of material data points (H1,M1,H2,M2,H3,M3,...)
@param np [in] number of material data points
@param p [in] lamination stacking factor
@param N [in] array of 3 numbers specifying cartesian coordinates of a vector normal to the lamination planes; if the pointer N is 0, the lamination planes are assumed to be perpendicular to the magnetization vector in the object to which the material is applied (the magnetization vector should be specified at the object creation)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatSatLamTab(int* n, double* MatData, int np, double p, double* N);

/** Creates a nonlinear anisotropic magnetic material. The magnetization vector component parallel to the easy axis is computed either by the formula: ms1*tanh(ksi1*(hpa-hc1)/ms1) + ms2*tanh(ksi2*(hpa-hc2)/ms2) + ms3*tanh(ksi3*(hpa-hc3)/ms3) + ksi0*(hpa-hc0), where hpa is the field strength vector component parallel to the easy axis, or by ksi0*hpa. The magnetization vector component perpendicular to the easy axis is computed either by the formula: ms1*tanh(ksi1*hpe/ms1) + ms2*tanh(ksi2*hpe/ms2) + ms3*tanh(ksi3*hpe/ms3) + ksi0*hpe, where hpe is the field strength vector component perpendicular to the easy axis, or by ksi0*hpe. At least one of the magnetization components should non-linearly depend on the field strength. The direction of the easy magnetisation axis is set up by the magnetization vector in the object to which the material is later applied.
@param mat [out] reference number of the material created
@param DataPar [in] flat array of constants defining the magnetic material behavior in the direction parallel to the easy magnetization axis. It can be {ksi1,ms1,hc1,ksi2,ms2,hc2,ksi3,ms3,hc3,ksi0,hc0} or {ksi0}. 
@param nDataPar [in] length of the array DataPar. Can be equal to 11 or 1, for the DataPar to be interpreted as {ksi1,ms1,hc1,ksi2,ms2,hc2,ksi3,ms3,hc3,ksi0,hc0} or {ksi0}.
@param DataPer [in] flat array of constants defining the magnetic material behavior in the direction perpendicular to the easy magnetization axis. It can be {ksi1,ms1,ksi2,ms2,ksi3,ms3,ksi0} or {ksi0}. 
@param nDataPer [in] length of the array DataPer. Can be equal to 7 or 1, for the DataPer to be interpreted as {ksi1,ms1,ksi2,ms2,ksi3,ms3,ksi0} or {ksi0}.
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadMatSatAniso(int* mat, double* DataPar, int nDataPar, double* DataPer, int nDataPer);

// RadRlxPre REMOVED (2026-01-31) - Use RadBuildMatrix instead
// The new API is: int handle = RadBuildMatrix(obj, image);
// where image is "+x", "-z", "+x-z", etc. for Image symmetry

/** Executes manual relaxation procedure for the interaction matrix intrc.
@param D [out] an array of four numbers specifying: [0] average absolute change in magnetization after previous iteration over all the objects participating in the relaxation, [1] maximum absolute value of magnetization over all the objects participating in the relaxation, [2] maximum absolute value of magnetic field strength over central points of all the objects participating in the relaxation, and [3] actual number of iterations done. The values [0]-[2] are those of last iteration.
@param n [out] length of array D
@param intrc [in] an integer number referencing the interaction matrix
@param meth [in] an integer number specifying the method of relaxation to be used
@param iter [in] an integer number specifying number of iterations to be made
@param rlxpar [in] a floating point number between 0 and 1 specifying the relaxation parameter
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadRlxMan(double* D, int* n, int intrc, int meth, int iter, double rlxpar);

/** Executes automatic relaxation procedure for the interaction matrix intrc.
The relaxation stops whenever the change of magnetization (averaged over all sub-elements) between two successive iterations is smaller than prec or the number of iterations is larger than iter.
@param D [out] an array of four numbers specifying: [0] average absolute change in magnetization after previous iteration over all the objects participating in the relaxation, [1] maximum absolute value of magnetization over all the objects participating in the relaxation, [2] maximum absolute value of magnetic field strength over central points of all the objects participating in the relaxation, and [3] actual number of iterations done. The values [0]-[2] are those of last iteration.
@param n [out] length of array D
@param intrc [in] an integer number referencing the interaction matrix
@param prec [in] a real number specifying an absolute precision value for magnetization (in Tesla), to be reached by the end of the relaxation
@param iter [in] maximum number of iterations permitted to reach the specified precision
@param meth [in] an integer number specifying the method of relaxation to be used (values 0, 3, 4, 5, 8, 9, 10; 0 means default method = 10 BiCGSTAB)
@param opt [in] pointer to an option string, which can be "ResetM->True" (default) or "ResetM->False"
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadRlxAuto(double* D, int* n, int intrc, double prec, int iter, int meth, const char* opt);

/** Updates external field data for the relaxation (to take into account e.g. modification of currents in coils, if any) without rebuilding the interaction matrix.
@param intrc [in] an integer number referencing the interaction object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadRlxUpdSrc(int intrc);

/** Builds an interaction matrix and performs a relaxation procedure. 
The relaxation stops whenever the change of magnetization (averaged over all sub-elements) between two successive iterations is smaller than prec or the number of iterations is larger than iter. The interaction matrix is deleted. 
@param D [out] an array of four numbers specifying: [0] average absolute change in magnetization after previous iteration over all the objects participating in the relaxation, [1] maximum absolute value of magnetization over all the objects participating in the relaxation, [2] maximum absolute value of magnetic field strength over central points of all the objects participating in the relaxation, and [3] actual number of iterations done. The values [0]-[2] are those of last iteration.
@param n [out] length of array D
@param obj [in] an integer number specifying the object to solve for magnetization
@param prec [in] a real number specifying an absolute precision value for magnetization (in Tesla), to be reached by the end of the relaxation
@param iter [in] maximum number of iterations permitted to reach the specified precision
@param meth [in] solver method: 0=LU, 1=BiCGSTAB, 2=HACApK
@param image [in] image symmetry string (e.g., "+x", "-z", "+x-z") or nullptr for no symmetry
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author P.E., O.C.
*/
EXP int CALL RadSolve(double* D, int* n, int obj, double prec, int iter, int meth, const char* image);

/** Builds an interaction matrix and performs a relaxation procedure with nonlinear method selection.
Similar to RadSolve but with additional nonlinear_method parameter for selecting convergence criterion.
@param D [out] an array of four numbers specifying: [0] residual, [1] max M, [2] max H, [3] iterations
@param n [out] length of array D
@param obj [in] an integer number specifying the object to solve for magnetization
@param prec [in] a real number specifying precision value for convergence
@param iter [in] maximum number of iterations
@param meth [in] linear solver method: 0=LU, 1=BiCGSTAB
@param nonl_method [in] nonlinear convergence method: 0=mucal1 (chi-change), 1=mucal2 (B-change/Newton)
@param image [in] image symmetry string (e.g., "+x", "-z", "+x-z") or nullptr for no symmetry
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author Radia Development Team
*/
EXP int CALL RadSolveNonl(double* D, int* n, int obj, double prec, int iter, int meth, int nonl_method, const char* image);

/** Builds interaction matrix for magnetostatic problem without solving.
Allows inspection of the matrix before solving. The matrix is cached for subsequent Solve() calls.
@param n [out] reference number of the interaction matrix (for GetInteractMatrix)
@param ElemKey [in] reference number of the object to build matrix for
@param image [in] image symmetry string (e.g., "+x", "-z", "+x-z") or nullptr for no symmetry
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@note This replaces the deprecated PreRelax + Image + BuildImageMatrix workflow
@author Radia Development Team
*/
EXP int CALL RadBuildMatrix(int* n, int ElemKey, const char* image);

/** Returns solve statistics from the last Solve() call.
@param dOut [out] array of statistics: [t_matrix_build, t_linear_solve, linear_iterations, nonl_iterations, openmp_enabled, openmp_max_threads]
@param nOut [out] number of statistics returned (0 if no solve has been performed)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author Radia Development Team
*/
EXP int CALL RadGetSolveStats(double* dOut, int* nOut);

/** Computes magnetic field created by the object obj at one or many points.
@param B [out] flat array of all computed values of the magnetic field components (should be allocated by calling function)
@param nB [out] total number of calculated magnetic field component values
@param obj [in] reference number of the magnetic field source object
@param id [in] string identifying magnetic field components to be computed
@param Coords [in] flat array of coordinates of all points where the field should be computed (x1,y1,z1,x2,y2,z2,...)
@param np [in] number of points where the magnetic field should be calculated (the length of the array Coords is equal to 3*np)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFld(double* B, int* nB, int obj, char* id, double* Coords, int np);

/** Sets general absolute accuracy levels for computation of magnetic field induction (prcB), vector potential (prcA), induction integrals along straight line (prcBint), field force (prcFrc), relativistic particle trajectory coordinates (prcTrjCrd) and angles (prcTrjAng).
@param n [out] dummy
@param prcB [in] absolute accuracy level for magnetic field induction [T]
@param prcA [in] absolute accuracy level for vector potential
@param prcBInt [in] absolute accuracy level for magnetic field induction integrals along straight line [T*mm]
@param prcFrc [in] absolute accuracy level for the force [N]
@param prcTrjCrd [in] absolute accuracy level for particle trajectory coordinate
@param prcTrjAng [in] absolute accuracy level for particle trajectory angle
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldCmpCrt(int* n, double prcB, double prcA, double prcBInt, double prcFrc, double prcTrjCrd, double prcTrjAng);

/** Sets general absolute accuracy levels for computation of magnetic field induction (PrcB), vector potential (PrcA), induction integral along straight line (PrcBInt), field force (PrcForce), torque (PrcTorque), energy (PrcEnergy); relativistic charged particle trajectory coordinates (PrcCoord) and angles (PrcAngle). The function works according to the mechanism of string options. The name(s) of the option(s) should be: PrcB, PrcA, PrcBInt, PrcForce, PrcTorque, PrcEnergy, PrcCoord, PrcAngle.
@param n [out] dummy
@param opt [in] pointer to an option string, which can be "PrcB->..." or "PrcA->..." or "PrcBInt->..." or "PrcForce->..." or "PrcTorque->..." or "PrcEnergy->..." or "PrcCoord->..." or "PrcAngle->...", where "..." should be replaced by the appropriate absolute accuracy level.
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldCmpPrc(int* n, char* opt);

/** Switches on or off the randomization of all the length values. The randomization magnitude can be set by the function radFldLenTol.
@param n [out] dummy
@param OnOrOff [in] string containing either "on" or "off"
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldLenRndSw(int* n, char* OnOrOff);

/** Sets absolute and relative randomization magnitudes for all the length values, including coordinates and dimensions of the objects producing magnetic field, and coordinates of points where the field is computed. Optimal values of the variables can be: RelVal=10^(-11), AbsVal=L*RelVal, ZeroVal=AbsVal, where L is the distance scale value (in mm) for the problem to be solved. Too small randomization magnitudes can result in run-time code errors.
@param n [out] dummy
@param AbsVal [in] absolute position/length randomization magnitude [mm]
@param RelVal [in] relative position/length randomization magnitude
@param ZeroVal [in] absolute zero tolerance [mm]
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldLenTol(int* n, double AbsVal, double RelVal, double ZeroVal);

/* RadFldEnr / RadFldEnrFrc / RadFldEnrTrq (energy-based force/torque/energy API) REMOVED (Phase C, 2026-04-16) */

/** Computes force of the field produced by the object obj into a shape defined by shape. shape can be the result of RadObjRecMag (parallelepiped) or RadFldFrcShpRtg (rectangular surface). This function uses the algorithm based on Maxwell tensor, which may not always provide high efficiency.
@param f [out] computed force component(s)
@param nf [out] number of force components computed
@param obj [in] reference number of the magnetic field source object
@param shape [in] reference number of the shape object
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldFrc(double* f, int* nf, int obj, int shape);

/** Computes "focusing potential" for trajectory of relativistic charged particle in magnetic field produced by the object obj. The integration is made from P1 to P2 with np equidistant points.
@param d [out] computed focusing potential value
@param obj [in] reference number of the magnetic field source object
@param P1 [in] array of 3 real numbers specifying cartesian coordinates of an edge point of the integration segment
@param P2 [in] array of 3 real numbers specifying cartesian coordinates of an edge point of the integration segment
@param np [in] number of points for the integration
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/

/** Computes magnetic field integral produced by the object obj along a straight line specified by points P1 and P2; depending on the InfOrFin variable value, the integral is infinite ("inf") or finite ("fin"), from P1 to P2; the field integral component is specified by the id input variable. The unit is T*mm.
@param f [out] computed field integral component(s)
@param nf [out] number of field integral components computed
@param obj [in] reference number of the object creating the magnetic field
@param InfOrFin [in] string specifying the type of field integral: finite ("fin") or infinite ("inf")
@param id [in] string identifying the field integral components to be computed (ibx|iby|ibz)
@param P1 [in] array of 3 real numbers specifying cartesian coordinates of a point on the integration line
@param P2 [in] array of 3 real numbers specifying cartesian coordinates of another point on the integration line
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldInt(double* f, int* nf, int obj, char* InfOrFin, char* id, double* P1, double* P2);

/** Computes magnetic field created by object obj in np equidistant points along a line segment from P1 to P2; the field component is specified by the id input variable.
@param B [out] computed field value(s)
@param nB [out] number of field values computed
@param obj [in] reference number of the object creating the magnetic field
@param id [in] string identifying magnetic field components to be computed
@param P1 [in] array of 3 real numbers specifying cartesian coordinates of an edge point of the line segment
@param P2 [in] array of 3 real numbers specifying cartesian coordinates of another edge point of the line segment
@param np [in] number of points where the magnetic field should be calculated
@param ArgOrNoArg [in] string specifying whether or not to output a longitudinal position for each point where the field is computed ("arg|noarg")
@param start [in] start value for the longitudinal position
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldLst(double* B, int* nB, int obj, char* id, double* P1, double* P2, int np, char* ArgOrNoArg, double start);



/** Computes a virtual "shim signature", i.e. variation of a given magnetic field component introduced by given displacement of magnetic field source object.
@param f [out] computed array of field component values
@param nf [out] number of values in the field component array computed
@param obj [in] reference number of the magnetic field source object
@param id [in] string identifying magnetic field component to be computed (can be e.g. "bx", "bz", "ix", "iz",...)
@param V [in] array of 3 real numbers specifying cartesian coordinates of the field source object displacement
@param P1 [in] array of 3 real numbers specifying cartesian coordinates of an edge point of the line segment along which the "shim signature" should be computed
@param P2 [in] array of 3 real numbers specifying cartesian coordinates of another edge point of the line segment
@param np [in] number of points where the magnetic field should be computed
@param Vi [in] array of 3 real numbers specifying cartesian coordinates of a vector defining the integration line (is taken into account only if id string specifies field integral component)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@authors: O.C.
*/
EXP int CALL RadFldShimSig(double* f, int* nf, int obj, char* id, double* V, double* P1, double* P2, int np, double* Vi);

/** Creates a rectangle with central point P and dimensions W (to be used for force computation via Maxwell tensor).
@param n [out] reference number of the object created
@param P [in] array of 3 real numbers specifying cartesian coordinates of the center point
@param W [in] array of 2 real numbers specifying dimensions of the rectangle
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadFldFrcShpRtg(int* n, double* P, double* W);

/** Deletes object obj.
@param n [out] dummy
@param obj [in] reference number of the object to be deleted
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadUtiDel(int* n, int obj);

/** Deletes all previously created objects.
@param n [out] dummy
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadUtiDelAll(int* n);

// RadUtiDmp / RadUtiDmpPrs REMOVED (Phase B1, 2026-04-15) -
// .rad save/load is no longer supported.

/** Sets interruption time quanta in seconds for platforms with no preemptive multitasking.
@param t [in] interruption time quanta [s]
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author P.E., O.C.
*/
EXP int CALL RadUtiIntrptTim(double* d, double t);

/** Returns data resulting from previous calculations in cases when the data size was not known 'a priori', e.g. after executing functions RadObjM, RadObjCenFld, ...
@param size [out] pointer to the resulting data (to be allocated in calling function)
@param typeData [in] string identifying type of the data: \"mad\" for multi-dim. array of double, \"mai\" for multi-dim. array of integer, \"bin\" for byte array, \"asc\" for ASCII string, \"d\" for double, \"i\" for integer
@param key [in] additional identifier of the data to be extracted, e.g. to ensure thread safety (not implemented yet)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/
EXP int CALL RadUtiDataGet(char* pcData, const char typeData[3], long key=0); //OC04102018
//EXP int CALL RadUtiDataGet(char* pcData, char typeData[3], long key=0); //OC27092018
//EXP int CALL RadUtiDataGet(double* pData, long key); //OC15092018

/** Identifies the version number of the Radia DLL.
@param d [out] version number
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author P.E., O.C.
*/ 
EXP int CALL RadUtiVer(double* d);

/** initializes or finalizes the Message Passing Inteface (MPI) for parallel calculations and returns list of basic MPI process parameters (in the case of initialization): rank of a process and total number of processes.
@param arPar [out] array of basic MPI parameters: rank of a process [0] and total number of processes [1]
@param OnOrOff [in] string containing either "on" or "off"
@param arData [in] array of data to be shared among different processes
@param pnData [in/out] pointer to length of array of data to be shared among different processes (for processes receiving data this param will be set upon function return)
@param pRankFrom [in/out] pointer to rank of process to take data (to be shered) from (the default value means take the from process with rank 0)
@param rankTo [in/out] pointer to rank of process(es) to send the data to (the default value means sharing the data anomg all processes)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author O.C.
*/ 
EXP int CALL RadUtiMPI(int* arPar, char* OnOrOff, double* arData=0, long* pnData=0, long* pRankFrom=0, long* pRankTo=0); //OC19032020
//EXP int CALL RadUtiMPI(int* arPar, char* OnOrOff);

EXP int CALL RadUtiYeldFuncSet(int (*pExtFunc)());

#ifdef RADIA_USE_HACAPK
/** Sets HACApK (H-matrix) parameters for BiCGSTAB solver with H-matrix acceleration.
These parameters control the H-matrix construction and compression.
Must be called before RadSolve with method=2 (BiCGSTAB+HACApK).
@param n [out] dummy output
@param eps [in] ACA+ compression tolerance (default: 1e-4, use 1e-8 for high accuracy)
@param leaf_size [in] minimum cluster size in elements (default: 32, ELF uses 10)
@param eta [in] admissibility parameter (default: 2.0)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author Radia Development Team
*/
EXP int CALL RadSetHACApKParams(int* n, double eps, int leaf_size, double eta);


/** Sets only the H-matrix ACA epsilon (tolerance) for HACApK solver.
ELF-compatible API: magic.set_hmatrix_epsilon(eps)
@param n [out] dummy output
@param eps [in] ACA+ compression tolerance (default: 1e-4, use 1e-8 for high accuracy)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author Radia Development Team
*/
EXP int CALL RadSetHMatrixEpsilon(int* n, double eps);

/** Gets HACApK (H-matrix) statistics after a solve with method=2.
Returns information about the H-matrix structure and performance.
@param dOut [out] array of 7 doubles: [0] n_lowrank, [1] n_dense, [2] max_rank,
                  [3] n_leaves, [4] n_dof, [5] compression_ratio, [6] build_time_sec
@param nOut [out] number of values written to dOut (7)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author Radia Development Team
*/
EXP int CALL RadGetHACApKStats(double* dOut, int* nOut);
#endif

/** Sets BiCGSTAB inner loop tolerance for iterative solvers (Method 1 and 2).
Default: 1e-4 (ELF-compatible). Lower values give higher accuracy but slower convergence.
@param n [out] dummy output
@param tol [in] BiCGSTAB tolerance (default: 1e-4)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author Radia Development Team
*/
EXP int CALL RadSetBiCGSTABTol(int* n, double tol);

/** Gets current BiCGSTAB inner loop tolerance.
@param tol [out] current BiCGSTAB tolerance
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author Radia Development Team
*/
EXP int CALL RadGetBiCGSTABTol(double* tol);

/** Select multipole-moment method-2 Krylov solver: 0=BiCGSTAB, 1=restarted GMRES. */
EXP int CALL RadSetMomentKrylovSolver(int* n, int solver);
EXP int CALL RadGetMomentKrylovSolver(int* solver);

/** Configure restarted GMRES restart length for multipole-moment method-2. */
EXP int CALL RadSetMomentGMRESRestart(int* n, int restart);
EXP int CALL RadGetMomentGMRESRestart(int* restart);

/** Configure safeguarded Anderson acceleration depth for multipole-moment method-2 (0=off, 1=depth-1). */
EXP int CALL RadSetMomentAndersonDepth(int* n, int depth);
EXP int CALL RadGetMomentAndersonDepth(int* depth);

/** Sets under-relaxation coefficient for nonlinear iteration.
@param n [out] dummy output (set to 1)
@param relax [in] relaxation coefficient (0.0 = full step, 0.0-1.0 = under-relaxation)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author Radia Development Team

Formula: chi_new = chi_new*(1-relax) + chi_old*relax
- relax=0.0: full Newton step (default, fastest convergence)
- relax>0.0: damped Newton step (improves stability for difficult cases)
*/
EXP int CALL RadSetRelaxParam(int* n, double relax);

/** Gets current under-relaxation coefficient.
@param relax [out] current relaxation coefficient
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
@author Radia Development Team
*/
EXP int CALL RadGetRelaxParam(double* relax);

/** Enables/disables Newton-Raphson nonlinear iteration.
When enabled, uses differential susceptibility chi_d = (dB/dH)/mu_0 - 1 for the system matrix
and adds a correction term to the RHS for quadratic convergence.
@param use_newton [in] 1 to enable Newton, 0 for Picard (default)
*/
/** Keep magnetization from previous Solve (for Newton->Picard workflow).
@param n [out] always 1
@param keep [in] 1 to keep, 0 to reset (default)
*/
EXP int CALL RadSetKeepMagnetization(int* n, int keep);
EXP int CALL RadGetKeepMagnetization(int* keep);

EXP int CALL RadSetNewtonMethod(int* n, int use_newton);

/** Gets current Newton method setting.
@param use_newton [out] 1 if Newton is enabled, 0 if Picard
*/
EXP int CALL RadGetNewtonMethod(int* use_newton);

/** Configure Newton-Raphson line search damping.
@param n [out] Number of elements modified (always 1)
@param enabled [in] 1 to enable damping, 0 to disable
@param max_iter [in] Max line search iterations (default: 5)
@param min_omega [in] Minimum omega threshold (default: 0.01)
*/
EXP int CALL RadSetNewtonDamping(int* n, int enabled, int max_iter, double min_omega);

/** Get Newton line search damping configuration.
@param enabled [out] 1 if damping enabled, 0 otherwise
@param max_iter [out] Max line search iterations
@param min_omega [out] Minimum omega threshold
*/
EXP int CALL RadGetNewtonDampingStats(int* enabled, int* max_iter, double* min_omega);

/** Enable/disable B-input Newton-Raphson for energy-based hysteresis.
@param n [out] always 1
@param enabled [in] 1 to enable B-input Newton, 0 to disable
*/
EXP int CALL RadSetBInputNewton(int* n, int enabled);

/** Get B-input Newton setting.
@param enabled [out] 1 if enabled, 0 otherwise
*/
EXP int CALL RadGetBInputNewton(int* enabled);

/** Enable/disable B-input Hantila polarization method for energy-based hysteresis.
@param n [out] always 1
@param enabled [in] 1 to enable, 0 to disable
*/
EXP int CALL RadSetBInputHantila(int* n, int enabled);

/** Get B-input Hantila setting.
@param enabled [out] 1 if enabled, 0 otherwise
*/
EXP int CALL RadGetBInputHantila(int* enabled);

/** Set Hantila polarization parameter alpha.
@param n [out] always 1
@param alpha [in] polarization parameter (0 = auto-compute)
*/
EXP int CALL RadSetHantilaAlpha(int* n, double alpha);

/** Get Hantila alpha parameter.
@param alpha [out] current alpha value
*/
EXP int CALL RadGetHantilaAlpha(double* alpha);

/** Set Hantila under-relaxation parameter.
@param n [out] always 1
@param relax [in] under-relaxation (0 = full step)
*/
EXP int CALL RadSetHantilaRelax(int* n, double relax);

/** Get Hantila relaxation parameter.
@param relax [out] current relax value
*/
EXP int CALL RadGetHantilaRelax(double* relax);

// RadSetIMASymmetry REMOVED (2026-01-31) - Use RadBuildMatrix(obj, image) instead
// The new unified API handles both interaction creation and Image symmetry
// Example: int handle = RadBuildMatrix(obj, "+x-z");

// RadBuildIMAMatrix REMOVED (2026-01-31) - Use RadBuildMatrix(obj, image) instead
// The new unified API builds the Image matrix in a single call

/** Classifies evaluation points as inside/near/far relative to mesh elements.
@param classification [out] array of classification (0=inside, 1=near, 2=far)
@param nearest_elem [out] array of nearest element indices
@param n_points [in] number of evaluation points
@param points [in] evaluation point coordinates (n_points * 3)
@param container_handle [in] Radia container handle
@param near_threshold [in] near zone multiplier (typically 3.0)
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
*/
EXP int CALL RadClassifyPoints(int* classification, int* nearest_elem, int n_points,
                               double* points, int container_handle, double near_threshold);

/** Computes B and H fields at multiple points (batch computation).
@param B_out [out] B field values (n_points * 3)
@param H_out [out] H field values (n_points * 3)
@param n_points [in] number of evaluation points
@param points [in] evaluation point coordinates (n_points * 3)
@param container_handle [in] Radia container handle
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
*/
EXP int CALL RadFldBatch(double* B_out, double* H_out, int n_points,
                         double* points, int container_handle);

/** Computes magnetic scalar potential at multiple points.
@param phi_out [out] scalar potential values (n_points)
@param n_points [in] number of evaluation points
@param points [in] evaluation point coordinates (n_points * 3)
@param container_handle [in] Radia container handle
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
*/
EXP int CALL RadFldPhi(double* phi_out, int n_points, double* points, int container_handle);

/** Computes magnetic vector potential at multiple points.
@param A_out [out] vector potential values (n_points * 3)
@param n_points [in] number of evaluation points
@param points [in] evaluation point coordinates (n_points * 3)
@param container_handle [in] Radia container handle
@return integer error code (0 : no error, >0 : error number, <0 : warning number)
*/
EXP int CALL RadFldA(double* A_out, int n_points, double* points, int container_handle);


// Replaced by Python-based PEEC topology solver (peec_topology.py, peec_coupled.py)
// and FastHenry .inp parser (fasthenry_parser.py).
// See CLAUDE.md for the new PEEC API reference.




#ifdef __cplusplus
}
#endif
