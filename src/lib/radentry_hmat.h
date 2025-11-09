/*-------------------------------------------------------------------------
*
* File name:      radentry_hmat.h
*
* Project:        RADIA
*
* Description:    H-matrix field evaluation API declarations
*
* Author(s):      Radia Development Team
*
* First release:  2025
*
-------------------------------------------------------------------------*/

#ifndef __RADENTRY_HMAT_H
#define __RADENTRY_HMAT_H


#ifdef __cplusplus
extern "C" {
#endif

/** Batch field evaluation with optional H-matrix acceleration
*
* Computes magnetic field at multiple observation points in batch mode.
* Significantly faster than calling RadFld repeatedly for many points.
*
* @param B [out] output array of field values (size: 3*np)
* @param nB [out] number of field components returned (should be 3*np)
* @param obj [in] object index
* @param id [in] field type: "bx|by|bz|hx|hy|hz|b|h" etc.
* @param Coords [in] observation points, flattened array [x1,y1,z1,x2,y2,z2,...]
* @param np [in] number of observation points
* @param use_hmatrix [in] 0=direct calculation, 1=H-matrix (if available)
* @return integer error code (0: no error, >0: error number, <0: warning number)
*/
EXP int CALL RadFldBatch(double* B, int* nB, int obj, char* id, double* Coords, int np, int use_hmatrix);

/** Enable or disable H-matrix acceleration for field evaluation
*
* Global setting that affects subsequent field calculations.
* Enables caching for non-linear iterations.
*
* @param enabled [in] 0=disable, 1=enable H-matrix acceleration
* @param tol [in] HACApK ACA tolerance (smaller = more accurate but slower) (default: 1e-6)
* @return integer error code (0: no error, >0: error number, <0: warning number)
*/
EXP int CALL RadSetHMatrixFieldEval(int enabled, double tol);

/** Clear H-matrix field evaluation cache
*
* Frees memory used by H-matrix cache.
* Should be called after geometry modifications.
*
* @return integer error code (0: no error, >0: error number, <0: warning number)
*/
EXP int CALL RadClearHMatrixCache(void);

/** Get H-matrix field evaluation statistics
*
* Returns information about H-matrix cache usage.
*
* @param stats [out] array of statistics [is_enabled, num_cached, total_memory_MB, ...]
* @param nstats [out] number of statistics returned
* @return integer error code (0: no error, >0: error number, <0: warning number)
*/
EXP int CALL RadGetHMatrixStats(double* stats, int* nstats);

/** Update magnetization without rebuilding H-matrix
*
* Fast update for non-linear relaxation: updates magnetic moments while keeping
* geometry (and H-matrix) fixed. Much faster than full rebuild.
*
* @param obj [in] object index (must have existing H-matrix)
* @return integer error code (0: no error, -1: error, -2: H-matrix not built)
*/
EXP int CALL RadUpdateHMatrixMagnetization(int obj);

#ifdef __cplusplus
}
#endif

#endif  // __RADENTRY_HMAT_H
