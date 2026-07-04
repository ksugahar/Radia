/*-------------------------------------------------------------------------
*
* File name:      rad_hacapk_field.h
*
* Project:        RADIA
*
* Description:    H-matrix-accelerated field evaluation (rad.Fld) via the
*                 EMBED-IN-SQUARE trick on the square-only RadHACApKBase.
*
*                 The field at N_obs observation points from N_src source
*                 elements is  field[obs] = sum_src G(obs,src) . M[src],  with
*                 G_ab(obs,src) = a-component of the field (B flux density, or A
*                 vector potential) at obs from src carrying unit magnetization
*                 in direction b (the field is LINEAR in M, so this is the exact
*                 response, NOT a dipole approximation).
*
*                 RadHACApKBase is SQUARE-ONLY (one coordinate set, one
*                 cluster tree, m_ndof x m_ndof).  An obs x src field operator
*                 is rectangular, so we EMBED it in a square SYMMETRIC operator
*                 over the COMBINED [obs; src] point set:
*
*                     A = [[ 0    K  ]     K[3o+a, 3s+b] = G_ab(obs_o, src_s)
*                          [ K^T  0  ]]
*
*                 DOF order (uniform 3 per element): obs elements 0..N_obs-1
*                 (field slots), then src elements 0..N_src-1 (magnetization
*                 slots).  With  x = [0 ; M_flat],  y = A x = [K M ; 0] reads
*                 out y[:3N_obs] = field (the K^T block sees only the zero obs
*                 slots).  A is SYMMETRIC on purpose: HACApK's cluster-tree
*                 build / fill is only correct on the symmetric path -- the
*                 one-sided [[0,K],[0,0]] mis-fills its sub-blocks under a fine
*                 multi-cluster tree.  This reuses
*                 the whole square HACApK cluster-tree / ACA / matvec pipeline
*                 (~4x the work of a true rectangular H-matrix, but still
*                 O(N log N) << the direct O(N_obs*N_src)).
*
*                 The field kernel is bit-consistent with rad.Fld: it calls
*                 the same radTg3d::B_genComp used by RadFieldUnified.  This is
*                 the far-field-dominated ACA regime (well-separated obs/src
*                 blocks compress to rank ~10-20); it does NOT reintroduce FMM
*                 (removed 2026-03-06) -- it is HACApK ACA on the field
*                 kernel's far blocks (a different, validated geometry class).
*
* First release:  2026-07
*
-------------------------------------------------------------------------*/

#ifndef __RAD_HACAPK_FIELD_H
#define __RAD_HACAPK_FIELD_H

#include "rad_hacapk.h"

#include <vector>
#include <string>

class radTg3dRelax;
class radTPolyhedron;

/**
 * RadHACApKFieldEval: embed-in-square H-matrix for rad.Fld-style field
 * evaluation.  Sources are the leaf magnetic elements of a Radia container
 * (resolved from a container handle); observation points are supplied as a
 * flat [x,y,z, ...] array.  field_type "b" (flux density B, Tesla) or "a" (vector potential A, T*m --
 * the analytic open-boundary source term for A-formulation eddy-current FEM coupling).  Uniform 3-DOF.
 *
 * Usage (Python, via the _FieldEvalHMatrix pybind):
 *   G  = _FieldEvalHMatrix(container_handle, obs_points, eps=...)
 *   x  = [0]*(3*G.n_obs()) + list(G.src_magnetization())
 *   y  = G.matvec(x)                       # O(N log N) H-matvec
 *   field = y[:3*G.n_obs()]                # B at the observation points
 *
 * Ground truth for validation is the direct rad.Fld(container,'b',obs);
 * G.entry(i,j) is the EXACT (uncompressed, eps-independent) kernel, so the
 * matvec-vs-entry error is a direct measure of the ACA compression accuracy
 * (mirrors validation_test/hacapk/test_compression_accuracy_vs_analytic.py).
 */
class RadHACApKFieldEval : public RadHACApKBase {
public:
    // container_handle : a Radia object (magnet / container of magnets).
    // obs_points       : flat [x0,y0,z0, x1,y1,z1, ...], length 3*N_obs.
    // field_type       : "b" (magnetic flux density, Tesla) or "a" (vector potential, T*m).
    //                    A is the source term for A-formulation eddy-current FEM coupling -- the
    //                    analytic open-boundary A supplied by the integral representation (HDiv gives
    //                    B and would need a curl-inverse to recover A).
    // image : an IMA (image method) string like "+x-z" ("" = none).  A reduced mirror-symmetry model
    //         (fundamental domain in the positive octant, mirror planes through the origin) reproduces
    //         the field of the full model: field += sum over the 2^P-1 image reflections.  Same
    //         convention as rad.Solve(image=...) (field PARALLEL to a mirror -> '+', PERPENDICULAR -> '-').
    RadHACApKFieldEval(int container_handle, std::vector<double> obs_points,
                       const std::string& field_type = "b", const std::string& image = "");
    ~RadHACApKFieldEval() override;

    int NObs() const { return m_nObs; }
    int NSrc() const { return m_nSrc; }
    bool IsAField() const { return m_isA; }   // true => vector potential A, false => flux density B
    int NImages() const { return (int)m_imgAxmask.size(); }   // number of IMA image reflections (0 if none)

    // Internal normalisation: the H-matrix stores entries scaled to O(1) (the raw B-response is
    // O(mu0*H) ~ 1e-13 for far obs, and HACApK's ACA stopping tolerance is effectively ABSOLUTE, so
    // sub-1e-8 raw entries never get enough rank).  entry()/matvec() are PHYSICAL (Tesla): the pybind
    // wrappers multiply by Scale().  physical_field = Scale() * matvec([0; M])[:3*n_obs].
    double Scale() const { return m_scale; }

    // Actual (solved / assigned) magnetization of every source element,
    // flat [Mx0,My0,Mz0, ...] length 3*N_src -- the src-slot half of the
    // embed matvec vector x = [0 (obs) ; M_flat].
    const std::vector<double>& SrcMagnetization() const { return m_srcM; }

    // The +N field-response entry: obs-row x src-col block = G_ab(obs,src);
    // 0 in the obs-obs / src-obs / src-src (embed-zero) blocks.  0-based
    // ORIGINAL DOF indices (RadHACApKCallback::ComputeEntry already un-lod's).
    double GetInteractionMatrixElement(int dof_i, int dof_j) const override;

protected:
    void ExtractCoordinates() override;
    void OnBeforeBuild() override;   // computes m_scale (sampled max |raw entry|) for O(1) normalisation
    void InitializeInvChi() override { m_inv_chi.assign(m_ndof, 0.0); }  // no 1/chi shift
    bool IsVariableDOF() const override { return false; }               // uniform 3-DOF
    int  GetUniformNFFC() const override { return 3; }

private:
    // 3x3 field-response of source s at observation point o (B or A per unit M).  Row-major G[a*3+b],
    // a = obs component, b = src component.  Reads the 3 fixed unit-M clones m_clones[3*s+b] READ-ONLY
    // (no elem->Magn mutation), so the HACApK fill runs fully PARALLEL under the caller's TaskManager
    // (the earlier mutate-under-mutex kernel serialised the whole build).
    void Compute3x3(int o, int s, double* G) const;

    int m_handle;
    int m_nObs;
    int m_nSrc;
    std::vector<double> m_obs;               // 3*N_obs observation coordinates
    std::vector<radTg3dRelax*> m_src;        // N_src source leaves (NOT owned)
    std::vector<double> m_srcM;              // 3*N_src actual magnetization
    // 3 OWNED unit-M clones per source (Magn = e_x/e_y/e_z fixed), flat [3*s + b].  radia's own safe
    // deep copy (radTPolyhedron copy ctor) -- read-only per-thread B_genComp, so no mutation race.
    std::vector<radTPolyhedron*> m_clones;
    bool m_isA;                              // true => A (vector potential, T*m); false => B (Tesla)
    double m_physScale;                      // per-field-type factor to physical units (B: mu0; A: 1)
    std::vector<int> m_imgAxmask;            // IMA: per-image axis bitmask (bit0=x, bit1=y, bit2=z)
    std::vector<double> m_imgSign;           // IMA: per-image field-eval sign (prod(plane signs) * (-1)^popcount)
    double m_scale;                          // O(1) normalisation (max sampled |raw entry|); physical = scale * stored
};

#endif // __RAD_HACAPK_FIELD_H
