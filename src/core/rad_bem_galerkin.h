/*-------------------------------------------------------------------------
*
* File name:      rad_bem_galerkin.h
*
* Project:        RADIA
*
* Description:    Fast C++ Galerkin BEM assembler for the Laplace single-
*                 layer (SL) and double-layer (DL) operators on flat or
*                 P2-curved triangular surface meshes.  Uses Sauter-Schwab
*                 Duffy quadrature for the singular pairs (identical /
*                 common-edge / common-vertex) and tensor-product Gauss
*                 for the regular pairs.  OpenMP parallelization over
*                 the (T_a, T_b) outer-pair loop.
*
* First release:  2026-04-30
*
* References:
*   - S. Sauter, C. Schwab, "Boundary Element Methods", Springer 2010,
*     Chapter 5.2 (Sauter-Schwab Duffy 4D quadrature)
*   - NGSolve.bem intrules_SauterSchwab.cpp (the SS reference-tri
*     formulas; we transcribe with attribution)
*
-------------------------------------------------------------------------*/

#ifndef __RAD_BEM_GALERKIN_H
#define __RAD_BEM_GALERKIN_H

#include <vector>
#include <cstdint>

namespace radia {
namespace bem {

/**
 * Build the Laplace SL and DL Galerkin matrices on a P1-H1 surface mesh
 * (geometry up to P2-curved) in one pass.
 *
 * Outputs row-major dense (n_v, n_v) matrices.  Convention matches
 * NGSolve.bem (verified bit-exact to 1e-10):
 *   SL_ij = int_T_a int_T_b  L_i^a(r) [1/(4*pi*|r-r'|)]      L_j^b(r') dS dS'
 *   DL_ij = int_T_a int_T_b  L_i^a(r) [(r-r')*n_y/(4*pi*|r-r'|^3)]
 *                                                            L_j^b(r') dS dS'
 *
 * Singular cases (vertex/edge/identical) are evaluated with
 * Sauter-Schwab Duffy 4D quadrature.  DL edge case applies a parity
 * sign correction for opposite-orientation tris (closed surface
 * convention).
 *
 * Args:
 *   verts        (n_v, 3) row-major: vertex coords
 *   tris         (n_t, 3) row-major: triangle corner-vertex indices
 *   p2_nodes     (n_t, 6, 3) row-major: per-tri P2 Lagrange node
 *                 positions in order [v0, v1, v2, mid01, mid12, mid20].
 *                 Pass corner mid-points (e.g., 0.5*(v0+v1)) for FLAT
 *                 triangles; pass curved-mesh nodes (extracted via
 *                 mesh.GetTrafo at the 6 P2 reference points) for
 *                 curved meshes.
 *   n_v, n_t     dimensions
 *   regular_quad_degree   triangle Gauss order (5 / 7 / >=8 for
 *                          Duffy-collapsed n*n GL).  Default: 11.
 *   singular_n_q          1D Gauss-Legendre order for SS sub-cubes.
 *                          n_q^4 nodes per sub-cube; 2/5/6 sub-cubes
 *                          for vertex/edge/identical.  Default: 8.
 *   n_threads             0 = use OpenMP default (typically all cores)
 *
 * Outputs:
 *   SL_out, DL_out: pre-allocated n_v * n_v doubles, row-major.
 *                    Will be filled with the computed matrices.
 */
void AssembleSLDL(
    const double* verts, int n_v,
    const int64_t* tris, int n_t,
    const double* p2_nodes,
    int regular_quad_degree,
    int singular_n_q,
    int n_threads,
    double* SL_out,
    double* DL_out);

/**
 * Build the Laplace SL and DL Galerkin matrices on a P2 H1 Lagrange
 * surface basis (6 DOFs per triangle: 3 vertices + 3 edge midpoints),
 * with curved-P2 geometry.  Output dense (n_dof, n_dof) matrices.
 *
 * Convention: the 6 Lagrange basis functions on T_ref are
 *   N_0 = (1-xi-eta) (2*(1-xi-eta) - 1)   (vertex v0)
 *   N_1 = xi (2*xi - 1)                    (vertex v1)
 *   N_2 = eta (2*eta - 1)                  (vertex v2)
 *   N_3 = 4 * (1-xi-eta) * xi              (mid-edge v0-v1)
 *   N_4 = 4 * xi * eta                     (mid-edge v1-v2)
 *   N_5 = 4 * (1-xi-eta) * eta             (mid-edge v2-v0)
 *
 * Singular pairs use the same Sauter-Schwab Duffy 4D rules as the P1
 * version; the only changes are the per-tri 6-hat evaluation and the
 * 6-element index permutation that maps mid-edge basis indices through
 * a corner permutation.
 *
 * Args:
 *   verts        (n_v, 3): row-major corner vertex coords
 *   tris         (n_t, 3): row-major triangle corner-vertex indices
 *   p2_nodes     (n_t, 6, 3): per-tri curved-geometry P2 node positions
 *                  in [v0, v1, v2, mid01, mid12, mid20] order
 *   dofs_per_tri (n_t, 6): row-major P2 Lagrange DOF index per local
 *                  basis function, in the same [v0..mid20] order.  This
 *                  is what extract_surface_p2_lagrange returns on the
 *                  Python side.
 *   n_v, n_t     mesh dimensions (n_v = number of corner vertices)
 *   n_dof        size of the output square matrices (= n_v + n_unique_edges)
 *   regular_quad_degree   triangle Gauss order, same semantics as the P1 entry
 *   singular_n_q          1D Gauss-Legendre order for SS sub-cubes
 *   n_threads             0 = OpenMP default
 *
 * Outputs:
 *   SL_out, DL_out: pre-allocated n_dof * n_dof doubles, row-major.
 */
void AssembleSLDL_P2(
    const double* verts, int n_v,
    const int64_t* tris, int n_t,
    const double* p2_nodes,
    const int64_t* dofs_per_tri,
    int n_dof,
    int regular_quad_degree,
    int singular_n_q,
    int n_threads,
    double* SL_out,
    double* DL_out);

}}  // namespace radia::bem

#endif  // __RAD_BEM_GALERKIN_H
