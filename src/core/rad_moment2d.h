//-------------------------------------------------------------------------
// rad_moment2d.h
//
// 2D planar (per-unit-length) collocation MMMM -- the Multipole Magnetic
// Moment Method on a 2D cross-section mesh of triangles / quadrilaterals.
// The 2D twin of the 3D moment kernel in rad_interaction.cpp: each element
// carries ONE uniform line-charge DOF per EDGE, and the constitutive law
// M = chi H is imposed on the field MOMENTS about the element centroid:
//   1 monopole (net charge / div B = 0) + 2 dipole (M(c)=chi H(c))
//   + (nEdge-3) quadrupole residual-eigenmode rows (chi grad H(c)).
// Triangle (3 edges): 1 monopole + 2 dipole, NO quadrupole (the 2D simplex,
// tet-analog).  Quad (4 edges): + 1 quadrupole row.
//
// Kernel: 2D Laplace G = -ln(r)/(2 pi).  The field and its gradient (Hessian)
// of a uniform line charge on a segment have closed forms (below).  The
// evaluation point is always an element centroid (interior -> never on an
// edge -> singularity-free), so no analytic self-term is needed.
//
// Reference (formulation validated in an independent numpy PoC):
//   disk demag 1/2, ellipse 2:1 -> 1/3, 2/3, chi-sweep chi/(1+chi/2).
//-------------------------------------------------------------------------
#ifndef __RAD_MOMENT2D_H
#define __RAD_MOMENT2D_H

namespace rad_moment2d {

// Solve the LINEAR 2D planar MMMM demag problem.
//
//   nElem            number of elements
//   voff[nElem+1]    element -> vertex-start offsets into vxy (element k has
//                    vertices voff[k]..voff[k+1); nEdge_k = voff[k+1]-voff[k])
//   vxy[2*nVert]     concatenated vertex coordinates (x,y), any winding
//                    (each element is re-oriented CCW internally)
//   chi[nElem]       per-element susceptibility (mu_r - 1)
//   Hext[2*nElem]    per-element applied field at the element centroid (Hx,Hy)
//   Mout[2*nElem]    OUT: per-element magnetization (Mx,My)
//
// Returns 0 on success, non-zero on a solver / geometry error.
int SolveLinear(int nElem, const int* voff, const double* vxy,
                const double* chi, const double* Hext, double* Mout);

// Solve the LINEAR problem for nRHS applied fields sharing the SAME geometry + chi (the moment
// matrix is Hext-independent, so it is assembled + LU-factored ONCE and back-substituted for all
// nRHS -- e.g. a rotation / angle sweep where only the applied field rotates).
//   HextMulti[nRHS*nElem*2]  RHS r, element k, comp c at [(r*nElem+k)*2+c]
//   MoutMulti[nRHS*nElem*2]  OUT: same layout, per-RHS per-element magnetization
int SolveMulti(int nElem, const int* voff, const double* vxy, const double* chi,
               int nRHS, const double* HextMulti, double* MoutMulti);

} // namespace rad_moment2d

#endif
