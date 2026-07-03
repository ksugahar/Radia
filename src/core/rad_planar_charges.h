//-------------------------------------------------------------------------
// rad_planar_charges.h
//
// SHARED 2D planar postprocessing: exterior field + Maxwell-stress torque from
// a point-charge cloud (2D Laplace kernel G = -ln(r)/(2 pi)).  Method-agnostic:
// BOTH the collocation MMMM (radia.mmmm2d) and the HDiv-VIM (radia.vim._vim2d)
// evaluate their exterior field / torque through this one routine -- each just
// supplies its charge cloud (Xq, Q).  The natural, method-independent cloud is
// the M.n equivalent bound charge of the solved per-element magnetization,
// sampled on the element edges (built in Python by radia.planar_charges).
//
// H(P) = (1/2 pi) sum_a Q_a (P - X_a) / |P - X_a|^2            (a 2D "point" charge
//        = a z-infinite line charge; Q_a already carries density x length x weight)
// Maxwell torque per unit length on a circle of radius Rc in AIR (B = mu0 H):
//   T = mu0 Rc^2 \oint H_r H_phi dphi,  with H = H_body(cloud) + H_ext(uniform).
//-------------------------------------------------------------------------
#ifndef __RAD_PLANAR_CHARGES_H
#define __RAD_PLANAR_CHARGES_H

namespace rad_planar_charges {

// H at nP observation points from a cloud of nq planar point charges.
//   Xq[2*nq]   charge positions (x,y)
//   Q[nq]      charge magnitudes (already density x length x quad weight)
//   P[2*nP]    observation points (x,y)
//   Hout[2*nP] OUT: field (Hx,Hy).  A point that coincides with a charge (r^2==0)
//              skips that charge (self term undefined for a point observation).
void Field(int nq, const double* Xq, const double* Q,
           int nP, const double* P, double* Hout);

// Maxwell-stress torque per unit length about (cx,cy) on a circle of radius Rc
// in air.  The TOTAL field on the circle is the cloud (body) field PLUS the
// UNIFORM applied field (hextx,hexty): a uniform field's self-torque integrates
// to zero, but its CROSS term with the body field is the reluctance/alignment
// torque (= mu0 A (M x H0)).  Pass hext=0 for the body-only self-torque.
//   T = mu0 Rc^2 (2 pi / n) sum_i H_r(phi_i) H_phi(phi_i),  H = H_body + H_ext.
double MaxwellTorqueCircle(int nq, const double* Xq, const double* Q,
                           double Rc, double cx, double cy, int n,
                           double hextx, double hexty);

} // namespace rad_planar_charges

#endif
