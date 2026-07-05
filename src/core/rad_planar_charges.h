//-------------------------------------------------------------------------
// rad_planar_charges.h
//
// SHARED 2D planar postprocessing: exterior field + Maxwell-stress torque from
// a point-charge cloud (2D Laplace kernel G = -ln(r)/(2 pi)).  The HDiv-VIM
// planar layer and dense planar helpers evaluate their exterior field / torque
// through this one routine.  The natural cloud is the M.n equivalent bound charge
// of the solved per-element magnetization,
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

// Out-of-plane vector potential A_z at nP points from the same cloud (for the reduced-FEM eddy /
// maglev coupling): A_z = (mu0 / 2 pi) sum_a Q_a atan2(P_y-X_y, P_x-X_x).  (dA_z/dy = mu0 H_x,
// -dA_z/dx = mu0 H_y.)  BRANCH-CUT CAVEAT: the atan2 form has a cut along the -x ray of every
// charge -- valid when the eval set sees the charges from ONE side (as in radia.vim._vim2d.Az_at);
// for points SURROUNDING the body use the polar-integrated single-valued construction.
void FieldAz(int nq, const double* Xq, const double* Q,
             int nP, const double* P, double* Azout);

// Maxwell-stress torque per unit length about (cx,cy) on a circle of radius Rc
// in air.  The TOTAL field on the circle is the cloud (body) field PLUS the
// UNIFORM applied field (hextx,hexty): a uniform field's self-torque integrates
// to zero, but its CROSS term with the body field is the reluctance/alignment
// torque (= mu0 A (M x H0)).  Pass hext=0 for the body-only self-torque.
//   T = mu0 Rc^2 (2 pi / n) sum_i H_r(phi_i) H_phi(phi_i),  H = H_body + H_ext.
double MaxwellTorqueCircle(int nq, const double* Xq, const double* Q,
                           double Rc, double cx, double cy, int n,
                           double hextx, double hexty);

// Maxwell-stress FORCE per unit length on a circle of radius Rc in air.  Fout[2] = (Fx, Fy):
//   F_i = mu0 Rc (2 pi / n) sum [ H_r H_i - 1/2 |H|^2 n_i ],  H = H_body(cloud) + H_ext(uniform).
// A UNIFORM applied field gives ZERO net force (force needs a field gradient / a second body); the
// useful case is a cloud CONCATENATING more than one body with the circle enclosing ONE of them
// (maglev / actuator: the force on body A in the field of body B).
void MaxwellForceCircle(int nq, const double* Xq, const double* Q,
                        double Rc, double cx, double cy, int n,
                        double hextx, double hexty, double* Fout);

} // namespace rad_planar_charges

#endif
