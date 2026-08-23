/* rad_hdiv_vim.h -- analytic charge-potential kernels for the BDM1 HDiv-VIM.
 *
 * Mesh topology, HDiv spaces, and charge maps are supplied by NGSolve.  This
 * header contains only the geometry kernels used by the production Coulomb
 * Gram and field reconstruction paths.
 */
#ifndef RAD_HDIV_VIM_H
#define RAD_HDIV_VIM_H

#include <vector>
#include <array>

namespace rad_hdiv {

typedef std::array<double, 3> Vec3;

//-------------------------------------------------------------------------
// M2: analytic charge-Gram potentials (the accurate Wilton surface + phi_tet volume integrals,
// for the C++ scalable ChargeGram H-matrix entry; reused from the dense Python reference).
//-------------------------------------------------------------------------

/* Exact INT_T 1/|r-r'| dA' over a flat triangle (V = 3 verts, [3][3]) at obs r -- the Wilton/Graglia
 * analytic triangle potential (pure 1/r integral, NO 1/4pi).  Matches the dense Python
 * radia.vim._core.tri_potential to ~1e-12 and rad_poly_analytical's I_scalar Wilton edge formula. */
double TriPotential(const double V[3][3], const double r[3]);

/* Newtonian potential INT_tet 1/|P-r'| dV' of a uniform tetrahedron (V = 4 verts, [4][3]) at P, via the
 * divergence theorem reusing TriPotential.  Matches radia.vim._core.phi_tet. */
double PhiTet(const double V[4][3], const double P[3]);

/* Field INT_T (r-r')/|r-r'|^3 dA' of a uniform (sigma=1) flat triangle = -grad TriPotential (the Wilton
 * triangle FIELD; out[3], NO 1/4pi).  Matches radia.vim.flat_triangle_charge_field. */
void TriField(const double V[3][3], const double r[3], double out[3]);

/* Field INT_tet (P-r')/|P-r'|^3 dV' of a uniform tetrahedron = -grad PhiTet, via the divergence theorem
 * (n_face*TriPotential + d_face*TriField).  out[3], NO 1/4pi.  Matches radia.vim.tet_self_volume_field*4pi. */
void TetField(const double V[4][3], const double P[3], double out[3]);

/* ---- degree-1/2 polynomial-charge field kernels (the order<=2 fast path; NO 1/4pi) ----
 * Ports of radia.vim._field, validated entry-by-entry vs the Python (machine precision). */
void TriMoment1(const double V[3][3], const double r[3], double out[3]);          /* INT_T r'/R dS' */
void TriMoment2(const double V[3][3], const double r[3], double out[3][3]);       /* INT_T r'(x)r'/R dS' */
void TetMoment1(const double V[4][3], const double r[3], double out[3]);          /* INT_V r'/R dV' */
/* Forward directional derivative of the same closed-form PhiTet/TetMoment1
 * path.  value/direction use moment ordering [1,x,y,z]. */
void TetPotentialMomentsDirectionalUpTo1(
    const double V[4][3], const double dV[4][3], const double r[3], const double dr[3],
    double value[4], double direction[4]);
double TetPotentialPolynomial(const double V[4][3], const double r[3],
                              const std::vector<std::array<int,3>>& exps,
                              const std::vector<double>& coeffs);                 /* SUM c_a INT_V x^a/R dV */
/* Selected reference-coordinate moments INT_V xi^a/R dV through total degree
 * 18.  `out[k]` follows `exps[k]`.  This is the shared hot kernel used by the
 * HCurl local-polynomial HACApK host blocks, so all polynomial charges on one
 * source tetrahedron reuse one degree ladder. */
void TetReferencePotentialMoments(const double V[4][3], const double r[3],
                                  const std::vector<std::array<int,3>>& exps,
                                  double* out);
void TetReferencePotentialMomentsDirectional(
    const double V[4][3], const double dV[4][3],
    const double r[3], const double dr[3],
    const std::vector<std::array<int,3>>& exps,
    double* value, double* direction);
void TetPotentialMomentsUpTo3(const double V[4][3], const double r[3],
                              double out[20]);                                    /* total-degree <= 3 moments */
void TetPotentialMomentsUpTo6(const double V[4][3], const double r[3],
                              double out[84]);                                    /* total-degree <= 6 moments */
void TetPotentialMomentsDirectionalUpTo3(const double V[4][3],const double dV[4][3],
    const double r[3],const double dr[3],double value[20],double direction[20]);
void TetPotentialMomentsDirectionalUpTo6(const double V[4][3],const double dV[4][3],
    const double r[3],const double dr[3],double value[84],double direction[84]);
/* Reduced HCurl-VIM vector-potential Gram on affine tetrahedra.
 *
 * Each vector mode is supplied as three polynomials in the tetrahedron's
 * reference coordinates (lambda1,lambda2,lambda3).  The source integral uses
 * analytic reference-monomial Newton potentials through total degree 18;
 * only the smooth outer integral uses the caller-supplied tetrahedron rule.
 * coefficients layout: [mode][cell][monomial][component].  The returned
 * row-major matrix includes 1/(4*pi), but not permeability.
 */
std::vector<double> TetHCurlReducedGram(
    const std::vector<double>& cell_verts,
    const std::vector<std::array<int,3>>& exponents,
    const std::vector<double>& coefficients,
    int n_modes,
    const std::vector<double>& ref_points,
    const std::vector<double>& ref_weights);
void TriPotentialMomentsUpTo4(const double V[3][3], const double r[3],
                              double out[35]);                                    /* total-degree <= 4 moments */
void TriPotentialMomentsUpTo2(const double V[3][3], const double r[3],
                              double out[10]);                                    /* total-degree <= 2 moments */
void TriPotentialMomentsDirectionalUpTo4(
    const double V[3][3],const double dV[3][3],const double r[3],const double dr[3],
    double value[35],double direction[35]);
void TriPotentialMomentsDirectionalUpTo2(
    const double V[3][3],const double dV[3][3],const double r[3],const double dr[3],
    double value[10],double direction[10]);
void TetVolFieldLinear(const double V[4][3], const double r[3], double rho0,
                       const double g[3], double out[3]);                          /* linear volume charge */
void TetVolFieldQuadratic(const double V[4][3], const double r[3], double rho0,
                          const double g[3], const double Q[3][3], double out[3]); /* quadratic volume charge */
void TetVolFieldCubic(const double V[4][3], const double r[3],
                      const double coefficient[20], double out[3]);              /* total-degree <= 3 volume charge */
void TetVolFieldCubicBasis(const double V[4][3], const double r[3],
                           double out[20][3]);                                   /* fields of physical monomials */
void LinTriField(const double V[3][3], const double r[3], double sigma0,
                 const double s[3], double out[3]);                                /* linear surface charge */
void QuadTriField(const double V[3][3], const double r[3], double sigma0,
                  const double s[3], const double S[3][3], double out[3]);         /* quadratic surface charge */
void QuadTriFieldBasis(const double V[3][3], const double r[3],
                       double out[10][3]);                                        /* fields of physical monomials */
void CubicTriField(const double V[3][3], const double r[3],
                   const double coefficient[20], double out[3]);                  /* total-degree <= 3 surface charge */
/* Forward geometry/coefficient derivative of the same exact affine TET/TRI
 * field kernels.  Both value and direction exclude 1/(4*pi).  These routines
 * differentiate the closed-form moment formulas; no finite-difference or
 * observation quadrature is used. */
void TetVolFieldLinearDirectional(
    const double V[4][3], const double dV[4][3],
    const double r[3], const double dr[3],
    double rho0, double drho0, const double g[3], const double dg[3],
    double value[3], double direction[3]);
void QuadTriFieldDirectional(
    const double V[3][3], const double dV[3][3],
    const double r[3], const double dr[3],
    double sigma0, double dsigma0,
    const double s[3], const double ds[3],
    const double S[3][3], const double dS[3][3],
    double value[3], double direction[3]);

/* Closest point on a flat triangle / in a tetrahedron to p -- the Duffy singularity origin x0 for the
 * order>=3 / curved singular-quadrature charge-potential path (PhiAtHO_Duffy).  out[3]. */
void ClosestPointTriangle(const double p[3], const double a[3], const double b[3], const double c[3],
                          double out[3]);
void ClosestPointTet(const double V[4][3], const double p[3], double out[3]);

/* CURVED (isoparametric P2) panel support.  ClosestRefTri: xi0 = argmin|X(xi)-p|^2 in reference coords
 * (out[2]).  CurvedTriPotential: INT_curvedtri xi^e0 eta^e1 / |p-X(xi)| dA_curved via the reference Duffy from
 * xi0, evaluating the curved P2 map X(xi) + curved area element at each point.  nodes = 6x3 P2 nodes (corners
 * 0,1,2 ; mid-edges 3=(0-1),4=(1-2),5=(2-0)); gl/gw = an nq-point Gauss-Legendre rule on [0,1]. */
void ClosestRefTri(const double nodes[6][3], const double p[3], double xi0[2]);
double CurvedTriPotential(const double nodes[6][3], int e0, int e1, const double p[3],
                          const double* gl, const double* gw, int nq,
                          bool include_measure = true);

/* CURVED (isoparametric P2) tetrahedron VOLUME-charge support.  nodes = 10x3 P2 nodes (corners 0,1,2,3 ;
 * mid-edges 4=(0-1),5=(1-2),6=(2-0),7=(0-3),8=(1-3),9=(2-3)).  CurvedTetPotential: INT_curvedtet
 * xi^e0 eta^e1 zeta^e2 / |p-X(xi)| dV_curved via the reference Duffy from xi0=ClosestRefTet. */
void ClosestRefTet(const double nodes[10][3], const double p[3], double xi0[3]);
double CurvedTetPotential(const double nodes[10][3], int e0, int e1, int e2, const double p[3],
                          const double* gl, const double* gw, int nq,
                          bool include_measure = true);

/* OUTER-quadrature helpers for the curved charge Gram: curved physical point X(xi) + curved MEASURE (area
 * element dA = |Xu x Xv| for a P2 face, volume element dV = |det dX/dxi| for a P2 tet) at a reference point. */
void CurvedTriMapMeasure(const double nodes[6][3], double xi, double eta, double X[3], double& dA);
void CurvedTetMapMeasure(const double nodes[10][3], double xi, double eta, double zeta, double X[3], double& dV);

/* Affine constant-volume-charge self energy and analytic directional shape
 * derivatives for one TET/HEX/WEDGE cell.  The scalar is
 *   (4*pi)^-1 INT_D INT_D |x-y|^-1 dx dy.
 * `cell_type` is 0=TET (4 nodes), 1=HEX (8 nodes), 2=WEDGE (6 nodes), with
 * the usual corner ordering. `node_velocities` is row-major [mode][node][3].
 * The implementation uses the existing exact PhiTet singular inner kernel
 * and the Hadamard boundary derivative; no finite difference and no FE basis
 * reconstruction occurs here.  Return layout is [value,dvalue_0,...]. */
std::vector<double> AffineCellSelfEnergyShapeDerivative(
    int cell_type, const std::vector<double>& nodes,
    const std::vector<double>& node_velocities, int n_modes);

} // namespace rad_hdiv
#endif
