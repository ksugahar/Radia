/* rad_hdiv_vim.cpp -- analytic charge-potential kernels for BDM1 HDiv-VIM. */
#include "rad_hdiv_vim.h"
#include <cmath>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <vector>

#include "rad_parallel.h"

namespace rad_hdiv {

static const double PI = 3.14159265358979323846;
static const double INV_FOUR_PI = 1.0 / (4.0 * PI);

// ---- analytic charge-Gram potentials (M2: the Wilton surface + phi_tet volume integrals) ----
// Small raw-double3 vector helpers (local; the file's Vec3 is std::array<double,3>).
static inline double v3dot(const double a[3], const double b[3]) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static inline void   v3cross(const double a[3], const double b[3], double o[3])
{ o[0]=a[1]*b[2]-a[2]*b[1]; o[1]=a[2]*b[0]-a[0]*b[2]; o[2]=a[0]*b[1]-a[1]*b[0]; }
static inline double v3nrm(const double a[3]) { return std::sqrt(v3dot(a,a)); }

// Exact INT_T 1/|r-r'| dA' over a flat triangle (V0,V1,V2) at obs r -- the Wilton/Graglia analytic
// triangle potential (Wilton-Rao-Glisson, IEEE TAP 32(3):276, 1984).  Pure 1/r integral (NO 1/4pi).
// Port of examples reference radia.vim._core.tri_potential (scalar form); same edge formula as
// rad_poly_analytical.cpp::RadScalarPotentialFromTriangleFaceGlobal's I_scalar.
double TriPotential(const double V[3][3], const double r[3])
{
    double e1[3], e2[3], n[3];
    for (int k=0;k<3;k++){ e1[k]=V[1][k]-V[0][k]; e2[k]=V[2][k]-V[0][k]; }
    v3cross(e1,e2,n); double nl=v3nrm(n); if (nl<1e-300) return 0.0;
    for (int k=0;k<3;k++) n[k]/=nl;
    double rmv0[3]; for (int k=0;k<3;k++) rmv0[k]=r[k]-V[0][k];
    double d = v3dot(rmv0,n);                                   // signed height above the plane
    double p[3]; for (int k=0;k<3;k++) p[k]=r[k]-d*n[k];         // projection onto the plane
    double ad = std::fabs(d);
    double I = 0.0;
    for (int i=0;i<3;i++){
        const double* a = V[i];
        const double* b = V[(i+1)%3];
        double lh[3]; for (int k=0;k<3;k++) lh[k]=b[k]-a[k];
        double ll=v3nrm(lh); if (ll<1e-300) continue; for (int k=0;k<3;k++) lh[k]/=ll;
        double uh[3]; v3cross(lh,n,uh);                         // in-plane unit normal to the edge
        double ap[3]; for (int k=0;k<3;k++) ap[k]=a[k]-p[k];
        double bp[3]; for (int k=0;k<3;k++) bp[k]=b[k]-p[k];
        double P0 = v3dot(ap,uh);
        double sm = v3dot(ap,lh), sp = v3dot(bp,lh);
        double ra[3], rb[3]; for (int k=0;k<3;k++){ ra[k]=r[k]-a[k]; rb[k]=r[k]-b[k]; }
        double Rm = v3nrm(ra), Rp = v3nrm(rb);
        double R0sq = P0*P0 + d*d;
        double dm = Rm+sm, dp = Rp+sp;
        double f = (dp>1e-300 && dm>1e-300) ? std::log(dp/dm) : 0.0;
        double beta = std::atan2(P0*sp, R0sq+ad*Rp) - std::atan2(P0*sm, R0sq+ad*Rm);
        I += P0*f - ad*beta;
    }
    return I;
}

// Exact field INT_T (r-r')/|r-r'|^3 dA' of a UNIFORM (sigma=1) flat triangle = -grad_r TriPotential
// (the Wilton/Graglia triangle FIELD; vector form reusing the SAME per-edge quantities as TriPotential).
// NO 1/4pi.  Validated entry-by-entry vs radia.vim.flat_triangle_charge_field (machine precision).
void TriField(const double V[3][3], const double r[3], double out[3])
{
    out[0]=out[1]=out[2]=0.0;
    double e1[3], e2[3], n[3];
    for (int k=0;k<3;k++){ e1[k]=V[1][k]-V[0][k]; e2[k]=V[2][k]-V[0][k]; }
    v3cross(e1,e2,n); double nl=v3nrm(n); if (nl<1e-300) return;
    for (int k=0;k<3;k++) n[k]/=nl;
    double rmv0[3]; for (int k=0;k<3;k++) rmv0[k]=r[k]-V[0][k];
    double d = v3dot(rmv0,n);                                   // signed height
    double p[3]; for (int k=0;k<3;k++) p[k]=r[k]-d*n[k];        // projection
    double ad = std::fabs(d);
    double omega = 0.0;
    for (int i=0;i<3;i++){
        const double* a = V[i];
        const double* b = V[(i+1)%3];
        double lh[3]; for (int k=0;k<3;k++) lh[k]=b[k]-a[k];
        double ll=v3nrm(lh); if (ll<1e-300) continue; for (int k=0;k<3;k++) lh[k]/=ll;
        double uh[3]; v3cross(lh,n,uh);                         // in-plane unit normal to the edge
        double ap[3]; for (int k=0;k<3;k++) ap[k]=a[k]-p[k];
        double bp[3]; for (int k=0;k<3;k++) bp[k]=b[k]-p[k];
        double P0 = v3dot(ap,uh);
        double sm = v3dot(ap,lh), sp = v3dot(bp,lh);
        double ra[3], rb[3]; for (int k=0;k<3;k++){ ra[k]=r[k]-a[k]; rb[k]=r[k]-b[k]; }
        double Rm = v3nrm(ra), Rp = v3nrm(rb);
        double R0sq = P0*P0 + d*d;
        double dm = Rm+sm, dp = Rp+sp;
        double f = (dp>1e-300 && dm>1e-300) ? std::log(dp/dm) : 0.0;
        double beta = std::atan2(P0*sp, R0sq+ad*Rp) - std::atan2(P0*sm, R0sq+ad*Rm);
        for (int k=0;k<3;k++) out[k] += f*uh[k];               // tangential (per-edge log term)
        omega += beta;
    }
    double sgn = (d>0)?1.0:((d<0)?-1.0:0.0);
    for (int k=0;k<3;k++) out[k] += sgn*omega*n[k];            // normal (solid-angle term)
}

// Newtonian potential INT_tet 1/|P-r'| dV' of a uniform tetrahedron (4 verts) at P, via the divergence
// theorem (nabla'^2 R = 2/R): INT_V 1/R dV = (1/2) sum_{4 faces} d_face * INT_face 1/R dA', reusing
// TriPotential.  Port of radia.vim._core.phi_tet.
double PhiTet(const double V[4][3], const double P[3])
{
    double cen[3]={0,0,0};
    for (int i=0;i<4;i++) for (int k=0;k<3;k++) cen[k]+=V[i][k]*0.25;
    static const int FACES[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    double tot=0.0;
    for (int fi=0;fi<4;fi++){
        double Fv[3][3];
        for (int j=0;j<3;j++) for (int k=0;k<3;k++) Fv[j][k]=V[FACES[fi][j]][k];
        double e1[3], e2[3], nrm[3];
        for (int k=0;k<3;k++){ e1[k]=Fv[1][k]-Fv[0][k]; e2[k]=Fv[2][k]-Fv[0][k]; }
        v3cross(e1,e2,nrm); double nl=v3nrm(nrm); if (nl<1e-300) continue;
        for (int k=0;k<3;k++) nrm[k]/=nl;
        double fc[3]={0,0,0}; for (int j=0;j<3;j++) for (int k=0;k<3;k++) fc[k]+=Fv[j][k]/3.0;
        double ov[3]; for (int k=0;k<3;k++) ov[k]=fc[k]-cen[k];
        if (v3dot(ov,nrm)<0) for (int k=0;k<3;k++) nrm[k]=-nrm[k];   // outward
        double fv0mP[3]; for (int k=0;k<3;k++) fv0mP[k]=Fv[0][k]-P[k];
        double dd = v3dot(fv0mP,nrm);
        tot += dd * TriPotential(Fv, P);
    }
    return 0.5*tot;
}

// Field INT_tet (P-r')/|P-r'|^3 dV' of a uniform tetrahedron = -grad_P PhiTet, via the divergence
// theorem: = 0.5 sum_faces [ n_face * TriPotential_face + d_face * TriField_face ].  Exact near AND far
// (no quadrature).  Validated vs radia.vim.tet_self_volume_field (the spherical ray-trace) * 4pi.
void TetField(const double V[4][3], const double P[3], double out[3])
{
    out[0]=out[1]=out[2]=0.0;
    double cen[3]={0,0,0};
    for (int i=0;i<4;i++) for (int k=0;k<3;k++) cen[k]+=V[i][k]*0.25;
    static const int FACES[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    for (int fi=0;fi<4;fi++){
        double Fv[3][3];
        for (int j=0;j<3;j++) for (int k=0;k<3;k++) Fv[j][k]=V[FACES[fi][j]][k];
        double e1[3], e2[3], nrm[3];
        for (int k=0;k<3;k++){ e1[k]=Fv[1][k]-Fv[0][k]; e2[k]=Fv[2][k]-Fv[0][k]; }
        v3cross(e1,e2,nrm); double nl=v3nrm(nrm); if (nl<1e-300) continue;
        for (int k=0;k<3;k++) nrm[k]/=nl;
        double fc[3]={0,0,0}; for (int j=0;j<3;j++) for (int k=0;k<3;k++) fc[k]+=Fv[j][k]/3.0;
        double ov[3]; for (int k=0;k<3;k++) ov[k]=fc[k]-cen[k];
        if (v3dot(ov,nrm)<0) for (int k=0;k<3;k++) nrm[k]=-nrm[k];   // outward
        double fv0mP[3]; for (int k=0;k<3;k++) fv0mP[k]=Fv[0][k]-P[k];
        double dd = v3dot(fv0mP,nrm);
        double tp = TriPotential(Fv, P);
        double tf[3]; TriField(Fv, P, tf);
        for (int k=0;k<3;k++) out[k] += nrm[k]*tp + dd*tf[k];
    }
    for (int k=0;k<3;k++) out[k] *= 0.5;
}

// ===== degree-1/2 polynomial-charge field kernels (the order<=2 fast path) =====
// Port of radia.vim._field {triangle_potential_moment / _moment2, tet_newtonian_moment,
// tet_volume_field_linear / _quadratic, linear_ / quadratic_triangle_charge_field}.  NO 1/4pi.
// Validated entry-by-entry vs the Python references (machine precision).

// in-plane OUTWARD edge normal m of edge (A,B) on a triangle with unit normal n and centroid cen
static inline void edge_outnormal(const double A[3], const double B[3], const double n[3],
                                  const double cen[3], double m[3])
{
    double t[3]; for (int k=0;k<3;k++) t[k]=B[k]-A[k];
    double L=v3nrm(t); for (int k=0;k<3;k++) t[k]/=L;
    v3cross(t,n,m);
    double mid[3]; for (int k=0;k<3;k++) mid[k]=0.5*(A[k]+B[k])-cen[k];
    if (v3dot(m,mid)<0) for (int k=0;k<3;k++) m[k]=-m[k];
}

// INT_{edge A->B} R dl  (= INT sqrt(u^2+d2) du), and (optionally) INT_{edge} R*(r'-r_p) dl into Rxi.
static double edge_R_dl(const double A[3], const double B[3], const double r[3])
{
    double t[3]; for (int k=0;k<3;k++) t[k]=B[k]-A[k];
    double L=v3nrm(t); if (L<1e-300) return 0.0;
    double th[3]; for (int k=0;k<3;k++) th[k]=t[k]/L;
    double w[3]; for (int k=0;k<3;k++) w[k]=r[k]-A[k];
    double l0=v3dot(w,th); double d2=v3dot(w,w)-l0*l0; if (d2<0) d2=0;
    double u1=-l0, u2=L-l0;
    double s1 = (d2<1e-300)? 0.5*u1*std::fabs(u1) : 0.5*(u1*std::sqrt(u1*u1+d2)+d2*std::asinh(u1/std::sqrt(d2)));
    double s2 = (d2<1e-300)? 0.5*u2*std::fabs(u2) : 0.5*(u2*std::sqrt(u2*u2+d2)+d2*std::asinh(u2/std::sqrt(d2)));
    return s2 - s1;
}

// INT_{edge A->B} R*(r'-r_p) dl = (A-r_p) INT R dl + that INT R l dl  (closed form).
static void edge_R_xi_dl(const double A[3], const double B[3], const double r[3],
                         const double r_p[3], double out[3])
{
    out[0]=out[1]=out[2]=0.0;
    double t[3]; for (int k=0;k<3;k++) t[k]=B[k]-A[k];
    double L=v3nrm(t); if (L<1e-300) return;
    double th[3]; for (int k=0;k<3;k++) th[k]=t[k]/L;
    double w[3]; for (int k=0;k<3;k++) w[k]=r[k]-A[k];
    double l0=v3dot(w,th); double d2=v3dot(w,w)-l0*l0; if (d2<0) d2=0;
    double u1=-l0, u2=L-l0;
    double Fsq1=(d2<1e-300)?0.5*u1*std::fabs(u1):0.5*(u1*std::sqrt(u1*u1+d2)+d2*std::asinh(u1/std::sqrt(d2)));
    double Fsq2=(d2<1e-300)?0.5*u2*std::fabs(u2):0.5*(u2*std::sqrt(u2*u2+d2)+d2*std::asinh(u2/std::sqrt(d2)));
    double gR=Fsq2-Fsq1;                                       // INT R dl
    double IRl=((u2*u2+d2)*std::sqrt(u2*u2+d2)-(u1*u1+d2)*std::sqrt(u1*u1+d2))/3.0 + l0*gR;  // INT R l dl
    for (int k=0;k<3;k++) out[k]=(A[k]-r_p[k])*gR + th[k]*IRl;
}

// INT_{edge A->B} (r-r')/R dl  (= G_e for the linear surface field), closed form.
static void edge_field_dl(const double A[3], const double B[3], const double r[3], double out[3])
{
    out[0]=out[1]=out[2]=0.0;
    double t[3]; for (int k=0;k<3;k++) t[k]=B[k]-A[k];
    double L=v3nrm(t); if (L<1e-300) return;
    double th[3]; for (int k=0;k<3;k++) th[k]=t[k]/L;
    double w[3]; for (int k=0;k<3;k++) w[k]=r[k]-A[k];
    double l0=v3dot(w,th); double d2=v3dot(w,w)-l0*l0; if (d2<0) d2=0; double d=std::sqrt(d2);
    double u1=-l0, u2=L-l0;
    double as1,as2;
    if (d<1e-300){ as1=(std::fabs(u1)>0)?(u1>0?1:-1)*std::log(2*std::fabs(u1)):0.0;
                   as2=(std::fabs(u2)>0)?(u2>0?1:-1)*std::log(2*std::fabs(u2)):0.0; }
    else { as1=std::asinh(u1/d); as2=std::asinh(u2/d); }
    double int_1R=as2-as1;                                     // INT 1/R dl
    double int_lR=(std::sqrt(u2*u2+d2)-std::sqrt(u1*u1+d2)) + l0*(as2-as1);   // INT l/R dl
    for (int k=0;k<3;k++) out[k]=(r[k]-A[k])*int_1R - th[k]*int_lR;
}

// M1 = INT_T r'/R dS'  (first moment, 3-vector) via the surface divergence theorem.
void TriMoment1(const double V[3][3], const double r[3], double out[3])
{
    double e1[3],e2[3],n[3];
    for (int k=0;k<3;k++){ e1[k]=V[1][k]-V[0][k]; e2[k]=V[2][k]-V[0][k]; }
    v3cross(e1,e2,n); double nl=v3nrm(n); for (int k=0;k<3;k++) n[k]/=nl;
    double h=0; { double d[3]; for(int k=0;k<3;k++) d[k]=r[k]-V[0][k]; h=v3dot(d,n); }
    double r_p[3]; for (int k=0;k<3;k++) r_p[k]=r[k]-h*n[k];
    double cen[3]={0,0,0}; for (int j=0;j<3;j++) for(int k=0;k<3;k++) cen[k]+=V[j][k]/3.0;
    double I0=TriPotential(V,r);
    out[0]=out[1]=out[2]=0.0;
    for (int i=0;i<3;i++){
        const double* A=V[i]; const double* B=V[(i+1)%3];
        double m[3]; edge_outnormal(A,B,n,cen,m);
        double gR=edge_R_dl(A,B,r);
        for (int k=0;k<3;k++) out[k]+=m[k]*gR;
    }
    for (int k=0;k<3;k++) out[k]+=r_p[k]*I0;
}

// M2 = INT_T r'(x)r'/R dS' (symmetric 3x3, row-major out[3][3]) via the Hessian-of-R^3 identity.
void TriMoment2(const double V[3][3], const double r[3], double out[3][3])
{
    double e1[3],e2[3],n[3];
    for (int k=0;k<3;k++){ e1[k]=V[1][k]-V[0][k]; e2[k]=V[2][k]-V[0][k]; }
    v3cross(e1,e2,n); double nl=v3nrm(n); for (int k=0;k<3;k++) n[k]/=nl;
    double h=0; { double d[3]; for(int k=0;k<3;k++) d[k]=r[k]-V[0][k]; h=v3dot(d,n); }
    double r_p[3]; for (int k=0;k<3;k++) r_p[k]=r[k]-h*n[k];
    double cen[3]={0,0,0}; for (int j=0;j<3;j++) for(int k=0;k<3;k++) cen[k]+=V[j][k]/3.0;
    double I0=TriPotential(V,r);
    double Mxi1[3]={0,0,0}, ohm[3][3]={{0,0,0},{0,0,0},{0,0,0}}, sum_m_dot=0.0;
    for (int i=0;i<3;i++){
        const double* A=V[i]; const double* B=V[(i+1)%3];
        double m[3]; edge_outnormal(A,B,n,cen,m);
        double ARxi[3]; edge_R_xi_dl(A,B,r,r_p,ARxi);
        double gR=edge_R_dl(A,B,r);
        for (int k=0;k<3;k++) Mxi1[k]+=m[k]*gR;
        for (int a=0;a<3;a++) for (int b=0;b<3;b++) ohm[a][b]+=ARxi[a]*m[b];
        sum_m_dot += v3dot(m,ARxi);
    }
    double intR=(sum_m_dot + h*h*I0)/3.0;                       // INT_T R dS'
    for (int a=0;a<3;a++) for (int b=0;b<3;b++){
        double Pproj = (a==b?1.0:0.0) - n[a]*n[b];
        double Mxi2 = ohm[a][b] - Pproj*intR;
        out[a][b] = r_p[a]*r_p[b]*I0 + r_p[a]*Mxi1[b] + Mxi1[a]*r_p[b] + Mxi2;
    }
}

// V1 = INT_V r'/R dV' over a tet (3-vector) = (1/3)[r PhiTet - SUM_f h_f M1_f].
void TetMoment1(const double V[4][3], const double r[3], double out[3])
{
    double cen[3]={0,0,0}; for (int i=0;i<4;i++) for(int k=0;k<3;k++) cen[k]+=V[i][k]*0.25;
    static const int FACES[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    double Phi=PhiTet(V,r);
    for (int k=0;k<3;k++) out[k]=r[k]*Phi;
    for (int fi=0;fi<4;fi++){
        double Fv[3][3]; for (int j=0;j<3;j++) for(int k=0;k<3;k++) Fv[j][k]=V[FACES[fi][j]][k];
        double e1[3],e2[3],nrm[3];
        for (int k=0;k<3;k++){ e1[k]=Fv[1][k]-Fv[0][k]; e2[k]=Fv[2][k]-Fv[0][k]; }
        v3cross(e1,e2,nrm); double nl=v3nrm(nrm); if (nl<1e-300) continue; for (int k=0;k<3;k++) nrm[k]/=nl;
        double fc[3]={0,0,0}; for (int j=0;j<3;j++) for(int k=0;k<3;k++) fc[k]+=Fv[j][k]/3.0;
        double ov[3]; for (int k=0;k<3;k++) ov[k]=fc[k]-cen[k];
        if (v3dot(ov,nrm)<0) for (int k=0;k<3;k++) nrm[k]=-nrm[k];
        double hf; { double d[3]; for(int k=0;k<3;k++) d[k]=r[k]-Fv[0][k]; hf=v3dot(d,nrm); }
        double M1[3]; TriMoment1(Fv,r,M1);
        for (int k=0;k<3;k++) out[k]-=hf*M1[k];
    }
    for (int k=0;k<3;k++) out[k]/=3.0;
}

namespace {

constexpr int POLY_MAX_DEG = 18;
constexpr int POLY_MAX_MOMENTS =
    (POLY_MAX_DEG + 1)*(POLY_MAX_DEG + 2)*(POLY_MAX_DEG + 3)/6;

static inline double small_comb(int n, int k)
{
    if (k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    double r = 1.0;
    for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;
    return r;
}

struct TriPolyEdge {
    double xiA[2];
    double t2[2];
    double m2[2];
    double L;
    double l0;
    double d2;
};

struct TriPolySetup {
    double n[3];
    double h;
    double rp[3];
    double e1[3];
    double e2[3];
    TriPolyEdge edges[3];
};

static bool tri_poly_setup(const double P[3][3], const double r[3], TriPolySetup& g)
{
    double a1[3], a2[3];
    for (int k = 0; k < 3; ++k) { a1[k] = P[1][k] - P[0][k]; a2[k] = P[2][k] - P[0][k]; }
    v3cross(a1, a2, g.n);
    double nl = v3nrm(g.n);
    if (nl < 1e-300) return false;
    for (int k = 0; k < 3; ++k) g.n[k] /= nl;
    double rv[3]; for (int k = 0; k < 3; ++k) rv[k] = r[k] - P[0][k];
    g.h = v3dot(rv, g.n);
    for (int k = 0; k < 3; ++k) g.rp[k] = r[k] - g.h * g.n[k];
    double e1n = v3nrm(a1);
    if (e1n < 1e-300) return false;
    for (int k = 0; k < 3; ++k) g.e1[k] = a1[k] / e1n;
    v3cross(g.n, g.e1, g.e2);
    double cen[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) cen[k] += P[i][k] / 3.0;
    for (int i = 0; i < 3; ++i) {
        const double* A = P[i];
        const double* B = P[(i + 1) % 3];
        TriPolyEdge& e = g.edges[i];
        double t[3]; for (int k = 0; k < 3; ++k) t[k] = B[k] - A[k];
        e.L = v3nrm(t);
        if (e.L < 1e-300) return false;
        double th[3]; for (int k = 0; k < 3; ++k) th[k] = t[k] / e.L;
        double m[3]; v3cross(th, g.n, m);
        double mid[3]; for (int k = 0; k < 3; ++k) mid[k] = 0.5 * (A[k] + B[k]) - cen[k];
        if (v3dot(m, mid) < 0.0) for (int k = 0; k < 3; ++k) m[k] = -m[k];
        e.m2[0] = v3dot(m, g.e1);
        e.m2[1] = v3dot(m, g.e2);
        double Amrp[3]; for (int k = 0; k < 3; ++k) Amrp[k] = A[k] - g.rp[k];
        e.xiA[0] = v3dot(Amrp, g.e1);
        e.xiA[1] = v3dot(Amrp, g.e2);
        e.t2[0] = v3dot(th, g.e1);
        e.t2[1] = v3dot(th, g.e2);
        double w[3]; for (int k = 0; k < 3; ++k) w[k] = r[k] - A[k];
        e.l0 = v3dot(w, th);
        e.d2 = v3dot(w, w) - e.l0 * e.l0;
        if (e.d2 < 0.0) e.d2 = 0.0;
    }
    return true;
}

static double edge_R_for_poly(const TriPolyEdge& e)
{
    const double u1 = -e.l0, u2 = e.L - e.l0, d2 = e.d2;
    auto F = [d2](double u) {
        if (d2 < 1e-300) return 0.5 * u * std::fabs(u);
        const double d = std::sqrt(d2);
        return 0.5 * (u * std::sqrt(u * u + d2) + d2 * std::asinh(u / d));
    };
    return F(u2) - F(u1);
}

static void edge_l_moments_poly(const TriPolyEdge& e, int nmax, double Jl[POLY_MAX_DEG + 3])
{
    nmax = std::min(nmax, POLY_MAX_DEG + 2);
    const double d2 = e.d2;
    const double d = std::sqrt(d2);
    const double u1 = -e.l0, u2 = e.L - e.l0;
    double W[POLY_MAX_DEG + 3] = {};
    auto asinh_safe = [d](double u) {
        if (d > 1e-300) return std::asinh(u / d);
        if (std::fabs(u) == 0.0) return 0.0;
        return (u > 0.0 ? 1.0 : -1.0) * std::log(2.0 * std::fabs(u));
    };
    W[0] = asinh_safe(u2) - asinh_safe(u1);
    if (nmax >= 1) W[1] = std::sqrt(u2 * u2 + d2) - std::sqrt(u1 * u1 + d2);
    for (int n = 2; n <= nmax; ++n) {
        const double term = (std::pow(u2, n - 1) * std::sqrt(u2 * u2 + d2)
                           - std::pow(u1, n - 1) * std::sqrt(u1 * u1 + d2)) / n;
        W[n] = term - ((n - 1.0) * d2 / n) * W[n - 2];
    }
    for (int n = 0; n <= nmax; ++n) {
        double s = 0.0;
        for (int i = 0; i <= n; ++i) s += small_comb(n, i) * std::pow(e.l0, n - i) * W[i];
        Jl[n] = s;
    }
}

static double edge_inplane_monomial_poly(const TriPolyEdge& e, int a, int b,
                                         const double Jl[POLY_MAX_DEG + 3])
{
    double poly[POLY_MAX_DEG + 3] = {};
    double tmp[POLY_MAX_DEG + 3] = {};
    int deg = 0;
    poly[0] = 1.0;
    auto mul_linear = [&](double c0, double c1) {
        std::fill(std::begin(tmp), std::end(tmp), 0.0);
        for (int i = 0; i <= deg; ++i) {
            tmp[i] += poly[i] * c0;
            tmp[i + 1] += poly[i] * c1;
        }
        ++deg;
        for (int i = 0; i <= deg; ++i) poly[i] = tmp[i];
    };
    for (int i = 0; i < a; ++i) mul_linear(e.xiA[0], e.t2[0]);
    for (int i = 0; i < b; ++i) mul_linear(e.xiA[1], e.t2[1]);
    double s = 0.0;
    for (int i = 0; i <= deg; ++i) s += poly[i] * Jl[i];
    return s;
}

static void triangle_inplane_A_moments(const double P[3][3], const double r[3], int degree,
                                       double A[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1],
                                       TriPolySetup& g)
{
    for (int i = 0; i <= POLY_MAX_DEG; ++i)
        for (int j = 0; j <= POLY_MAX_DEG; ++j) A[i][j] = 0.0;
    degree = std::min(degree, POLY_MAX_DEG);
    if (!tri_poly_setup(P, r, g)) return;
    double Jl[3][POLY_MAX_DEG + 3] = {};
    for (int i = 0; i < 3; ++i) edge_l_moments_poly(g.edges[i], degree + 2, Jl[i]);
    A[0][0] = TriPotential(P, r);
    if (degree >= 1) {
        double A1[2] = {0.0, 0.0};
        for (int i = 0; i < 3; ++i) {
            const double er = edge_R_for_poly(g.edges[i]);
            A1[0] += g.edges[i].m2[0] * er;
            A1[1] += g.edges[i].m2[1] * er;
        }
        A[1][0] = A1[0];
        A[0][1] = A1[1];
    }
    auto Eedge = [&](int j, int p, int q) {
        double s = 0.0;
        for (int i = 0; i < 3; ++i)
            s += g.edges[i].m2[j] * edge_inplane_monomial_poly(g.edges[i], p, q, Jl[i]);
        return s;
    };
    auto Eneg1 = [&](int a, int b) {
        return Eedge(0, a + 1, b) + Eedge(1, a, b + 1);
    };
    for (int k = 2; k <= degree; ++k) {
        for (int a = k; a >= 0; --a) {
            const int b = k - a;
            double h2B;
            if (a >= 1) h2B = g.h * g.h * ((a - 1) * (a >= 2 ? A[a - 2][b] : 0.0) - Eedge(0, a - 1, b));
            else        h2B = g.h * g.h * ((b - 1) * (b >= 2 ? A[a][b - 2] : 0.0) - Eedge(1, a, b - 1));
            A[a][b] = (Eneg1(a, b) - h2B) / (k + 1.0);
        }
    }
}

static void poly2_mul_linear(double poly[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1], int& deg,
                             double c0, double c1, double c2)
{
    double tmp[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1] = {};
    for (int a = 0; a <= deg; ++a) {
        for (int b = 0; b <= deg - a; ++b) {
            const double v = poly[a][b];
            tmp[a][b] += v * c0;
            tmp[a + 1][b] += v * c1;
            tmp[a][b + 1] += v * c2;
        }
    }
    ++deg;
    for (int a = 0; a <= deg; ++a)
        for (int b = 0; b <= deg - a; ++b) poly[a][b] = tmp[a][b];
}

static double SurfacePotentialMonomial(const double P[3][3], const double r[3],
                                       const int alpha[3], int degree)
{
    double A[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1];
    TriPolySetup g;
    triangle_inplane_A_moments(P, r, degree, A, g);
    double poly[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1] = {};
    int deg = 0;
    poly[0][0] = 1.0;
    for (int coord = 0; coord < 3; ++coord)
        for (int p = 0; p < alpha[coord]; ++p)
            poly2_mul_linear(poly, deg, g.rp[coord], g.e1[coord], g.e2[coord]);
    double s = 0.0;
    for (int a = 0; a <= deg; ++a)
        for (int b = 0; b <= deg - a; ++b) s += poly[a][b] * A[a][b];
    return s;
}

static int PotentialMomentIndex(int ax, int ay, int az)
{
    const int degree = ax + ay + az;
    int idx = 0;
    for (int d = 0; d < degree; ++d) idx += (d + 1)*(d + 2)/2;
    for (int x = 0; x < ax; ++x) idx += degree - x + 1;
    return idx + ay;
}

static void SurfacePotentialMomentsUpTo(const double P[3][3], const double r[3],
                                        int degree, double* out)
{
    double A[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1];
    TriPolySetup g;
    triangle_inplane_A_moments(P, r, degree, A, g);
    int idx = 0;
    for (int total = 0; total <= degree; ++total)
        for (int ax = 0; ax <= total; ++ax)
            for (int ay = 0; ay <= total - ax; ++ay) {
                const int alpha[3] = {ax, ay, total - ax - ay};
                double poly[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1] = {};
                int poly_degree = 0;
                poly[0][0] = 1.0;
                for (int coord = 0; coord < 3; ++coord)
                    for (int p = 0; p < alpha[coord]; ++p)
                        poly2_mul_linear(poly, poly_degree, g.rp[coord], g.e1[coord], g.e2[coord]);
                double s = 0.0;
                for (int a = 0; a <= poly_degree; ++a)
                    for (int b = 0; b <= poly_degree - a; ++b) s += poly[a][b]*A[a][b];
                out[idx++] = s;
            }
}

static void TetPotentialMomentsUpTo(const double V[4][3], const double r[3],
                                    int degree, double* out)
{
    static const int FACES[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    double face_moments[4][POLY_MAX_MOMENTS] = {};
    double h[4] = {};
    double cen[3] = {0,0,0};
    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 3; ++k) cen[k] += 0.25*V[i][k];
    for (int fi = 0; fi < 4; ++fi) {
        double Fv[3][3];
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) Fv[j][k] = V[FACES[fi][j]][k];
        double e1[3], e2[3], nrm[3];
        for (int k = 0; k < 3; ++k) { e1[k] = Fv[1][k]-Fv[0][k]; e2[k] = Fv[2][k]-Fv[0][k]; }
        v3cross(e1, e2, nrm);
        const double nl = v3nrm(nrm);
        if (nl < 1e-300) continue;
        for (int k = 0; k < 3; ++k) nrm[k] /= nl;
        double fc[3] = {0,0,0};
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) fc[k] += Fv[j][k]/3.0;
        double outward[3];
        for (int k = 0; k < 3; ++k) outward[k] = fc[k]-cen[k];
        if (v3dot(outward, nrm) < 0.0)
            for (int k = 0; k < 3; ++k) nrm[k] = -nrm[k];
        double rmf[3];
        for (int k = 0; k < 3; ++k) rmf[k] = r[k]-Fv[0][k];
        h[fi] = v3dot(rmf, nrm);
        SurfacePotentialMomentsUpTo(Fv, r, degree, face_moments[fi]);
    }
    out[0] = PhiTet(V, r);
    for (int total = 1; total <= degree; ++total)
        for (int ax = 0; ax <= total; ++ax)
            for (int ay = 0; ay <= total - ax; ++ay) {
                const int az = total - ax - ay;
                const int idx = PotentialMomentIndex(ax, ay, az);
                double s = 0.0;
                for (int fi = 0; fi < 4; ++fi) s -= h[fi]*face_moments[fi][idx];
                const int alpha[3] = {ax, ay, az};
                for (int k = 0; k < 3; ++k)
                    if (alpha[k] > 0) {
                        int lower[3] = {ax, ay, az};
                        --lower[k];
                        s += r[k]*alpha[k]*out[PotentialMomentIndex(lower[0], lower[1], lower[2])];
                    }
                out[idx] = s/(total + 2.0);
            }
}

static bool TetReferenceInverse(const double V[4][3], double invJ[3][3], double& detJ)
{
    const double J[3][3] = {
        {V[1][0]-V[0][0], V[2][0]-V[0][0], V[3][0]-V[0][0]},
        {V[1][1]-V[0][1], V[2][1]-V[0][1], V[3][1]-V[0][1]},
        {V[1][2]-V[0][2], V[2][2]-V[0][2], V[3][2]-V[0][2]}
    };
    detJ = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
         - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
         + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
    if (std::fabs(detJ) < 1e-300) return false;
    const double id = 1.0/detJ;
    invJ[0][0] =  (J[1][1]*J[2][2]-J[1][2]*J[2][1])*id;
    invJ[0][1] = -(J[0][1]*J[2][2]-J[0][2]*J[2][1])*id;
    invJ[0][2] =  (J[0][1]*J[1][2]-J[0][2]*J[1][1])*id;
    invJ[1][0] = -(J[1][0]*J[2][2]-J[1][2]*J[2][0])*id;
    invJ[1][1] =  (J[0][0]*J[2][2]-J[0][2]*J[2][0])*id;
    invJ[1][2] = -(J[0][0]*J[1][2]-J[0][2]*J[1][0])*id;
    invJ[2][0] =  (J[1][0]*J[2][1]-J[1][1]*J[2][0])*id;
    invJ[2][1] = -(J[0][0]*J[2][1]-J[0][1]*J[2][0])*id;
    invJ[2][2] =  (J[0][0]*J[1][1]-J[0][1]*J[1][0])*id;
    return true;
}

static void SurfaceReferencePotentialMomentsUpTo(
    const double P[3][3], const double r[3], int degree,
    const double v0[3], const double invJ[3][3], double* out)
{
    double A[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1];
    TriPolySetup g;
    triangle_inplane_A_moments(P, r, degree, A, g);
    double rpmv0[3];
    for (int k = 0; k < 3; ++k) rpmv0[k] = g.rp[k] - v0[k];
    int idx = 0;
    for (int total = 0; total <= degree; ++total)
        for (int ax = 0; ax <= total; ++ax)
            for (int ay = 0; ay <= total - ax; ++ay) {
                const int alpha[3] = {ax, ay, total - ax - ay};
                double poly[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1] = {};
                int poly_degree = 0;
                poly[0][0] = 1.0;
                for (int coord = 0; coord < 3; ++coord) {
                    const double c0 = invJ[coord][0]*rpmv0[0]
                                    + invJ[coord][1]*rpmv0[1]
                                    + invJ[coord][2]*rpmv0[2];
                    const double c1 = invJ[coord][0]*g.e1[0]
                                    + invJ[coord][1]*g.e1[1]
                                    + invJ[coord][2]*g.e1[2];
                    const double c2 = invJ[coord][0]*g.e2[0]
                                    + invJ[coord][1]*g.e2[1]
                                    + invJ[coord][2]*g.e2[2];
                    for (int p = 0; p < alpha[coord]; ++p)
                        poly2_mul_linear(poly, poly_degree, c0, c1, c2);
                }
                double value = 0.0;
                for (int a = 0; a <= poly_degree; ++a)
                    for (int b = 0; b <= poly_degree - a; ++b)
                        value += poly[a][b]*A[a][b];
                out[idx++] = value;
            }
}

/* Stable analogue of TetPotentialMomentsUpTo in the source tetrahedron's
 * reference coordinates.  Keeping lambda in its natural O(1) frame avoids
 * the severe cancellation caused by expanding fifth-order lambda monomials
 * into global x/y/z powers on millimetre-scale cells. */
static void TetReferencePotentialMomentsUpTo(
    const double V[4][3], const double r[3], int degree, double* out)
{
    static const int FACES[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    double invJ[3][3], detJ = 0.0;
    if (!TetReferenceInverse(V, invJ, detJ)) {
        const int count = (degree + 1)*(degree + 2)*(degree + 3)/6;
        std::fill(out, out + count, 0.0);
        return;
    }
    double xi_r[3] = {};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            xi_r[i] += invJ[i][k]*(r[k]-V[0][k]);

    double face_moments[4][POLY_MAX_MOMENTS] = {};
    double h[4] = {};
    double cen[3] = {0,0,0};
    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 3; ++k) cen[k] += 0.25*V[i][k];
    for (int fi = 0; fi < 4; ++fi) {
        double Fv[3][3];
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) Fv[j][k] = V[FACES[fi][j]][k];
        double e1[3], e2[3], nrm[3];
        for (int k = 0; k < 3; ++k) {
            e1[k] = Fv[1][k]-Fv[0][k];
            e2[k] = Fv[2][k]-Fv[0][k];
        }
        v3cross(e1, e2, nrm);
        const double nl = v3nrm(nrm);
        if (nl < 1e-300) continue;
        for (int k = 0; k < 3; ++k) nrm[k] /= nl;
        double fc[3] = {0,0,0};
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) fc[k] += Fv[j][k]/3.0;
        double outward[3];
        for (int k = 0; k < 3; ++k) outward[k] = fc[k]-cen[k];
        if (v3dot(outward, nrm) < 0.0)
            for (int k = 0; k < 3; ++k) nrm[k] = -nrm[k];
        double rmf[3];
        for (int k = 0; k < 3; ++k) rmf[k] = r[k]-Fv[0][k];
        h[fi] = v3dot(rmf, nrm);
        SurfaceReferencePotentialMomentsUpTo(
            Fv, r, degree, V[0], invJ, face_moments[fi]);
    }

    out[0] = PhiTet(V, r);
    for (int total = 1; total <= degree; ++total)
        for (int ax = 0; ax <= total; ++ax)
            for (int ay = 0; ay <= total - ax; ++ay) {
                const int az = total - ax - ay;
                const int idx = PotentialMomentIndex(ax, ay, az);
                double value = 0.0;
                for (int fi = 0; fi < 4; ++fi)
                    value -= h[fi]*face_moments[fi][idx];
                const int alpha[3] = {ax, ay, az};
                for (int k = 0; k < 3; ++k) {
                    if (alpha[k] <= 0) continue;
                    int lower[3] = {ax, ay, az};
                    --lower[k];
                    value += xi_r[k]*alpha[k]
                           * out[PotentialMomentIndex(lower[0], lower[1], lower[2])];
                }
                out[idx] = value/(total + 2.0);
            }
}

struct TetMomentMemo {
    bool seen[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1][POLY_MAX_DEG + 1] = {};
    double val[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1][POLY_MAX_DEG + 1] = {};
};

static double TetPotentialMomentRec(const double V[4][3], const double r[3], const int alpha[3],
                                    int degree, TetMomentMemo& memo)
{
    const int d = alpha[0] + alpha[1] + alpha[2];
    if (d == 0) return PhiTet(V, r);
    bool& seen = memo.seen[alpha[0]][alpha[1]][alpha[2]];
    double& cached = memo.val[alpha[0]][alpha[1]][alpha[2]];
    if (seen) return cached;
    double cen[3] = {0, 0, 0};
    for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) cen[k] += V[i][k] * 0.25;
    static const int FACES[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    double s = 0.0;
    for (int fi = 0; fi < 4; ++fi) {
        double Fv[3][3];
        for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) Fv[j][k] = V[FACES[fi][j]][k];
        double e1[3], e2[3], nrm[3];
        for (int k = 0; k < 3; ++k) { e1[k] = Fv[1][k] - Fv[0][k]; e2[k] = Fv[2][k] - Fv[0][k]; }
        v3cross(e1, e2, nrm);
        double nl = v3nrm(nrm);
        if (nl < 1e-300) continue;
        for (int k = 0; k < 3; ++k) nrm[k] /= nl;
        double fc[3] = {0,0,0};
        for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) fc[k] += Fv[j][k] / 3.0;
        double ov[3]; for (int k = 0; k < 3; ++k) ov[k] = fc[k] - cen[k];
        if (v3dot(ov, nrm) < 0.0) for (int k = 0; k < 3; ++k) nrm[k] = -nrm[k];
        double rmf[3]; for (int k = 0; k < 3; ++k) rmf[k] = r[k] - Fv[0][k];
        const double h = v3dot(rmf, nrm);
        s -= h * SurfacePotentialMonomial(Fv, r, alpha, degree);
    }
    for (int i = 0; i < 3; ++i) {
        if (alpha[i] <= 0) continue;
        int am[3] = {alpha[0], alpha[1], alpha[2]};
        am[i] -= 1;
        s += r[i] * alpha[i] * TetPotentialMomentRec(V, r, am, degree, memo);
    }
    cached = s / (d + 2.0);
    seen = true;
    return cached;
}

} // namespace

double TetPotentialPolynomial(const double V[4][3], const double r[3],
                              const std::vector<std::array<int,3>>& exps,
                              const std::vector<double>& coeffs)
{
    const size_t n = std::min(exps.size(), coeffs.size());
    int degree = 0;
    for (size_t i = 0; i < n; ++i)
        degree = std::max(degree, exps[i][0] + exps[i][1] + exps[i][2]);
    if (degree > POLY_MAX_DEG)
        throw std::runtime_error("TetPotentialPolynomial: degree exceeds POLY_MAX_DEG");
    TetMomentMemo memo;
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) {
        int a[3] = {exps[i][0], exps[i][1], exps[i][2]};
        if (a[0] < 0 || a[1] < 0 || a[2] < 0 || a[0] + a[1] + a[2] > POLY_MAX_DEG)
            throw std::runtime_error("TetPotentialPolynomial: invalid exponent");
        s += coeffs[i] * TetPotentialMomentRec(V, r, a, degree, memo);
    }
    return s;
}

void TetPotentialMomentsUpTo3(const double V[4][3], const double r[3], double out[20])
{
    TetPotentialMomentsUpTo(V, r, 3, out);
}

void TetPotentialMomentsUpTo6(const double V[4][3], const double r[3], double out[84])
{
    TetPotentialMomentsUpTo(V, r, 6, out);
}

std::vector<double> TetHCurlReducedGram(
    const std::vector<double>& cell_verts,
    const std::vector<std::array<int,3>>& exponents,
    const std::vector<double>& coefficients,
    int n_modes,
    const std::vector<double>& ref_points,
    const std::vector<double>& ref_weights)
{
    if (cell_verts.empty() || cell_verts.size() % 12 != 0)
        throw std::invalid_argument("TetHCurlReducedGram: cell_verts must have shape (n_cell,4,3)");
    if (n_modes <= 0)
        throw std::invalid_argument("TetHCurlReducedGram: n_modes must be positive");
    if (exponents.empty())
        throw std::invalid_argument("TetHCurlReducedGram: exponents must not be empty");
    if (ref_points.empty() || ref_points.size() % 3 != 0
        || ref_weights.size() != ref_points.size()/3)
        throw std::invalid_argument("TetHCurlReducedGram: invalid outer tetrahedron rule");
    const int n_cells = static_cast<int>(cell_verts.size()/12);
    const int n_mono = static_cast<int>(exponents.size());
    const size_t expected = static_cast<size_t>(n_modes)*n_cells*n_mono*3;
    if (coefficients.size() != expected)
        throw std::invalid_argument(
            "TetHCurlReducedGram: coefficients must have shape (n_mode,n_cell,n_mono,3)");
    int degree = 0;
    std::vector<int> moment_index(static_cast<size_t>(n_mono));
    for (int m = 0; m < n_mono; ++m) {
        const auto& e = exponents[static_cast<size_t>(m)];
        if (e[0] < 0 || e[1] < 0 || e[2] < 0)
            throw std::invalid_argument("TetHCurlReducedGram: exponents must be non-negative");
        const int d = e[0] + e[1] + e[2];
        if (d > POLY_MAX_DEG)
            throw std::invalid_argument("TetHCurlReducedGram: polynomial degree exceeds 18");
        degree = std::max(degree, d);
        moment_index[static_cast<size_t>(m)] = PotentialMomentIndex(e[0], e[1], e[2]);
    }
    for (double value : cell_verts)
        if (!std::isfinite(value))
            throw std::invalid_argument("TetHCurlReducedGram: cell_verts contains non-finite values");
    for (double value : coefficients)
        if (!std::isfinite(value))
            throw std::invalid_argument("TetHCurlReducedGram: coefficients contains non-finite values");
    for (double value : ref_points)
        if (!std::isfinite(value))
            throw std::invalid_argument("TetHCurlReducedGram: ref_points contains non-finite values");
    for (double value : ref_weights)
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument("TetHCurlReducedGram: ref_weights must be finite and positive");

    auto coeff_at = [&](int mode, int cell, int mono, int component) -> double {
        const size_t index = (((static_cast<size_t>(mode)*n_cells + cell)*n_mono + mono)*3
                              + component);
        return coefficients[index];
    };
    auto vertex_at = [&](int cell, int vertex, int component) -> double {
        return cell_verts[(static_cast<size_t>(cell)*4 + vertex)*3 + component];
    };

    std::vector<double> cell_blocks(
        static_cast<size_t>(n_cells)*n_modes*n_modes, 0.0);
    ngcore::ParallelFor(ngcore::IntRange(n_cells), [&](size_t target_index) {
        const int target = static_cast<int>(target_index);
        double Vt[4][3];
        for (int a = 0; a < 4; ++a)
            for (int k = 0; k < 3; ++k) Vt[a][k] = vertex_at(target, a, k);
        double invJ[3][3], detJ = 0.0;
        if (!TetReferenceInverse(Vt, invJ, detJ))
            throw std::runtime_error("TetHCurlReducedGram: degenerate target tetrahedron");
        const double abs_detJ = std::fabs(detJ);
        std::vector<double> target_values(static_cast<size_t>(n_modes)*3);
        std::vector<double> source_potential(static_cast<size_t>(n_modes)*3);
        std::vector<double> monomial_values(static_cast<size_t>(n_mono));
        double moments[POLY_MAX_MOMENTS] = {};

        for (size_t q = 0; q < ref_weights.size(); ++q) {
            const double xi[3] = {
                ref_points[3*q], ref_points[3*q + 1], ref_points[3*q + 2]
            };
            double point[3];
            for (int k = 0; k < 3; ++k) {
                point[k] = Vt[0][k]
                         + (Vt[1][k]-Vt[0][k])*xi[0]
                         + (Vt[2][k]-Vt[0][k])*xi[1]
                         + (Vt[3][k]-Vt[0][k])*xi[2];
            }
            for (int m = 0; m < n_mono; ++m) {
                const auto& e = exponents[static_cast<size_t>(m)];
                monomial_values[static_cast<size_t>(m)] =
                    std::pow(xi[0], e[0])*std::pow(xi[1], e[1])*std::pow(xi[2], e[2]);
            }
            std::fill(target_values.begin(), target_values.end(), 0.0);
            for (int mode = 0; mode < n_modes; ++mode)
                for (int m = 0; m < n_mono; ++m)
                    for (int k = 0; k < 3; ++k)
                        target_values[static_cast<size_t>(mode)*3 + k] +=
                            coeff_at(mode, target, m, k)*monomial_values[static_cast<size_t>(m)];

            std::fill(source_potential.begin(), source_potential.end(), 0.0);
            for (int source = 0; source < n_cells; ++source) {
                double Vs[4][3];
                for (int a = 0; a < 4; ++a)
                    for (int k = 0; k < 3; ++k) Vs[a][k] = vertex_at(source, a, k);
                std::fill(std::begin(moments), std::end(moments), 0.0);
                TetReferencePotentialMomentsUpTo(Vs, point, degree, moments);
                for (int mode = 0; mode < n_modes; ++mode)
                    for (int m = 0; m < n_mono; ++m) {
                        const double potential = moments[moment_index[static_cast<size_t>(m)]];
                        for (int k = 0; k < 3; ++k)
                            source_potential[static_cast<size_t>(mode)*3 + k] +=
                                coeff_at(mode, source, m, k)*potential;
                    }
            }

            const double weight = ref_weights[q]*abs_detJ;
            double* block = cell_blocks.data()
                          + static_cast<size_t>(target)*n_modes*n_modes;
            for (int i = 0; i < n_modes; ++i)
                for (int j = 0; j < n_modes; ++j) {
                    double dot = 0.0;
                    for (int k = 0; k < 3; ++k)
                        dot += target_values[static_cast<size_t>(i)*3 + k]
                             * source_potential[static_cast<size_t>(j)*3 + k];
                    block[static_cast<size_t>(i)*n_modes + j] += weight*dot;
                }
        }
    });

    std::vector<double> result(static_cast<size_t>(n_modes)*n_modes, 0.0);
    for (int target = 0; target < n_cells; ++target) {
        const double* block = cell_blocks.data()
                            + static_cast<size_t>(target)*n_modes*n_modes;
        for (int i = 0; i < n_modes; ++i)
            for (int j = 0; j < n_modes; ++j)
                result[static_cast<size_t>(i)*n_modes + j] +=
                    block[static_cast<size_t>(i)*n_modes + j];
    }
    for (int i = 0; i < n_modes; ++i)
        for (int j = i; j < n_modes; ++j) {
            const double value = 0.5*INV_FOUR_PI*(
                result[static_cast<size_t>(i)*n_modes + j]
              + result[static_cast<size_t>(j)*n_modes + i]);
            result[static_cast<size_t>(i)*n_modes + j] = value;
            result[static_cast<size_t>(j)*n_modes + i] = value;
        }
    return result;
}

void TriPotentialMomentsUpTo4(const double V[3][3], const double r[3], double out[35])
{
    SurfacePotentialMomentsUpTo(V, r, 4, out);
}

void TriPotentialMomentsUpTo2(const double V[3][3], const double r[3], double out[10])
{
    SurfacePotentialMomentsUpTo(V, r, 2, out);
}

// INT_V (rho0 + g.r')(r-r')/R^3 dV' (linear volume charge) = SUM_f n_f[rho0 I0_f + g.M1_f] - g PhiTet.
void TetVolFieldLinear(const double V[4][3], const double r[3], double rho0, const double g[3], double out[3])
{
    out[0]=out[1]=out[2]=0.0;
    double cen[3]={0,0,0}; for (int i=0;i<4;i++) for(int k=0;k<3;k++) cen[k]+=V[i][k]*0.25;
    static const int FACES[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    for (int fi=0;fi<4;fi++){
        double Fv[3][3]; for (int j=0;j<3;j++) for(int k=0;k<3;k++) Fv[j][k]=V[FACES[fi][j]][k];
        double e1[3],e2[3],nrm[3];
        for (int k=0;k<3;k++){ e1[k]=Fv[1][k]-Fv[0][k]; e2[k]=Fv[2][k]-Fv[0][k]; }
        v3cross(e1,e2,nrm); double nl=v3nrm(nrm); if (nl<1e-300) continue; for (int k=0;k<3;k++) nrm[k]/=nl;
        double fc[3]={0,0,0}; for (int j=0;j<3;j++) for(int k=0;k<3;k++) fc[k]+=Fv[j][k]/3.0;
        double ov[3]; for (int k=0;k<3;k++) ov[k]=fc[k]-cen[k];
        if (v3dot(ov,nrm)<0) for (int k=0;k<3;k++) nrm[k]=-nrm[k];
        double I0=TriPotential(Fv,r); double M1[3]; TriMoment1(Fv,r,M1);
        double w=rho0*I0 + v3dot(g,M1);
        for (int k=0;k<3;k++) out[k]+=nrm[k]*w;
    }
    double Phi=PhiTet(V,r);
    for (int k=0;k<3;k++) out[k]-=g[k]*Phi;
}

// INT_V (rho0 + g.r' + r'^T Q r')(r-r')/R^3 dV' (Q symmetric, row-major Q[3][3])
// = SUM_f n_f[rho0 I0_f + g.M1_f + Q:M2_f] - (g PhiTet + 2 Q.V1).
void TetVolFieldQuadratic(const double V[4][3], const double r[3], double rho0,
                          const double g[3], const double Q[3][3], double out[3])
{
    out[0]=out[1]=out[2]=0.0;
    double cen[3]={0,0,0}; for (int i=0;i<4;i++) for(int k=0;k<3;k++) cen[k]+=V[i][k]*0.25;
    static const int FACES[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    for (int fi=0;fi<4;fi++){
        double Fv[3][3]; for (int j=0;j<3;j++) for(int k=0;k<3;k++) Fv[j][k]=V[FACES[fi][j]][k];
        double e1[3],e2[3],nrm[3];
        for (int k=0;k<3;k++){ e1[k]=Fv[1][k]-Fv[0][k]; e2[k]=Fv[2][k]-Fv[0][k]; }
        v3cross(e1,e2,nrm); double nl=v3nrm(nrm); if (nl<1e-300) continue; for (int k=0;k<3;k++) nrm[k]/=nl;
        double fc[3]={0,0,0}; for (int j=0;j<3;j++) for(int k=0;k<3;k++) fc[k]+=Fv[j][k]/3.0;
        double ov[3]; for (int k=0;k<3;k++) ov[k]=fc[k]-cen[k];
        if (v3dot(ov,nrm)<0) for (int k=0;k<3;k++) nrm[k]=-nrm[k];
        double I0=TriPotential(Fv,r); double M1[3]; TriMoment1(Fv,r,M1);
        double M2[3][3]; TriMoment2(Fv,r,M2);
        double QM2=0.0; for (int a=0;a<3;a++) for (int b=0;b<3;b++) QM2+=Q[a][b]*M2[a][b];
        double w=rho0*I0 + v3dot(g,M1) + QM2;
        for (int k=0;k<3;k++) out[k]+=nrm[k]*w;
    }
    double Phi=PhiTet(V,r); double V1[3]; TetMoment1(V,r,V1);
    for (int k=0;k<3;k++){
        double QV1=Q[k][0]*V1[0]+Q[k][1]*V1[1]+Q[k][2]*V1[2];
        out[k]-=(g[k]*Phi + 2.0*QV1);
    }
}

// INT_T (sigma0 + s.r')(r-r')/R^3 dS' (linear surface charge)
// = (sigma0 + s.r_p) F_const - SUM_e (s.m_e) G_e - I0 s_par.
void LinTriField(const double V[3][3], const double r[3], double sigma0, const double s[3], double out[3])
{
    double e1[3],e2[3],n[3];
    for (int k=0;k<3;k++){ e1[k]=V[1][k]-V[0][k]; e2[k]=V[2][k]-V[0][k]; }
    v3cross(e1,e2,n); double nl=v3nrm(n); for (int k=0;k<3;k++) n[k]/=nl;
    double h=0; { double d[3]; for(int k=0;k<3;k++) d[k]=r[k]-V[0][k]; h=v3dot(d,n); }
    double r_p[3]; for (int k=0;k<3;k++) r_p[k]=r[k]-h*n[k];
    double cen[3]={0,0,0}; for (int j=0;j<3;j++) for(int k=0;k<3;k++) cen[k]+=V[j][k]/3.0;
    double Fc[3]; TriField(V,r,Fc);
    double I0=TriPotential(V,r);
    double sn=v3dot(s,n); double s_par[3]; for (int k=0;k<3;k++) s_par[k]=s[k]-sn*n[k];
    double coef=sigma0+v3dot(s,r_p);
    for (int k=0;k<3;k++) out[k]=coef*Fc[k] - I0*s_par[k];
    for (int i=0;i<3;i++){
        const double* A=V[i]; const double* B=V[(i+1)%3];
        double m[3]; edge_outnormal(A,B,n,cen,m);
        double G[3]; edge_field_dl(A,B,r,G);
        double sm=v3dot(s,m);
        for (int k=0;k<3;k++) out[k]-=sm*G[k];
    }
}

// edge 1D monomial moments Jl[0..2] = INT_edge l^n/R dl, plus in-plane edge geometry (for QuadTriField).
static void edge_Jmoments(const double A[3], const double B[3], const double r[3],
                          double Jl[3], double& L, double xiA[2], double t2[2],
                          const double e1[3], const double e2[3], const double r_p[3])
{
    double t[3]; for (int k=0;k<3;k++) t[k]=B[k]-A[k];
    L=v3nrm(t); double th[3]; for (int k=0;k<3;k++) th[k]=t[k]/L;
    double w[3]; for (int k=0;k<3;k++) w[k]=r[k]-A[k];
    double l0=v3dot(w,th); double d2=v3dot(w,w)-l0*l0; if (d2<0) d2=0; double d=std::sqrt(d2);
    double u1=-l0, u2=L-l0;
    double as1,as2;
    if (d<1e-300){ as1=(std::fabs(u1)>0)?(u1>0?1:-1)*std::log(2*std::fabs(u1)):0.0;
                   as2=(std::fabs(u2)>0)?(u2>0?1:-1)*std::log(2*std::fabs(u2)):0.0; }
    else { as1=std::asinh(u1/d); as2=std::asinh(u2/d); }
    double W0=as2-as1;
    double W1=std::sqrt(u2*u2+d2)-std::sqrt(u1*u1+d2);
    double W2=(0.5*(u2*std::sqrt(u2*u2+d2)-d2*as2)) - (0.5*(u1*std::sqrt(u1*u1+d2)-d2*as1));
    Jl[0]=W0; Jl[1]=W1+l0*W0; Jl[2]=W2+2*l0*W1+l0*l0*W0;        // l = u + l0
    double Av[3]; for (int k=0;k<3;k++) Av[k]=A[k]-r_p[k];
    xiA[0]=v3dot(Av,e1); xiA[1]=v3dot(Av,e2);
    t2[0]=v3dot(th,e1);  t2[1]=v3dot(th,e2);
}

// INT_T (sigma0 + s.r' + r'^T S r')(r-r')/R^3 dS' (S symmetric) via the in-plane/normal split.
void QuadTriField(const double V[3][3], const double r[3], double sigma0,
                  const double s[3], const double S[3][3], double out[3])
{
    double e1u[3],e2u[3],n[3];
    for (int k=0;k<3;k++){ e1u[k]=V[1][k]-V[0][k]; e2u[k]=V[2][k]-V[0][k]; }
    v3cross(e1u,e2u,n); double nl=v3nrm(n); for (int k=0;k<3;k++) n[k]/=nl;
    double e1[3]; double e1l=v3nrm(e1u); for (int k=0;k<3;k++) e1[k]=e1u[k]/e1l;
    double e2[3]; v3cross(n,e1,e2);
    double h=0; { double d[3]; for(int k=0;k<3;k++) d[k]=r[k]-V[0][k]; h=v3dot(d,n); }
    double r_p[3]; for (int k=0;k<3;k++) r_p[k]=r[k]-h*n[k];
    double cen[3]={0,0,0}; for (int j=0;j<3;j++) for(int k=0;k<3;k++) cen[k]+=V[j][k]/3.0;
    double Fc[3]; TriField(V,r,Fc);
    double I0=TriPotential(V,r);
    double M1[3]; TriMoment1(V,r,M1);
    double J3_0=v3dot(n,Fc)/h;                                 // INT_T 1/R^3
    double intxi1[3]={0,0,0};                                  // INT_T xi/R^3
    double xixi[3][3];                                         // INT_T xi(x)xi/R^3 = P I0 - SUM (Gxi)(x)m
    for (int a=0;a<3;a++) for (int b=0;b<3;b++) xixi[a][b]=((a==b?1.0:0.0)-n[a]*n[b])*I0;
    double inplane[3]={0,0,0};
    for (int i=0;i<3;i++){
        const double* A=V[i]; const double* B=V[(i+1)%3];
        double m[3]; edge_outnormal(A,B,n,cen,m);
        double Jl[3], L, xiA[2], t2[2];
        edge_Jmoments(A,B,r,Jl,L,xiA,t2,e1,e2,r_p);
        double m1=v3dot(m,e1), m2=v3dot(m,e2);
        // INT_edge xi/R dl = (A-r_p) J0 + th J1 ; in (e1,e2): Gxi = xiA*J0 + t2*J1
        double Gxi[2]={ xiA[0]*Jl[0]+t2[0]*Jl[1], xiA[1]*Jl[0]+t2[1]*Jl[1] };
        double Gxi3[3]; for (int k=0;k<3;k++) Gxi3[k]=Gxi[0]*e1[k]+Gxi[1]*e2[k];
        for (int k=0;k<3;k++) intxi1[k]-=m[k]*Jl[0];
        for (int a=0;a<3;a++) for (int b=0;b<3;b++) xixi[a][b]-=Gxi3[a]*m[b];
        // in-plane: m_e * INT_edge sigma/R dl ; sigma = sigma0 + s.r' + r'^T S r' as a poly c0+c1 l+c2 l^2
        // along r'(l) = A + l*th  (th = unit edge tangent)
        double sA=v3dot(s,A);
        double tt[3]; for (int k=0;k<3;k++) tt[k]=B[k]-A[k]; double Lt=v3nrm(tt); for(int k=0;k<3;k++) tt[k]/=Lt;
        double AStA=0, AStt=0, ttStt=0;
        for (int a=0;a<3;a++) for (int b=0;b<3;b++){ AStA+=A[a]*S[a][b]*A[b]; AStt+=A[a]*S[a][b]*tt[b]; ttStt+=tt[a]*S[a][b]*tt[b]; }
        double c0=sigma0+sA+AStA;
        double c1=v3dot(s,tt)+2.0*AStt;
        double c2=ttStt;
        double esig=c0*Jl[0]+c1*Jl[1]+c2*Jl[2];
        for (int k=0;k<3;k++) inplane[k]+=m[k]*esig;
    }
    double J3_1[3]; for (int k=0;k<3;k++) J3_1[k]=intxi1[k]+r_p[k]*J3_0;          // INT_T r'/R^3
    double J3_2[3][3];                                                            // INT_T r'(x)r'/R^3
    for (int a=0;a<3;a++) for (int b=0;b<3;b++)
        J3_2[a][b]=xixi[a][b] + r_p[a]*intxi1[b] + intxi1[a]*r_p[b] + r_p[a]*r_p[b]*J3_0;
    // in-plane: - (P s I0 + 2 P S M1)
    double sn=v3dot(s,n); double Psn[3]; for (int k=0;k<3;k++) Psn[k]=s[k]-sn*n[k];
    double SM1[3]; for (int a=0;a<3;a++){ SM1[a]=S[a][0]*M1[0]+S[a][1]*M1[1]+S[a][2]*M1[2]; }
    double SM1n=v3dot(SM1,n); double PSM1[3]; for (int k=0;k<3;k++) PSM1[k]=SM1[k]-SM1n*n[k];
    for (int k=0;k<3;k++) inplane[k]-=(Psn[k]*I0 + 2.0*PSM1[k]);
    // normal: h n [sigma0 J3_0 + s.J3_1 + S:J3_2]
    double SJ2=0; for (int a=0;a<3;a++) for (int b=0;b<3;b++) SJ2+=S[a][b]*J3_2[a][b];
    double nrmscal=h*(sigma0*J3_0 + v3dot(s,J3_1) + SJ2);
    for (int k=0;k<3;k++) out[k]=inplane[k]+n[k]*nrmscal;
}

// ---- closest-point helpers for the Duffy singularity origin (x0) ----------------------------------------
// Ericson, Real-Time Collision Detection: closest point on a triangle to p (Voronoi-region method).
void ClosestPointTriangle(const double p[3], const double a[3], const double b[3], const double c[3],
                          double out[3])
{
    auto dot = [](const double u[3], const double v[3]) { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; };
    double ab[3], ac[3], ap[3];
    for (int k = 0; k < 3; ++k) { ab[k] = b[k]-a[k]; ac[k] = c[k]-a[k]; ap[k] = p[k]-a[k]; }
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) { for (int k=0;k<3;++k) out[k]=a[k]; return; }
    double bp[3]; for (int k=0;k<3;++k) bp[k]=p[k]-b[k];
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) { for (int k=0;k<3;++k) out[k]=b[k]; return; }
    const double vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) { const double v=d1/(d1-d3); for(int k=0;k<3;++k) out[k]=a[k]+v*ab[k]; return; }
    double cp[3]; for (int k=0;k<3;++k) cp[k]=p[k]-c[k];
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) { for (int k=0;k<3;++k) out[k]=c[k]; return; }
    const double vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) { const double w=d2/(d2-d6); for(int k=0;k<3;++k) out[k]=a[k]+w*ac[k]; return; }
    const double va = d3*d6 - d5*d4;
    if (va <= 0 && (d4-d3) >= 0 && (d5-d6) >= 0) {
        const double w=(d4-d3)/((d4-d3)+(d5-d6)); for(int k=0;k<3;++k) out[k]=b[k]+w*(c[k]-b[k]); return;
    }
    const double denom = 1.0/(va+vb+vc), v = vb*denom, w = vc*denom;
    for (int k=0;k<3;++k) out[k] = a[k] + ab[k]*v + ac[k]*w;
}

// Closest point of a tetrahedron to p: p itself if inside (barycentric >= 0), else the min over the 4 faces.
void ClosestPointTet(const double V[4][3], const double p[3], double out[3])
{
    // barycentric via the affine inverse of [V1-V0, V2-V0, V3-V0]
    double E[3][3];
    for (int k = 0; k < 3; ++k) { E[k][0]=V[1][k]-V[0][k]; E[k][1]=V[2][k]-V[0][k]; E[k][2]=V[3][k]-V[0][k]; }
    const double det = E[0][0]*(E[1][1]*E[2][2]-E[1][2]*E[2][1])
                     - E[0][1]*(E[1][0]*E[2][2]-E[1][2]*E[2][0])
                     + E[0][2]*(E[1][0]*E[2][1]-E[1][1]*E[2][0]);
    if (std::fabs(det) > 1e-300) {
        const double d[3] = {p[0]-V[0][0], p[1]-V[0][1], p[2]-V[0][2]};
        const double inv = 1.0/det;
        const double l0 = inv*( (E[1][1]*E[2][2]-E[1][2]*E[2][1])*d[0] + (E[0][2]*E[2][1]-E[0][1]*E[2][2])*d[1] + (E[0][1]*E[1][2]-E[0][2]*E[1][1])*d[2] );
        const double l1 = inv*( (E[1][2]*E[2][0]-E[1][0]*E[2][2])*d[0] + (E[0][0]*E[2][2]-E[0][2]*E[2][0])*d[1] + (E[0][2]*E[1][0]-E[0][0]*E[1][2])*d[2] );
        const double l2 = inv*( (E[1][0]*E[2][1]-E[1][1]*E[2][0])*d[0] + (E[0][1]*E[2][0]-E[0][0]*E[2][1])*d[1] + (E[0][0]*E[1][1]-E[0][1]*E[1][0])*d[2] );
        if (l0 >= -1e-12 && l1 >= -1e-12 && l2 >= -1e-12 && (l0+l1+l2) <= 1.0+1e-12) {
            out[0]=p[0]; out[1]=p[1]; out[2]=p[2]; return;
        }
    }
    const int F[4][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
    double best[3]; double bd = 1e300;
    for (int f = 0; f < 4; ++f) {
        double q[3]; ClosestPointTriangle(p, V[F[f][0]], V[F[f][1]], V[F[f][2]], q);
        const double dx=p[0]-q[0], dy=p[1]-q[1], dz=p[2]-q[2]; const double dd=dx*dx+dy*dy+dz*dz;
        if (dd < bd) { bd=dd; best[0]=q[0]; best[1]=q[1]; best[2]=q[2]; }
    }
    out[0]=best[0]; out[1]=best[1]; out[2]=best[2];
}

// ---- CURVED (isoparametric) panel support: P2 (6-node) triangle geometry + curved-panel Duffy ----------
static inline double _ipow(double b, int e) { double r = 1.0; for (int i = 0; i < e; ++i) r *= b; return r; }
// P2 triangle shape functions + derivatives on the reference triangle {L1=1-xi-eta, L2=xi, L3=eta}; node
// order: 0,1,2 corners, 3=mid(0-1), 4=mid(1-2), 5=mid(2-0).  (Matches the validated curved_duffy.py.)
static void P2TriShape(double xi, double eta, double N[6], double dNxi[6], double dNeta[6])
{
    const double L1 = 1.0 - xi - eta, L2 = xi, L3 = eta;
    N[0]=L1*(2*L1-1); N[1]=L2*(2*L2-1); N[2]=L3*(2*L3-1); N[3]=4*L1*L2; N[4]=4*L2*L3; N[5]=4*L3*L1;
    dNxi[0]=1-4*L1; dNxi[1]=4*L2-1; dNxi[2]=0.0;      dNxi[3]=4*(L1-L2); dNxi[4]=4*L3;  dNxi[5]=-4*L3;
    dNeta[0]=1-4*L1; dNeta[1]=0.0;  dNeta[2]=4*L3-1;  dNeta[3]=-4*L2;    dNeta[4]=4*L2; dNeta[5]=4*(L1-L3);
}

// X(xi,eta) and the two tangents dX/dxi, dX/deta for a P2 curved triangle (nodes = 6 x 3).
static void CurvedTriEval(const double nodes[6][3], double xi, double eta,
                          double X[3], double Xu[3], double Xv[3])
{
    double N[6], dNxi[6], dNeta[6]; P2TriShape(xi, eta, N, dNxi, dNeta);
    for (int k = 0; k < 3; ++k) { X[k]=0; Xu[k]=0; Xv[k]=0; }
    for (int i = 0; i < 6; ++i) for (int k = 0; k < 3; ++k) {
        X[k]  += N[i]    * nodes[i][k];
        Xu[k] += dNxi[i] * nodes[i][k];
        Xv[k] += dNeta[i]* nodes[i][k];
    }
}

// xi0 = argmin |X(xi)-p|^2 over the reference triangle (Gauss-Newton from a coarse scan, clamped to the tri).
void ClosestRefTri(const double nodes[6][3], const double p[3], double xi0[2])
{
    auto clamp = [](double& a, double& b) {
        if (a < 0) a = 0; if (b < 0) b = 0;
        if (a + b > 1.0) { const double s = a + b; a /= s; b /= s; }
    };
    // coarse scan over a grid for a robust starting point
    double bx = 1.0/3, by = 1.0/3, bd = 1e300;
    for (int i = 0; i <= 6; ++i) for (int j = 0; j <= 6 - i; ++j) {
        double a = i/6.0, b = j/6.0, X[3], Xu[3], Xv[3]; CurvedTriEval(nodes, a, b, X, Xu, Xv);
        const double dd = (X[0]-p[0])*(X[0]-p[0])+(X[1]-p[1])*(X[1]-p[1])+(X[2]-p[2])*(X[2]-p[2]);
        if (dd < bd) { bd = dd; bx = a; by = b; }
    }
    double a = bx, b = by;
    for (int it = 0; it < 30; ++it) {
        double X[3], Xu[3], Xv[3]; CurvedTriEval(nodes, a, b, X, Xu, Xv);
        const double g[3] = {X[0]-p[0], X[1]-p[1], X[2]-p[2]};
        // Gauss-Newton: H = J^T J (2x2), rhs = J^T g (2), J = [Xu, Xv]
        const double h00 = Xu[0]*Xu[0]+Xu[1]*Xu[1]+Xu[2]*Xu[2];
        const double h01 = Xu[0]*Xv[0]+Xu[1]*Xv[1]+Xu[2]*Xv[2];
        const double h11 = Xv[0]*Xv[0]+Xv[1]*Xv[1]+Xv[2]*Xv[2];
        const double r0 = Xu[0]*g[0]+Xu[1]*g[1]+Xu[2]*g[2];
        const double r1 = Xv[0]*g[0]+Xv[1]*g[1]+Xv[2]*g[2];
        const double det = h00*h11 - h01*h01;
        if (std::fabs(det) < 1e-300) break;
        const double da = -( h11*r0 - h01*r1) / det;
        const double db = -(-h01*r0 + h00*r1) / det;
        double na = a + da, nb = b + db; clamp(na, nb);
        if (std::fabs(na-a) + std::fabs(nb-b) < 1e-13) { a = na; b = nb; break; }
        a = na; b = nb;
    }
    xi0[0] = a; xi0[1] = b;
}

// CURVED-panel surface-charge inner potential  INT_curvedtri  xi^e0 eta^e1 / |p - X(xi,eta)|  dA_curved
// via the reference Duffy from xi0 = closest reference point (3 SIGNED 2D reference sub-triangles), evaluating
// the curved map X(xi) and curved area element J=|Xu x Xv| at each reference Duffy point.  gl/gw = an nq-point
// Gauss-Legendre rule on [0,1].  Validated vs the Python prototype (curved_duffy.py) to ~1e-5..1e-7.
double CurvedTriPotential(const double nodes[6][3], int e0, int e1, const double p[3],
                          const double* gl, const double* gw, int nq, bool include_measure)
{
    double xi0[2]; ClosestRefTri(nodes, p, xi0);
    static const double C[3][2] = {{0,0},{1,0},{0,1}};        // reference-triangle corners
    double acc = 0.0;
    for (int k = 0; k < 3; ++k) {
        const double* A = C[k]; const double* B = C[(k+1)%3];
        const double e1x = A[0]-xi0[0], e1y = A[1]-xi0[1];
        const double e2x = B[0]-xi0[0], e2y = B[1]-xi0[1];
        const double sgn2 = e1x*e2y - e1y*e2x;                // signed 2*area of the 2D reference sub-triangle
        for (int a = 0; a < nq; ++a) {
            const double u = gl[a];
            for (int b = 0; b < nq; ++b) {
                const double v = gl[b];
                const double xi  = xi0[0] + u*e1x + u*v*(e2x - e1x);
                const double eta = xi0[1] + u*e1y + u*v*(e2y - e1y);
                double X[3], Xu[3], Xv[3]; CurvedTriEval(nodes, xi, eta, X, Xu, Xv);
                const double cr[3] = {Xu[1]*Xv[2]-Xu[2]*Xv[1], Xu[2]*Xv[0]-Xu[0]*Xv[2], Xu[0]*Xv[1]-Xu[1]*Xv[0]};
                const double J = std::sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);   // curved area element
                const double dx=p[0]-X[0], dy=p[1]-X[1], dz=p[2]-X[2];
                const double r = std::sqrt(dx*dx+dy*dy+dz*dz);
                if (r < 1e-300) continue;
                const double measure = include_measure ? J : 1.0;
                acc += gw[a]*gw[b]*(u*sgn2)*measure*_ipow(xi, e0)*_ipow(eta, e1)/r;
            }
        }
    }
    return acc;
}

// ---- CURVED P2 (10-node) tetrahedron geometry + curved-panel Duffy (VOLUME charge) ---------------------
// node order: 0,1,2,3 corners ; 4=mid(0-1),5=mid(1-2),6=mid(2-0),7=mid(0-3),8=mid(1-3),9=mid(2-3).
static void P2TetShape(double xi, double eta, double zeta, double N[10], double dN[10][3])
{
    const double L1=1-xi-eta-zeta, L2=xi, L3=eta, L4=zeta;
    N[0]=L1*(2*L1-1); N[1]=L2*(2*L2-1); N[2]=L3*(2*L3-1); N[3]=L4*(2*L4-1);
    N[4]=4*L1*L2; N[5]=4*L2*L3; N[6]=4*L3*L1; N[7]=4*L1*L4; N[8]=4*L2*L4; N[9]=4*L3*L4;
    const double dL1[3]={-1,-1,-1}, dL2[3]={1,0,0}, dL3[3]={0,1,0}, dL4[3]={0,0,1};
    for (int k = 0; k < 3; ++k) {
        dN[0][k]=(4*L1-1)*dL1[k]; dN[1][k]=(4*L2-1)*dL2[k]; dN[2][k]=(4*L3-1)*dL3[k]; dN[3][k]=(4*L4-1)*dL4[k];
        dN[4][k]=4*(dL1[k]*L2+L1*dL2[k]); dN[5][k]=4*(dL2[k]*L3+L2*dL3[k]); dN[6][k]=4*(dL3[k]*L1+L3*dL1[k]);
        dN[7][k]=4*(dL1[k]*L4+L1*dL4[k]); dN[8][k]=4*(dL2[k]*L4+L2*dL4[k]); dN[9][k]=4*(dL3[k]*L4+L3*dL4[k]);
    }
}

static void CurvedTetEval(const double nodes[10][3], double xi, double eta, double zeta,
                          double X[3], double Jac[3][3])     // Jac[k][c] = dX_k/dxi_c
{
    double N[10], dN[10][3]; P2TetShape(xi, eta, zeta, N, dN);
    for (int k = 0; k < 3; ++k) { X[k]=0; for (int c=0;c<3;++c) Jac[k][c]=0; }
    for (int i = 0; i < 10; ++i) for (int k = 0; k < 3; ++k) {
        X[k] += N[i]*nodes[i][k];
        for (int c = 0; c < 3; ++c) Jac[k][c] += dN[i][c]*nodes[i][k];
    }
}

static double det3(const double M[3][3])
{
    return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
         - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
         + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
}

// xi0 = argmin |X(xi)-p|^2 over the reference tet {xi,eta,zeta>=0, xi+eta+zeta<=1} (Gauss-Newton + clamp).
void ClosestRefTet(const double nodes[10][3], const double p[3], double xi0[3])
{
    auto clamp = [](double z[3]) {
        for (int k=0;k<3;++k) if (z[k] < 0) z[k] = 0;
        const double s = z[0]+z[1]+z[2]; if (s > 1.0) { z[0]/=s; z[1]/=s; z[2]/=s; }
    };
    double best[3] = {0.25,0.25,0.25}, bd = 1e300;
    for (int i=0;i<=4;++i) for (int j=0;j<=4-i;++j) for (int k=0;k<=4-i-j;++k) {
        double z[3]={i/4.0,j/4.0,k/4.0}, X[3], Jac[3][3]; CurvedTetEval(nodes,z[0],z[1],z[2],X,Jac);
        const double dd=(X[0]-p[0])*(X[0]-p[0])+(X[1]-p[1])*(X[1]-p[1])+(X[2]-p[2])*(X[2]-p[2]);
        if (dd<bd){bd=dd;best[0]=z[0];best[1]=z[1];best[2]=z[2];}
    }
    double z[3]={best[0],best[1],best[2]};
    for (int it=0; it<30; ++it) {
        double X[3], Jac[3][3]; CurvedTetEval(nodes,z[0],z[1],z[2],X,Jac);
        const double g[3]={X[0]-p[0],X[1]-p[1],X[2]-p[2]};
        double H[3][3], rhs[3];
        for (int a=0;a<3;++a){ rhs[a]=0; for (int k=0;k<3;++k) rhs[a]+=Jac[k][a]*g[k];
            for (int b=0;b<3;++b){ H[a][b]=0; for (int k=0;k<3;++k) H[a][b]+=Jac[k][a]*Jac[k][b]; } }
        const double dH = det3(H); if (std::fabs(dH) < 1e-300) break;
        double d[3];
        for (int c=0;c<3;++c){ double Hc[3][3]; for(int a=0;a<3;++a)for(int b=0;b<3;++b)Hc[a][b]=(b==c)?rhs[a]:H[a][b]; d[c]=det3(Hc)/dH; }
        double nz[3]={z[0]-d[0],z[1]-d[1],z[2]-d[2]}; clamp(nz);
        const double mv=std::fabs(nz[0]-z[0])+std::fabs(nz[1]-z[1])+std::fabs(nz[2]-z[2]);
        z[0]=nz[0]; z[1]=nz[1]; z[2]=nz[2];
        if (mv < 1e-13) break;
    }
    xi0[0]=z[0]; xi0[1]=z[1]; xi0[2]=z[2];
}

// CURVED-panel VOLUME-charge inner potential  INT_curvedtet xi^e0 eta^e1 zeta^e2 / |p - X(xi)|  dV_curved
// via the reference Duffy from xi0 (4 SIGNED 3D reference sub-tets), evaluating X(xi) + the curved volume
// element Jv=|det dX/dxi| per point.  The REFERENCE tet is always +oriented and Jv=|det| -> NO host-sign
// correction (unlike the flat physical tet Duffy).  gl/gw = an nq-pt Gauss-Legendre rule on [0,1].
double CurvedTetPotential(const double nodes[10][3], int e0, int e1, int e2, const double p[3],
                          const double* gl, const double* gw, int nq, bool include_measure)
{
    double xi0[3]; ClosestRefTet(nodes, p, xi0);
    static const double C[4][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
    static const int FC[4][3] = {{1,2,3},{0,3,2},{0,1,3},{2,1,0}};
    double acc = 0.0;
    for (int f = 0; f < 4; ++f) {
        // A product Duffy rule resolves the sharp angular variation near a
        // curved boundary much better than a fixed-degree simplex cubature.
        // Average its three cyclic face orderings: the nested w rule already
        // makes the other two vertices exchange-symmetric, so this is fully
        // permutation invariant without evaluating all six permutations.
        for (int lead = 0; lead < 3; ++lead) {
            const double* b1=C[FC[f][lead]];
            const double* b2=C[FC[f][(lead+1)%3]];
            const double* b3=C[FC[f][(lead+2)%3]];
            double d1[3],d2[3],d3[3],e21[3],e32[3];
            for (int k=0;k<3;++k){ d1[k]=b1[k]-xi0[k]; d2[k]=b2[k]-xi0[k]; d3[k]=b3[k]-xi0[k];
                                   e21[k]=b2[k]-b1[k]; e32[k]=b3[k]-b2[k]; }
            const double cr[3]={d2[1]*d3[2]-d2[2]*d3[1], d2[2]*d3[0]-d2[0]*d3[2],
                                d2[0]*d3[1]-d2[1]*d3[0]};
            const double D=d1[0]*cr[0]+d1[1]*cr[1]+d1[2]*cr[2];
            if (std::fabs(D) < 1e-300) continue;
            for (int a=0;a<nq;++a){ const double u=gl[a];
                for (int b=0;b<nq;++b){ const double v=gl[b];
                    for (int c=0;c<nq;++c){ const double w=gl[c];
                        double z[3]; for (int k=0;k<3;++k) z[k]=xi0[k]+u*(d1[k]+v*(e21[k]+w*e32[k]));
                        double X[3], Jac[3][3]; CurvedTetEval(nodes, z[0], z[1], z[2], X, Jac);
                        const double Jv = std::fabs(det3(Jac));
                        const double dx=p[0]-X[0], dy=p[1]-X[1], dz=p[2]-X[2];
                        const double r=std::sqrt(dx*dx+dy*dy+dz*dz);
                        if (r<1e-300) continue;
                        const double measure = include_measure ? Jv : 1.0;
                        acc += (gw[a]*gw[b]*gw[c]/3.0)*(u*u*v*D)*measure
                             *_ipow(z[0],e0)*_ipow(z[1],e1)*_ipow(z[2],e2)/r;
                    }}}
        }
    }
    return acc;
}

// ---- OUTER-quadrature helpers for the curved charge Gram: curved physical point X(xi) + curved MEASURE at a
//      reference point (area element for a face, volume element for a cell).  Reuse CurvedTri/TetEval. --------
void CurvedTriMapMeasure(const double nodes[6][3], double xi, double eta, double X[3], double& dA)
{
    double Xu[3], Xv[3]; CurvedTriEval(nodes, xi, eta, X, Xu, Xv);
    const double cr[3] = {Xu[1]*Xv[2]-Xu[2]*Xv[1], Xu[2]*Xv[0]-Xu[0]*Xv[2], Xu[0]*Xv[1]-Xu[1]*Xv[0]};
    dA = std::sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
}

void CurvedTetMapMeasure(const double nodes[10][3], double xi, double eta, double zeta, double X[3], double& dV)
{
    double Jac[3][3]; CurvedTetEval(nodes, xi, eta, zeta, X, Jac);
    dV = std::fabs(det3(Jac));
}

} // namespace rad_hdiv
