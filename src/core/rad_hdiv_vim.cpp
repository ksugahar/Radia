/* rad_hdiv_vim.cpp -- structured-hex RT0 topology + charge map B + dense Coulomb Gram G +
 * symmetric demag operator N = B^T G B.  See rad_hdiv_vim.h.  Hand-enumerated (no NGSolve),
 * validated against the NGSolve prototype golden (3x3x3 -> ndof=108, n_loop=28). */
#include "rad_hdiv_vim.h"
#include <cmath>

namespace rad_hdiv {

static const double PI = 3.14159265358979323846;
static const double INV_FOUR_PI = 1.0 / (4.0 * PI);
static const double C_CUBE = 1.88231;   // <1/r>_unitcube  (cube self-energy constant)
static const double C_SQ   = 2.97321;   // <1/r>_unitsquare (square self-energy constant)

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

// ---- trilinear hex geometry (NGSolve vertex order v0(000)..v7(110) in (x,y,z) local) ----
static void hex_shape(double xi, double eta, double ze, double N[8], double dN[8][3])
{
    double x0=1-xi, x1=xi, y0=1-eta, y1=eta, z0=1-ze, z1=ze;
    N[0]=x0*y0*z0; N[1]=x0*y0*z1; N[2]=x0*y1*z1; N[3]=x0*y1*z0;
    N[4]=x1*y0*z0; N[5]=x1*y0*z1; N[6]=x1*y1*z1; N[7]=x1*y1*z0;
    dN[0][0]=-y0*z0; dN[1][0]=-y0*z1; dN[2][0]=-y1*z1; dN[3][0]=-y1*z0;
    dN[4][0]= y0*z0; dN[5][0]= y0*z1; dN[6][0]= y1*z1; dN[7][0]= y1*z0;
    dN[0][1]=-x0*z0; dN[1][1]=-x0*z1; dN[2][1]= x0*z1; dN[3][1]= x0*z0;
    dN[4][1]=-x1*z0; dN[5][1]=-x1*z1; dN[6][1]= x1*z1; dN[7][1]= x1*z0;
    dN[0][2]=-x0*y0; dN[1][2]= x0*y0; dN[2][2]= x0*y1; dN[3][2]=-x0*y1;
    dN[4][2]=-x1*y0; dN[5][2]= x1*y0; dN[6][2]= x1*y1; dN[7][2]=-x1*y1;
}
static Vec3 hex_map(const std::array<Vec3,8>& V, double xi, double eta, double ze)
{
    double N[8], dN[8][3]; hex_shape(xi, eta, ze, N, dN);
    Vec3 p{0,0,0};
    for (int i=0;i<8;i++) for (int a=0;a<3;a++) p[a]+=N[i]*V[i][a];
    return p;
}
static double hex_detJ(const std::array<Vec3,8>& V, double xi, double eta, double ze)
{
    double N[8], dN[8][3]; hex_shape(xi, eta, ze, N, dN);
    double J[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for (int i=0;i<8;i++) for (int a=0;a<3;a++) for (int b=0;b<3;b++) J[a][b]+=V[i][a]*dN[i][b];
    return J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
         - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
         + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
}
static const double GP2[2] = {0.5 - 0.5/1.7320508075688772, 0.5 + 0.5/1.7320508075688772};
static double hex_volume(const std::array<Vec3,8>& V)
{
    double vol = 0.0;
    for (int a=0;a<2;a++) for (int b=0;b<2;b++) for (int c=0;c<2;c++)
        vol += 0.125 * std::fabs(hex_detJ(V, GP2[a], GP2[b], GP2[c]));
    return vol;
}
// ---- bilinear quad geometry (corner order v0(0,0) v1(0,1) v2(1,1) v3(1,0)) ----
static Vec3 quad_map(const Vec3 W[4], double u, double v)
{
    double a0=(1-u)*(1-v), a1=(1-u)*v, a2=u*v, a3=u*(1-v);
    Vec3 p; for (int a=0;a<3;a++) p[a]=a0*W[0][a]+a1*W[1][a]+a2*W[2][a]+a3*W[3][a];
    return p;
}
static double quad_area(const Vec3 W[4])
{
    double A=0.0; const double e=1e-6;
    for (int a=0;a<2;a++) for (int b=0;b<2;b++) {
        Vec3 xu, xv, pu0=quad_map(W,GP2[a]-e,GP2[b]), pu1=quad_map(W,GP2[a]+e,GP2[b]);
        Vec3 pv0=quad_map(W,GP2[a],GP2[b]-e), pv1=quad_map(W,GP2[a],GP2[b]+e);
        for (int d=0;d<3;d++){ xu[d]=(pu1[d]-pu0[d])/(2*e); xv[d]=(pv1[d]-pv0[d])/(2*e); }
        Vec3 cr{ xu[1]*xv[2]-xu[2]*xv[1], xu[2]*xv[0]-xu[0]*xv[2], xu[0]*xv[1]-xu[1]*xv[0] };
        A += 0.25 * std::sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
    }
    return A;
}

Mesh BuildStructuredRT0(int nx, int ny, int nz, double h, double distort)
{
    Mesh m; m.nx = nx; m.ny = ny; m.nz = nz; m.n_cell = nx*ny*nz;
    auto cell_id = [=](int i, int j, int k) { return (i*ny + j)*nz + k; };
    const double L = nx * h, d = distort;
    auto node = [=](int i, int j, int k) -> Vec3 {
        double x=i*h, y=j*h, z=k*h;
        if (d == 0.0) return {x, y, z};
        double sy = std::sin(PI*y/L);
        return { x + d*sy*z, y + 0.83*d*x*z, z + 0.67*d*x*sy };   // smooth node displacement
    };
    auto mean = [](const Vec3* w, int n) -> Vec3 {
        Vec3 c{0,0,0}; for (int t=0;t<n;t++) for (int a=0;a<3;a++) c[a]+=w[t][a]/n; return c; };

    // cells: 8 corner vertices (NGSolve order), centroid, volume
    m.cell_verts.resize(m.n_cell); m.cell_c.resize(m.n_cell); m.cell_V.resize(m.n_cell);
    static const int CV[8][3] = {{0,0,0},{0,0,1},{0,1,1},{0,1,0},{1,0,0},{1,0,1},{1,1,1},{1,1,0}};
    for (int i=0;i<nx;i++) for (int j=0;j<ny;j++) for (int k=0;k<nz;k++) {
        int c = cell_id(i,j,k);
        for (int t=0;t<8;t++) m.cell_verts[c][t] = node(i+CV[t][0], j+CV[t][1], k+CV[t][2]);
        m.cell_c[c] = mean(m.cell_verts[c].data(), 8);
        m.cell_V[c] = hex_volume(m.cell_verts[c]);
    }
    // faces: 4 corner vertices (around the quad), centroid, area; topological lo/hi/bnd by axis
    auto add_face = [&](int ax, int i, int j, int k, int lo, int hi) {
        Face f; f.ax = ax; f.lo = lo; f.hi = hi; f.bnd = (lo < 0 || hi < 0);
        int u0=(ax+1)%3, u1=(ax+2)%3;
        static const int Q[4][2] = {{0,0},{0,1},{1,1},{1,0}};
        for (int t=0;t<4;t++) {
            int off[3]={0,0,0}; off[u0]=Q[t][0]; off[u1]=Q[t][1];
            f.v[t] = node(i+off[0], j+off[1], k+off[2]);
        }
        f.c = mean(f.v, 4); f.area = quad_area(f.v);
        m.faces.push_back(f);
    };
    for (int i=0;i<=nx;i++) for (int j=0;j<ny;j++) for (int k=0;k<nz;k++)
        add_face(0, i, j, k, (i>0)?cell_id(i-1,j,k):-1, (i<nx)?cell_id(i,j,k):-1);
    for (int i=0;i<nx;i++) for (int j=0;j<=ny;j++) for (int k=0;k<nz;k++)
        add_face(1, i, j, k, (j>0)?cell_id(i,j-1,k):-1, (j<ny)?cell_id(i,j,k):-1);
    for (int i=0;i<nx;i++) for (int j=0;j<ny;j++) for (int k=0;k<=nz;k++)
        add_face(2, i, j, k, (k>0)?cell_id(i,j,k-1):-1, (k<nz)?cell_id(i,j,k):-1);
    return m;
}

// boundary-face -> its sigma charge row (n_cell + running index); -1 if not boundary.
static std::vector<int> bnd_charge_rows(const Mesh& m, int& n_bnd)
{
    std::vector<int> row(m.n_face(), -1);
    int r = 0;
    for (int f = 0; f < m.n_face(); ++f)
        if (m.faces[f].bnd) row[f] = m.n_cell + (r++);
    n_bnd = r;
    return row;
}

void AssembleChargeMap(const Mesh& m, std::vector<double>& B, int& n_charge, int& n_bnd)
{
    std::vector<int> brow = bnd_charge_rows(m, n_bnd);
    n_charge = m.n_cell + n_bnd;
    const int nf = m.n_face();
    B.assign((size_t)n_charge * nf, 0.0);
    for (int f = 0; f < nf; ++f) {
        const Face& fc = m.faces[f];
        // cell on the LO side: this face is that cell's HI face -> global normal points OUT -> +1
        // div contribution; rho = -div => -(+1) = -1, per unit volume.
        if (fc.lo >= 0) B[(size_t)fc.lo * nf + f] += -(+1.0) / m.cell_V[fc.lo];
        // cell on the HI side: this face is that cell's LO face -> normal points IN -> -1 div;
        // rho = -(-1) = +1, per unit volume.
        if (fc.hi >= 0) B[(size_t)fc.hi * nf + f] += -(-1.0) / m.cell_V[fc.hi];
        // boundary face: sigma = M . n_OUTWARD.  Global face normal is +axis; outward (out of the
        // domain) is +global if the cell sits on the LO side (domain HIGH boundary), -global if on
        // the HI side (domain LOW boundary).  Using the global normal for all boundary faces flips
        // sigma on the low boundary -> a spurious monopole surface charge -> unphysical demag
        // factors (>1).  (Symmetry + loop-nullity do NOT catch this; the physics test does.)
        if (fc.bnd) {
            double out_sign = (fc.lo >= 0) ? 1.0 : -1.0;
            B[(size_t)brow[f] * nf + f] += out_sign / fc.area;
        }
    }
}

void BuildChargeMapCSC(const Mesh& m, ChargeMapCSC& csc)
{
    int n_bnd; std::vector<int> brow = bnd_charge_rows(m, n_bnd);
    csc.n_charge = m.n_cell + n_bnd;
    const int nf = m.n_face();
    csc.rows.assign(nf, std::array<int, 2>{ -1, -1 });
    csc.coef.assign(nf, std::array<double, 2>{ 0.0, 0.0 });
    for (int f = 0; f < nf; ++f) {
        const Face& fc = m.faces[f];
        int k = 0;   // <=2: interior face -> (lo, hi); boundary face -> (its one cell, sigma)
        if (fc.lo >= 0) { csc.rows[f][k] = fc.lo;     csc.coef[f][k] = -(+1.0) / m.cell_V[fc.lo]; ++k; }
        if (fc.hi >= 0) { csc.rows[f][k] = fc.hi;     csc.coef[f][k] = -(-1.0) / m.cell_V[fc.hi]; ++k; }
        if (fc.bnd) {
            double out_sign = (fc.lo >= 0) ? 1.0 : -1.0;
            csc.rows[f][k] = brow[f];                 csc.coef[f][k] = out_sign / fc.area;        ++k;
        }
    }
}

void BuildChargeQuad(const Mesh& m, int nsub, ChargeQuad& q)
{
    int n_bnd; std::vector<int> brow = bnd_charge_rows(m, n_bnd);
    q.n_cell   = m.n_cell;
    q.n_charge = m.n_cell + n_bnd;
    // centroids + measures: volume cells then boundary faces (same row order as AssembleChargeMap)
    q.cent.assign(q.n_charge, Vec3{0,0,0});
    q.meas.assign(q.n_charge, 0.0);
    for (int c = 0; c < m.n_cell; ++c) { q.cent[c] = m.cell_c[c]; q.meas[c] = m.cell_V[c]; }
    for (int f = 0; f < m.n_face(); ++f)
        if (m.faces[f].bnd) { int r = brow[f]; q.cent[r] = m.faces[f].c; q.meas[r] = m.faces[f].area; }

    q.sp.clear(); q.sw.clear();
    if (nsub <= 0) return;   // centroid-monopole mode: no sub-points

    // ---- accurate: sub-point quadrature on the ACTUAL geometry (trilinear hex / bilinear quad).
    // Volume cell -> nsub^3 trilinear sub-points (|detJ| sub-weights); boundary face -> nsub^2
    // bilinear sub-points (|x_u x x_v| sub-weights).  Per-sub-point weights (not uniform) -> the
    // sub-point cloud fills the actual DISTORTED cell, making the self-energy + G geometry-exact. --
    q.sp.assign(q.n_charge, std::vector<Vec3>());
    q.sw.assign(q.n_charge, std::vector<double>());
    for (int c = 0; c < m.n_cell; ++c) {
        double inv = 1.0 / (double)(nsub*nsub*nsub);
        const std::array<Vec3,8>& V = m.cell_verts[c];
        q.sp[c].reserve((size_t)nsub*nsub*nsub); q.sw[c].reserve((size_t)nsub*nsub*nsub);
        for (int i = 0; i < nsub; ++i)
            for (int j = 0; j < nsub; ++j)
                for (int k = 0; k < nsub; ++k) {
                    double xi=(i+0.5)/nsub, eta=(j+0.5)/nsub, ze=(k+0.5)/nsub;
                    q.sp[c].push_back(hex_map(V, xi, eta, ze));
                    q.sw[c].push_back(std::fabs(hex_detJ(V, xi, eta, ze)) * inv);
                }
    }
    const double e = 1e-6;
    for (int f = 0; f < m.n_face(); ++f) {
        if (!m.faces[f].bnd) continue;
        const Face& fc = m.faces[f]; int r = brow[f];
        double inv = 1.0 / (double)(nsub*nsub);
        q.sp[r].reserve((size_t)nsub*nsub); q.sw[r].reserve((size_t)nsub*nsub);
        for (int i = 0; i < nsub; ++i)
            for (int j = 0; j < nsub; ++j) {
                double u=(i+0.5)/nsub, v=(j+0.5)/nsub;
                q.sp[r].push_back(quad_map(fc.v, u, v));
                Vec3 pu0=quad_map(fc.v,u-e,v), pu1=quad_map(fc.v,u+e,v);
                Vec3 pv0=quad_map(fc.v,u,v-e), pv1=quad_map(fc.v,u,v+e);
                Vec3 xu, xv; for (int d=0;d<3;d++){ xu[d]=(pu1[d]-pu0[d])/(2*e); xv[d]=(pv1[d]-pv0[d])/(2*e); }
                Vec3 cr{ xu[1]*xv[2]-xu[2]*xv[1], xu[2]*xv[0]-xu[0]*xv[2], xu[0]*xv[1]-xu[1]*xv[0] };
                q.sw[r].push_back(std::sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]) * inv);
            }
    }
}

double CoulombGramEntry(const ChargeQuad& q, int a, int b)
{
    const bool accurate = !q.sp.empty();
    if (!accurate) {
        // ---- centroid-monopole off-diagonal + cube/square self-energy ----
        if (a == b)
            return (a < q.n_cell ? C_CUBE * std::pow(q.meas[a], 5.0/3.0)
                                 : C_SQ   * std::pow(q.meas[a], 1.5)) * INV_FOUR_PI;
        double dx = q.cent[a][0]-q.cent[b][0], dy = q.cent[a][1]-q.cent[b][1], dz = q.cent[a][2]-q.cent[b][2];
        return q.meas[a]*q.meas[b] * INV_FOUR_PI / std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    // ---- accurate: sub-point quadrature ----
    const std::vector<Vec3>&   pa = q.sp[a]; const std::vector<double>& wa = q.sw[a];
    if (a == b) {
        bool a_vol = (a < q.n_cell);
        double diag = 0.0;
        for (size_t p = 0; p < pa.size(); ++p) {
            for (size_t r = 0; r < pa.size(); ++r) {
                if (p == r) continue;
                double dx = pa[p][0]-pa[r][0], dy = pa[p][1]-pa[r][1], dz = pa[p][2]-pa[r][2];
                diag += wa[p]*wa[r] * INV_FOUR_PI / std::sqrt(dx*dx + dy*dy + dz*dz);
            }
            diag += (a_vol ? C_CUBE * std::pow(wa[p], 5.0/3.0) : C_SQ * std::pow(wa[p], 1.5)) * INV_FOUR_PI;
        }
        return diag;
    }
    const std::vector<Vec3>&   pb = q.sp[b]; const std::vector<double>& wb = q.sw[b];
    double g = 0.0;
    for (size_t p = 0; p < pa.size(); ++p)
        for (size_t r = 0; r < pb.size(); ++r) {
            double dx = pa[p][0]-pb[r][0], dy = pa[p][1]-pb[r][1], dz = pa[p][2]-pb[r][2];
            g += wa[p]*wb[r] * INV_FOUR_PI / std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    return g;
}

void AssembleCoulombGram(const Mesh& m, std::vector<double>& G, int& n_charge, int nsub)
{
    ChargeQuad q; BuildChargeQuad(m, nsub, q);
    n_charge = q.n_charge;
    G.assign((size_t)n_charge * n_charge, 0.0);
    for (int a = 0; a < n_charge; ++a) {
        G[(size_t)a*n_charge + a] = CoulombGramEntry(q, a, a);
        for (int b = a+1; b < n_charge; ++b) {
            double g = CoulombGramEntry(q, a, b);
            G[(size_t)a*n_charge + b] = g;
            G[(size_t)b*n_charge + a] = g;
        }
    }
}

void AssembleN(const Mesh& m, std::vector<double>& N, int nsub)
{
    std::vector<double> B, G;
    int n_charge, n_bnd, n_charge_g;
    AssembleChargeMap(m, B, n_charge, n_bnd);
    AssembleCoulombGram(m, G, n_charge_g, nsub);
    const int nf = m.n_face();
    // N = B^T G B, row-major (nf x nf).  GB (n_charge x nf) first, then B^T (GB).
    std::vector<double> GB((size_t)n_charge * nf, 0.0);
    for (int a = 0; a < n_charge; ++a)
        for (int c = 0; c < n_charge; ++c) {
            double gac = G[(size_t)a*n_charge + c];
            if (gac == 0.0) continue;
            const double* Brow = &B[(size_t)c * nf];
            double* GBrow = &GB[(size_t)a * nf];
            for (int f = 0; f < nf; ++f) GBrow[f] += gac * Brow[f];
        }
    N.assign((size_t)nf * nf, 0.0);
    for (int a = 0; a < n_charge; ++a) {
        const double* Brow  = &B[(size_t)a * nf];
        const double* GBrow = &GB[(size_t)a * nf];
        for (int i = 0; i < nf; ++i) {
            double bi = Brow[i];
            if (bi == 0.0) continue;
            double* Nrow = &N[(size_t)i * nf];
            for (int j = 0; j < nf; ++j) Nrow[j] += bi * GBrow[j];
        }
    }
}

void AssembleMass(const Mesh& m, std::vector<double>& M_mass)
{
    const int nf = m.n_face();
    M_mass.assign((size_t)nf * nf, 0.0);
    // per cell, per axis: the (lo_face, hi_face) pair gets the 2x2 block (1/h)[[1/3,1/6],[1/6,1/3]].
    // Find each cell's lo/hi face per axis from the face table (a face is its lo-cell's HI face and
    // its hi-cell's LO face).
    std::vector<std::array<int, 2>> cell_axis_faces((size_t)m.n_cell * 3, {-1, -1});
    auto CAF = [&](int c, int ax) -> std::array<int, 2>& { return cell_axis_faces[(size_t)c * 3 + ax]; };
    for (int f = 0; f < nf; ++f) {
        const Face& fc = m.faces[f];
        if (fc.lo >= 0) CAF(fc.lo, fc.ax)[1] = f;   // cell on LO side -> face is its HI face
        if (fc.hi >= 0) CAF(fc.hi, fc.ax)[0] = f;   // cell on HI side -> face is its LO face
    }
    for (int c = 0; c < m.n_cell; ++c) {
        // cell size h from its volume (structured cubic cell): h = V^{1/3}
        double h = std::cbrt(m.cell_V[c]);
        double d = (1.0 / h) * (1.0 / 3.0), o = (1.0 / h) * (1.0 / 6.0);
        for (int ax = 0; ax < 3; ++ax) {
            std::array<int, 2> idx = CAF(c, ax);
            int lo = idx[0], hi = idx[1];
            M_mass[(size_t)lo * nf + lo] += d;
            M_mass[(size_t)hi * nf + hi] += d;
            M_mass[(size_t)lo * nf + hi] += o;
            M_mass[(size_t)hi * nf + lo] += o;
        }
    }
}

void BuildMassCOO(const Mesh& m, std::vector<int>& I, std::vector<int>& J,
                  std::vector<double>& V, std::vector<double>& diag)
{
    const int nf = m.n_face();
    I.clear(); J.clear(); V.clear();
    diag.assign(nf, 0.0);
    std::vector<std::array<int, 2>> cell_axis_faces((size_t)m.n_cell * 3, {-1, -1});
    auto CAF = [&](int c, int ax) -> std::array<int, 2>& { return cell_axis_faces[(size_t)c * 3 + ax]; };
    for (int f = 0; f < nf; ++f) {
        const Face& fc = m.faces[f];
        if (fc.lo >= 0) CAF(fc.lo, fc.ax)[1] = f;
        if (fc.hi >= 0) CAF(fc.hi, fc.ax)[0] = f;
    }
    auto push = [&](int i, int j, double v) { I.push_back(i); J.push_back(j); V.push_back(v); };
    for (int c = 0; c < m.n_cell; ++c) {
        double h = std::cbrt(m.cell_V[c]);
        double d = (1.0 / h) * (1.0 / 3.0), o = (1.0 / h) * (1.0 / 6.0);
        for (int ax = 0; ax < 3; ++ax) {
            std::array<int, 2> idx = CAF(c, ax);
            int lo = idx[0], hi = idx[1];
            push(lo, lo, d); diag[lo] += d;
            push(hi, hi, d); diag[hi] += d;
            push(lo, hi, o);
            push(hi, lo, o);
        }
    }
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
                          const double* gl, const double* gw, int nq)
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
                acc += gw[a]*gw[b]*(u*sgn2)*J*_ipow(xi, e0)*_ipow(eta, e1)/r;
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
                          const double* gl, const double* gw, int nq)
{
    double xi0[3]; ClosestRefTet(nodes, p, xi0);
    static const double C[4][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
    static const int FC[4][3] = {{1,2,3},{0,3,2},{0,1,3},{2,1,0}};
    double acc = 0.0;
    for (int f = 0; f < 4; ++f) {
        const double* b1=C[FC[f][0]]; const double* b2=C[FC[f][1]]; const double* b3=C[FC[f][2]];
        double d1[3],d2[3],d3[3],e21[3],e32[3];
        for (int k=0;k<3;++k){ d1[k]=b1[k]-xi0[k]; d2[k]=b2[k]-xi0[k]; d3[k]=b3[k]-xi0[k];
                               e21[k]=b2[k]-b1[k]; e32[k]=b3[k]-b2[k]; }
        const double cr[3]={d2[1]*d3[2]-d2[2]*d3[1], d2[2]*d3[0]-d2[0]*d3[2], d2[0]*d3[1]-d2[1]*d3[0]};
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
                    acc += gw[a]*gw[b]*gw[c]*(u*u*v*D)*Jv*_ipow(z[0],e0)*_ipow(z[1],e1)*_ipow(z[2],e2)/r;
                }}}
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
