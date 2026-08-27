/* rad_hdiv_vim.cpp -- analytic charge-potential kernels for BDM1 HDiv-VIM. */
#include "rad_hdiv_vim.h"
#include <cmath>
#include <array>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <map>
#include <vector>

#include "rad_parallel.h"

namespace rad_hdiv {

std::vector<double> AffineCellSelfEnergyShapeDerivative(
    int cell_type, const std::vector<double>& nodes,
    const std::vector<double>& node_velocities, int n_modes)
{
    const int nn = cell_type == 0 ? 4 : (cell_type == 1 ? 8 : (cell_type == 2 ? 6 : 0));
    if (!nn) throw std::invalid_argument("cell_type must be 0 (TET), 1 (HEX), or 2 (WEDGE)");
    if ((int)nodes.size() != 3*nn)
        throw std::invalid_argument("nodes has the wrong size for cell_type");
    if (n_modes < 0 || (int)node_velocities.size() != n_modes*3*nn)
        throw std::invalid_argument("node_velocities must have shape (n_modes,nodes,3)");
    for (double x : nodes) if (!std::isfinite(x))
        throw std::invalid_argument("nodes must be finite");
    for (double x : node_velocities) if (!std::isfinite(x))
        throw std::invalid_argument("node_velocities must be finite");

    std::vector<std::array<int,4>> tets;
    if (cell_type == 0) tets = {{0,1,2,3}};
    else if (cell_type == 1) tets = {
        {0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}};
    else tets = {{0,1,2,3},{1,2,4,3},{2,4,5,3}};

    auto tet_vertices = [&](const std::array<int,4>& ti, double V[4][3]) {
        for (int i=0;i<4;++i) for (int k=0;k<3;++k) V[i][k]=nodes[3*ti[i]+k];
    };
    auto det6 = [](const double V[4][3]) {
        const double ax=V[1][0]-V[0][0], ay=V[1][1]-V[0][1], az=V[1][2]-V[0][2];
        const double bx=V[2][0]-V[0][0], by=V[2][1]-V[0][1], bz=V[2][2]-V[0][2];
        const double cx=V[3][0]-V[0][0], cy=V[3][1]-V[0][1], cz=V[3][2]-V[0][2];
        return ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)+az*(bx*cy-by*cx);
    };
    for (const auto& ti:tets) { double V[4][3]; tet_vertices(ti,V); if (std::fabs(det6(V)) < 1e-18)
        throw std::invalid_argument("degenerate affine cell decomposition"); }

    // Same 4^3 smooth outer rule used by the production affine TET ChargeGram.
    static const double gx[4]={0.06943184420297371,0.33000947820757187,0.66999052179242813,0.93056815579702629};
    static const double gw[4]={0.17392742256872693,0.32607257743127307,0.32607257743127307,0.17392742256872693};
    constexpr double inv4pi=0.07957747154594766788;
    double energy=0.0;
    for (const auto& tt:tets) {
        double VT[4][3]; tet_vertices(tt,VT); const double jac=std::fabs(det6(VT));
        for(int ia=0;ia<4;++ia) for(int ib=0;ib<4;++ib) for(int ic=0;ic<4;++ic) {
            const double a=gx[ia], b=gx[ib], c=gx[ic];
            const double l1=a, l2=b*(1-a), l3=c*(1-a)*(1-b);
            double p[3]; for(int k=0;k<3;++k) p[k]=VT[0][k]+l1*(VT[1][k]-VT[0][k])+l2*(VT[2][k]-VT[0][k])+l3*(VT[3][k]-VT[0][k]);
            double phi=0.0; for(const auto& ts:tets) { double VS[4][3]; tet_vertices(ts,VS); phi+=PhiTet(VS,p); }
            energy += gw[ia]*gw[ib]*gw[ic]*(1-a)*(1-a)*(1-b)*jac*phi;
        }
    }

    // Extract boundary triangles: a sorted face occurring once is on dD.
    struct Face { std::array<int,3> oriented; int count=0; };
    std::map<std::array<int,3>,Face> faces;
    static const int lf[4][3]={{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
    for(const auto& t:tets) for(const auto& f:lf) {
        std::array<int,3> o={t[f[0]],t[f[1]],t[f[2]]}, key=o;
        std::sort(key.begin(),key.end()); auto& rec=faces[key]; if(rec.count++==0) rec.oriented=o;
    }
    double center[3]={0,0,0}; for(int i=0;i<nn;++i) for(int k=0;k<3;++k) center[k]+=nodes[3*i+k]/nn;
    std::vector<double> out((size_t)n_modes+1,0.0); out[0]=energy*inv4pi;
    std::vector<double> mean_velocity((size_t)n_modes*3,0.0);
    for(int m=0;m<n_modes;++m) for(int i=0;i<nn;++i) for(int k=0;k<3;++k)
        mean_velocity[3*m+k]+=node_velocities[(m*nn+i)*3+k]/nn;
    static const double dun[7][4]={{1./3,1./3,1./3,.225},{.0597158717,.4701420641,.4701420641,.1323941527},{.4701420641,.0597158717,.4701420641,.1323941527},{.4701420641,.4701420641,.0597158717,.1323941527},{.7974269853,.1012865073,.1012865073,.1259391805},{.1012865073,.7974269853,.1012865073,.1259391805},{.1012865073,.1012865073,.7974269853,.1259391805}};
    for(const auto& kv:faces) if(kv.second.count==1) {
        auto f=kv.second.oriented; double A[3],B[3],C[3]; for(int k=0;k<3;++k){A[k]=nodes[3*f[0]+k];B[k]=nodes[3*f[1]+k];C[k]=nodes[3*f[2]+k];}
        double e1[3]={B[0]-A[0],B[1]-A[1],B[2]-A[2]},e2[3]={C[0]-A[0],C[1]-A[1],C[2]-A[2]};
        double av[3]={.5*(e1[1]*e2[2]-e1[2]*e2[1]),.5*(e1[2]*e2[0]-e1[0]*e2[2]),.5*(e1[0]*e2[1]-e1[1]*e2[0])};
        double fc[3]={(A[0]+B[0]+C[0])/3,(A[1]+B[1]+C[1])/3,(A[2]+B[2]+C[2])/3};
        if(av[0]*(fc[0]-center[0])+av[1]*(fc[1]-center[1])+av[2]*(fc[2]-center[2])<0) for(double& z:av) z=-z;
        for(const auto& q:dun) {
            double p[3]; for(int k=0;k<3;++k)p[k]=q[0]*A[k]+q[1]*B[k]+q[2]*C[k];
            double phi=0; for(const auto& ts:tets){double V[4][3];tet_vertices(ts,V);phi+=PhiTet(V,p);}
            for(int m=0;m<n_modes;++m){double vn=0;for(int k=0;k<3;++k){double v=q[0]*node_velocities[(m*nn+f[0])*3+k]+q[1]*node_velocities[(m*nn+f[1])*3+k]+q[2]*node_velocities[(m*nn+f[2])*3+k]-mean_velocity[3*m+k];vn+=v*av[k];}out[m+1]+=2*inv4pi*q[3]*vn*phi;}
        }
    }
    for(int m=0;m<n_modes;++m) {
        bool translation=true;
        for(int i=1;i<nn&&translation;++i) for(int k=0;k<3;++k)
            if(node_velocities[(m*nn+i)*3+k]!=node_velocities[(m*nn)*3+k]) { translation=false; break; }
        if(translation) out[m+1]=0.0;
    }
    return out;
}

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
    out[0]=out[1]=out[2]=0.0;
    double e1[3],e2[3],n[3];
    for (int k=0;k<3;k++){ e1[k]=V[1][k]-V[0][k]; e2[k]=V[2][k]-V[0][k]; }
    v3cross(e1,e2,n); double nl=v3nrm(n);
    // Degenerate (zero-area) face: TriPotential/TriField already return the
    // zero limit here, so match them instead of dividing by zero and pushing
    // NaN into every moment that consumes this face.
    if (nl<1e-300) return;
    for (int k=0;k<3;k++) n[k]/=nl;
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
    for (int a=0;a<3;a++) for (int b=0;b<3;b++) out[a][b]=0.0;
    double e1[3],e2[3],n[3];
    for (int k=0;k<3;k++){ e1[k]=V[1][k]-V[0][k]; e2[k]=V[2][k]-V[0][k]; }
    v3cross(e1,e2,n); double nl=v3nrm(n);
    // Degenerate (zero-area) face: TriPotential/TriField already return the
    // zero limit here, so match them instead of dividing by zero and pushing
    // NaN into every moment that consumes this face.
    if (nl<1e-300) return;
    for (int k=0;k<3;k++) n[k]/=nl;
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

// Forward directional scalar used by the affine-face analytic moment path.
// Branches are selected by value; away from degenerate triangles this gives the
// exact directional derivative of the same closed-form expression as double.
struct ValueDirectional {
    double value = 0.0;
    double direction = 0.0;
    ValueDirectional() = default;
    ValueDirectional(double v, double d=0.0) : value(v), direction(d) {}
};
static inline ValueDirectional operator+(ValueDirectional a,ValueDirectional b){return {a.value+b.value,a.direction+b.direction};}
static inline ValueDirectional operator-(ValueDirectional a,ValueDirectional b){return {a.value-b.value,a.direction-b.direction};}
static inline ValueDirectional operator-(ValueDirectional a){return {-a.value,-a.direction};}
static inline ValueDirectional operator*(ValueDirectional a,ValueDirectional b){return {a.value*b.value,a.direction*b.value+a.value*b.direction};}
static inline ValueDirectional operator/(ValueDirectional a,ValueDirectional b){return {a.value/b.value,(a.direction*b.value-a.value*b.direction)/(b.value*b.value)};}
static inline ValueDirectional vd_sqrt(ValueDirectional a){
    if(a.value<=1e-300)return {std::sqrt(std::max(0.0,a.value)),0.0};
    const double s=std::sqrt(a.value);return {s,a.direction/(2.0*s)};
}
static inline ValueDirectional vd_asinh(ValueDirectional a){return {std::asinh(a.value),a.direction/std::sqrt(1.0+a.value*a.value)};}
static inline ValueDirectional vd_abs(ValueDirectional a){return {std::fabs(a.value),(a.value<0.0?-1.0:1.0)*a.direction};}
static inline ValueDirectional vd_log(ValueDirectional a){return {std::log(a.value),a.direction/a.value};}
static inline ValueDirectional vd_atan2(ValueDirectional y, ValueDirectional x){
    const double den=x.value*x.value+y.value*y.value;
    if(den<=1e-300)return {std::atan2(y.value,x.value),0.0};
    return {std::atan2(y.value,x.value),(x.value*y.direction-y.value*x.direction)/den};
}
static inline ValueDirectional vd_pow(ValueDirectional a,int n){if(n==0)return {1.0,0.0};const double p=std::pow(a.value,n);return {p,n*std::pow(a.value,n-1)*a.direction};}

static inline ValueDirectional vd_dot(const ValueDirectional a[3],const ValueDirectional b[3])
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static inline void vd_cross(const ValueDirectional a[3],const ValueDirectional b[3],ValueDirectional o[3])
{ o[0]=a[1]*b[2]-a[2]*b[1];o[1]=a[2]*b[0]-a[0]*b[2];o[2]=a[0]*b[1]-a[1]*b[0]; }
static inline ValueDirectional vd_norm(const ValueDirectional a[3]) { return vd_sqrt(vd_dot(a,a)); }

static ValueDirectional TriPotentialDirectionalValue(
    const ValueDirectional V[3][3],const ValueDirectional r[3])
{
    ValueDirectional e1[3],e2[3],n[3];
    for(int k=0;k<3;++k){e1[k]=V[1][k]-V[0][k];e2[k]=V[2][k]-V[0][k];}
    vd_cross(e1,e2,n);const ValueDirectional nl=vd_norm(n);if(nl.value<1e-300)return {};
    for(auto& x:n)x=x/nl;
    ValueDirectional rmv0[3];for(int k=0;k<3;++k)rmv0[k]=r[k]-V[0][k];
    const ValueDirectional d=vd_dot(rmv0,n);
    ValueDirectional p[3];for(int k=0;k<3;++k)p[k]=r[k]-d*n[k];
    const ValueDirectional ad=vd_abs(d);ValueDirectional I;
    for(int i=0;i<3;++i){
        const ValueDirectional* a=V[i];const ValueDirectional* b=V[(i+1)%3];
        ValueDirectional lh[3];for(int k=0;k<3;++k)lh[k]=b[k]-a[k];
        const ValueDirectional ll=vd_norm(lh);if(ll.value<1e-300)continue;for(auto& x:lh)x=x/ll;
        ValueDirectional uh[3];vd_cross(lh,n,uh);
        ValueDirectional ap[3],bp[3],ra[3],rb[3];
        for(int k=0;k<3;++k){ap[k]=a[k]-p[k];bp[k]=b[k]-p[k];ra[k]=r[k]-a[k];rb[k]=r[k]-b[k];}
        const ValueDirectional P0=vd_dot(ap,uh),sm=vd_dot(ap,lh),sp=vd_dot(bp,lh);
        const ValueDirectional Rm=vd_norm(ra),Rp=vd_norm(rb),R0sq=P0*P0+d*d;
        const ValueDirectional dm=Rm+sm,dp=Rp+sp;
        const ValueDirectional f=(dp.value>1e-300&&dm.value>1e-300)?vd_log(dp/dm):ValueDirectional{};
        const ValueDirectional beta=vd_atan2(P0*sp,R0sq+ad*Rp)-vd_atan2(P0*sm,R0sq+ad*Rm);
        I=I+P0*f-ad*beta;
    }
    return I;
}

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

struct TriPolyEdgeDirectional { ValueDirectional xiA[2],t2[2],m2[2],L,l0,d2; };
struct TriPolySetupDirectional {
    ValueDirectional n[3],h,rp[3],e1[3],e2[3];
    TriPolyEdgeDirectional edges[3];
};

static bool tri_poly_setup_directional(const ValueDirectional P[3][3],
                                       const ValueDirectional r[3],TriPolySetupDirectional& g)
{
    ValueDirectional a1[3],a2[3];
    for(int k=0;k<3;++k){a1[k]=P[1][k]-P[0][k];a2[k]=P[2][k]-P[0][k];}
    vd_cross(a1,a2,g.n);const ValueDirectional nl=vd_norm(g.n);if(nl.value<1e-300)return false;
    for(auto& x:g.n)x=x/nl;
    ValueDirectional rv[3];for(int k=0;k<3;++k)rv[k]=r[k]-P[0][k];
    g.h=vd_dot(rv,g.n);for(int k=0;k<3;++k)g.rp[k]=r[k]-g.h*g.n[k];
    const ValueDirectional e1n=vd_norm(a1);if(e1n.value<1e-300)return false;
    for(int k=0;k<3;++k)g.e1[k]=a1[k]/e1n;vd_cross(g.n,g.e1,g.e2);
    ValueDirectional cen[3];for(int i=0;i<3;++i)for(int k=0;k<3;++k)cen[k]=cen[k]+P[i][k]/3.0;
    for(int i=0;i<3;++i){
        const ValueDirectional* A=P[i];const ValueDirectional* B=P[(i+1)%3];auto& e=g.edges[i];
        ValueDirectional t[3];for(int k=0;k<3;++k)t[k]=B[k]-A[k];e.L=vd_norm(t);if(e.L.value<1e-300)return false;
        ValueDirectional th[3];for(int k=0;k<3;++k)th[k]=t[k]/e.L;
        ValueDirectional m[3];vd_cross(th,g.n,m);
        ValueDirectional mid[3];for(int k=0;k<3;++k)mid[k]=(A[k]+B[k])/2.0-cen[k];
        if(vd_dot(m,mid).value<0.0)for(auto& x:m)x=-x;
        e.m2[0]=vd_dot(m,g.e1);e.m2[1]=vd_dot(m,g.e2);
        ValueDirectional Amrp[3];for(int k=0;k<3;++k)Amrp[k]=A[k]-g.rp[k];
        e.xiA[0]=vd_dot(Amrp,g.e1);e.xiA[1]=vd_dot(Amrp,g.e2);
        e.t2[0]=vd_dot(th,g.e1);e.t2[1]=vd_dot(th,g.e2);
        ValueDirectional w[3];for(int k=0;k<3;++k)w[k]=r[k]-A[k];
        e.l0=vd_dot(w,th);e.d2=vd_dot(w,w)-e.l0*e.l0;
        if(e.d2.value<0.0)e.d2={0.0,0.0};
    }
    return true;
}

static ValueDirectional edge_R_for_poly_directional(const TriPolyEdgeDirectional& e)
{
    const ValueDirectional u1=-e.l0,u2=e.L-e.l0,d2=e.d2;
    auto F=[&](ValueDirectional u){
        if(d2.value<1e-300)return 0.5*u*vd_abs(u);
        const ValueDirectional d=vd_sqrt(d2);
        return 0.5*(u*vd_sqrt(u*u+d2)+d2*vd_asinh(u/d));
    };
    return F(u2)-F(u1);
}

static void edge_l_moments_poly_directional(const TriPolyEdgeDirectional& e,int nmax,
                                             ValueDirectional Jl[POLY_MAX_DEG+3])
{
    nmax=std::min(nmax,POLY_MAX_DEG+2);const auto d2=e.d2,d=vd_sqrt(d2),u1=-e.l0,u2=e.L-e.l0;
    ValueDirectional W[POLY_MAX_DEG+3]{};
    auto asinh_safe=[&](ValueDirectional u){
        if(d.value>1e-300)return vd_asinh(u/d);
        if(std::fabs(u.value)==0.0)return ValueDirectional{};
        return (u.value>0.0?1.0:-1.0)*vd_log(2.0*vd_abs(u));
    };
    W[0]=asinh_safe(u2)-asinh_safe(u1);
    if(nmax>=1)W[1]=vd_sqrt(u2*u2+d2)-vd_sqrt(u1*u1+d2);
    for(int n=2;n<=nmax;++n){
        const auto term=(vd_pow(u2,n-1)*vd_sqrt(u2*u2+d2)-vd_pow(u1,n-1)*vd_sqrt(u1*u1+d2))/double(n);
        W[n]=term-((n-1.0)/n)*d2*W[n-2];
    }
    for(int n=0;n<=nmax;++n){ValueDirectional s;for(int i=0;i<=n;++i)s=s+small_comb(n,i)*vd_pow(e.l0,n-i)*W[i];Jl[n]=s;}
}

static ValueDirectional edge_inplane_monomial_poly_directional(
    const TriPolyEdgeDirectional& e,int a,int b,const ValueDirectional Jl[POLY_MAX_DEG+3])
{
    ValueDirectional poly[POLY_MAX_DEG+3]{},tmp[POLY_MAX_DEG+3]{};int deg=0;poly[0]={1.0,0.0};
    auto mul=[&](ValueDirectional c0,ValueDirectional c1){std::fill(std::begin(tmp),std::end(tmp),ValueDirectional{});
        for(int i=0;i<=deg;++i){tmp[i]=tmp[i]+poly[i]*c0;tmp[i+1]=tmp[i+1]+poly[i]*c1;}++deg;for(int i=0;i<=deg;++i)poly[i]=tmp[i];};
    for(int i=0;i<a;++i)mul(e.xiA[0],e.t2[0]);for(int i=0;i<b;++i)mul(e.xiA[1],e.t2[1]);
    ValueDirectional s;for(int i=0;i<=deg;++i)s=s+poly[i]*Jl[i];return s;
}

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

static void triangle_inplane_A_moments_directional(
    const ValueDirectional P[3][3],const ValueDirectional r[3],int degree,
    ValueDirectional A[POLY_MAX_DEG+1][POLY_MAX_DEG+1],TriPolySetupDirectional& g)
{
    for(int i=0;i<=POLY_MAX_DEG;++i)for(int j=0;j<=POLY_MAX_DEG;++j)A[i][j]={};
    degree=std::min(degree,POLY_MAX_DEG);if(!tri_poly_setup_directional(P,r,g))return;
    ValueDirectional Jl[3][POLY_MAX_DEG+3]{};for(int i=0;i<3;++i)edge_l_moments_poly_directional(g.edges[i],degree+2,Jl[i]);
    A[0][0]=TriPotentialDirectionalValue(P,r);
    if(degree>=1){ValueDirectional A1[2]{};for(int i=0;i<3;++i){const auto er=edge_R_for_poly_directional(g.edges[i]);A1[0]=A1[0]+g.edges[i].m2[0]*er;A1[1]=A1[1]+g.edges[i].m2[1]*er;}A[1][0]=A1[0];A[0][1]=A1[1];}
    auto Eedge=[&](int j,int p,int q){ValueDirectional s;for(int i=0;i<3;++i)s=s+g.edges[i].m2[j]*edge_inplane_monomial_poly_directional(g.edges[i],p,q,Jl[i]);return s;};
    auto Eneg1=[&](int a,int b){return Eedge(0,a+1,b)+Eedge(1,a,b+1);};
    for(int k=2;k<=degree;++k)for(int a=k;a>=0;--a){const int b=k-a;ValueDirectional h2B;
        if(a>=1)h2B=g.h*g.h*((a-1.0)*(a>=2?A[a-2][b]:ValueDirectional{})-Eedge(0,a-1,b));
        else h2B=g.h*g.h*((b-1.0)*(b>=2?A[a][b-2]:ValueDirectional{})-Eedge(1,a,b-1));
        A[a][b]=(Eneg1(a,b)-h2B)/(k+1.0);
    }
}

static void poly2_mul_linear_directional(
    ValueDirectional poly[POLY_MAX_DEG+1][POLY_MAX_DEG+1],int& deg,
    ValueDirectional c0,ValueDirectional c1,ValueDirectional c2)
{
    ValueDirectional tmp[POLY_MAX_DEG+1][POLY_MAX_DEG+1]{};
    for(int a=0;a<=deg;++a)for(int b=0;b<=deg-a;++b){const auto v=poly[a][b];tmp[a][b]=tmp[a][b]+v*c0;tmp[a+1][b]=tmp[a+1][b]+v*c1;tmp[a][b+1]=tmp[a][b+1]+v*c2;}
    ++deg;for(int a=0;a<=deg;++a)for(int b=0;b<=deg-a;++b)poly[a][b]=tmp[a][b];
}

static void SurfacePotentialMomentsUpToDirectional(
    const double P[3][3],const double dP[3][3],const double r[3],const double dr[3],int degree,double* value,double* direction)
{
    ValueDirectional Pd[3][3],rd[3];for(int i=0;i<3;++i)for(int k=0;k<3;++k)Pd[i][k]={P[i][k],dP[i][k]};for(int k=0;k<3;++k)rd[k]={r[k],dr[k]};
    ValueDirectional A[POLY_MAX_DEG+1][POLY_MAX_DEG+1];TriPolySetupDirectional g;
    triangle_inplane_A_moments_directional(Pd,rd,degree,A,g);int idx=0;
    for(int total=0;total<=degree;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){
        const int alpha[3]={ax,ay,total-ax-ay};ValueDirectional poly[POLY_MAX_DEG+1][POLY_MAX_DEG+1]{};int pd=0;poly[0][0]={1.0,0.0};
        for(int coord=0;coord<3;++coord)for(int q=0;q<alpha[coord];++q)poly2_mul_linear_directional(poly,pd,g.rp[coord],g.e1[coord],g.e2[coord]);
        ValueDirectional s;for(int a=0;a<=pd;++a)for(int b=0;b<=pd-a;++b)s=s+poly[a][b]*A[a][b];value[idx]=s.value;direction[idx]=s.direction;++idx;
    }
}

static int PotentialMomentIndex(int ax,int ay,int az);

static void TetPotentialMomentsUpToDirectional(
    const double V[4][3],const double dV[4][3],const double r[3],const double dr[3],
    int degree,double* value,double* direction)
{
    static const int FACES[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    ValueDirectional Vd[4][3],rd[3],cen[3];
    for(int i=0;i<4;++i)for(int k=0;k<3;++k){Vd[i][k]={V[i][k],dV[i][k]};cen[k]=cen[k]+Vd[i][k]/4.0;}
    for(int k=0;k<3;++k)rd[k]={r[k],dr[k]};
    ValueDirectional fm[4][POLY_MAX_MOMENTS]{},h[4]{},phi;
    for(int fi=0;fi<4;++fi){double F[3][3],dF[3][3],fv[POLY_MAX_MOMENTS]{},fd[POLY_MAX_MOMENTS]{};ValueDirectional Fd[3][3];
        for(int j=0;j<3;++j)for(int k=0;k<3;++k){const int vi=FACES[fi][j];F[j][k]=V[vi][k];dF[j][k]=dV[vi][k];Fd[j][k]=Vd[vi][k];}
        SurfacePotentialMomentsUpToDirectional(F,dF,r,dr,degree,fv,fd);const int count=(degree+1)*(degree+2)*(degree+3)/6;for(int q=0;q<count;++q)fm[fi][q]={fv[q],fd[q]};
        ValueDirectional e1[3],e2[3],n[3],fc[3],outward[3],rmf[3];for(int k=0;k<3;++k){e1[k]=Fd[1][k]-Fd[0][k];e2[k]=Fd[2][k]-Fd[0][k];fc[k]=(Fd[0][k]+Fd[1][k]+Fd[2][k])/3.0;outward[k]=fc[k]-cen[k];rmf[k]=rd[k]-Fd[0][k];}
        vd_cross(e1,e2,n);const auto nl=vd_norm(n);if(nl.value<1e-300)continue;for(auto& x:n)x=x/nl;if(vd_dot(outward,n).value<0)for(auto& x:n)x=-x;
        h[fi]=vd_dot(rmf,n);phi=phi+0.5*(-h[fi])*fm[fi][0];
    }
    ValueDirectional out[POLY_MAX_MOMENTS]{};out[0]=phi;
    for(int total=1;total<=degree;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){const int az=total-ax-ay,idx=PotentialMomentIndex(ax,ay,az);ValueDirectional s;
        for(int fi=0;fi<4;++fi)s=s-h[fi]*fm[fi][idx];const int alpha[3]={ax,ay,az};for(int k=0;k<3;++k)if(alpha[k]>0){int lo[3]={ax,ay,az};--lo[k];s=s+rd[k]*double(alpha[k])*out[PotentialMomentIndex(lo[0],lo[1],lo[2])];}out[idx]=s/(total+2.0);}
    const int count=(degree+1)*(degree+2)*(degree+3)/6;for(int i=0;i<count;++i){value[i]=out[i].value;direction[i]=out[i].direction;}
}

// Zero only the a+b <= degree triangle of a dense (POLY_MAX_DEG+1)^2 scratch
// block.  The whole block is 2.9 kB, but these scratches sit in the innermost
// loop of the H-matrix charge-Gram fill and the production charge orders touch
// only the first few rows (degree 2 uses 6 of 361 entries), so zeroing the
// full block costs far more than the arithmetic it protects.
static inline void poly2_zero_triangle(
    double poly[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1], int degree)
{
    for (int a = 0; a <= degree; ++a)
        std::fill(poly[a], poly[a] + (degree - a + 1), 0.0);
}

static void poly2_mul_linear(double poly[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1], int& deg,
                             double c0, double c1, double c2)
{
    double tmp[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1];
    // The accumulation below touches a+b <= deg+1 and uses +=, so exactly that
    // triangle must start at zero.
    poly2_zero_triangle(tmp, deg + 1);
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
    double poly[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1];
    poly2_zero_triangle(poly, alpha[0] + alpha[1] + alpha[2]);
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
                double poly[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1];
                poly2_zero_triangle(poly, total);
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
    // A degenerate face keeps h[fi] = 0 but its row is still read below, so the
    // used prefix must start zeroed -- only the used prefix, not all 1330
    // entries (degree 2 needs 10).
    const int n_moments = (degree + 1)*(degree + 2)*(degree + 3)/6;
    double face_moments[4][POLY_MAX_MOMENTS];
    for (int fi = 0; fi < 4; ++fi)
        std::fill(face_moments[fi], face_moments[fi] + n_moments, 0.0);
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
                double poly[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1];
                poly2_zero_triangle(poly, total);
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

    // Same degenerate-face contract as TetPotentialMomentsUpTo: zero the used
    // prefix of every row, not the full 1330-entry block.
    const int n_moments = (degree + 1)*(degree + 2)*(degree + 3)/6;
    double face_moments[4][POLY_MAX_MOMENTS];
    for (int fi = 0; fi < 4; ++fi)
        std::fill(face_moments[fi], face_moments[fi] + n_moments, 0.0);
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

// `val` is read only where `seen` is true, and `seen` is set only after the
// matching `val` is written, so `val` needs no initialization at all.  Zeroing
// the full pair cost ~62 kB of memset per TetPotentialPolynomial call, and only
// the [0..degree]^3 corner of `seen` is ever addressed.
struct TetMomentMemo {
    bool seen[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1][POLY_MAX_DEG + 1];
    double val[POLY_MAX_DEG + 1][POLY_MAX_DEG + 1][POLY_MAX_DEG + 1];

    explicit TetMomentMemo(int degree)
    {
        const int span = std::min(std::max(degree, 0), POLY_MAX_DEG) + 1;
        for (int i = 0; i < span; ++i)
            for (int j = 0; j < span; ++j)
                std::fill(seen[i][j], seen[i][j] + span, false);
    }
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
    TetMomentMemo memo(degree);
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) {
        int a[3] = {exps[i][0], exps[i][1], exps[i][2]};
        if (a[0] < 0 || a[1] < 0 || a[2] < 0 || a[0] + a[1] + a[2] > POLY_MAX_DEG)
            throw std::runtime_error("TetPotentialPolynomial: invalid exponent");
        s += coeffs[i] * TetPotentialMomentRec(V, r, a, degree, memo);
    }
    return s;
}

void TetReferencePotentialMoments(const double V[4][3], const double r[3],
                                  const std::vector<std::array<int,3>>& exps,
                                  double* out)
{
    if (!out) throw std::invalid_argument("TetReferencePotentialMoments: out is null");
    int degree = 0;
    for (const auto& e : exps) {
        if (e[0] < 0 || e[1] < 0 || e[2] < 0)
            throw std::invalid_argument("TetReferencePotentialMoments: negative exponent");
        degree = std::max(degree, e[0] + e[1] + e[2]);
    }
    if (degree > POLY_MAX_DEG)
        throw std::invalid_argument("TetReferencePotentialMoments: degree exceeds 18");
    // The callee writes (or zeroes) exactly this prefix on both of its paths.
    const int n_moments = (degree + 1)*(degree + 2)*(degree + 3)/6;
    double moments[POLY_MAX_MOMENTS];
    std::fill(moments, moments + n_moments, 0.0);
    TetReferencePotentialMomentsUpTo(V, r, degree, moments);
    for (size_t i = 0; i < exps.size(); ++i) {
        const auto& e = exps[i];
        out[i] = moments[PotentialMomentIndex(e[0], e[1], e[2])];
    }
}

void TetPotentialMomentsUpTo3(const double V[4][3], const double r[3], double out[20])
{
    TetPotentialMomentsUpTo(V, r, 3, out);
}

void TetPotentialMomentsUpTo6(const double V[4][3], const double r[3], double out[84])
{
    TetPotentialMomentsUpTo(V, r, 6, out);
}

void TetReferencePotentialMomentsDirectional(
    const double V[4][3], const double dV[4][3],
    const double r[3], const double dr[3],
    const std::vector<std::array<int,3>>& exps,
    double* value, double* direction)
{
    if (!value || !direction)
        throw std::invalid_argument("TetReferencePotentialMomentsDirectional: null output");
    int degree = 0;
    for (const auto& e : exps) {
        if (e[0] < 0 || e[1] < 0 || e[2] < 0)
            throw std::invalid_argument("TetReferencePotentialMomentsDirectional: negative exponent");
        degree = std::max(degree, e[0] + e[1] + e[2]);
    }
    if (degree > POLY_MAX_DEG)
        throw std::invalid_argument("TetReferencePotentialMomentsDirectional: degree exceeds 18");

    ValueDirectional J[3][3];
    for (int row=0; row<3; ++row) for (int col=0; col<3; ++col)
        J[row][col] = {V[col+1][row]-V[0][row],
                       dV[col+1][row]-dV[0][row]};
    const ValueDirectional det =
        J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
      - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
      + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
    if (std::fabs(det.value) < 1e-300) {
        std::fill(value, value+exps.size(), 0.0);
        std::fill(direction, direction+exps.size(), 0.0);
        return;
    }
    ValueDirectional inv[3][3];
    inv[0][0]=(J[1][1]*J[2][2]-J[1][2]*J[2][1])/det;
    inv[0][1]=-(J[0][1]*J[2][2]-J[0][2]*J[2][1])/det;
    inv[0][2]=(J[0][1]*J[1][2]-J[0][2]*J[1][1])/det;
    inv[1][0]=-(J[1][0]*J[2][2]-J[1][2]*J[2][0])/det;
    inv[1][1]=(J[0][0]*J[2][2]-J[0][2]*J[2][0])/det;
    inv[1][2]=-(J[0][0]*J[1][2]-J[0][2]*J[1][0])/det;
    inv[2][0]=(J[1][0]*J[2][1]-J[1][1]*J[2][0])/det;
    inv[2][1]=-(J[0][0]*J[2][1]-J[0][1]*J[2][0])/det;
    inv[2][2]=(J[0][0]*J[1][1]-J[0][1]*J[1][0])/det;

    ValueDirectional rd[3], v0[3], xiR[3];
    for (int k=0;k<3;++k) { rd[k]={r[k],dr[k]}; v0[k]={V[0][k],dV[0][k]}; }
    for (int i=0;i<3;++i) for (int k=0;k<3;++k) xiR[i]=xiR[i]+inv[i][k]*(rd[k]-v0[k]);

    static const int FACES[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    ValueDirectional face[4][POLY_MAX_MOMENTS]{}, h[4]{}, cen[3]{}, phi;
    for(int i=0;i<4;++i)for(int k=0;k<3;++k)cen[k]=cen[k]+ValueDirectional(V[i][k],dV[i][k])/4.0;
    for(int fi=0;fi<4;++fi){
        ValueDirectional P[3][3],e1[3],e2[3],n[3],fc[3],rmf[3];
        for(int a=0;a<3;++a)for(int k=0;k<3;++k){int vi=FACES[fi][a];P[a][k]={V[vi][k],dV[vi][k]};fc[k]=fc[k]+P[a][k]/3.0;}
        for(int k=0;k<3;++k){e1[k]=P[1][k]-P[0][k];e2[k]=P[2][k]-P[0][k];rmf[k]=rd[k]-P[0][k];}
        vd_cross(e1,e2,n);auto nl=vd_norm(n);if(nl.value<1e-300)continue;for(auto&x:n)x=x/nl;
        ValueDirectional outward[3];for(int k=0;k<3;++k)outward[k]=fc[k]-cen[k];if(vd_dot(outward,n).value<0)for(auto&x:n)x=-x;
        h[fi]=vd_dot(rmf,n);
        ValueDirectional A[POLY_MAX_DEG+1][POLY_MAX_DEG+1];TriPolySetupDirectional g;
        triangle_inplane_A_moments_directional(P,rd,degree,A,g);
        ValueDirectional rpmv0[3];for(int k=0;k<3;++k)rpmv0[k]=g.rp[k]-v0[k];
        int idx=0;
        for(int total=0;total<=degree;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){
            const int alpha[3]={ax,ay,total-ax-ay};ValueDirectional poly[POLY_MAX_DEG+1][POLY_MAX_DEG+1]{};int pd=0;poly[0][0]={1,0};
            for(int coord=0;coord<3;++coord){ValueDirectional c0,c1,c2;for(int k=0;k<3;++k){c0=c0+inv[coord][k]*rpmv0[k];c1=c1+inv[coord][k]*g.e1[k];c2=c2+inv[coord][k]*g.e2[k];}for(int q=0;q<alpha[coord];++q)poly2_mul_linear_directional(poly,pd,c0,c1,c2);}
            ValueDirectional s;for(int a=0;a<=pd;++a)for(int b=0;b<=pd-a;++b)s=s+poly[a][b]*A[a][b];face[fi][idx++]=s;
        }
        phi=phi+0.5*(-h[fi])*face[fi][0];
    }
    ValueDirectional out[POLY_MAX_MOMENTS]{};out[0]=phi;
    for(int total=1;total<=degree;++total)for(int ax=0;ax<=total;++ax)for(int ay=0;ay<=total-ax;++ay){int az=total-ax-ay,idx=PotentialMomentIndex(ax,ay,az);ValueDirectional s;for(int fi=0;fi<4;++fi)s=s-h[fi]*face[fi][idx];const int alpha[3]={ax,ay,az};for(int k=0;k<3;++k)if(alpha[k]){int lo[3]={ax,ay,az};--lo[k];s=s+xiR[k]*double(alpha[k])*out[PotentialMomentIndex(lo[0],lo[1],lo[2])];}out[idx]=s/(total+2.0);}
    for(size_t i=0;i<exps.size();++i){const auto&e=exps[i];const auto&x=out[PotentialMomentIndex(e[0],e[1],e[2])];value[i]=x.value;direction[i]=x.direction;}
}
void TetPotentialMomentsDirectionalUpTo3(const double V[4][3],const double dV[4][3],
    const double r[3],const double dr[3],double value[20],double direction[20])
{ TetPotentialMomentsUpToDirectional(V,dV,r,dr,3,value,direction); }
void TetPotentialMomentsDirectionalUpTo6(const double V[4][3],const double dV[4][3],
    const double r[3],const double dr[3],double value[84],double direction[84])
{ TetPotentialMomentsUpToDirectional(V,dV,r,dr,6,value,direction); }

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
        // Only the degree-dependent prefix is written and read; the full
        // 1330-entry block is re-zeroed once per (quadrature point, source).
        const int n_moments = (degree + 1)*(degree + 2)*(degree + 3)/6;
        double moments[POLY_MAX_MOMENTS];

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
                std::fill(moments, moments + n_moments, 0.0);
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

void TriPotentialMomentsDirectionalUpTo4(
    const double V[3][3],const double dV[3][3],const double r[3],const double dr[3],
    double value[35],double direction[35])
{
    SurfacePotentialMomentsUpToDirectional(V,dV,r,dr,4,value,direction);
}

void TriPotentialMomentsDirectionalUpTo2(
    const double V[3][3],const double dV[3][3],const double r[3],const double dr[3],
    double value[10],double direction[10])
{
    SurfacePotentialMomentsUpToDirectional(V,dV,r,dr,2,value,direction);
}

void TetPotentialMomentsDirectionalUpTo1(
    const double V[4][3],const double dV[4][3],const double r[3],const double dr[3],
    double value[4],double direction[4])
{
    // Keep degree one on the same positive-Newtonian-potential recurrence used
    // by all higher orders.  The retired special case returned -PhiTet and
    // assumed an x,y,z storage order where PotentialMomentIndex is z,y,x;
    // callers then needed a second sign reversal.  One physical convention is
    // now shared by the Gram derivative and field derivative paths.
    TetPotentialMomentsUpToDirectional(V,dV,r,dr,1,value,direction);
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

// General total-degree <= 3 physical polynomial density.  coefficient uses
// PotentialMomentIndex(ax,ay,az), i.e. increasing total degree and then ax,ay.
// Integration by parts gives
//   int_V rho grad'(1/R) = int_dV rho n/R - int_V grad(rho)/R.
// The existing analytic triangle/tetrahedron potential moments therefore close
// the cubic field without target or source quadrature.
void TetVolFieldCubic(const double V[4][3], const double r[3],
                      const double coefficient[20], double out[3])
{
    out[0]=out[1]=out[2]=0.0;
    double cen[3]={0,0,0};
    for(int i=0;i<4;++i)for(int k=0;k<3;++k)cen[k]+=0.25*V[i][k];
    static const int FACES[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    for(int fi=0;fi<4;++fi){
        double face[3][3];
        for(int j=0;j<3;++j)for(int k=0;k<3;++k)
            face[j][k]=V[FACES[fi][j]][k];
        double e1[3],e2[3],normal[3];
        for(int k=0;k<3;++k){e1[k]=face[1][k]-face[0][k];e2[k]=face[2][k]-face[0][k];}
        v3cross(e1,e2,normal);const double norm=v3nrm(normal);
        if(norm<1e-300)continue;
        for(double& value:normal)value/=norm;
        double fc[3]={0,0,0};
        for(int j=0;j<3;++j)for(int k=0;k<3;++k)fc[k]+=face[j][k]/3.0;
        double outward[3];for(int k=0;k<3;++k)outward[k]=fc[k]-cen[k];
        if(v3dot(outward,normal)<0.0)for(double& value:normal)value=-value;
        double moments[35]={};TriPotentialMomentsUpTo4(face,r,moments);
        double weighted=0.0;
        for(int index=0;index<20;++index)weighted+=coefficient[index]*moments[index];
        for(int k=0;k<3;++k)out[k]+=normal[k]*weighted;
    }
    double volume_moments[20]={};TetPotentialMomentsUpTo3(V,r,volume_moments);
    for(int total=1;total<=3;++total)
        for(int ax=0;ax<=total;++ax)
            for(int ay=0;ay<=total-ax;++ay){
                const int az=total-ax-ay;
                const int alpha[3]={ax,ay,az};
                const double value=coefficient[PotentialMomentIndex(ax,ay,az)];
                if(value==0.0)continue;
                for(int k=0;k<3;++k)if(alpha[k]>0){
                    int lower[3]={ax,ay,az};--lower[k];
                    out[k]-=value*alpha[k]*volume_moments[
                        PotentialMomentIndex(lower[0],lower[1],lower[2])];
                }
            }
}

// All 20 physical-monomial fields in one geometry pass.  This is the batched
// form used by affine-HEX observation rows: the expensive triangle/tetrahedron
// moments depend on (V,r), not on the individual charge coefficient.  The
// ordering is PotentialMomentIndex(ax,ay,az), identical to coefficient[].
void TetVolFieldCubicBasis(const double V[4][3], const double r[3],
                           double out[20][3])
{
    for(int index=0;index<20;++index)
        for(int k=0;k<3;++k)out[index][k]=0.0;
    double cen[3]={0,0,0};
    for(int i=0;i<4;++i)for(int k=0;k<3;++k)cen[k]+=0.25*V[i][k];
    static const int FACES[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    for(int fi=0;fi<4;++fi){
        double face[3][3];
        for(int j=0;j<3;++j)for(int k=0;k<3;++k)
            face[j][k]=V[FACES[fi][j]][k];
        double e1[3],e2[3],normal[3];
        for(int k=0;k<3;++k){e1[k]=face[1][k]-face[0][k];e2[k]=face[2][k]-face[0][k];}
        v3cross(e1,e2,normal);const double norm=v3nrm(normal);
        if(norm<1e-300)continue;
        for(double& value:normal)value/=norm;
        double fc[3]={0,0,0};
        for(int j=0;j<3;++j)for(int k=0;k<3;++k)fc[k]+=face[j][k]/3.0;
        double outward[3];for(int k=0;k<3;++k)outward[k]=fc[k]-cen[k];
        if(v3dot(outward,normal)<0.0)for(double& value:normal)value=-value;
        double moments[35]={};TriPotentialMomentsUpTo4(face,r,moments);
        for(int index=0;index<20;++index)
            for(int k=0;k<3;++k)out[index][k]+=normal[k]*moments[index];
    }
    double volume_moments[20]={};TetPotentialMomentsUpTo3(V,r,volume_moments);
    for(int total=1;total<=3;++total)
        for(int ax=0;ax<=total;++ax)
            for(int ay=0;ay<=total-ax;++ay){
                const int az=total-ax-ay;
                const int alpha[3]={ax,ay,az};
                const int index=PotentialMomentIndex(ax,ay,az);
                for(int k=0;k<3;++k)if(alpha[k]>0){
                    int lower[3]={ax,ay,az};--lower[k];
                    out[index][k]-=alpha[k]*volume_moments[
                        PotentialMomentIndex(lower[0],lower[1],lower[2])];
                }
            }
}

// INT_T (sigma0 + s.r')(r-r')/R^3 dS' (linear surface charge)
// = (sigma0 + s.r_p) F_const - SUM_e (s.m_e) G_e - I0 s_par.
void LinTriField(const double V[3][3], const double r[3], double sigma0, const double s[3], double out[3])
{
    out[0]=out[1]=out[2]=0.0;
    double e1[3],e2[3],n[3];
    for (int k=0;k<3;k++){ e1[k]=V[1][k]-V[0][k]; e2[k]=V[2][k]-V[0][k]; }
    v3cross(e1,e2,n); double nl=v3nrm(n);
    // Degenerate (zero-area) face: TriPotential/TriField already return the
    // zero limit here, so match them instead of dividing by zero and pushing
    // NaN into every moment that consumes this face.
    if (nl<1e-300) return;
    for (int k=0;k<3;k++) n[k]/=nl;
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

// Exact coplanar exterior limit of INT_T 1/|r-r'|^3 dS'.  The quadratic
// surface-charge field normally obtains this moment as (n.TriField)/h.  That
// quotient is 0/0 when an observation lies in the panel plane but outside the
// triangle, even though the physical field is smooth there.  In two
// dimensions, div_{r'}((r'-p)/|r'-p|^3)=-1/|r'-p|^3, so the finite limit is a
// sum of elementary edge integrals.  Points on/in the triangle deliberately
// stay on the ordinary one-sided path because their surface field is singular
// or discontinuous.
static bool coplanar_outside_inverse_cube_moment(
    const double V[3][3], const double r[3], double& moment)
{
    double edge01[3],edge02[3],normal[3];
    for(int k=0;k<3;++k){edge01[k]=V[1][k]-V[0][k];edge02[k]=V[2][k]-V[0][k];}
    v3cross(edge01,edge02,normal);
    const double normal_length=v3nrm(normal);
    if(normal_length<1e-300)return false;
    for(double& value:normal)value/=normal_length;
    double offset[3];for(int k=0;k<3;++k)offset[k]=r[k]-V[0][k];
    const double height=v3dot(offset,normal);
    double scale=0.0;
    for(int vertex=0;vertex<3;++vertex)
        for(int k=0;k<3;++k)
            scale=std::max(scale,std::fabs(V[vertex][k]));
    for(int k=0;k<3;++k)scale=std::max(scale,std::fabs(r[k]));
    for(int edge=0;edge<3;++edge){
        double side[3];for(int k=0;k<3;++k)side[k]=V[(edge+1)%3][k]-V[edge][k];
        scale=std::max(scale,v3nrm(side));
    }
    if(scale<1e-300)return false;
    const double tolerance=32.0*std::numeric_limits<double>::epsilon()*scale;
    if(std::fabs(height)>tolerance)return false;
    double projection[3];for(int k=0;k<3;++k)projection[k]=r[k]-height*normal[k];
    double closest[3];ClosestPointTriangle(projection,V[0],V[1],V[2],closest);
    double separation[3];for(int k=0;k<3;++k)separation[k]=projection[k]-closest[k];
    if(v3nrm(separation)<=tolerance)return false;

    double centroid[3]={0.0,0.0,0.0};
    for(int vertex=0;vertex<3;++vertex)for(int k=0;k<3;++k)
        centroid[k]+=V[vertex][k]/3.0;
    double total=0.0;
    for(int edge=0;edge<3;++edge){
        const double* A=V[edge];const double* B=V[(edge+1)%3];
        double tangent[3];for(int k=0;k<3;++k)tangent[k]=B[k]-A[k];
        const double length=v3nrm(tangent);if(length<1e-300)return false;
        for(double& value:tangent)value/=length;
        double outward[3];edge_outnormal(A,B,normal,centroid,outward);
        double q[3];for(int k=0;k<3;++k)q[k]=A[k]-projection[k];
        const double perpendicular=v3dot(q,outward);
        const double u1=v3dot(q,tangent),u2=u1+length;
        if(std::fabs(perpendicular)<=tolerance){
            // An exterior point can lie on an edge-line extension.  Its edge
            // contribution has the finite zero limit when both endpoints are
            // on the same side; crossing the segment would be singular.
            if(u1*u2<=0.0)return false;
            continue;
        }
        const double perpendicular2=perpendicular*perpendicular;
        const double R1=std::hypot(perpendicular,u1);
        const double R2=std::hypot(perpendicular,u2);
        const double edge_integral=(u2/(perpendicular2*R2)
                                   -u1/(perpendicular2*R1));
        total-=perpendicular*edge_integral;
    }
    if(!std::isfinite(total)||total<=0.0)return false;
    moment=total;
    return true;
}

// INT_T (sigma0 + s.r' + r'^T S r')(r-r')/R^3 dS' (S symmetric) via the in-plane/normal split.
void QuadTriField(const double V[3][3], const double r[3], double sigma0,
                  const double s[3], const double S[3][3], double out[3])
{
    double basis[10][3];
    QuadTriFieldBasis(V,r,basis);
    double coefficient[10]={};
    coefficient[0]=sigma0;
    for(int axis=0;axis<3;++axis)
        coefficient[PotentialMomentIndex(axis==0,axis==1,axis==2)]=s[axis];
    for(int a=0;a<3;++a)for(int b=a;b<3;++b){
        int exponent[3]={0,0,0};++exponent[a];++exponent[b];
        coefficient[PotentialMomentIndex(exponent[0],exponent[1],exponent[2])]=
            a==b?S[a][a]:S[a][b]+S[b][a];
    }
    for(int k=0;k<3;++k){
        out[k]=0.0;
        for(int index=0;index<10;++index)out[k]+=coefficient[index]*basis[index][k];
    }
}

// Exact field of a total-degree <= 3 physical surface polynomial. If
// U(r)=int_T sigma(r')/|r-r'| dS', then the magnetic-charge field is
// -grad_r U. The potential-moment directional recurrence already carries
// the analytic target derivative, so three unit target directions close the
// cubic field without source or observation quadrature.
void CubicTriField(const double V[3][3], const double r[3],
                   const double coefficient[20], double out[3])
{
    const double zero_vertices[3][3] = {};
    for(int axis=0;axis<3;++axis){
        double direction_vector[3] = {};
        direction_vector[axis] = 1.0;
        double moments[35] = {}, direction[35] = {};
        TriPotentialMomentsDirectionalUpTo4(
            V,zero_vertices,r,direction_vector,moments,direction);
        out[axis] = 0.0;
        for(int index=0;index<20;++index)
            out[axis] -= coefficient[index]*direction[index];
    }
}

void QuadTriFieldBasis(const double V[3][3], const double r[3],
                       double out[10][3])
{
    for(int index=0;index<10;++index)
        for(int k=0;k<3;++k)out[index][k]=0.0;
    double e1u[3],e2u[3],n[3];
    for (int k=0;k<3;k++){ e1u[k]=V[1][k]-V[0][k]; e2u[k]=V[2][k]-V[0][k]; }
    v3cross(e1u,e2u,n); double nl=v3nrm(n);
    double e1l=v3nrm(e1u);
    // Degenerate (zero-area) face: TriPotential/TriField already return the
    // zero limit here, so match them instead of dividing by zero and pushing
    // NaN into every moment that consumes this face.
    if (nl<1e-300 || e1l<1e-300) return;
    for (int k=0;k<3;k++) n[k]/=nl;
    double e1[3]; for (int k=0;k<3;k++) e1[k]=e1u[k]/e1l;
    double e2[3]; v3cross(n,e1,e2);
    double h=0; { double d[3]; for(int k=0;k<3;k++) d[k]=r[k]-V[0][k]; h=v3dot(d,n); }
    double r_p[3]; for (int k=0;k<3;k++) r_p[k]=r[k]-h*n[k];
    double cen[3]={0,0,0}; for (int j=0;j<3;j++) for(int k=0;k<3;k++) cen[k]+=V[j][k]/3.0;
    double Fc[3]; TriField(V,r,Fc);
    double I0=TriPotential(V,r);
    double M1[3]; TriMoment1(V,r,M1);
    double J3_0=0.0;
    const bool coplanar_outside=
        coplanar_outside_inverse_cube_moment(V,r,J3_0);
    if(!coplanar_outside)J3_0=v3dot(n,Fc)/h;                    // INT_T 1/R^3
    double intxi1[3]={0,0,0};                                  // INT_T xi/R^3
    double xixi[3][3];                                         // INT_T xi(x)xi/R^3 = P I0 - SUM (Gxi)(x)m
    for (int a=0;a<3;a++) for (int b=0;b<3;b++) xixi[a][b]=((a==b?1.0:0.0)-n[a]*n[b])*I0;
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
        // In-plane edge terms for every physical monomial.  Along
        // r'(l)=A+l*tt they are polynomials through degree two.
        double tt[3]; for (int k=0;k<3;k++) tt[k]=B[k]-A[k]; double Lt=v3nrm(tt); for(int k=0;k<3;k++) tt[k]/=Lt;
        double edge_basis[10]={};
        edge_basis[0]=Jl[0];
        for(int axis=0;axis<3;++axis)
            edge_basis[PotentialMomentIndex(axis==0,axis==1,axis==2)]=
                A[axis]*Jl[0]+tt[axis]*Jl[1];
        for(int a=0;a<3;++a)for(int b=a;b<3;++b){
            int exponent[3]={0,0,0};++exponent[a];++exponent[b];
            const int index=PotentialMomentIndex(exponent[0],exponent[1],exponent[2]);
            edge_basis[index]=A[a]*A[b]*Jl[0]
                +(A[a]*tt[b]+A[b]*tt[a])*Jl[1]+tt[a]*tt[b]*Jl[2];
        }
        for(int index=0;index<10;++index)
            for(int k=0;k<3;++k)out[index][k]+=m[k]*edge_basis[index];
    }
    double J3_1[3]; for (int k=0;k<3;k++) J3_1[k]=intxi1[k]+r_p[k]*J3_0;          // INT_T r'/R^3
    double J3_2[3][3];                                                            // INT_T r'(x)r'/R^3
    for (int a=0;a<3;a++) for (int b=0;b<3;b++)
        J3_2[a][b]=xixi[a][b] + r_p[a]*intxi1[b] + intxi1[a]*r_p[b] + r_p[a]*r_p[b]*J3_0;
    // Interior derivative and normal terms, again one column per physical
    // monomial coefficient.
    for(int axis=0;axis<3;++axis){
        const int index=PotentialMomentIndex(axis==0,axis==1,axis==2);
        for(int k=0;k<3;++k){
            const double projected=(k==axis?1.0:0.0)-n[k]*n[axis];
            out[index][k]-=projected*I0;
            out[index][k]+=n[k]*h*J3_1[axis];
        }
    }
    for(int a=0;a<3;++a)for(int b=a;b<3;++b){
        int exponent[3]={0,0,0};++exponent[a];++exponent[b];
        const int index=PotentialMomentIndex(exponent[0],exponent[1],exponent[2]);
        for(int k=0;k<3;++k){
            double derivative=0.0;
            if(a==b)
                derivative=2.0*((k==a?1.0:0.0)-n[k]*n[a])*M1[a];
            else
                derivative=((k==a?1.0:0.0)-n[k]*n[a])*M1[b]
                          +((k==b?1.0:0.0)-n[k]*n[b])*M1[a];
            out[index][k]-=derivative;
            const double j2=a==b?J3_2[a][a]
                :0.5*(J3_2[a][b]+J3_2[b][a]);
            out[index][k]+=n[k]*h*j2;
        }
    }
    for(int k=0;k<3;++k)out[0][k]+=n[k]*h*J3_0;
}

// Forward-mode copy of the exact Wilton/Graglia triangle field.  Keeping the
// same branch decisions as TriField makes the derivative the tangent of the
// production closed form, including near-plane targets.
static void TriFieldDirectionalValue(
    const ValueDirectional V[3][3], const ValueDirectional r[3],
    ValueDirectional out[3])
{
    ValueDirectional e1[3],e2[3],n[3];
    for(int k=0;k<3;++k){e1[k]=V[1][k]-V[0][k];e2[k]=V[2][k]-V[0][k];}
    vd_cross(e1,e2,n);const auto nl=vd_norm(n);if(nl.value<1e-300)return;
    for(auto& x:n)x=x/nl;
    ValueDirectional rmv0[3];for(int k=0;k<3;++k)rmv0[k]=r[k]-V[0][k];
    const auto d=vd_dot(rmv0,n);
    ValueDirectional p[3];for(int k=0;k<3;++k)p[k]=r[k]-d*n[k];
    const auto ad=vd_abs(d);ValueDirectional omega;
    for(int i=0;i<3;++i){
        const auto* a=V[i];const auto* b=V[(i+1)%3];
        ValueDirectional lh[3];for(int k=0;k<3;++k)lh[k]=b[k]-a[k];
        const auto ll=vd_norm(lh);if(ll.value<1e-300)continue;
        for(auto& x:lh)x=x/ll;
        ValueDirectional uh[3];vd_cross(lh,n,uh);
        ValueDirectional ap[3],bp[3],ra[3],rb[3];
        for(int k=0;k<3;++k){
            ap[k]=a[k]-p[k];bp[k]=b[k]-p[k];
            ra[k]=r[k]-a[k];rb[k]=r[k]-b[k];
        }
        const auto P0=vd_dot(ap,uh),sm=vd_dot(ap,lh),sp=vd_dot(bp,lh);
        const auto Rm=vd_norm(ra),Rp=vd_norm(rb),R0sq=P0*P0+d*d;
        const auto dm=Rm+sm,dp=Rp+sp;
        const auto f=(dp.value>1e-300&&dm.value>1e-300)
            ?vd_log(dp/dm):ValueDirectional{};
        const auto beta=vd_atan2(P0*sp,R0sq+ad*Rp)
                       -vd_atan2(P0*sm,R0sq+ad*Rm);
        for(int k=0;k<3;++k)out[k]=out[k]+f*uh[k];
        omega=omega+beta;
    }
    const double sign=d.value>0.0?1.0:(d.value<0.0?-1.0:0.0);
    for(int k=0;k<3;++k)out[k]=out[k]+sign*omega*n[k];
}

void TetVolFieldLinearDirectional(
    const double V[4][3], const double dV[4][3],
    const double r[3], const double dr[3],
    double rho0, double drho0, const double g[3], const double dg[3],
    double value[3], double direction[3])
{
    static const int faces[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    ValueDirectional vertices[4][3],target[3],centroid[3],gradient[3];
    for(int i=0;i<4;++i)for(int k=0;k<3;++k){
        vertices[i][k]={V[i][k],dV[i][k]};
        centroid[k]=centroid[k]+vertices[i][k]/4.0;
    }
    for(int k=0;k<3;++k){target[k]={r[k],dr[k]};gradient[k]={g[k],dg[k]};}
    const ValueDirectional density0{rho0,drho0};
    ValueDirectional out[3];
    for(int fi=0;fi<4;++fi){
        double face[3][3],dface[3][3],moments[4],dmoments[4];
        ValueDirectional face_dual[3][3],e1[3],e2[3],normal[3],center[3];
        for(int j=0;j<3;++j)for(int k=0;k<3;++k){
            const int vi=faces[fi][j];
            face[j][k]=V[vi][k];dface[j][k]=dV[vi][k];
            face_dual[j][k]=vertices[vi][k];
            center[k]=center[k]+face_dual[j][k]/3.0;
        }
        for(int k=0;k<3;++k){
            e1[k]=face_dual[1][k]-face_dual[0][k];
            e2[k]=face_dual[2][k]-face_dual[0][k];
        }
        vd_cross(e1,e2,normal);const auto length=vd_norm(normal);
        if(length.value<1e-300)continue;
        for(auto& x:normal)x=x/length;
        ValueDirectional outward[3];
        for(int k=0;k<3;++k)outward[k]=center[k]-centroid[k];
        if(vd_dot(outward,normal).value<0.0)for(auto& x:normal)x=-x;
        SurfacePotentialMomentsUpToDirectional(
            face,dface,r,dr,1,moments,dmoments);
        ValueDirectional weighted=density0*ValueDirectional{moments[0],dmoments[0]};
        for(int k=0;k<3;++k){
            const int index=PotentialMomentIndex(k==0?1:0,k==1?1:0,k==2?1:0);
            weighted=weighted+gradient[k]*ValueDirectional{moments[index],dmoments[index]};
        }
        for(int k=0;k<3;++k)out[k]=out[k]+normal[k]*weighted;
    }
    double tet_moments[4],dtet_moments[4];
    TetPotentialMomentsDirectionalUpTo1(V,dV,r,dr,tet_moments,dtet_moments);
    const ValueDirectional potential{tet_moments[0],dtet_moments[0]};
    for(int k=0;k<3;++k){
        out[k]=out[k]-gradient[k]*potential;
        value[k]=out[k].value;direction[k]=out[k].direction;
    }
}

void QuadTriFieldDirectional(
    const double V[3][3], const double dV[3][3],
    const double r[3], const double dr[3],
    double sigma0, double dsigma0,
    const double s[3], const double ds[3],
    const double S[3][3], const double dS[3][3],
    double value[3], double direction[3])
{
    ValueDirectional P[3][3],target[3],slope[3],hessian[3][3];
    for(int i=0;i<3;++i)for(int k=0;k<3;++k)P[i][k]={V[i][k],dV[i][k]};
    for(int k=0;k<3;++k){target[k]={r[k],dr[k]};slope[k]={s[k],ds[k]};}
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)hessian[i][j]={S[i][j],dS[i][j]};
    const ValueDirectional constant{sigma0,dsigma0};
    TriPolySetupDirectional geometry;
    if(!tri_poly_setup_directional(P,target,geometry)){
        for(int k=0;k<3;++k)value[k]=direction[k]=0.0;
        return;
    }
    ValueDirectional field0[3];TriFieldDirectionalValue(P,target,field0);
    double moments[4],dmoments[4];
    SurfacePotentialMomentsUpToDirectional(V,dV,r,dr,1,moments,dmoments);
    const ValueDirectional I0{moments[0],dmoments[0]};
    ValueDirectional M1[3];
    for(int k=0;k<3;++k){
        const int index=PotentialMomentIndex(k==0?1:0,k==1?1:0,k==2?1:0);
        M1[k]={moments[index],dmoments[index]};
    }
    double coplanar_J30=0.0;
    const bool coplanar_outside=
        coplanar_outside_inverse_cube_moment(V,r,coplanar_J30);
    const auto J30=(coplanar_outside
        ?ValueDirectional{coplanar_J30,0.0}
        :vd_dot(geometry.n,field0)/geometry.h);
    ValueDirectional intxi1[3],xixi[3][3],inplane[3];
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        xixi[a][b]=((a==b?1.0:0.0)-geometry.n[a]*geometry.n[b])*I0;
    for(int edge=0;edge<3;++edge){
        const auto* A=P[edge];const auto* B=P[(edge+1)%3];
        const auto& eg=geometry.edges[edge];
        ValueDirectional tangent[3],normal[3];
        for(int k=0;k<3;++k)tangent[k]=(B[k]-A[k])/eg.L;
        for(int k=0;k<3;++k)
            normal[k]=eg.m2[0]*geometry.e1[k]+eg.m2[1]*geometry.e2[k];
        ValueDirectional Jl[POLY_MAX_DEG+3]{};
        edge_l_moments_poly_directional(eg,2,Jl);
        const ValueDirectional Gxi[2]={
            eg.xiA[0]*Jl[0]+eg.t2[0]*Jl[1],
            eg.xiA[1]*Jl[0]+eg.t2[1]*Jl[1]};
        ValueDirectional Gxi3[3];
        for(int k=0;k<3;++k)
            Gxi3[k]=Gxi[0]*geometry.e1[k]+Gxi[1]*geometry.e2[k];
        for(int k=0;k<3;++k)intxi1[k]=intxi1[k]-normal[k]*Jl[0];
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            xixi[a][b]=xixi[a][b]-Gxi3[a]*normal[b];
        ValueDirectional sA,AStA,AStt,ttStt;
        for(int a=0;a<3;++a){
            sA=sA+slope[a]*A[a];
            for(int b=0;b<3;++b){
                AStA=AStA+A[a]*hessian[a][b]*A[b];
                AStt=AStt+A[a]*hessian[a][b]*tangent[b];
                ttStt=ttStt+tangent[a]*hessian[a][b]*tangent[b];
            }
        }
        const auto c0=constant+sA+AStA;
        auto c1=vd_dot(slope,tangent)+2.0*AStt;
        const auto esig=c0*Jl[0]+c1*Jl[1]+ttStt*Jl[2];
        for(int k=0;k<3;++k)inplane[k]=inplane[k]+normal[k]*esig;
    }
    ValueDirectional J31[3],J32[3][3];
    for(int k=0;k<3;++k)J31[k]=intxi1[k]+geometry.rp[k]*J30;
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        J32[a][b]=xixi[a][b]+geometry.rp[a]*intxi1[b]
                   +intxi1[a]*geometry.rp[b]+geometry.rp[a]*geometry.rp[b]*J30;
    const auto sn=vd_dot(slope,geometry.n);
    ValueDirectional projected_slope[3],SM1[3];
    for(int k=0;k<3;++k)projected_slope[k]=slope[k]-sn*geometry.n[k];
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        SM1[a]=SM1[a]+hessian[a][b]*M1[b];
    const auto SM1n=vd_dot(SM1,geometry.n);
    for(int k=0;k<3;++k)
        inplane[k]=inplane[k]-(projected_slope[k]*I0
                              +2.0*(SM1[k]-SM1n*geometry.n[k]));
    ValueDirectional normal_scale=constant*J30+vd_dot(slope,J31);
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        normal_scale=normal_scale+hessian[a][b]*J32[a][b];
    if(coplanar_outside){
        // At h=0 the normal field is zero, while its directional derivative is
        // h' times the finite weighted inverse-cube moment.  The derivative of
        // that moment is multiplied by h and therefore vanishes in this limit.
        normal_scale={geometry.h.value*normal_scale.value,
                      geometry.h.direction*normal_scale.value};
    }else normal_scale=geometry.h*normal_scale;
    for(int k=0;k<3;++k){
        const auto out=inplane[k]+geometry.n[k]*normal_scale;
        value[k]=out.value;direction[k]=out.direction;
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
