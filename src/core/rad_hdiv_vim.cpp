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
// Port of examples reference radia.hdiv_vim._core.tri_potential (scalar form); same edge formula as
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
// NO 1/4pi.  Validated entry-by-entry vs radia.hdiv_vim.flat_triangle_charge_field (machine precision).
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
// TriPotential.  Port of radia.hdiv_vim._core.phi_tet.
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
// (no quadrature).  Validated vs radia.hdiv_vim.tet_self_volume_field (the spherical ray-trace) * 4pi.
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

} // namespace rad_hdiv
