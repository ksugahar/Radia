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

Mesh BuildStructuredRT0(int nx, int ny, int nz, double h)
{
    Mesh m; m.nx = nx; m.ny = ny; m.nz = nz; m.n_cell = nx*ny*nz;
    auto cell_id = [=](int i, int j, int k) { return (i*ny + j)*nz + k; };
    auto node = [=](int i, int j, int k) -> Vec3 { return {i*h, j*h, k*h}; };

    // x-faces: i in [0,nx], j in [0,ny), k in [0,nz); normal +x
    for (int i = 0; i <= nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k) {
                Vec3 nd = node(i, j, k);
                Face f; f.ax = 0; f.area = h*h;
                f.c = {nd[0], nd[1] + 0.5*h, nd[2] + 0.5*h};
                f.lo = (i > 0)  ? cell_id(i-1, j, k) : -1;
                f.hi = (i < nx) ? cell_id(i,   j, k) : -1;
                f.bnd = (i == 0 || i == nx);
                m.faces.push_back(f);
            }
    // y-faces: i in [0,nx), j in [0,ny], k in [0,nz); normal +y
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j <= ny; ++j)
            for (int k = 0; k < nz; ++k) {
                Vec3 nd = node(i, j, k);
                Face f; f.ax = 1; f.area = h*h;
                f.c = {nd[0] + 0.5*h, nd[1], nd[2] + 0.5*h};
                f.lo = (j > 0)  ? cell_id(i, j-1, k) : -1;
                f.hi = (j < ny) ? cell_id(i, j,   k) : -1;
                f.bnd = (j == 0 || j == ny);
                m.faces.push_back(f);
            }
    // z-faces: i in [0,nx), j in [0,ny), k in [0,nz]; normal +z
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k <= nz; ++k) {
                Vec3 nd = node(i, j, k);
                Face f; f.ax = 2; f.area = h*h;
                f.c = {nd[0] + 0.5*h, nd[1] + 0.5*h, nd[2]};
                f.lo = (k > 0)  ? cell_id(i, j, k-1) : -1;
                f.hi = (k < nz) ? cell_id(i, j, k)   : -1;
                f.bnd = (k == 0 || k == nz);
                m.faces.push_back(f);
            }
    m.cell_c.resize(m.n_cell); m.cell_V.assign(m.n_cell, h*h*h);
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
                m.cell_c[cell_id(i, j, k)] = {(i+0.5)*h, (j+0.5)*h, (k+0.5)*h};
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

void AssembleCoulombGram(const Mesh& m, std::vector<double>& G, int& n_charge)
{
    int n_bnd; std::vector<int> brow = bnd_charge_rows(m, n_bnd);
    n_charge = m.n_cell + n_bnd;
    // charge-cell centroids + measures: volume cells then boundary faces (same order as B's rows)
    std::vector<Vec3>   cent(n_charge);
    std::vector<double> meas(n_charge);
    for (int c = 0; c < m.n_cell; ++c) { cent[c] = m.cell_c[c]; meas[c] = m.cell_V[c]; }
    for (int f = 0; f < m.n_face(); ++f)
        if (m.faces[f].bnd) { int r = brow[f]; cent[r] = m.faces[f].c; meas[r] = m.faces[f].area; }

    G.assign((size_t)n_charge * n_charge, 0.0);
    for (int a = 0; a < n_charge; ++a) {
        for (int b = 0; b < n_charge; ++b) {
            if (a == b) continue;
            double dx = cent[a][0]-cent[b][0], dy = cent[a][1]-cent[b][1], dz = cent[a][2]-cent[b][2];
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            G[(size_t)a*n_charge + b] = meas[a]*meas[b] * INV_FOUR_PI / r;   // centroid monopole
        }
        // diagonal self-energy (placeholder): cube for volume cells, square for boundary faces
        double self;
        if (a < m.n_cell) self = C_CUBE * std::pow(meas[a], 5.0/3.0) * INV_FOUR_PI;
        else              self = C_SQ   * std::pow(meas[a], 1.5)     * INV_FOUR_PI;
        G[(size_t)a*n_charge + a] = self;
    }
}

void AssembleN(const Mesh& m, std::vector<double>& N)
{
    std::vector<double> B, G;
    int n_charge, n_bnd, n_charge_g;
    AssembleChargeMap(m, B, n_charge, n_bnd);
    AssembleCoulombGram(m, G, n_charge_g);
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

} // namespace rad_hdiv
