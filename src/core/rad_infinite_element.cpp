/* rad_infinite_element.cpp -- see rad_infinite_element.h.
 *
 * C++ port of the Python prototype acts7_32/7_33 (examples/kelvin_transformation/DtN_spectrum):
 * the orthogonal nodal radial decay basis + the tensor exterior energy (R1 (x) Mtil + R0 (x) Ktil)
 * + static condensation of the radial DOFs -> the DtN surface stiffness S_Gamma.  Pure C++ + MKL.
 */
#include "rad_infinite_element.h"

#include <cmath>
#include <vector>
#include <mkl.h>   // cblas_dgemm + LAPACKE_dgesv (row-major), matching rad_cln.cpp

namespace rad_ie {

namespace {
const double PI = 3.14159265358979323846;

// Legendre polynomials L_0..L_P at x (on [-1,1]); fills out[0..P].
void legendre_all(int P, double x, std::vector<double>& out) {
    out.resize(P + 1);
    out[0] = 1.0;
    if (P >= 1) out[1] = x;
    for (int k = 2; k <= P; ++k)
        out[k] = ((2.0 * k - 1.0) * x * out[k - 1] - (k - 1.0) * out[k - 2]) / k;
}
}  // namespace

// Gauss-Legendre nodes/weights on [0,1] (Newton iteration on the Legendre roots).
void GaussLegendre01(int n, std::vector<double>& x, std::vector<double>& w) {
    x.assign(n, 0.0);
    w.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double xi = std::cos(PI * (i + 0.75) / (n + 0.5));   // initial guess for root i
        double dpn = 1.0;
        for (int it = 0; it < 100; ++it) {
            // evaluate P_n and P_n' at xi by the recurrence
            double p0 = 1.0, p1 = xi;
            for (int k = 2; k <= n; ++k) {
                double p2 = ((2.0 * k - 1.0) * xi * p1 - (k - 1.0) * p0) / k;
                p0 = p1; p1 = p2;
            }
            dpn = n * (xi * p1 - p0) / (xi * xi - 1.0);       // P_n'(xi)
            double dx = -p1 / dpn;
            xi += dx;
            if (std::fabs(dx) < 1e-15) break;
        }
        double wi = 2.0 / ((1.0 - xi * xi) * dpn * dpn);      // GL weight on [-1,1]
        x[i] = 0.5 * (xi + 1.0);                             // map to [0,1]
        w[i] = 0.5 * wi;
    }
}

// Evaluate the nodal radial basis N_k and its derivative N_k' at parameter t in (0,1].
//   N_1 = t (vertex/trace), N_k = (L_k(2t-1) - L_{k-2}(2t-1))/(2k-1)  (bubble, k>=2)
//   N_1' = 1,               N_k' = L_{k-1}(2t-1) * 2
static void radial_basis_at(int P, double t, std::vector<double>& N, std::vector<double>& Np) {
    N.assign(P, 0.0);
    Np.assign(P, 0.0);
    N[0] = t;
    Np[0] = 1.0;
    std::vector<double> L;
    legendre_all(P, 2.0 * t - 1.0, L);
    for (int k = 2; k <= P; ++k) {
        N[k - 1] = (L[k] - L[k - 2]) / (2.0 * k - 1.0);
        Np[k - 1] = L[k - 1] * 2.0;
    }
}

void RadialOperators(int P, double a, std::vector<double>& R1, std::vector<double>& R0,
                     std::vector<double>& g, int nq) {
    std::vector<double> t, w;
    GaussLegendre01(nq, t, w);
    R1.assign(P * P, 0.0);
    R0.assign(P * P, 0.0);
    std::vector<double> N, Np;
    for (int q = 0; q < nq; ++q) {
        radial_basis_at(P, t[q], N, Np);
        const double inv_t2 = 1.0 / (t[q] * t[q]);
        for (int k = 0; k < P; ++k) {
            const double wk_Np = w[q] * Np[k];
            const double wk_N = w[q] * N[k] * inv_t2;
            for (int l = 0; l < P; ++l) {
                R1[k * P + l] += wk_Np * Np[l];     // a * int N_k' N_l' dt
                R0[k * P + l] += wk_N * N[l];        // a * int N_k N_l / t^2 dt
            }
        }
    }
    for (int i = 0; i < P * P; ++i) { R1[i] *= a; R0[i] *= a; }
    // trace g_k = N_k(t=1): N_1(1)=1, bubbles 0  -> g = e_1
    radial_basis_at(P, 1.0, N, Np);
    g.assign(N.begin(), N.end());
}

int DtNSurfaceOperator(const double* Mtil, const double* Ktil, int N, int P, double a,
                       std::vector<double>& S, int nq) {
    std::vector<double> R1, R0, g;
    RadialOperators(P, a, R1, R0, g, nq);

    // block (k,l) of the P-level tensor operator: A^{kl} = R1_kl*Mtil + R0_kl*Ktil  (N x N)
    auto fill_block = [&](int k, int l, double* dst, int ld, int roff, int coff) {
        const double r1 = R1[k * P + l], r0 = R0[k * P + l];
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                dst[(roff + i) * ld + (coff + j)] = r1 * Mtil[i * N + j] + r0 * Ktil[i * N + j];
    };

    // A00 (trace-trace, N x N) -> seed S
    S.assign((size_t)N * N, 0.0);
    fill_block(0, 0, S.data(), N, 0, 0);
    if (P == 1) return 0;   // no bubbles to condense

    const int nb = (P - 1) * N;                       // bubble-block dimension
    // A0b (N x nb, trace-vs-bubbles) and Abb (nb x nb, bubble-bubble), row-major
    std::vector<double> A0b((size_t)N * nb, 0.0);
    std::vector<double> Abb((size_t)nb * nb, 0.0);
    for (int bl = 0; bl < P - 1; ++bl)                // bubble column block (level l = bl+1)
        fill_block(0, bl + 1, A0b.data(), nb, 0, bl * N);
    for (int bk = 0; bk < P - 1; ++bk)                // bubble row block (level k = bk+1)
        for (int bl = 0; bl < P - 1; ++bl)
            fill_block(bk + 1, bl + 1, Abb.data(), nb, bk * N, bl * N);

    // X = Abb^{-1} * A0b^T   (nb x N).  Build RHS = A0b^T (nb x N), then dgesv overwrites it with X.
    std::vector<double> X((size_t)nb * N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int c = 0; c < nb; ++c)
            X[(size_t)c * N + i] = A0b[(size_t)i * nb + c];   // (A0b^T)[c,i] = A0b[i,c]
    std::vector<MKL_INT> ipiv(nb, 0);
    MKL_INT info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (MKL_INT)nb, (MKL_INT)N,
                                 Abb.data(), (MKL_INT)nb, ipiv.data(), X.data(), (MKL_INT)N);
    if (info != 0) return (int)info;

    // S = A00 - A0b * X   (N x N)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, nb,
                -1.0, A0b.data(), nb, X.data(), N, 1.0, S.data(), N);

    // symmetrize (kill round-off asymmetry)
    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j) {
            double s = 0.5 * (S[(size_t)i * N + j] + S[(size_t)j * N + i]);
            S[(size_t)i * N + j] = s;
            S[(size_t)j * N + i] = s;
        }
    return 0;
}

}  // namespace rad_ie
