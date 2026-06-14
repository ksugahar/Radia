// sauter_schwab_hex.hpp — singular 6D quadrature for self-coincident
// hex × hex Newton-kernel integrals.
//
// Computes K[i,j] = ∫∫_{[0,1]^6} φ_i(r)·diag(a²,b²,c²)·φ_j(r')
//                                   / |diag(a,b,c)·(r-r')| dr dr'
// for the cuboid eddy-current Volume Integral Method. The 1/|r-r'|
// singularity along the diagonal is removed by spherical coordinates on
// u = r - r': in du = ρ²dρdω, the ρ² volume factor cancels the 1/ρ
// kernel, leaving a regular ρ·F(ρω) integrand.
//
// Substitution chain:
//   r = c + u/2,   r' = c - u/2,    dr dr' = du dc
//   Domain in (c, u): u ∈ [-1, 1]^3, c ∈ ∏_i [|u_i|/2, 1 - |u_i|/2]
//   u = ρω (spherical),  du = ρ² dρ dω,  ρ ∈ [0, 1/max|ω_i|]
//
//   K = ∫_S² dω (1/|Λω|) ∫_0^{ρ_max(ω)} dρ ρ · F_ij(ρω)
//   F_ij(u) = ∫_{c-box(u)} (φ_i(c+u/2) · diag(a²,b²,c²) · φ_j(c-u/2)) dc
//
// Header-only, type-templated. Quadrature: tensor-product Gauss-Legendre
// on (θ, φ, ρ, c_x, c_y, c_z) — 6 nested loops. Per-entry cost scales as
// N_θ · N_φ · N_ρ · N_c³.
//
// Phase F-3a: this file. Single-entry computation, double precision validated
// against Mathematica reference. No parallelism, no compression.
// Phase F-4 will add OpenMP-parallel full K-matrix assembly.

#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "hdiv_divfree_hex_basis.hpp"

namespace radia_vim {

// ---------------------------------------------------------------------------
// Compile-time pi (long-double precision; sufficient for double quadrature,
// will be replaced with boost::math::constants::pi<T>() when QUAD/MPFR enabled).
// ---------------------------------------------------------------------------
template <typename T>
inline constexpr T radia_vim_pi() {
    return T(3.141592653589793238462643383279502884L);
}

// ---------------------------------------------------------------------------
// Gauss-Legendre nodes/weights on [-1, 1] via Newton iteration on Legendre
// polynomial roots. Templated on scalar T.
// ---------------------------------------------------------------------------
template <typename T>
inline void gauss_legendre_m1_1(int n,
                                std::vector<T>& nodes,
                                std::vector<T>& weights) {
    using std::abs;
    using std::cos;

    if (n < 1) throw std::invalid_argument("gauss_legendre_m1_1: n must be >= 1");

    nodes.assign(n, T(0));
    weights.assign(n, T(0));

    const T pi = radia_vim_pi<T>();
    // Tolerance: ~ epsilon for the type. For double, 1e-15 is safe.
    const T tol = T(1e-15);
    const int max_iter = 100;

    for (int k = 1; k <= n; ++k) {
        // Initial guess (Tricomi/Petras): x ~ cos(π(k - 1/4)/(n + 1/2))
        T x = cos(pi * (T(k) - T(1)/T(4)) / (T(n) + T(1)/T(2)));
        T x_prev = x + T(1);

        T pn = T(0), pn_prev_inner = T(0);

        for (int iter = 0; iter < max_iter; ++iter) {
            // Compute P_n(x) and P_{n-1}(x) by recurrence.
            T p_prev = T(1);   // P_0
            T p      = x;      // P_1
            for (int m = 1; m < n; ++m) {
                T p_next = ((T(2*m + 1) * x * p - T(m) * p_prev) / T(m + 1));
                p_prev = p;
                p = p_next;
            }
            pn = p;                  // P_n(x)
            pn_prev_inner = p_prev;  // P_{n-1}(x)

            // P_n'(x) = n / (1 - x²) · (P_{n-1}(x) - x · P_n(x))
            T pn_deriv = T(n) / (T(1) - x*x) * (pn_prev_inner - x * pn);
            x_prev = x;
            x = x - pn / pn_deriv;
            if (abs(x - x_prev) <= tol * (abs(x) + tol)) break;
        }

        nodes[k - 1] = x;

        // Weight: 2 / ((1 - x²) · (P_n'(x))²)
        T pn_deriv = T(n) / (T(1) - x*x) * (pn_prev_inner - x * pn);
        weights[k - 1] = T(2) / ((T(1) - x*x) * pn_deriv * pn_deriv);
    }

    // Sort ascending.
    std::vector<int> idx(n);
    for (int k = 0; k < n; ++k) idx[k] = k;
    std::sort(idx.begin(), idx.end(),
              [&](int p, int q) { return nodes[p] < nodes[q]; });
    std::vector<T> ns(n), ws(n);
    for (int k = 0; k < n; ++k) {
        ns[k] = nodes[idx[k]];
        ws[k] = weights[idx[k]];
    }
    nodes   = std::move(ns);
    weights = std::move(ws);
}


// ---------------------------------------------------------------------------
// Single K-matrix entry via spherical Duffy on the unit reference cube.
// Returns the bare integral (Mathematica multiplies by mu_0/(4π) afterwards).
//
// Defaults n_θ=16, n_φ=32, n_ρ=16, n_c=8 → ~16·32·16·8³ ≈ 2.7M evals/entry.
// At ~50 flops/eval (basis evaluations dominate) this is ~135 MFlops, so
// ~30 ms per entry on a modern CPU.
// ---------------------------------------------------------------------------
template <typename T>
T compute_K_entry_spherical_duffy(
    const HDivDivFreeHexBasis<T>& basis,
    int i, int j,
    T a, T b, T c,
    int n_omega_theta = 16,
    int n_omega_phi   = 32,
    int n_rho         = 16,
    int n_c           = 8) {

    using std::sin;
    using std::cos;
    using std::sqrt;
    using std::abs;

    const T pi = radia_vim_pi<T>();

    std::vector<T> gl_t_n, gl_t_w;
    std::vector<T> gl_p_n, gl_p_w;
    std::vector<T> gl_r_n, gl_r_w;
    std::vector<T> gl_c_n, gl_c_w;

    gauss_legendre_m1_1(n_omega_theta, gl_t_n, gl_t_w);
    gauss_legendre_m1_1(n_omega_phi,   gl_p_n, gl_p_w);
    gauss_legendre_m1_1(n_rho,         gl_r_n, gl_r_w);
    gauss_legendre_m1_1(n_c,           gl_c_n, gl_c_w);

    T result = T(0);

    for (int it = 0; it < n_omega_theta; ++it) {
        const T theta    = (gl_t_n[it] + T(1)) * pi / T(2);    // [0, π]
        const T jac_th   = pi / T(2);
        const T s_th     = sin(theta);
        const T c_th     = cos(theta);

        for (int ip = 0; ip < n_omega_phi; ++ip) {
            const T phi_ang = (gl_p_n[ip] + T(1)) * pi;          // [0, 2π]
            const T jac_ph  = pi;
            const T c_ph    = cos(phi_ang);
            const T s_ph    = sin(phi_ang);

            const T omega_x = s_th * c_ph;
            const T omega_y = s_th * s_ph;
            const T omega_z = c_th;

            T omega_max = abs(omega_x);
            if (abs(omega_y) > omega_max) omega_max = abs(omega_y);
            if (abs(omega_z) > omega_max) omega_max = abs(omega_z);
            if (omega_max < T(1e-30)) continue;
            const T rho_max = T(1) / omega_max;

            const T Lomega_norm = sqrt(a*a*omega_x*omega_x +
                                        b*b*omega_y*omega_y +
                                        c*c*omega_z*omega_z);

            const T sphere_w = gl_t_w[it] * gl_p_w[ip] *
                                jac_th * jac_ph * s_th / Lomega_norm;

            for (int ir = 0; ir < n_rho; ++ir) {
                const T rho_val = (gl_r_n[ir] + T(1)) * rho_max / T(2);
                const T jac_rh  = rho_max / T(2);

                const T u_x = rho_val * omega_x;
                const T u_y = rho_val * omega_y;
                const T u_z = rho_val * omega_z;

                const T cx_lo = abs(u_x) / T(2);
                const T cx_hi = T(1) - abs(u_x) / T(2);
                const T cy_lo = abs(u_y) / T(2);
                const T cy_hi = T(1) - abs(u_y) / T(2);
                const T cz_lo = abs(u_z) / T(2);
                const T cz_hi = T(1) - abs(u_z) / T(2);
                const T cx_jac = (cx_hi - cx_lo) / T(2);
                const T cy_jac = (cy_hi - cy_lo) / T(2);
                const T cz_jac = (cz_hi - cz_lo) / T(2);
                if (cx_jac <= T(0) || cy_jac <= T(0) || cz_jac <= T(0)) continue;

                T c_inner = T(0);

                for (int icx = 0; icx < n_c; ++icx) {
                    const T cx = cx_lo + (gl_c_n[icx] + T(1)) * cx_jac;
                    const T wcx = gl_c_w[icx];
                    for (int icy = 0; icy < n_c; ++icy) {
                        const T cy = cy_lo + (gl_c_n[icy] + T(1)) * cy_jac;
                        const T wcy = gl_c_w[icy];
                        for (int icz = 0; icz < n_c; ++icz) {
                            const T cz = cz_lo + (gl_c_n[icz] + T(1)) * cz_jac;
                            const T wcz = gl_c_w[icz];

                            const T rx  = cx + u_x / T(2);
                            const T ry  = cy + u_y / T(2);
                            const T rz  = cz + u_z / T(2);
                            const T rpx = cx - u_x / T(2);
                            const T rpy = cy - u_y / T(2);
                            const T rpz = cz - u_z / T(2);

                            const Vec3<T> phi_i = basis.evaluate(i, rx, ry, rz);
                            const Vec3<T> phi_j = basis.evaluate(j, rpx, rpy, rpz);

                            const T dot = phi_i.x * phi_j.x * a*a +
                                          phi_i.y * phi_j.y * b*b +
                                          phi_i.z * phi_j.z * c*c;

                            c_inner += wcx * wcy * wcz * dot;
                        }
                    }
                }
                c_inner *= cx_jac * cy_jac * cz_jac;

                // ρ-integrand: ρ² (vol) · 1/(ρ |Λω|) (kernel) · c_inner = ρ · c_inner / |Λω|
                // Sphere prefactor (sin θ / |Λω|) absorbed into sphere_w.
                result += sphere_w * gl_r_w[ir] * jac_rh * rho_val * c_inner;
            }
        }
    }

    return result;
}

}  // namespace radia_vim
