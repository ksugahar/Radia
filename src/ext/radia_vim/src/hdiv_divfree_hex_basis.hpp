// hdiv_divfree_hex_basis.hpp — HDiv div-free interior basis on the unit hex.
//
// Mathematica reference: ShapesHDivHexDivFree in
//   ngsolve_validation/hex_bem_hdiv_order3_overnight.wls
//
// The basis is a Schöberl-Zaglmayr-style interior basis on the reference cube
// [0,1]^3 (ξ = 2x-1, η = 2y-1, ζ = 2z-1 in [-1,1]^3) with the property that
// every basis function satisfies the conductor BC J·n = 0 on all 6 faces and
// has zero divergence inside. Total DOFs = 2p^3 + 3p^2  where p is the order.
//
// 5 sub-blocks, parameterised by integer indices (i, j, k) with i,j,k ∈ [0..p-1]:
//   Block A (p^3 shapes): w = -4 P_{j+1}(η),   Cross[grad(L_{i+2}(ξ) L_{k+2}(ζ)), e_y]
//   Block B (p^3 shapes): w =  4 P_{i+1}(ξ),   Cross[grad(L_{j+2}(η) L_{k+2}(ζ)), e_x]
//   Block C (p^2 shapes):                       Cross[grad(L_{j+2}(η) L_{k+2}(ζ)), grad(2ξ)]
//   Block D (p^2 shapes):                       Cross[grad(L_{i+2}(ξ) L_{k+2}(ζ)), grad(2η)]
//   Block E (p^2 shapes):                       Cross[grad(L_{i+2}(ξ) L_{j+2}(η)), grad(2ζ)]
//
// L_n(ξ) is the Mathematica `IntegratedLegendre` polynomial; it vanishes at
// ξ = ±1 for n ≥ 2, which is what enforces the J·n = 0 BC. P_n(ξ) is the
// classical Legendre polynomial.
//
// Header-only, type-templated. Boost.Multiprecision instantiations work
// transparently if the corresponding RADIA_VIM_WITH_QUAD / RADIA_VIM_WITH_MPFR
// CMake flags are on.

#pragma once

#include <vector>
#include <array>
#include <stdexcept>
#include <string>

namespace radia_vim {

template <typename T>
struct Vec3 {
    T x{}, y{}, z{};
};

enum class HexBasisBlock : int {
    A = 0, B = 1, C = 2, D = 3, E = 4
};

struct HexBasisDofKey {
    HexBasisBlock block;
    int i, j, k;            // tensor index; some are unused depending on block
};

template <typename T>
class HDivDivFreeHexBasis {
public:
    explicit HDivDivFreeHexBasis(int order) : order_(order) {
        if (order < 1) throw std::invalid_argument(
            "HDivDivFreeHexBasis: order must be >= 1");
        build_index_map();
    }

    int order()  const { return order_; }
    int n_dofs() const { return static_cast<int>(index_map_.size()); }

    // The (block, i, j, k) tuple describing DOF idx.
    HexBasisDofKey dof_key(int idx) const { return index_map_.at(idx); }

    // Evaluate basis function `idx` at reference point (x, y, z) in [0,1]^3.
    Vec3<T> evaluate(int idx, T x, T y, T z) const;

private:
    int order_;
    std::vector<HexBasisDofKey> index_map_;

    void build_index_map();

    // Integrated-Legendre polynomials on ξ ∈ [-1,1].
    //   out[m] = L_{m+2}(ξ)         for m = 0 .. p-1
    //   d_out[m] = dL_{m+2}/dξ      (derivative w.r.t. ξ, NOT w.r.t. x)
    void integrated_legendre_with_deriv(T xi,
                                        std::vector<T>& out,
                                        std::vector<T>& d_out) const;

    // Classical Legendre polynomials on ξ ∈ [-1,1].
    //   out[m] = P_m(ξ)             for m = 0 .. p   (length p+1)
    void legendre(T xi, std::vector<T>& out) const;
};


// =============================================================================
// Implementation
// =============================================================================

template <typename T>
void HDivDivFreeHexBasis<T>::build_index_map() {
    index_map_.clear();
    const int p = order_;
    index_map_.reserve(2 * p * p * p + 3 * p * p);
    // Block A: i, j, k in [0..p-1]^3
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        for (int k = 0; k < p; ++k)
          index_map_.push_back({HexBasisBlock::A, i, j, k});
    // Block B: i, j, k in [0..p-1]^3
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        for (int k = 0; k < p; ++k)
          index_map_.push_back({HexBasisBlock::B, i, j, k});
    // Block C: j, k in [0..p-1]^2 (i unused — set to 0)
    for (int j = 0; j < p; ++j)
      for (int k = 0; k < p; ++k)
        index_map_.push_back({HexBasisBlock::C, 0, j, k});
    // Block D: i, k in [0..p-1]^2
    for (int i = 0; i < p; ++i)
      for (int k = 0; k < p; ++k)
        index_map_.push_back({HexBasisBlock::D, i, 0, k});
    // Block E: i, j in [0..p-1]^2
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        index_map_.push_back({HexBasisBlock::E, i, j, 0});
}


template <typename T>
void HDivDivFreeHexBasis<T>::integrated_legendre_with_deriv(
    T xi, std::vector<T>& L, std::vector<T>& Lp) const {
    const int p = order_;
    L.assign(p, T(0));
    Lp.assign(p, T(0));
    if (p < 1) return;

    // Mirror Mathematica IntegratedLegendre[p+1, xi]. The recursion appends
    // result for j = 2..(p+1), giving p elements.
    //   p3 = 0, p2 = -1, p1 = xi   (initial state, NOT appended)
    //   for each j: shift down (p3 ← p2, p2 ← p1), then
    //     p1 ← ((2j-3)/j) xi p2  -  ((j-3)/j) p3
    //
    // Derivative recursion: differentiate the assignment for p1.
    //   p1' ← ((2j-3)/j) (p2_new + xi p2'_new)  -  ((j-3)/j) p3'_new
    T p3   = T(0), p2   = T(-1), p1   = xi;
    T p3p  = T(0), p2p  = T(0),  p1p  = T(1);

    for (int j = 2; j <= p + 1; ++j) {
        T a = T(2 * j - 3) / T(j);
        T b = T(j - 3) / T(j);

        T p3_new  = p2;
        T p2_new  = p1;
        T p3p_new = p2p;
        T p2p_new = p1p;

        T new_p1  = a * xi * p2_new - b * p3_new;
        T new_p1p = a * (p2_new + xi * p2p_new) - b * p3p_new;

        L [j - 2] = new_p1;
        Lp[j - 2] = new_p1p;

        p3  = p3_new;  p2  = p2_new;  p1  = new_p1;
        p3p = p3p_new; p2p = p2p_new; p1p = new_p1p;
    }
}


template <typename T>
void HDivDivFreeHexBasis<T>::legendre(T xi, std::vector<T>& P) const {
    const int p = order_;
    P.assign(p + 1, T(0));
    P[0] = T(1);
    if (p >= 1) P[1] = xi;
    // P_k = (2 - 1/k) xi P_{k-1} + (1/k - 1) P_{k-2}
    //     = ((2k-1)/k) xi P_{k-1} - ((k-1)/k) P_{k-2}    (standard form)
    for (int k = 2; k <= p; ++k) {
        T inv_k = T(1) / T(k);
        P[k] = (T(2) - inv_k) * xi * P[k-1] + (inv_k - T(1)) * P[k-2];
    }
}


template <typename T>
Vec3<T> HDivDivFreeHexBasis<T>::evaluate(int idx, T x, T y, T z) const {
    const HexBasisDofKey& key = index_map_.at(idx);

    // Map [0,1] -> [-1,1]
    const T xi   = T(2) * x - T(1);
    const T eta  = T(2) * y - T(1);
    const T zeta = T(2) * z - T(1);

    std::vector<T> LXi, LXi_d, LEta, LEta_d, LZeta, LZeta_d;
    std::vector<T> PXi, PEta, PZeta;
    integrated_legendre_with_deriv(xi,   LXi,   LXi_d);
    integrated_legendre_with_deriv(eta,  LEta,  LEta_d);
    integrated_legendre_with_deriv(zeta, LZeta, LZeta_d);
    legendre(xi,   PXi);
    legendre(eta,  PEta);
    legendre(zeta, PZeta);

    // --- Derivation note ----------------------------------------------------
    // For a function u(ξ, η, ζ), grad_x(u) = (∂u/∂x, ∂u/∂y, ∂u/∂z)
    //   = (2 ∂u/∂ξ, 2 ∂u/∂η, 2 ∂u/∂ζ)  via chain rule. The factor of 2
    // appears N times when evaluating Cross[grad(u), grad(v)]; for v ∈ {x,y,z}
    // we have grad(v) ∈ {e_x, e_y, e_z} so Cross gives a single ∂u/∂* term;
    // for v = 2 ξ etc. grad(v) carries a factor of 4. The constants below
    // (8 = 4·2 and -4·2) are the products that arise after all chain rules.
    //
    // Mathematica array indexing in the .wls file:
    //   LXi[[i+1]]  ↔  C++ LXi[i]   = L_{i+2}(ξ)        for i=0..p-1
    //   PEta[[j+2]] ↔  C++ PEta[j+1] = P_{j+1}(η)        for j=0..p-1
    //   PXi[[i+2]]  ↔  C++ PXi[i+1]  = P_{i+1}(ξ)        for i=0..p-1

    const T eight = T(8);
    Vec3<T> S{};
    const int i = key.i, j = key.j, k = key.k;

    switch (key.block) {
      case HexBasisBlock::A:
        // w = -4 P_{j+1}(η);  Cross[grad(L_{i+2}(ξ) L_{k+2}(ζ)), e_y]
        // = (-∂u/∂z, 0, ∂u/∂x).  After w·... and chain rule:
        S.x =  eight * LXi  [i] * PEta[j + 1] * LZeta_d[k];
        S.y =  T(0);
        S.z = -eight * LXi_d[i] * PEta[j + 1] * LZeta  [k];
        break;
      case HexBasisBlock::B:
        // w =  4 P_{i+1}(ξ);  Cross[grad(L_{j+2}(η) L_{k+2}(ζ)), e_x]
        S.x =  T(0);
        S.y =  eight * PXi[i + 1] * LEta  [j] * LZeta_d[k];
        S.z = -eight * PXi[i + 1] * LEta_d[j] * LZeta  [k];
        break;
      case HexBasisBlock::C:
        // Cross[grad(L_{j+2}(η) L_{k+2}(ζ)), grad(2ξ) = (4,0,0)]
        S.x =  T(0);
        S.y =  eight * LEta  [j] * LZeta_d[k];
        S.z = -eight * LEta_d[j] * LZeta  [k];
        break;
      case HexBasisBlock::D:
        // Cross[grad(L_{i+2}(ξ) L_{k+2}(ζ)), grad(2η) = (0,4,0)]
        S.x = -eight * LXi  [i] * LZeta_d[k];
        S.y =  T(0);
        S.z =  eight * LXi_d[i] * LZeta  [k];
        break;
      case HexBasisBlock::E:
        // Cross[grad(L_{i+2}(ξ) L_{j+2}(η)), grad(2ζ) = (0,0,4)]
        S.x =  eight * LXi  [i] * LEta_d[j];
        S.y = -eight * LXi_d[i] * LEta  [j];
        S.z =  T(0);
        break;
    }
    return S;
}

}  // namespace radia_vim
