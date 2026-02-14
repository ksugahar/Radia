# PEEC Panel: Analytical Edge Integration (FastImp Approach)

**Date**: 2026-02-13
**Status**: Implementation Plan (Based on FastImp Source Code Analysis)

---

## Key Insight from FastImp

FastImp does **NOT** use Gaussian quadrature for panel-panel potential coefficients.

Instead, it uses **analytical edge-based integration** (Hess-Smith and Newman methods).

**Why**:
- ✅ Exact handling of singularities (r → 0)
- ✅ Higher accuracy for near-field interactions
- ✅ Computationally efficient

---

## Three Methods for Panel-Panel Interaction

### Method 1: Analytical Diagonal (Self-Interaction)

For panel interacting with itself (i = j):

```
P_ii (double layer) = 2π  (analytical)
```

**No integration needed** - this is an exact result from potential theory.

### Method 2: Newman Method (Near-Field)

For panels that are close but not identical:

**Condition**: `distance(panel_i, panel_j) < k * characteristic_size`

**Approach**: Vertex-based analytical formulas

**Reference**: J. N. Newman, "Distributions of sources and normal dipoles over a quadrilateral panel," 1986.

### Method 3: Hess-Smith Method (Far-Field)

For well-separated panels:

**Approach**: Integration along panel edges using analytical formulas

**Reference**: J. L. Hess and A. M. O. Smith, "Calculation of potential flow about arbitrary bodies," 1967.

---

## Implementation Plan for Radia

### Phase 1: Triangle Self-Potential (Analytical)

**Formula** (from Wilton et al., 1984):

For a flat triangular panel with vertices v₀, v₁, v₂:

```
P_self = (1 / 4πε₀) * Σ_edges [l_edge * ln(...)]
```

Where integration is performed analytically along each edge.

**Code Structure**:

```cpp
double PEECMatrixBuilder::SelfPotentialTriangle(const PEECPanel& triangle) const {
    // Analytical formula based on edge lengths and vertex positions
    // NO numerical integration needed

    double sum = 0;
    for (int edge = 0; edge < 3; ++edge) {
        TVector3d v0 = triangle.vertices[edge];
        TVector3d v1 = triangle.vertices[(edge + 1) % 3];
        TVector3d v2 = triangle.vertices[(edge + 2) % 3];

        // Edge vector
        TVector3d edge_vec = v1 - v0;
        double l_edge = edge_vec.Abs();

        // Projection of v2 onto edge
        TVector3d r0 = v2 - v0;
        double t1 = dot(r0, edge_vec) / (l_edge * l_edge);
        double h = cross(r0, edge_vec).Abs() / l_edge;  // Height

        // Analytical edge integral (from Wilton)
        double R0 = r0.Abs();
        double R1 = (v2 - v1).Abs();
        double ln_term = std::log((R0 + R1 + l_edge) / (R0 + R1 - l_edge));

        sum += l_edge * ln_term;
    }

    return sum / (4.0 * RadConst::PI * PEEC_EPS_0);
}
```

### Phase 2: Triangle-Triangle Mutual Potential

**Far-Field** (distance > 3 * panel_size):

Use **centroid approximation**:
```
P_ij ≈ (Area_i * Area_j) / (4πε₀ * |r_i - r_j|)
```

**Near-Field** (distance < 3 * panel_size):

Use **edge-based integration** (Hess-Smith):
```
P_ij = (1 / 4πε₀) * Σ_edges_i Σ_edges_j [analytical_edge_integral(e_i, e_j)]
```

**Code Structure**:

```cpp
double PEECMatrixBuilder::MutualPotentialTriangle(const PEECPanel& tri_i,
                                                  const PEECPanel& tri_j) const {
    // Distance check
    double dist = distance(tri_i.centroid, tri_j.centroid);
    double char_size = std::sqrt(std::max(tri_i.area, tri_j.area));

    if (dist > 3.0 * char_size) {
        // Far-field: Centroid approximation
        return (tri_i.area * tri_j.area) /
               (4.0 * RadConst::PI * PEEC_EPS_0 * dist);
    } else {
        // Near-field: Analytical edge integration (Hess-Smith)
        return HessSmithIntegration(tri_i, tri_j);
    }
}

double PEECMatrixBuilder::HessSmithIntegration(const PEECPanel& tri_i,
                                               const PEECPanel& tri_j) const {
    // Integrate 1/R over both panels using edge-based formulas
    // Reference: Hess & Smith, 1967

    double sum = 0;
    for (int edge_i = 0; edge_i < 3; ++edge_i) {
        for (int edge_j = 0; edge_j < 3; ++edge_j) {
            // Analytical formula for edge-edge interaction
            sum += EdgeEdgeIntegral(tri_i, edge_i, tri_j, edge_j);
        }
    }

    return sum / (4.0 * RadConst::PI * PEEC_EPS_0);
}
```

### Phase 3: Quadrilateral Panels

**Split into two triangles** and sum contributions:

```cpp
double PEECMatrixBuilder::SelfPotentialQuad(const PEECPanel& quad) const {
    // Split quad into two triangles
    PEECPanel tri1 = MakeTriangle(quad.vertices[0], quad.vertices[1], quad.vertices[2]);
    PEECPanel tri2 = MakeTriangle(quad.vertices[0], quad.vertices[2], quad.vertices[3]);

    // Sum self-potentials
    return SelfPotentialTriangle(tri1) + SelfPotentialTriangle(tri2);
}
```

---

## Comparison: Analytical vs Gaussian Quadrature

| Aspect | Analytical Edge Integration | Gaussian Quadrature |
|--------|----------------------------|---------------------|
| **Self-interaction** | Exact | Diverges (singularity) |
| **Near-field** | High accuracy | Poor (near-singular) |
| **Far-field** | Exact or centroid approx | Accurate |
| **Computation** | Edge-based (9 integrals for tri-tri) | Point-based (N² evaluations) |
| **Implementation** | Complex (analytical formulas) | Simple (quadrature rules) |

**Recommendation**: Use **analytical edge integration** (FastImp approach) for production.

---

## References

1. **J. N. Newman**, "Distributions of sources and normal dipoles over a quadrilateral panel," J. Engineering Mathematics, vol. 20, pp. 113-126, 1986.

2. **J. L. Hess and A. M. O. Smith**, "Calculation of potential flow about arbitrary bodies," Progress in Aerospace Sciences, vol. 8, pp. 1-138, 1967.

3. **D. R. Wilton, S. M. Rao, A. W. Glisson, D. H. Schaubert, O. M. Al-Bundak, and C. M. Butler**, "Potential integrals for uniform and linear source distributions on polygonal and polyhedral domains," IEEE Trans. Antennas and Propagation, vol. 32, no. 3, pp. 276-281, Mar. 1984.

4. **R. D. Graglia**, "On the numerical integration of the linear shape functions times the 3-D Green's function or its gradient on a plane triangle," IEEE Trans. Antennas and Propagation, vol. 41, no. 10, pp. 1448-1455, Oct. 1993.

5. **FastImp source code**: https://github.com/ediloren/FastImp
   - Key files: `calcpForOneOverR.h`, `element.cc`, `formulation.cc`

---

## Implementation Priority

1. **HIGH**: Triangle self-potential (analytical formula from Wilton et al.)
2. **HIGH**: Triangle-triangle mutual potential (centroid approx for far-field)
3. MEDIUM: Hess-Smith edge integration for near-field
4. MEDIUM: Quadrilateral support (split into triangles)
5. LOW: Newman method for very close panels

---

**Status**: Ready for implementation
**Next Action**: Implement triangle self-potential using Wilton analytical formula
