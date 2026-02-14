# True 2D Panel Implementation for PEEC

**Date**: 2026-02-13
**Status**: Implementation Plan

---

## Current Limitation

**Point Approximation** (current):
```cpp
// Panel = Point with area
struct PEECNode {
    TVector3d position;  // Point position
    double area;         // Associated area
};

// Self-potential: Disk approximation
P_ii = 1 / (4π ε₀ √(π A))

// Mutual potential: Point-to-point distance
P_ij = 1 / (4π ε₀ r_ij)
```

**Problem**:
- Inaccurate for close panels (distance < panel size)
- Breaks down for large panels
- Not FastImp-compatible

---

## True 2D Panel Approach

**Surface Patch** (target):
```cpp
struct PEECPanel {
    std::vector<TVector3d> vertices;  // 3 (tri) or 4 (quad) vertices
    TVector3d center;                 // Panel centroid
    TVector3d normal;                 // Outward normal
    double area;                      // Panel area

    // Gauss quadrature points (precomputed)
    std::vector<TVector3d> quad_points;
    std::vector<double> quad_weights;
};

// Self-potential: Surface double integral
P_ii = ∫∫_panel ∫∫_panel' (1/4πε₀) dS dS' / |r - r'|

// Mutual potential: Double surface integral
P_ij = ∫∫_panel_i ∫∫_panel_j (1/4πε₀) dS_i dS_j / |r_i - r_j|
```

---

## Implementation Plan

### Phase 1: Data Structure

Add `PEECPanel` structure alongside existing `PEECNode`:

```cpp
// rad_peec_matrices.h

struct PEECPanel {
    std::vector<TVector3d> vertices;  // 3 or 4 vertices
    TVector3d center;                 // Centroid
    TVector3d normal;                 // Outward normal
    double area;                      // Panel area

    enum Type { Triangle, Quadrilateral } type;

    // Gauss quadrature (4-point for triangles, 4-point for quads)
    std::vector<TVector3d> quad_points;
    std::vector<double> quad_weights;

    PEECPanel();
    void ComputeQuadrature();  // Precompute integration points
};
```

### Phase 2: Quadrature Rules

**Triangle (4-point Gauss)**:
```
Points (barycentric coordinates):
  (1/3, 1/3, 1/3) - weight = -27/48
  (3/5, 1/5, 1/5) - weight = 25/48
  (1/5, 3/5, 1/5) - weight = 25/48
  (1/5, 1/5, 3/5) - weight = 25/48
```

**Quadrilateral (2x2 Gauss)**:
```
ξ = ±1/√3
Points in local coords:
  (-1/√3, -1/√3) - weight = 1
  (+1/√3, -1/√3) - weight = 1
  (+1/√3, +1/√3) - weight = 1
  (-1/√3, +1/√3) - weight = 1
```

### Phase 3: Surface Integrals

```cpp
// rad_peec_matrices.cpp

double PEECMatrixBuilder::SelfPotentialPanel(const PEECPanel& panel) const {
    // Double surface integral (singularity treatment needed)
    // Use analytical formula for flat panels

    double sum = 0;
    for (size_t i = 0; i < panel.quad_points.size(); ++i) {
        for (size_t j = 0; j < panel.quad_points.size(); ++j) {
            if (i == j) continue;  // Skip singular point

            TVector3d r1 = panel.quad_points[i];
            TVector3d r2 = panel.quad_points[j];
            double dist = distance(r1, r2);

            double w1 = panel.quad_weights[i];
            double w2 = panel.quad_weights[j];

            sum += w1 * w2 / dist;
        }
    }

    // Singularity correction (analytical for flat panels)
    double singularity_term = ComputeSingularityTerm(panel);

    return (sum + singularity_term) / (4.0 * RadConst::PI * PEEC_EPS_0);
}

double PEECMatrixBuilder::MutualPotentialPanel(const PEECPanel& panel_i,
                                               const PEECPanel& panel_j) const {
    // Double surface integral (well-separated panels)
    double sum = 0;

    for (size_t i = 0; i < panel_i.quad_points.size(); ++i) {
        for (size_t j = 0; j < panel_j.quad_points.size(); ++j) {
            TVector3d r1 = panel_i.quad_points[i];
            TVector3d r2 = panel_j.quad_points[j];
            double dist = distance(r1, r2);

            double w1 = panel_i.quad_weights[i];
            double w2 = panel_j.quad_weights[j];

            sum += w1 * w2 / dist;
        }
    }

    return sum / (4.0 * RadConst::PI * PEEC_EPS_0);
}
```

### Phase 4: Singularity Treatment

For self-potential, the integrand has a singularity at r = r'. Use analytical formulas:

**Triangle** (Wilton et al., 1984):
```
P_self(triangle) = (Area / 4πε₀) * f(geometry)
```

**Quadrilateral**:
Split into 2 triangles and sum.

### Phase 5: API Extension

```cpp
// rad_peec_matrices.h

class PEECMatrixBuilder {
public:
    // ... existing methods ...

    /**
     * @brief Add 2D surface panel (Star element)
     */
    void AddPanel(const PEECPanel& panel);

    /**
     * @brief Add panels from surface mesh
     */
    void AddPanelsFromSurfaceMesh(const std::vector<std::vector<TVector3d>>& triangles,
                                  const std::vector<std::vector<TVector3d>>& quads);

private:
    std::vector<PEECPanel> panels_;  // 2D panels (alternative to point nodes)
    bool use_true_panels_;           // Toggle between point/panel mode

    double SelfPotentialPanel(const PEECPanel& panel) const;
    double MutualPotentialPanel(const PEECPanel& p1, const PEECPanel& p2) const;
    double ComputeSingularityTerm(const PEECPanel& panel) const;
};
```

---

## Python API

```python
from peec_matrices import PEECBuilder

builder = PEECBuilder()

# Add filaments (1D)
for edge in filament_edges:
    builder.create_wire(p0, p1, width, height, 1, sigma)

# Add panels (2D) - NEW API
for triangle in panel_triangles:
    vertices = [coords[n] for n in triangle]
    builder.add_panel_triangle(vertices)

for quad in panel_quads:
    vertices = [coords[n] for n in quad]
    builder.add_panel_quad(vertices)

# Build matrices with true 2D panel integration
L, R, P, M_LS = builder.build(use_true_panels=True)
```

---

## Testing Plan

1. **Unit test**: Single triangle self-potential vs analytical
2. **Unit test**: Two triangles mutual potential vs point approximation
3. **Validation**: Parallel plate capacitor (known analytical result)
4. **Benchmark**: Compare point vs panel for circular coil

---

## Priority

- **HIGH**: Self-potential and mutual potential for triangles
- MEDIUM: Quadrilateral support (split into triangles as workaround)
- LOW: Adaptive quadrature for high-accuracy cases

---

## References

1. R. D. Graglia, "On the numerical integration of the linear shape functions times the 3-D Green's function or its gradient on a plane triangle," IEEE TAP, 1993.
2. D. R. Wilton et al., "Potential integrals for uniform and linear source distributions on polygonal and polyhedral domains," IEEE TAP, 1984.
3. Z. Zhu et al., "Algorithms in FastImp," IEEE TCAD, 2005.

---

**Status**: Implementation plan ready
**Next Action**: Implement `PEECPanel` structure and triangle self-potential
