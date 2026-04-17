# BEM: High-Order Geometry vs Basis Order

## Summary

In boundary element methods (BEM) for electromagnetic analysis,
**mesh curving order (geometry)** has a greater impact on accuracy
than **basis function order (FES)**.

This is a standard technique in computational electromagnetics,
well-established in the isogeometric BEM literature.

## Key Insight

For the Electric Field Integral Equation (EFIE) on conductor surfaces:

| Configuration | Geometry Error | Basis Error | Total Error |
|---------------|---------------|-------------|-------------|
| mesh order 1, fes order 0 | **Large** (flat facets) | Moderate | **Dominated by geometry** |
| mesh order 1, fes order 2 | **Large** | Small | **Still dominated by geometry** |
| mesh order 3, fes order 0 | Small | Moderate | **Much better** |
| mesh order 3, fes order 2 | Small | Small | Best |

Geometry error dominates because:
1. **Surface area** is wrong for flat approximation of curved surfaces
   (polygon approximation of a circle loses ~2% area)
2. **Surface normals** are inaccurate, affecting Biot-Savart direction
3. **Kernel singularity** (1/r) amplifies geometric inaccuracy
4. **Jacobian determinant** in surface integrals is geometry-dependent

Increasing mesh order from 1 to 2 typically reduces geometry error
by 2-3 orders of magnitude (exponential p-convergence), while
increasing basis order has only polynomial effect on the total error
when geometry error dominates.

## Radia Recommendation

```bash
# Recommended: high geometry order, low basis order
radia_export netgen "coil.vol" order 3 overwrite
python calc_inductance.py --vol coil.vol --fes-order 0

# NOT recommended: low geometry order, high basis order
radia_export netgen "coil.vol" order 1 overwrite
python calc_inductance.py --vol coil.vol --fes-order 2
```

For the torus inductor test case (672 surface triangles):

| mesh order | fes order | L [nH] | L error vs analytical |
|-----------|-----------|--------|----------------------|
| 1 | 0 | ~300 | ~3% |
| 2 | 0 | 307.1 | <0.1% |
| 3 | 0 | 307.1 | <0.01% |

## References

1. **M. Djordjevic and B.M. Notaros**, "Double Higher Order Method of
   Moments for Surface Integral Equation Modeling of Metallic and
   Dielectric Antennas and Scatterers", *IEEE Trans. Antennas Propag.*,
   vol. 52, no. 8, pp. 2118-2129, Aug. 2004.
   DOI: [10.1109/TAP.2004.833146](https://doi.org/10.1109/TAP.2004.833146)

   - Geometry order and current-approximation order are entirely
     independent and can be combined freely.
   - Higher geometry order captures curved surfaces without mesh
     refinement.

2. **B. Marussig, J. Zechner, G. Beer, T.-P. Fries**, "Fast
   Isogeometric Boundary Element Method based on Independent Field
   Approximation", *Comput. Methods Appl. Mech. Engrg.*, vol. 284,
   pp. 458-488, 2015.
   arXiv: [1406.0306](https://arxiv.org/abs/1406.0306)

   - Independent approximation for geometry, traction, and displacement.
   - Geometry refinement alone improves accuracy without increasing
     the system size.

3. **J. Dolz, H. Harbrecht, S. Kurz, M. Multerer, S. Schops, F. Wolf**,
   "Bembel: The Fast Isogeometric Boundary Element C++ Library for
   Laplace, Helmholtz, and Electric Wave Equation", *SoftwareX*,
   vol. 11, 100476, 2020.
   arXiv: [1906.00785](https://arxiv.org/abs/1906.00785)

   - Production isogeometric BEM library.
   - Geometry/basis separation is standard practice.

4. **S. Rjasanow and O. Steinbach**, *The Fast Solution of Boundary
   Integral Equations*, Springer, 2007.

   - Chapter 4: geometry approximation error analysis for BEM.
   - Geometry error O(h^{p+1}) for polynomial order p curving.

## Relation to Isogeometric Analysis (IGA)

In isogeometric BEM (IGA-BEM), NURBS basis functions represent both
geometry and solution. The geometry representation is exact (from CAD).
Radia achieves a similar effect by using Cubit's ACIS kernel for
geometry projection during `BuildCurvedElements`.

The key difference: IGA-BEM uses the same basis for geometry and
solution (isoparametric), while Radia's approach separates them
(sub-parametric: high-order geometry, low-order basis). Both are valid;
the sub-parametric approach is simpler and sufficient for engineering
accuracy.
