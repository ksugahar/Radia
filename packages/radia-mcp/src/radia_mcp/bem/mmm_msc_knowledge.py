"""Magnetic Moment Method (MMM) and Magnetic Surface Charge (MSC) — Radia's core."""

MMM_OVERVIEW = r"""
# Magnetic Moment Method (MMM) — Radia's CORE method

★ This is the foundation of `radia` C++.  Not "BEM" in the Galerkin RWG
sense, but a closely related integral equation method.

## Formulation

For magnetostatics in a volume V_iron with magnetization M(r):

```
H(r) = H_ext(r) + N · M(r)
```

where N is the demagnetization tensor: for source M(r'), H_demag(r) is
given by the **volume integral**

```
H_demag(r) = -(1/4π) ∇ ∫_V ∇' · [M(r')/|r-r'|] dV'
           = -(1/4π) ∫_V [3(M·r̂)r̂ - M] / |r-r'|³  dV'   (dipole form)
```

For linear material M = χ·H, this gives a closed system:

```
M(r) = χ(r) · [H_ext(r) + N · M(r)]
(I - χ·N) M = χ · H_ext
```

Solve for M, then compute H, B everywhere.

## Discretization (constant M per element)

Discretize V_iron into N tetrahedra T_e, each with constant M_e ∈ R³.
- DOFs: 3N (Mx_e, My_e, Mz_e for each element)
- Interaction matrix N_{ef}: H_demag at center of element e due to
  unit M at element f
- Each N_{ef} is a 3×3 dense block

System size: 3N × 3N dense.

## Why MMM and not Galerkin BEM (ngsolve.bem)?

| Feature | MMM (Radia) | Galerkin RWG-BEM |
|---------|-------------|-------------------|
| DOFs | 3 per tet | ~3 per edge (RWG) |
| Nonlinear material | natural (M = χ(H)·H) | hard (need iteration) |
| Mesh requirement | volume tets | surface triangles only |
| Asymptotic accuracy | O(h²) for smooth | O(h^p) for HO basis |
| Memory (dense) | (3N)² | (3N)² |
| H-matrix compression | YES (HACApK) | YES (ACA) |

For ACCELERATOR MAGNETS (Radia's heritage):
- Permanent magnets + nonlinear soft iron yokes
- MMM is natural: it directly solves for M
- Galerkin BEM would require iteration for nonlinear materials

## Code

```python
import radia as rad

# Hex magnet with magnetization
m = rad.ObjHexahedron(verts, [Mx, My, Mz])  # M in A/m

# Soft iron with BH curve
iron = rad.ObjHexahedron(iron_verts, [0, 0, 0])  # zero initial M
rad.MatApl(iron, rad.MatSatIsoTab(BH_DATA))

# Solve nonlinear system
assembly = rad.ObjCnt([m, iron])
result = rad.Solve(assembly, 1e-4, 1000, 2)  # method=2 = HACApK
B = rad.Fld(assembly, 'b', [0, 0, 0.05])
```

See `radia_mcp.radia_ngsolve.radia_usage` for the full code recipe.

## References (this folder)

[LOCAL] 02_mmm_surface_charge/radia/A_3D_magnetostatics_computer_code_insertion_devices.pdf
        — Chubar-Elleaume-Chavanne 1998 (Radia paper)
[LOCAL] 02_mmm_surface_charge/Efficient_IE_Method_3D_Magnetostatic.pdf
[LOCAL] 02_mmm_surface_charge/Magnetostatic_IE_Greens_theorems.pdf
[LOCAL] 02_mmm_surface_charge/nonlinear/*.pdf — nonlinear MMM
[LOCAL] 02_mmm_surface_charge/nonlinear/LargeScale_Fast_Nonlinear_MMM_ACA.pdf
        — MMM + HACApK ACA acceleration
"""


MSC_HEXAHEDRA = r"""
# Magnetic Surface Charge (MSC) — for hexahedra and wedges

For hexahedral mesh (6 faces per element), Radia uses MSC instead of MMM.

## Formulation

The bound surface charge σ_bound = M · n̂ on each face is the unknown.
Use solid-angle integration (EIEM2 convention) for the integral
equation:

```
∫_face_i G(r, r') σ(r') dS' = -H_n(r_i)
```

with evaluation point at the **half-center** (CLAUDE.md POLICY):

```
EvalPt = 0.5 * (FaceCenter[i] + ElementCenter)
```

## DOFs per element

| Element | DOFs |
|---------|------|
| Hexahedron (6 faces) | 6 (one σ per face) |
| Wedge (5 faces) | 5 |
| Tet (4 faces) | 4 (but MMM is preferred → 3 DOFs instead) |

## Mixed elements

Radia supports mixed hex+wedge+tet meshes with variable DOF offsets:
- `m_elemDOF`, `m_elemDOFOffset`, `m_totalDOF` arrays

Cross-block interaction matrices: 3×3 (tet-tet), 5×5 (wedge), 6×6 (hex),
and 3-5-6 cross blocks.

## Why MSC for hexahedra?

For hex elements, surface charge gives a **better-conditioned** system
than volume MMM:
- Hex face integrals have closed-form solid-angle formulas
- Lower DOF count for similar accuracy
- Less prone to spurious internal currents

For tetrahedra, MMM is preferred:
- 3 DOFs vs 4 (saves one)
- More natural for nonlinear iron

## Code (single hex)

```python
verts = [(0,0,0), (1,0,0), (1,1,0), (0,1,0),
         (0,0,1), (1,0,1), (1,1,1), (0,1,1)]
hex = rad.ObjHexahedron(verts, [Mx, My, Mz])  # M in A/m
```

For mesh-imported hexahedra, use `radia.netgen_mesh_to_radia()`.

## References

[LOCAL] 02_mmm_surface_charge/3D_Surface_Charge_Model_Permeability_GeneralTheory.pdf
[LOCAL] 02_mmm_surface_charge/3D_Magnetostatic_SingleLayer_BIE_DifferenceField.pdf
[LOCAL] 02_mmm_surface_charge/Solution_Dense_Linear_IE_Formulations.pdf

## Code reference

- `src/core/rad_polyhedron.cpp` — element dispatch
- `src/core/rad_poly_analytical.cpp` — triangle/quad integration
- `src/core/rad_interaction.cpp` — interaction matrix
- See `radia_mcp.radia_ngsolve.radia_usage`
"""


RNA_MMM_COUPLING = r"""
# Reluctance Network + MMM Coupling (transformer modeling)

For transformer / reactor modeling, the lab has explored coupling
MMM (full Maxwell on iron core) with Reluctance Network Analysis
(circuit-level) for efficiency.

## References (this folder)

[LOCAL] 02_mmm_surface_charge/rna_mmm/RNA_MMM_Mixed_Transformer.pdf
[LOCAL] 02_mmm_surface_charge/rna_mmm/RNA_Simplified_MMM_CT.pdf
[LOCAL] 02_mmm_surface_charge/rna_mmm/PEEC_MMM_Coupling.pdf
        — PEEC (winding) + MMM (core) coupled

## When to use

| Problem | Pure MMM | RNA+MMM hybrid |
|---------|----------|----------------|
| Full transformer (geom-detailed) | works (slow) | faster |
| Lumped equivalent for SPICE | hard | natural |
| Saturation in core | natural | natural |
| Leakage flux estimation | accurate | approximate |

## Status in lab

Documented in this folder; not yet a production path in
`radia.sparsesolv_ngsolve` or radia panels.  Production transformer
modeling uses FEM (NGSolve) + PEEC for windings.
"""


def get_mmm_msc_knowledge(topic: str = "mmm_overview") -> str:
    """Dispatch MMM/MSC topics.

    Topics:
        mmm_overview    - Magnetic Moment Method foundation (DEFAULT)
        msc_hexahedra   - Magnetic Surface Charge for hex elements
        rna_mmm         - RNA + MMM coupling (transformer modeling)
        all             - Everything
    """
    topic = topic.lower().strip()
    if topic in ("mmm_overview", "mmm", "overview"):
        return MMM_OVERVIEW
    if topic in ("msc_hexahedra", "msc", "hex"):
        return MSC_HEXAHEDRA
    if topic in ("rna_mmm", "rna", "transformer"):
        return RNA_MMM_COUPLING
    if topic == "all":
        return "\n\n".join([MMM_OVERVIEW, MSC_HEXAHEDRA, RNA_MMM_COUPLING])
    return (f"Unknown topic '{topic}'. Available: mmm_overview, "
            "msc_hexahedra, rna_mmm, all.")
