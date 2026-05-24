"""Reluctance Network Analysis (RNA) / Magnetic Equivalent Circuit (MEC)."""

OVERVIEW = r"""
# Reluctance Network Analysis (RNA) / Magnetic Equivalent Circuit (MEC)

Convert a magnetic circuit into an electrical-circuit-like network with:
- **Reluctance R = l / (μ A)** (analogous to electrical resistance)
- **MMF F = N I** (analogous to voltage source)
- **Flux Φ = F / R** (analogous to current)

Solve via:
- Nodal analysis (flux at nodes)
- Mesh analysis (loop fluxes)

## When to use RNA vs FEM

| Scenario | RNA | FEM |
|----------|-----|-----|
| Quick concept design | ★ best | overkill |
| Lumped equivalent for SPICE coupling | ★ natural | hard |
| Detailed flux distribution | rough | ★ best |
| Saturation in core | works (lookup) | works (BH curve) |
| Hysteresis (MEC + dynamic) | ★ lab focus | works (Play/Energy) |
| Leakage flux | hard | natural |
| Optimization (sensitivity) | fast | slow |

## Nodal vs Mesh-based MEC

| Method | Pros | Cons |
|--------|------|------|
| **Nodal** (Φ unknown at nodes) | Sparse matrix, KCL-like | Need cohomology cuts for multiply-connected |
| **Mesh** (Φ_loop unknown) | Cleaner for multiply-connected | More unknowns, KVL-like |

Reference comparison: "A Comparison of Nodal- and Mesh-Based Magnetic
Equivalent Circuit" (in W:/.../13_RNA/external/).

## Lab core papers (lab_sugahara/)

田中 / 羽根 lineage:
- 1609 マグ研_田中
- 1702 マグ研_田中
- 1801 磁気学会_羽根_採録
- 1803 マグ研_羽根
- 1902 マグ研_羽根
- 修士論文_本審査_久田
- 博士論文_本審査_吉田
- AC_DC_Bridgeless_Flyback

## Cross-references

- `radia_mcp.bem.mmm_msc.rna_mmm` — RNA + MMM coupling
- `radia_mcp.magnetic_materials.hysteresis_models.lab_core` — Play model
"""


RNA_MMM_COUPLING = r"""
# RNA + MMM coupling (transformer / reactor)

For transformer / reactor analysis, hybrid RNA + MMM combines:
- **RNA** for the gross magnetic circuit (legs, yoke, gaps)
- **MMM** (full 3D Radia) for local-detail field evaluation around
  windings and air-gap regions

## Workflow

1. Build RNA model:
   ```
   R_leg = l_leg / (μ_iron A_leg)
   R_gap = g / (μ_0 A_gap)
   F = N I (per winding)
   Solve KCL → flux per branch
   ```

2. Use RNA flux as boundary condition for MMM:
   ```
   Local 3D Radia model around critical air gap
   Boundary: Φ from RNA
   Compute local B distribution, force, loss
   ```

3. Couple back to RNA if iron saturates:
   - Use B from MMM to update μ_iron in RNA
   - Iterate until convergence

## Reference (lab)

W:/.../13_RNA/external/:
- A Comparison of Nodal- and Mesh-Based Magnetic Equivalent Circuit
- Dynamic Hysteresis Modeling for Magnetic Circuit Analysis by Inco...
- Topology Optimization of Magnetic Actuator based on Reluctance Network

W:/.../11_BEM_モーメント法/02_mmm_surface_charge/rna_mmm/:
- RNA_MMM_Mixed_Transformer (2 versions)
- RNA_Simplified_MMM_CT (2 versions)
- PEEC_MMM_Coupling

## Cross-references

- `radia_mcp.bem.mmm_msc.rna_mmm` — full detail on the coupling
- `radia_mcp.peec` — PEEC for winding inductance + RNA core
"""




def get_knowledge(topic: str = "overview") -> str:
    """Dispatch RNA / MEC topics.

    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("rna_mmm_coupling", "rna_mmm", "transformer"):
        return RNA_MMM_COUPLING
    if topic == "all":
        return "\n\n".join([OVERVIEW, RNA_MMM_COUPLING])
    return (f"Unknown topic '{topic}'. Available: overview, rna_mmm_coupling, "
            "dynamic_hysteresis, all.")
