---
name: cohomology-cuts
description: >-
  Compute topological "cuts" (cohomology / homology generators) for a
  multiply-connected magnetostatic domain, then use them in NGSolve. Use this
  whenever a magnetic SCALAR potential (total or reduced phi / T-Omega) needs
  cuts so phi stays single-valued around a coil hole, when you need the
  net-current loops / coil topology (Betti number b1, e.g. torus b1=2) for
  stream-function coil design, or when a question mentions cohomology cuts,
  thick/thin cuts, Ampere loops around conductors, or "h_k" curl-free fields.
  PROJECT DIRECTION is a GMSH-FREE path (native FEEC de Rham cohomology +
  analytic generators on parametric winding surfaces); the Gmsh homology
  solver (Pellikka 2013) is kept only as an independent reference/validator
  (b1 confirmation, cut cross-check). Mesh-source agnostic (Netgen tet / Cubit
  hex), in-memory transfer (no .msh file). Wraps src/radia/cohomology_cut.py
  and radia_mcp.streamfunction.
---

# Cohomology cuts for magnetostatics (NGSolve)

The "cuts" = the **cohomology generators** `h_k` of a multiply-connected region:
curl-free `H(curl)` fields with unit circulation around each coil hole. They let a
magnetic SCALAR potential carry coil current — `H = -grad(phi) + sum_k I_k h_k` — with
no Biot-Savart source field, no manual cut surface, no Gibbs artifacts at the coil.

> **PROJECT DIRECTION (user 2026-06-12): compute the cuts GMSH-FREE.** The production
> path is native — analytic generators on parametric winding surfaces, and the discrete
> FEEC / de Rham cohomology on a general NGSolve mesh (the feec_vim direction). The Gmsh
> homology solver (Pellikka 2013, later sections) is retained **only as an independent
> reference / validator** — confirm `b1`, cross-check the cut representatives — not the
> production engine.

The current engine `src/radia/cohomology_cut.py` (`CohomologyCutSolver`) still calls Gmsh;
deep theory is in `radia_mcp.streamfunction` (`streamfunction_knowledge.py`). This skill
is the discoverable recipe.

## Why cuts (the topology)

A magnetic *scalar* potential `phi` with `H = -grad(phi)` only exists where the
domain is curl-free. Around a current-carrying hole (a coil), `oint H.dl = NI != 0`,
so `phi` is **multivalued** — you must either cut the domain or add topological
basis fields. The total-scalar-potential formulation does the latter:

    H = -grad(phi) + sum_k I_k h_k ,   I_k = N_k * I_coil  (ampere-turns)

where each `h_k` is a **cohomology generator**: a curl-free `H(curl)` field with

    curl(h_k) = 0            (irrotational)
    oint_{loop_j} h_k . dl = delta_jk   (unit circulation around independent loop j)

The number of generators = the first Betti number `b1` of the domain = number of
independent coils / holes. (A torus surface has `b1 = 2`: one toroidal + one
poloidal loop — the two stream-function net-current generators.)

## Preferred path: GMSH-FREE generators

1. **Analytic generators (parametric surfaces).** When the winding surface has a known
   parametrization (torus, helical band, VMEC plasma boundary), the H1 generators are
   closed-form: on a torus, `grad(zeta)` (net-toroidal) and `grad(theta)` (net-poloidal)
   — the two net-current loops, `b1 = 2`. No topology computation, no Gmsh. This is the
   production path for stream-function coil design (see `radia_mcp.streamfunction`).

2. **Native discrete de Rham cohomology (general meshes).** The 1-cohomology is
   `ker(curl) / im(grad)` in the discrete `H(curl)`. Build the discrete operators in
   NGSolve — `grad: H1 -> HCurl` and `curl: HCurl -> HDiv` as sparse matrices — then
   `dim H^1 = b1 = nullity(curl) - rank(grad)`. Get representatives either as harmonic
   1-forms (lowest modes of the Hodge-Laplacian restricted to `ker(curl)` orthogonal to
   `im(grad)`) or by a spanning-tree / cotree edge-cut construction on the mesh edge graph
   (Kotiuga / Webb-Forghani — pure combinatorics, no Gmsh). This is the feec_vim direction
   and runs natively on any Netgen tet or Cubit hex mesh.

3. **Validate against Gmsh (oracle).** Use the Gmsh homology solver (below) as the
   independent reference: it must return the same `b1` and homologous cut representatives.
   Keep it as a regression check, not in the production path.

## Gmsh homology solver — reference / validator (NOT the production path)

Gmsh's homology solver finds representative chains of the (relative) (co)homology
basis and stores **each generator as a physical group** in the mesh.

    import gmsh
    # domain = 3D physical group of the FIELD region (air+iron), EXCLUDING the coil hole
    # subdomain (relative-to) = list of 2D physical groups, or [] for absolute
    gmsh.model.mesh.addHomologyRequest("Cohomology", [domain_tag], [], [1])  # thick cuts, H^1
    # ... (queue more requests if needed) ...
    dimTags = gmsh.model.mesh.computeHomology()   # runs all queued requests

- `addHomologyRequest(type, domainTags, subdomainTags, dims)`:
  - `type`: `"Homology"` or `"Cohomology"`.
  - `domainTags` / `subdomainTags`: **physical-group tags** (relative (co)homology of
    domain modulo subdomain). `[]` subdomain = absolute.
  - `dims`: list of chain dimensions, e.g. `[1]` for the 1-(co)homology = the loop cuts
    (use `[0,1,2,3]` to get the full Betti sequence; `dims=[1]` is what coils need).
- `computeHomology()` returns the dimTags of the created physical groups (one per
  generator). Read the chains back with
  `gmsh.model.getEntitiesForPhysicalGroup(1, tag)` + `gmsh.model.mesh.getElements(1, ent)`.
- **Thin vs thick cuts**: a `"Homology"` request modulo the non-terminal boundary gives
  *thin cuts* (1-chains / curves); a `"Cohomology"` request modulo the terminals gives
  *thick cuts* (the curl-free cochain you integrate against — this is what the scalar
  potential wants). For the standard coil problem use **`"Cohomology", dims=[1]`**.
- A *simply-connected* domain yields **0 generators** — `computeHomology()` returns
  empty / the solver raises. That is the signal you forgot to exclude the coil hole.

## Mesh transfer — in-memory via the API, NO `.msh` file

**Do not round-trip through a `.msh`/`.vol` file.** A file transfer discards node
identity, which forces a fragile **KDTree re-match (1e-8 tol)** of Gmsh nodes to NGSolve
vertices and pins you to `Mesh.MshFileVersion 2.2`. Move the mesh **in memory through the
Gmsh Python API**, recording the exact `gmsh_node_tag -> point` map as you build it — then
the cohomology edge chains map to NGSolve edges with **zero tolerance** (the map is the
identity, no KDTree). Two file-free directions:

1. **Gmsh meshed it -> build the Netgen mesh in memory.** Read Gmsh's mesh with
   `gmsh.model.mesh.getNodes()` / `getElements(dim)` and add to a `netgen.meshing.Mesh`,
   capturing `g2n = gmsh_tag -> PointId` while adding the points:
   ```python
   from netgen.meshing import Mesh as NMesh, MeshPoint, Pnt, Element3D
   import ngsolve
   nm = NMesh(dim=3); nm.SetMaterial(1, "domain")
   ntags, xyz, _ = gmsh.model.mesh.getNodes()
   g2n = {int(t): nm.Add(MeshPoint(Pnt(xyz[3*i], xyz[3*i+1], xyz[3*i+2])))
          for i, t in enumerate(ntags)}
   etypes, _, conn = gmsh.model.mesh.getElements(3)          # type 4 = 4-node tet
   for k, et in enumerate(etypes):
       if et != 4: continue
       c = conn[k]
       for j in range(len(c) // 4):
           nm.Add(Element3D(1, [g2n[int(c[4*j+m])] for m in range(4)]))
   # + Element2D/FaceDescriptor per 2D physical group for the boundary BCs
   mesh = ngsolve.Mesh(nm)
   ```
   `g2n` is exact, so mapping the cut edges to `HCurl(order=0)` DOFs needs no KDTree.

2. **You already have a Netgen/Cubit mesh -> push it into Gmsh (no file).** Inject the
   mesh through the API instead of `gmsh.open(...)` on a written `.msh`:
   `gmsh.model.addDiscreteEntity(3, tag)`,
   `gmsh.model.mesh.addNodes(3, tag, nodeTags, coords)`,
   `gmsh.model.mesh.addElements(3, tag, [4], [elemTags], [tetNodeTags])`,
   `gmsh.model.addPhysicalGroup(3, [tag], pg)`, then `addHomologyRequest(...)` +
   `computeHomology()`. The cuts return in **your own node numbering** (the tags you
   pushed) -> the map to NGSolve DOFs is the identity. Cubit hex: push 8-node hexes
   (Gmsh element type 5) — cohomology works on cell complexes, not just simplices.

> Watch element orientation: Gmsh and Netgen can disagree on tet node order — check for
> inverted elements (negative volume) after import and swap two nodes if needed.

> Lab policy (gmsh server): Gmsh is the lightweight **topology / post-processing** engine;
> **mesh generation is Netgen or Cubit**. Cohomology is a topology post-step —
> mesh-source-agnostic, and the mesh should move **by API, not by file**.

## Use it via the existing solver (preferred)

    import gmsh
    from radia.cohomology_cut import CohomologyCutSolver

    gmsh.initialize()
    # ... build geometry: 3D physical group "domain" (air+iron, NOT the coil hole),
    #     2D physical group "outer" for the Dirichlet boundary ...
    solver = CohomologyCutSolver()
    n_coils = solver.setup_from_gmsh(domain_physical_name="domain",
                                     boundary_physical_name="outer")
    # n_coils == b1 == number of cohomology generators (assert it matches your coils)
    solver.solve([100.0], {"iron": 1000.0})       # NI per generator, mu_r per material
    B = solver.get_B()                            # NGSolve CF (Tesla, dim 3)
    B_hdiv = solver.project_to_hdiv(order=2)      # div-conforming B
    gmsh.finalize()

- Open boundary: `solver.solve(..., kelvin_region="shell", kelvin_radius=R,
  kelvin_center=(0,0,0))` applies the `(R/r)^2` Kelvin-transform weight in the exterior
  shell (pairs with the radia Kelvin-transformation knowledge).
- Nonlinear iron: `solver.solve_nonlinear(NI_list, bh_data, iron_domain="iron")`.
- The generators are mapped to **`HCurl(mesh, order=0)`** edge DOFs (`h_gf.vec[edge_nr]
  = +/-1` along the chain); `get_cohomology_basis()` returns the list of `h_k`
  GridFunctions. NOTE: the current `_transfer_mesh_to_ngsolve` / `_build_vertex_map`
  still uses the **`.msh`->`.vol`->KDTree** path — replace it with the in-memory transfer
  above (exact `g2n` map) to drop the file + the tolerance.

## Stream-function coil design link

For a winding *surface* (not a volume), the same `b1` counts the surface cohomology
generators = the **net-current** degrees of freedom. On a torus `b1 = 2`: the two
generators are the net-toroidal and net-poloidal currents (`grad(zeta)`, `grad(theta)`
are the analytic generators on the standard torus; a general shaped surface uses the
Gmsh cohomology engine). A single-valued stream function `psi` carries **zero** net
current through each hole — the net current is the *secular / multivalued* term, the
cohomology part. See `radia_mcp.streamfunction` (`addHomologyRequest('Cohomology',
[pg], [], [1]); computeHomology() -> 2` confirms `b1==2`).

## Pitfalls

- **Exclude the conductor.** The domain physical group must be the field region with
  the coil *hole removed*, or it is simply-connected and you get 0 cuts.
- `dims=[1]` for coil loops; using `[0,1,2,3]` also returns 0-, 2-, 3-generators you
  usually do not want as sources.
- **Cohomology (thick cuts), not Homology**, for the scalar-potential source term.
- **No `.msh` round-trip.** Move the mesh in memory (see above); the file path forces
  `MshFileVersion 2.2` + a KDTree node re-match and is the main fragility to avoid.
- Generators are lowest-order (`HCurl` order 0) edge cochains — the circulation is
  exactly `+/-1` per edge along the cut; do not "Set" them from a CF.
- `n_coils` returned by `setup_from_gmsh` must equal your physical coil count; a
  mismatch means a hole is missing or a terminal group is mis-tagged.

## References (papers / open tools only)

- M. Pellikka, S. Suuriniemi, L. Kettunen, C. Geuzaine, "Homology and cohomology
  computation in finite element modeling," *SIAM J. Sci. Comput.* **35**(5):1195-1214 (2013).
  (The algorithm Gmsh's `computeHomology` implements.)
- P. R. Kotiuga, theory of cuts for the magnetic scalar potential (multiply-connected
  magnetoquasistatics).
- Z. Ren (2002), T-Omega dual / thick-cut formulation.
- Gmsh tutorial **t14** ("Homology and cohomology computation").
