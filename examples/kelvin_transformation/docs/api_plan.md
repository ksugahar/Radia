# Kelvin Transformation Helper API: Design Plan

A staged plan for refactoring the
`examples/kelvin_transformation/A-formulation/` scripts into a reusable
helper API in `src/radia/kelvin_*.py`. Now that the 1-form A and 2-form
B pullbacks have been derived and unit-tested
(`pullback_derivation_3D.md`, `tests/test_kelvin_source.py`), the math
foundation is solid enough to design a general-purpose API on top.

## 0. Goals

- Make the Sugahara two-sphere Kelvin pattern a **one-line geometry
  helper** instead of 50 lines of OCC + periodic-BC bookkeeping.
- Provide a **single source of truth** for the Kelvin nu modulation,
  the A_s pullback, and the L extraction.
- Allow **arbitrary external sources** (Radia ObjArcCur, ObjRecCur,
  ObjRaceTrk, CoilBuilder filaments, ...) via a uniform interface.
- Make the validation harness
  (`validate_radia_HB_kelvin.py`) an exemplar that any new source
  type can use to self-check inductance against a Radia reference.

## 1. Layered architecture

```
+---------------------------------------------------------------+
|  L4  validation harness                                       |
|       compare_against_radia_self_inductance(...)              |
+---------------------------------------------------------------+
|  L3  FEM drivers                                              |
|       solve_full_A_kelvin(geo, J_source_cf, ...)              |
|       solve_reduced_A_kelvin(geo, A_s_factory, ...)           |
+---------------------------------------------------------------+
|  L2  field / material CFs (mesh-aware)                        |
|       make_kelvin_nu_cf(mesh, R_K, offset, ...)               |
|       make_kelvin_aware_A_s_cf(mesh, A_phys_factory, ...)     |
|       make_kelvin_aware_J_source_cf(...) -- if needed         |
+---------------------------------------------------------------+
|  L1  geometry helpers (OCC + periodic BC)                     |
|       add_kelvin_exterior_domain(inner_shape, offset, R, ...) |
|       split_into_air_kelvin_materials(...)                    |
+---------------------------------------------------------------+
|  L0  math primitives  (DONE 2026-04-15)                       |
|       kelvin_map_3d, kelvin_pullback_vector,                  |
|       kelvin_pullback_B_pseudovector,                         |
|       biot_savart_A_at_points, ...                            |
+---------------------------------------------------------------+
```

## 2. Layer-by-layer specification

### L1 -- geometry helpers (`src/radia/kelvin_geometry.py`, new)

```python
def add_kelvin_exterior_domain(
        inner_shape,            # OCC TopoDS_Shape: physical inner part
        offset,                 # (3,) Kelvin exterior sphere center
        R_K,                    # Kelvin sphere radius
        inner_kelvin_face_name="kelvin_int",
        outer_kelvin_face_name="kelvin_ext",
        kelvin_mat="kelvin",
        gnd_vertex_name="GND",
        outer_maxh_factor=2.0):
    """Glue inner_shape with an offset Kelvin sphere, identify the
    inner sphere boundary with the outer Kelvin sphere boundary via a
    periodic BC, attach a GND vertex at the Kelvin center.

    Returns the assembled OCC compound and a small dict with the named
    sub-shapes (for the caller to feed to GenerateMesh)."""
```

Encapsulates the boilerplate currently in
`Coil_3D_A_HCurl_with_Kelvin.py` lines 70-120. A single call replaces
~50 lines of Glue + Identify + GND-vertex logic.

### L2 -- mesh-aware coefficient functions (`src/radia/kelvin_material.py`,
              `src/radia/kelvin_source.py` extended)

```python
def make_kelvin_nu_cf(
        mesh, R_K, offset, nu_0,
        kelvin_mats=("kelvin",),
        convention="sugahara",        # or "energy_invariant" (rho'/R)^2
        ):
    """Return an NGSolve CF nu(r) modulated in Kelvin material.

    convention == 'sugahara':         nu = (R/rho')^2 * nu_0
                  (matches Coil_3D_A_HCurl_with_Kelvin.py code,
                   gives correct L for full-A volume-J FEM)
    convention == 'energy_invariant': nu = (rho'/R)^2 * nu_0
                  (matches the energy-invariance derivation, which
                   uses the validated A/B pullback formulas; this
                   yields different numerical L until the open
                   question in pullback_derivation_3D.md sec 8 is
                   resolved.)
    """

def make_kelvin_aware_A_s_cf(
        mesh, A_phys_factory, R_K, offset,
        kelvin_mats=("kelvin",)):
    """A_phys_factory(x_cf, y_cf, z_cf) -> vector CF for A in PHYSICAL
    coordinates. This helper returns a vector CF that evaluates:

      - in non-Kelvin material: A_phys_factory(x, y, z)
      - in 'kelvin' material:   A_phys_factory(kel_x, kel_y, kel_z)
                                 then 1-form pullback applied:
                                   A_comp = (R/rho')^2 * H * A_phys

    The pullback is the formula derived and unit-tested 2026-04-15
    (kelvin_pullback_vector). The Householder reflection is included.

    Per-component MaterialCF switching is built-in.
    """
```

### L3 -- FEM drivers (`src/radia/kelvin_solver.py`, new)

```python
def solve_full_A_kelvin(
        mesh, J_source_cf,
        R_K, offset, nu_0=NU_0,
        order=1, dirichlet_bbnd="GND",
        gauge_eps=1e-8,
        convention="sugahara"):
    """Full-A 3D HCurl FEM with Sugahara Kelvin convention.

    Solves curl(nu_kelvin curl A) = J in computational coords. Returns
    (gfu, fes, nu_cf) for downstream postprocessing.
    """

def solve_reduced_A_kelvin(
        mesh, A_s_cf,
        R_K, offset, nu_0=NU_0,
        order=1, dirichlet_bbnd="GND",
        gauge_eps=1e-8,
        convention="sugahara"):
    """Reduced-A 3D HCurl FEM with external A_s source.

    Source term: -int (nu_kelvin - nu_0) curl(A_s_cf) . curl(v) dV
                  over Kelvin material only (zero elsewhere).
    """
```

### L4 -- validation harness (`src/radia/kelvin_validate.py`, new)

```python
def compare_against_radia_self_inductance(
        radia_obj_handle,                 # rad.ObjArcCur / RecCur / ...
        coil_label="coil",                # mesh material containing coil
        kelvin_offset=(0.15, 0, 0),
        kelvin_R=0.06,
        outer_maxh=12e-3,
        coil_maxh=3e-3,
        order=1,
        I_total=1.0,
        verbose=True):
    """Build inner sphere + Kelvin exterior, FEM solve volume-J, report

      L_FEM      (from int 0.5 nu |curl A|^2 dV)
      L_Radia    (from rad J.A volume integral via fine quadrature)
      diff       (rel error)

    Acts as the canonical end-to-end validation that the Kelvin
    convention + FEM setup reproduces a Radia analytical reference.
    """
```

## 3. Key design decisions

### 3.1 nu convention choice

Until the open question (pullback_derivation_3D.md sec 8) is resolved,
expose `convention="sugahara"` (matches working code) as default with
`"energy_invariant"` as opt-in for derivation-driven research.

### 3.2 Material naming convention

Standardize on the substrings:
- `"air"` for inner physical air
- `"coil"` for current-carrying region (volume-J source case)
- `"kelvin"` for the Kelvin exterior domain
- `"GND"` for the Dirichlet vertex at infinity image

Substring-match (case-insensitive) so `"kelvin_outer"`,
`"air_inner"`, etc. also work.

### 3.3 1-form pullback default

Always use the full Phase 2 pullback (with Householder reflection),
not the scalar approximation. Empirically the scalar form gave a
debugging detour 2026-04-15.

### 3.4 Compatibility with existing examples

The new helpers MUST allow rewriting:
- `Coil_3D_A_HCurl_with_Kelvin.py` (full-A volume-J baseline)
- `validate_radia_HB_kelvin.py` (reduced-A external A_s)
- `Coil_3D_A_HCurl_PEEC_source.py` (filament A_s when filament dev
  matures)

into <50 lines each, using the layered API.

## 4. Milestones

### M1 -- L1 + L2 helpers (~ 1 day)
- `src/radia/kelvin_geometry.py` with `add_kelvin_exterior_domain`
- `src/radia/kelvin_material.py` with `make_kelvin_nu_cf`
- Extend `src/radia/kelvin_source.py` with
  `make_kelvin_aware_A_s_cf`
- Smoke test: rebuild `Coil_3D_A_HCurl_with_Kelvin.py` baseline
  geometry via the helpers, verify mesh statistics match.

### M2 -- L3 FEM drivers (~ 1 day)
- `src/radia/kelvin_solver.py` with the two solve functions
- Refactor `Coil_3D_A_HCurl_with_Kelvin.py` to use the driver:
  call must reproduce L = 89.44 nH within 0.1%.
- Refactor `validate_radia_HB_kelvin.py` similarly.

### M3 -- L4 validation harness (~ 1 day)
- `src/radia/kelvin_validate.py` with
  `compare_against_radia_self_inductance`
- Standalone smoke script that runs the comparison on
  `rad.ObjArcCur` (square section), `rad.ObjRaceTrk`, and
  `rad.ObjRecCur`. Each must achieve a documented FEM/Radia
  agreement target (e.g. <2% on the medium mesh).

### M4 -- documentation + porting (~ 1-2 days)
- Add a tutorial `examples/kelvin_transformation/docs/quickstart.md`
- Port one IH example to use the new helpers end-to-end
- Update MCP `kelvin_knowledge.py` with the new API surface

## 5. Open questions to resolve before finalizing

1. **nu_kelvin convention**: pure-math derivation says
   `(rho'/R)^2 * nu_0`; production code uses `(R/rho')^2 * nu_0`;
   both give correct L for full-A volume-J FEM. Until reconciled,
   the API exposes both as `convention=...`.
2. **2D axisymmetric variant**: the H-formulation
   `Coil_A_formulation_with_Kelvin.py` uses Z-offset Kelvin
   (different topology). Should the helpers handle 2D too? Probably
   yes, with a `dim=2` switch or a separate `kelvin_geometry_2d.py`.
3. **Multiple Kelvin centers**: for very large coils where one
   Kelvin sphere is too small, would multiple stitched Kelvin
   spheres be useful? Out of scope for M1-M4.

## 6. Order of operations (suggested)

1. Resolve open question #1 above by re-deriving with paper [6] in
   hand, OR running a controlled numerical experiment that pins down
   which `nu_kelvin` is actually consistent with the validated
   pullback formulas.
2. Implement M1 (geometry + material).
3. Implement M2 (FEM drivers), refactor the existing two scripts as
   the canonical regression tests.
4. Implement M3 (validation harness), use it in CI.
5. M4 docs + porting.

The math foundation (L0) is already in place and unit-tested.
