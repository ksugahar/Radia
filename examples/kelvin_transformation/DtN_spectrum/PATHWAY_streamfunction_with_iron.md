# New pathway: stream-function coil design with magnetic material via a Kelvin-FEM material-aware DtN matrix

> **SCOPE (2026-06-15): this is TRACK B — a SEPARATE paper from the DtN+Kelvin core.**
> Decision: keep two distinct tracks so the core is not diluted by the application.
> - **Track A — DtN + Kelvin (core; the SA / Hachinohe paper):** the DtN-spectrum datasheet, the sparse
>   Kelvin open-boundary, the Sommerfeld isomorphism/surrogate, the directly-assembled material-aware DtN
>   matrix and what it IS (FEM-condensed, not BEM). Demos d…dd + x/y/z/aa/cc.
> - **Track B — stream-function coil design with iron (this document; a separate paper):** uses the
>   Track-A operator as the *material-aware design kernel*. Demos ee/ff + a future general-iron design.
> The two share machinery but are written up independently.

*(Consolidation of the stream-function track, 2026-06-15. Honest novelty status at the bottom — a targeted check is running; phrase any claim as "to our knowledge".)*

## The idea (one line)
The **transfer / DtN matrix** that a stream-function (surface-current / current-potential) coil design
inverts — the linear map `psi -> field` — is, **with magnetic material present**, the system's
*material Green operator*. Generate it **sparsely and Green-function-free** as the Schur complement of a
**Kelvin-transformed FEM** (which carries arbitrary `mu(x)` in the inverted exterior). Then coil design
is the same clean linear inverse as in free space, but with the correct **material-aware** kernel.

## Why it was hard (the user's observation: "流れ関数法は磁性体があると楽じゃない")
- Free space: `psi -> field` kernel = **Biot-Savart** (analytic, easy). Design = a clean linear solve.
- With iron (yoke / shield / core): total field = coil field **+ iron reaction**; the kernel becomes the
  **material Green operator**. For planar/cylindrical iron that is the (hard) layered/**Sommerfeld** Green
  function; for **arbitrary** iron there is **no closed-form Green function at all** (a volume integral
  equation revives the dense volume unknown). So the clean "psi x kernel" structure is lost.

## The mechanism (what is new)
1. Mesh **iron + coil surface + open exterior** once; the Kelvin inversion makes the unbounded exterior a
   bounded sparse SPD volume (no Green function, infinity exact).
2. **Condense (Schur)** onto {coil surface, target} -> a small **dense, material-aware transfer/DtN matrix
   M** (the operator BEM would need a Green function to build).
3. **Design = invert M**: `psi = M^+ B_target`. The matrix is the deliverable (this is precisely the
   "operator is the deliverable" case — you do NOT just solve one field; the inverse design consumes M).

## Evidence chain (all committed, verified)
| demo | establishes |
|---|---|
| `demo_t` | FEM-Kelvin carries arbitrary exterior material (layered shell matches analytic ~1e-4) |
| `demo_v` | the material-loaded exterior DtN/Green **matrix** is assembled directly (Schur), spectrum = analytic |
| `demo_w`,`demo_bb` | the matrix for **arbitrary geometry** (cube O_h split; non-layered on-axis inclusion C∞v split) — where **no Sommerfeld Green function exists** |
| `demo_x`,`demo_y`,`demo_z`,`demo_aa` | Kelvin-FEM is the Sommerfeld operator (static isomorphism; multilayer kernel; works DC->wave; low-freq eddy-current = Bannister complex image) |
| `demo_cc` | it is still **FEM** (condensed substructure / SBFEM), not BEM — the Green-function criterion |
| `demo_dd` | **when** to form the matrix: only when the operator is the deliverable (not to solve one field) |
| `demo_ee` | coil + iron shield: free-space kernel off by up to **~16x**; material-aware matches ~1e-4 |
| `demo_ff` | **design inverts M**: with the material-aware M the target is hit (2e-16); the free-space-designed coil misses by **77%** in the iron system |

## To turn into a contribution (next steps)
1. **General (arbitrary iron) coil-design demo**: non-concentric / non-spherical iron -> a genuinely
   coupled dense M; design a real surface-current `psi` (not just modal amplitudes); forward-verify the
   designed coil in an independent full Kelvin-FEM solve (and show free-space design fails).
2. **Benchmark** M-build (sparse Kelvin-FEM Schur) vs the dense layered-Green / FE-BEM baseline
   (conditioning, sparsity, FE-coupling) — the selling point is "sparse, material-aware, no Green fn".
3. **Manuscript**: position as a *formulation* contribution; cite the free-space transformed-FE prior art
   (Brunotte 1992, Meeker 2013), the stream-function/target-field design lineage, and the author's own
   Sugahara 2022 (uniform specimen) as the foundation extended.

## Honest novelty status
- The broad survey found the **transformed-FE-open-boundary** and **stream-function-coil-design**
  literatures **never meet**; the closest "material in the Kelvin exterior" is the author's **own
  Sugahara 2022** (uniform specimen, eddy-current — NOT stream-function design). So the specific
  combination here is **plausibly original to Sugahara**.
- A **targeted novelty check is running** (stream-function x Kelvin/transformed-FE DtN x magnetic
  material). Coverage gaps remain (IEEE/grey-lit/Japanese proceedings; broad-survey confidence ~0.78).
- **Phrasing for the paper: "to our knowledge / we are not aware of ..."** — not a bare "world-first"
  until the targeted check + a Sugahara-2022 forward-citation-graph + 和文 grey-lit scan are done.
