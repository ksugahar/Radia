"""
Bibliography of the source PDFs distilled into this MCP server.

Each entry includes: file path, language, length (pages), and which
topics of this server it informs.
"""

BIBLIOGRAPHY = """
# Source bibliography for radia_mcp.differential_forms

The knowledge in this MCP server was distilled from the following six
canonical sources.  Each entry includes the original file path, language,
length, and a one-line orientation.

## 1. A. Bossavit, "Computational Electromagnetism"
(Variational Formulations, Complementarity, Edge Elements),
Academic Press, 1998, ISBN 0-12-118710-1.

Path:
  public-safe curated corpus
  Computational Electromagnetism\\
  {01Preface,02Contents,1..9,A,B,C}.pdf

Pages: ~353 (chapter-split English PDFs)
Language: English
Relevance: PRIMARY SOURCE. Defines Whitney elements, Maxwell's house,
  weak operators, tree-cotree, complementarity.  Chapter 5 (Whitney
  Elements) is the most-cited single chapter in computational EM.

Informs: basics, whitney, de_rham, homology, maxwell (Ch.1).

## 2. D.N. Arnold, R.S. Falk, R. Winther, "Finite element exterior
calculus, homological techniques, and applications", Acta Numerica
15 (2006), pp. 1-155. doi:10.1017/S0962492906210018.

Path:
  public-safe curated corpus
  コホモロジー\\Finite element exterior calculus, homological
  techniques, and applications.pdf

Pages: 155
Language: English
Relevance: THE modern unifying reference.  Introduces  P_r Λ^k  and
  P_r^- Λ^k  families, characterizes them via Koszul complex, proves
  stability of mixed FEM via cochain projections.

Informs: feec, basics (§2 exterior algebra), de_rham, whitney
  (via §5 finite-element differential forms).

## 3. 微分形式.pdf  (lecture notes / monograph, Japanese)

Path:
  public-safe curated corpus
  コホモロジー\\微分形式.pdf
  (Moved 2026-05-20 from public-safe curated corpus to
   consolidate cohomology / forms literature in one place.  Sibling
   papers in the same folder: 微分形式とその外微分.pdf, 微分形式とその積分.pdf,
   微分幾何講義.pdf, 位相解析の基礎.pdf.)

Pages: 13
Language: Japanese (mathematical)
Relevance: Self-contained Japanese-language introduction.  Builds
  differential forms from tangent vectors and tensor products with
  Dirac bra-ket notation.  §6-7 derive Maxwell's equations and the
  Hodge operator in Minkowski space.

Informs: basics (entire), maxwell (§6 Maxwell-as-forms, §7 Hodge).

## 4. A. Kameari (亀有 昭久),
"辺要素を巡って — 「張力」って？—",
SA-11-013 / RM-11-013, 電気学会静止器・回転機合同研究会, 2011.

Path:
  public-safe curated corpus
  コホモロジー\\(2011) 亀有、辺要素を巡ってー「張力」って？－、SA-11-013,
  RM-11-013.pdf
  (Moved 2026-05-20 from .../Bossavit/ to consolidate cohomology
   literature; see コホモロジー folder for sibling papers.)

Pages: 6
Language: Japanese
Relevance: Historical retrospective.  Recounts how edge elements
  emerged from COMPUMAG 1987 (Albanese-Rubinacci, Bossavit),
  Whitney's 1957 formula  φ_i dφ_j − φ_j dφ_i, Kameari's
  symmetric 2nd-order element, the Tonti diagram (fig.6 with primal
  and dual mesh).

Informs: whitney (history + 2nd-order), maxwell (twisted vs
  straight, Kameari diagram), de_rham (overview).

## 5. H. Whitney, "Geometric Integration Theory",
Princeton University Press, 1957;
republished by Dover Publications, 2005.

(Not directly read for this server, but referenced via Bossavit Ch.5
and Kameari 2011.  The Whitney form formula in IV.27 is the
historical origin of the W^k basis.)

Informs: whitney (historical foundation), homology (Whitney
discretization).

## 6. D.J. Cartano, "Exterior Calculus: Wedge Products, Manifolds
and Differential Forms", Internet Publications Ltd, New York, 2023.
ISBN 9798398260243.

Path:
  public-safe curated corpus
  コホモロジー\\Exterior Calculs.pdf

Pages: 474
Language: English
Relevance: Self-published exhaustive textbook on exterior algebra,
  wedge products, differentiable manifolds, and differential forms.
  Scanned PDF — OCR text extraction is unreliable (sample: "VNedge
  Products", "DavidJ.Cartano"), so this reference was NOT directly
  distilled.  Listed as a beginner-friendly general textbook for
  readers who want a longer, slower exposition than 微分形式.pdf
  provides.

Informs: bibliography only (referenced as background textbook).

## 7. 電磁気学とベクトル解析  (吉田 善章 著, 谷島賢二 編),
共立出版, "数学と物理の交差点 / Crossroads of Mathematics and
Physics 2".

Path:
  public-safe curated corpus
  洋書・専門\\電磁気学とベクトル解析.pdf

Pages: 315
Language: Japanese
Relevance: Background reference for the vector-calculus-to-
  differential-forms transition.  NOT directly distilled (the PDF
  is a scanned image — OCR was unreliable).  Listed for completeness.

Informs: bibliography only (referenced as background).

## 8. 五十嵐 一, 亀有 昭久, 加川 幸雄, 西口 磯春, A. Bossavit,
"新しい計算電磁気学 [基礎と数理]", 培風館, 2003.

Path:
  public-safe curated corpus
  新しい計算電磁気学.pdf

Pages: 307
Language: Japanese
Relevance: HIGH.  Chapter 5 was written specifically for this book
  by Dr. Bossavit himself and contains post-1998 (post-Computational-
  Electromagnetism-book) research advances:
  - 5.1 Affine space, piecewise smooth manifolds, inner/outer
        orientation, chains, boundary operator, metric
  - 5.2 Maxwell in differential-form language: integration,
        Stokes, magnetic field as 2-form, Faraday + Ampère, Hodge
  - 5.3 Discretization: primal/dual mesh, "discrete kit"
  - 5.4 FEM: consistency, stability, Whitney forms, HIGHER-ORDER
        Whitney forms (5.4.5), and **Whitney forms for non-simplices
        (5.4.6)** — extends Whitney to hex/prism/polyhedra

  Chapters 1-4 (by Igarashi, Kameari, Kagawa, Nishiguchi) cover
  Japanese-language exposition of edge-element FEM, gauge fixing,
  equivalent circuit representations, electro-mechanical coupling,
  and Section 4.5 on differential forms (push-forward, pull-back,
  geometric meaning, Maxwell in form language).

  Source PDF is EasyOCR'd, text is largely readable (minor kanji
  confusions like 媒→嫌, 電→雷 in some places; section numbers and
  equations preserved).

Informs: whitney (higher-order, non-simplex), de_rham (orientation
  + chains), maxwell (Japanese expression of forms).

## 9. A. Bossavit, L. Kettunen,
"Yee-Like Schemes on Staggered Cellular Grids: A Synthesis Between
FIT and FEM Approaches", IEEE Trans. Magn. 36 (4), 2000, pp. 861-867.

Path:
  public-safe curated corpus
  コホモロジー\\Yee Like Schemes on Staggered Cellular Grids A
  Synthesis Between DIT and FEM Approaches.pdf

Pages: 7
Language: English
Relevance: Convergence analysis for staggered cellular Yee-like
  schemes via Lax theorem (consistency + stability → convergence).
  Distinguishes "network equations" (topological, exact, e.g.
  Faraday on each face) from "network constitutive laws" (metric,
  approximate, e.g. discrete Hodge ⋆_h).  Shows the two main
  constructions of ⋆_h:
    (a) orthogonal primal-dual grid + diagonal-Hodge (FDTD / FIT)
    (b) simplicial mesh + Galerkin mass matrix (Whitney FEM)
  and PROVES both converge.  Companion to Bossavit Ch.5; together
  they form the canonical "differential forms in computational EM"
  reference set.  Distilled into:
    - de_rham (consistency error formulation)
    - whitney (mass matrix Hodge = Galerkin-Hodge)
    - feec    (Lax theorem framework)

## 10. L. Codecasa, R. Specogna, F. Trevisan,
"A new set of basis functions for the discrete geometric approach",
Journal of Computational Physics 229 (2010), pp. 7401-7410.

Path:
  public-safe curated corpus
  コホモロジー\\A new set of basis functions for the discrete
  geometric approach.pdf

Pages: 10
Language: English
Relevance: MODERN extension of Whitney to general polyhedra.
  Introduces 4 new vector basis function families (edge/face,
  primal/dual) on polyhedral primal grids subdivided barycentrically.
  Basis functions are **piecewise uniform** on subregions (NOT
  polynomial like Whitney), but satisfy:
    Property 1: ∫_{r_j} v^r_i · dr = δ_{ij}   (duality)
    Property 2: exact uniform-field interpolation
    Property 3: ∫_v v^r_i dv = ~r_i           (consistency)
  Energetic approach gives symmetric, positive-definite, consistent
  constitutive matrices via  M_{ij} = ∫_v v_i · m v_j dv.
  Cited by `mathematica_recipes('hex_dga')` for symbolic recipe.
  Directly relevant to Radia C++ MSC (Magnetic Surface Charge)
  hexahedral mass-matrix construction.

Informs: whitney ("Whitney for non-simplices"),
  feec (energetic method, alternative to Galerkin),
  mathematica_recipes (hex_dga).

## 11. F. Henrotte (with A. Bossavit's framework),
"Handbook for the computation of electromagnetic forces in a
continuous medium", ICS Newsletter (post-2004), 7 pp.

Path:
  public-safe curated corpus
  Bossavit\\Handbook for the computation of electromagnetic forces
  in a continuous medium.pdf

Pages: 7
Language: English
Relevance: THE unified differential-form treatment of EM forces.
  In 7 pages, derives the Maxwell stress tensor σ_em = b⊗h̃ -
  {h̃·b - ρΨ}I from first principles (energy balance with material
  manifold + placement map), shows ALL classical force methods
  (Maxwell stress, eggshell, Arkkio, Coulomb, Kameari's local
  force, Ren-Razek) emerge as special cases.  Cites Kameari 1993
  as Ref [12].  Used as the primary source for the `forces`
  knowledge module.

Informs: forces (entire), maxwell (mechanical coupling).

## 12. A. Bossavit, "Forces inside a magnet", ICS Newsletter
11(1):4-12, March 2004.
Also: A. Bossavit, "Forces in magnetostatics and their computation",
JAP 67(9):5812-5814, 1990;  A. Bossavit, "Edge element computation
of the force field in deformable bodies", IEEE Trans. Mag. 28(2):
1263-1266, 1992;  A. Bossavit, "On local computation of the
electromagnetic force field in deformable bodies", IJAEM 2:333-343,
1992.

Path:
  (Not directly stored in public-safe curated corpus  Cited via Bossavit-Henrotte 2004
  Handbook reference list.)

Language: English
Relevance: Bossavit's own decades-long investigation of the force
  problem, from which the Bossavit-Henrotte unified framework
  emerged.  Forces inside a magnet (2004) is the most accessible
  expository entry.

Informs: forces (background motivation).

## 13. A. Kameari (亀有),
"Local force calculation in 3D FEM with edge elements",
International Journal of Applied Electromagnetics in Materials
(IJAEM), Vol. 3, pp. 231-240, 1993.

Path:
  (Not directly stored in public-safe curated corpus — bibliographic reference only.
  Cited as Ref [12] in Bossavit-Henrotte 2004 Handbook.)

Language: English
Relevance: The canonical "Kameari's nodal force method" for
3D edge-element FEM.  Generalizes Coulomb's nodal-force formula
(which uses nodal shape functions) to the case where the magnetic
field is discretized on edges (Whitney 1-forms).  Method (7) in
the unified table in `forces('force_methods')`.

Used in: Radia / NGSolve / FEMTET / and most modern edge-element
based magnetic FEM codes.

Informs: forces.force_methods (Kameari's local force method).

## 14. A. Kameari, "有限要素法電磁気解析における種々な手法のアイデア
(無限要素、三角形高次辺要素と高周波固有値問題)", SA-97-6, RM-97-65,
IEEJ Static Apparatus & Rotating Machinery joint study group,
August 1997.

Path:
  public-safe curated corpus
  (1997) 亀有、有限要素法電磁気解析における種々な手法のアイデア,
  SA-97-6_RM-97-65.pdf

Pages: 6
Language: Japanese
Relevance: Kameari proposes (a) infinite-element method for open
  boundary problems, (b) higher-order triangular/tetrahedral edge
  elements including the **"Kikko" (亀甲/tortoise-shell) symmetric
  element** which is invariant under 3-fold (triangular) and
  4-fold (tetrahedral) rotations.  Companion to SA-11-013 in
  Kameari's series of expository papers.  Scanned PDF, OCR
  marginal but figures + key equations recoverable.

Informs: whitney (Kikko symmetric high-order element).

## 15. F. Henrotte, H. Vande Sande, G. Deliége, K. Hameyer,
"Electromagnetic Force Density in a Ferromagnetic Material",
IEEE Trans. Magn. 40(2):553-556, March 2004.

Path:
  public-safe curated corpus
  Electromagnetic Force Density in a Ferromagnetic Material.pdf

Pages: 4
Language: English
Relevance: Derivation of Maxwell stress tensor for 6 different
ferromagnetic material models (linear, PM-as-1-form, PM-as-2-form,
polycrystalline, monocrystal with Weiss domains, magnetostrictive).
Shows that magnetically-equivalent models give DIFFERENT Maxwell
stress tensors — the source of long-standing controversy about
EM forces in iron.  Key reference for `forces('henrotte_ferromagnetic')`.

## 16. TEAM Workshop Problems (force benchmarks)

Path:
  public-safe curated corpus)\\
    problem20_3-D Static Force Problem.pdf           (13 pp)
    problem23_Forces in Permanent Magnets.pdf         (8 pp)
    problem-33b_V3_Experimental Validation of Electric Local
      Force Formulations.pdf                          (6 pp)

Language: English
Relevance: Standard force-computation benchmarks for FEM code
validation.  Problem 20 = saturated steel C-yoke with coil.
Problem 23 = SmCo / NdFeB magnet + coil with measured force data.
Problem 33b = gelatin dielectric in electric field with measured
deformation.  All three are usable as regression tests for any
σ_em implementation.  Distilled into `forces('validation')`.

## 17. T. Iino, Y. Okamoto, A. Ahagon, Y. Kida, K. Semba, T. Yamada
SA-20-003 (magnetostatic) / SA-20-034 (eddy current),
IEEJ Static Apparatus & Rotating Machinery joint study group, 2020.

Path:
  public-safe curated corpus
    SA-20-003_ニ段階誤差補正による有限要素法に基づく電磁力解析の計算
    精度改善に関する検討.pdf  (12 pp)
    SA-20-034_誤差補正による時間領域渦電流解析から得られる電磁力の
    計算精度改善に関する検討.pdf  (6 pp)

Language: Japanese
Relevance: Two-stage error correction technique for FEM-based
electromagnetic force calculation.  Subtracts spurious single-
source contributions to recover only the genuine cross-coupling
force.  Magnetostatic version (SA-20-003) extended to time-domain
eddy current (SA-20-034).  Distilled into `forces('error_correction')`.

## 18. S. Niikura, A. Kameari (亀有),
"Analysis of Eddy Current and Force in Conductors with Motion",
IEEE Trans. Magn. 28(2):1450-1453, March 1992.

Path:
  public-safe curated corpus
  Analysis_of_eddy_current_and_force_in_conductors_with_motion.pdf

Pages: 4
Language: English
Relevance: Kameari's clever moving-source trick: keep conductor
stationary, move source magnetic field via Biot-Savart re-evaluation.
Avoids the unstable v×B advective term.  Two-code cross-validation
(EDDY3D edge-element FEM + EDDYCUFF circuit method).  Distilled
into `forces('kameari_eddy_current')`.

## 19. R. Hiptmair, F. Krämer, J. Ostrowski,
"A Robust Maxwell Formulation for All Frequencies",
IEEE Trans. Magn. 44(6):682-685, June 2008.

Path:
  public-safe curated corpus
  A Robust Maxwell Formulation for All Frequencies.pdf

Pages: 4
Language: English
Relevance: Generating-system approach to fix the low-frequency
instability of the A-V Maxwell formulation.  Introduces redundant
variable q on non-conducting domain with Gauss's law as redundant
equation; system becomes singular but iteratively solvable at all
frequencies including ω = 0.  Distilled into `feec('hiptmair_robust')`.

## 20. P. Gangl, U. Langer, A. Laurain, H. Meftahi, K. Sturm,
"Shape optimization of an electric motor subject to nonlinear
magnetostatics", arXiv:1501.04752 [math.OC], January 2015.

Plus the Gangl PhD thesis:
- "Shape and Topology Optimization Subject to 3D Nonlinear
   Magnetostatics, Part I (Sensitivity Analysis)" — 31 pp
- "Shape and Topology Optimization Subject to 3D Nonlinear
   Magnetostatics, Part II (Numerics)" — 51 pp

Path:
  public-safe curated corpus
  Shape and Topology Optimization subject to 3D Nonlinear
  Magnetostatics Part 1 sensitivity analysis.pdf
  Part II Numerics.pdf
  Shape optimization of an electric motor subject to nonlinear
  magnetostatics.pdf

Language: English
Relevance: Modern rigorous framework for shape optimization with
NONLINEAR PDE constraints.  Sturm's minimax Lagrangian sidesteps
the need to compute material derivatives of the state variable.
IPM motor cogging-torque minimization case study.  Foundation for
the new subpackage `radia_mcp.topology_optimization`.

## Cross-reference

| Topic in this server | Primary source | Supplementary |
|----------------------|----------------|----------------|
| basics               | 微分形式.pdf §1-5 | Arnold-Falk-Winther §2 |
| maxwell              | 微分形式.pdf §6-7 + Bossavit Ch.1 | Kameari fig.6, 新しい計算電磁気学 4.5/5.2 |
| whitney              | Bossavit Ch.5 §5.2 | Kameari §3-5, Whitney 1957 IV.27, 新しい計算電磁気学 5.4 |
| de_rham              | Bossavit Ch.5 §5.1 (Maxwell's house) | Arnold-Falk-Winther §2.2-2.3, Bossavit-Kettunen 2000 |
| homology             | Bossavit Ch.5 §5.3 | Kameari §1 fig.2, 新しい計算電磁気学 5.1.4 |
| feec                 | Arnold-Falk-Winther 2006 §1-7 | Bossavit-Kettunen 2000 (Lax theorem), Codecasa 2010 (DGA energetic) |
| forces               | Bossavit-Henrotte 2004 Handbook | Kameari 1993 local force, Kameari 2011 SA-11-013 shear question, 新しい計算電磁気学 Ch.4 |
| mathematica_recipes  | (composed)         | All — pairs with `radia_mcp.mathematica` |

## Related Radia knowledge servers
- `radia_mcp.radia_ngsolve` — concrete FEM/BEM examples using Whitney
  elements via NGSolve (H1, HCurl, HDiv).
- `radia_mcp.electromagnet` — accelerator magnet workflow using the
  Whitney complex for omega-reduced formulations.
- `radia_mcp.ih` — induction heating with edge / face element coupling.
"""


def get_bibliography_documentation(topic: str = "all") -> str:
    """Return the full bibliography. The ``topic`` argument is
    accepted for API uniformity but has only the value ``"all"``.
    """
    return BIBLIOGRAPHY
