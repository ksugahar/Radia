"""
Differential-Forms MCP Server (radia_mcp.differential_forms)

Provides knowledge tools for differential geometry in computational
electromagnetism:
- Basic differential-form machinery (tangent spaces, k-forms,
  wedge product, exterior derivative, Stokes theorem)
- Maxwell's equations in differential-form language (Hodge star,
  twisted vs straight forms, Kameari diagram)
- Whitney elements (the discrete realization of the de Rham complex
  on a simplicial mesh)
- de Rham complex and Sobolev spaces of forms (H^1, H(curl), H(div),
  Bossavit's "Maxwell's house")
- Homology (chains, boundary operator, Betti numbers, tree-cotree
  decomposition, Euler-Poincaré)
- Finite Element Exterior Calculus (Arnold-Falk-Winther 2006:
  P_r Λ^k and P_r^- Λ^k families, Koszul complex, cochain projections,
  Hodge Laplacian)
- Bibliography of source PDFs

Usage:
    mcp-server-differential-forms              # Start MCP server (stdio)
    mcp-server-differential-forms --selftest   # Run self-test

This server is read-only: every tool returns structured documentation
strings.  No solver, no symbolic computation.  Companion knowledge:
- `radia_mcp.radia_ngsolve` — concrete FEM examples
- `radia_mcp.electromagnet` — accelerator-magnet workflow
- `radia_mcp.ih` — induction heating with edge/face element coupling
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .basics_knowledge import get_basics_documentation
from .maxwell_knowledge import get_maxwell_documentation
from .whitney_knowledge import get_whitney_documentation
from .de_rham_knowledge import get_de_rham_documentation
from .homology_knowledge import get_homology_documentation
from .feec_knowledge import get_feec_documentation
from .bibliography_knowledge import get_bibliography_documentation
from .mathematica_recipes_knowledge import get_mathematica_recipes_documentation
from .forces_knowledge import get_forces_documentation
from .em_force_ngsolve_recipe_knowledge import get_em_force_ngsolve_recipe
from .em_force_extras_knowledge import get_em_force_extras

mcp = FastMCP("mcp-server-differential-forms")


# ============================================================
# MCP Tools
# ============================================================

@mcp.tool()
def differential_forms_basics(topic: str = "all") -> str:
    """
    Basic differential-form machinery.

    Topics:
      "all"                 - Concatenate everything below
      "tangent_vectors"     - Tangent space, dual (cotangent) space,
                              bra-ket notation, gradient as a 1-form
      "wedge_product"       - Tensor product, alternating tensor,
                              wedge product (∧), dimension counting
      "exterior_derivative" - d : Λ^k → Λ^{k+1}, d² = 0, Poincaré lemma,
                              grad/rot/div as the same operator in 3-D
      "stokes"              - Stokes' theorem ∫_σ dω = ∫_∂σ ω and its
                              classical instances (FTC, Stokes/Kelvin,
                              divergence theorem)
    """
    return get_basics_documentation(topic)


@mcp.tool()
def differential_forms_maxwell(topic: str = "all") -> str:
    """
    Maxwell's equations in differential-form language.

    Topics:
      "all"
      "maxwell"   - Each EM quantity (e, h, b, d, j, ρ) as a specific
                    form type, Faraday/Ampère/continuity as form-d
                    equations, 4-D F and G formulation
      "hodge"     - The Hodge star operator ⋆: Λ^k → Λ^{n-k},
                    constitutive relations d=εe, b=μh, j=σe as ⋆
                    operations, role of the metric
      "twisted"   - Inner vs outer orientation, straight vs twisted
                    (density/pseudo) forms, the Kameari/Tonti diagram
                    primal/dual mesh layout
    """
    return get_maxwell_documentation(topic)


@mcp.tool()
def differential_forms_whitney(topic: str = "all") -> str:
    """
    Whitney elements: discrete differential forms on a simplicial mesh.

    Topics:
      "all"
      "overview"   - Four basis families w_n, w_e, w_f, w_T; DoF = integrals
                     over node/edge/face/cell; Whitney complex diagram
      "edge"       - Edge elements (1-forms): formula  w_e = w_m ∇w_n
                     - w_n ∇w_m, tangential continuity, Nedelec hierarchy,
                     tree-cotree gauge
      "face"       - Face elements (2-forms): formula, normal continuity,
                     RWG/Raviart-Thomas connection, mass matrix
      "properties" - DR=0, RG=0 (exact), Betti numbers, mass matrices,
                     why mixed nodal+edge formulations have spurious modes
    """
    return get_whitney_documentation(topic)


@mcp.tool()
def differential_forms_de_rham(topic: str = "all") -> str:
    """
    de Rham complex, Sobolev spaces of forms, Maxwell's house.

    Topics:
      "all"
      "complex"       - The de Rham complex with d² = 0, de Rham cohomology,
                        Betti numbers, why the kernel of rot ≠ image of
                        grad in non-contractible domains
      "weak"          - Weak operators (grad/rot/div), Sobolev spaces
                        H^1 (= L²_grad), H(curl) (= IL²_rot), H(div)
                        (= IL²_div), traces, why FEM needs the weak forms
      "maxwell_house" - Bossavit's 4-pillar diagram: primal/dual columns,
                        constitutive ⋆ as horizontal arrows, time
                        derivative ∂_t as vertical link
    """
    return get_de_rham_documentation(topic)


@mcp.tool()
def differential_forms_homology(topic: str = "all") -> str:
    """
    Chain complex, homology, Betti numbers, tree-cotree gauge.

    Topics:
      "all"
      "chain_complex"   - p-chains, boundary operator ∂, ∂² = 0,
                          cycles vs boundaries, homology groups H_p,
                          de Rham duality with cohomology
      "tree_cotree"     - Spanning tree of edges (= maximal acyclic
                          subgraph), cotree, gauge fixing via setting
                          a_e = 0 on tree, belted-tree for non-
                          contractible domains
      "euler_poincare"  - χ = N − E + F − T = Σ (−1)^p b_p formula,
                          Betti-number table for common topologies
                          (ball, torus, hollow shell), cut-detection
    """
    return get_homology_documentation(topic)


@mcp.tool()
def differential_forms_feec(topic: str = "all") -> str:
    """
    Finite Element Exterior Calculus (Arnold-Falk-Winther 2006).

    Topics:
      "all"
      "overview"            - The two families P_r Λ^k and P_r^- Λ^k,
                              classical equivalents (Lagrange, RT, BDM,
                              Nedelec I/II), Whitney = P_1^- Λ^k
      "cochain_projection"  - Commuting diagram d Π_h = Π_h d, why this
                              gives stable mixed FEM, smoothed projections
                              (Schöberl, Christiansen)
      "koszul"              - Koszul complex (algebraic engine), κ d + d κ
                              = (r+k) id, characterization of P_r^- Λ^k,
                              affine invariance and finite design space
      "hodge_laplacian"     - Abstract Hodge Laplacian L = dδ + δd, mixed
                              variational FEM template, connection to
                              auxiliary-space preconditioners (Compact HX
                              / Hiptmair-Xu) used in Radia
    """
    return get_feec_documentation(topic)


@mcp.tool()
def differential_forms_forces(topic: str = "all") -> str:
    """
    Electromagnetic forces in differential-form language.

    The "question of forces in electromagnetism", long controversial
    in vector-calculus textbooks, becomes a single energy-balance
    identity (Bossavit-Henrotte 2004) where ALL classical force
    formulas (Maxwell stress, eggshell, Arkkio, Coulomb, Kameari's
    nodal force, Lorentz) emerge as special cases of one tensor
    σ_em.

    Topics:
      "all"
      "overview"             - Why forces need differential forms
      "mechanical_framework" - Material manifold M + Euclidean E +
                               placement map p_t : M → E (the conceptual
                               move that makes the theory clean)
      "material_derivative"  - Lie derivative formulas for 0/1/2/3-forms;
                               the algebraic backbone of every method
      "maxwell_stress"       - σ_em = b ⊗ h̃ - {h̃·b - ρΨ}I derivation
                               from energy balance; extends naturally to
                               magnetostrictive / anisotropic / hysteretic
                               materials
      "force_methods"        - Unified table of 7 force-computation
                               methods (force density, Maxwell stress,
                               eggshell, Arkkio, Coulomb-nodal,
                               Ren-Razek, **Kameari's local force
                               method for 3D edge-element FEM**)
                               — all derived from σ_em
      "kameari_stress"       - Kameari SA-11-013 (2011) open question
                               on differential-form representation of
                               shear stress / shear strain
    """
    return get_forces_documentation(topic)


@mcp.tool()
def differential_forms_em_force_recipe(topic: str = "all") -> str:
    """
    PRACTICAL NGSolve recipe for EM force computation (2026-05-22).

    Goes beyond `differential_forms_forces` (which catalogs the 7
    theoretical methods) and prescribes:
        - WHICH method to use in which scenario (decision tree)
        - HOW to set it up in NGSolve for full p-convergence
        - WHAT pitfalls destroy accuracy (image artifacts, BND trace,
          mu jumps, ...)
        - HOW to validate the result (Newton-3, path independence)

    Captures hard-won findings from implementing `calc_em_force.py`.
    Achieved 0.5% Newton-3 violation on iron+coil (mu_r=5000) benchmark
    -- 10x better than Pile-Devillers-Le Besnerais 2018 'production' 5%
    -- via Mertens-Hameyer 2007 v_order matching + far-field outer box.

    Topics:
        "method_choice"         - Decision tree: which method when
        "high_order_recipe"     - Mertens-Hameyer 2007 v_order match
        "common_pitfalls"       - 6 traps (image artifact, BND trace,
                                   mu-jump L2 projection, ...)
        "validation_protocol"   - Newton-3 + path-independence + multi-
                                   method consistency checks
        "implementation_status" - calc_em_force.py status: which methods
                                   are production / legacy / TODO
        "full_example"          - Complete working code template
        "all"                   - Everything

    Cross-reference: `differential_forms_forces` (theory),
    `motor_henrotte_lineage` (Henrotte 1993-2018 papers).
    """
    return get_em_force_ngsolve_recipe(topic)


@mcp.tool()
def differential_forms_em_force_extras(topic: str = "all") -> str:
    """
    Additional EM force topics beyond the 7-method catalog (2026-05-22).

    Complements `differential_forms_forces` (7 classical methods) and
    `differential_forms_em_force_recipe` (NGSolve practical recipe)
    with the following advanced/specialized topics:

    Topics:
        "permanent_magnet_force"  - PM force formulations: bound surface
                                     current K = M x n, Kelvin
                                     (M . grad) B, and equivalent
                                     magnetic scalar charge methods.
        "energy_method_classical" - Constant flux vs constant current,
                                     Legendre duality W <-> W',
                                     when to use which, FE recipes.
        "sensitivity_vwp"         - Proper shape-derivative VWP
                                     (Pironneau-Allaire) -- the
                                     mathematically correct replacement
                                     for IfPos virtual translation.
        "lorentz_canonical"       - Lorentz F = integral J x B dV as a
                                     first-class method (not just VWP
                                     equivalent); decomposable per-
                                     conductor, immune to BND trace.
        "meissner_force"          - Type-I superconductor mu_r -> 0,
                                     image method, levitation, Bean
                                     model TODO note.
        "all"                     - Everything

    Cross-reference: `differential_forms_em_force_recipe` for the
    practical NGSolve implementation status of each method.
    """
    return get_em_force_extras(topic)


@mcp.tool()
def differential_forms_bibliography(topic: str = "all") -> str:
    """
    Bibliography of the source PDFs distilled into this server.

    Returns the full list of references with file paths, languages,
    and which topics each source informs.
    """
    return get_bibliography_documentation(topic)


@mcp.tool()
def differential_forms_mathematica_recipes(topic: str = "all") -> str:
    """
    Wolfram Language recipes for symbolic verification, pairing
    `radia_mcp.differential_forms` (theory) with `radia_mcp.mathematica`
    (Wolfram subprocess bridge).

    Workflow: read theory here, copy the Wolfram code block from the
    recipe, send it to `mathematica_evaluate` via the mathematica MCP
    server, and read the symbolic confirmation.

    Topics:
      "all"
      "dsquared"      - Verify d² = 0  (rot grad = 0, div rot = 0)
      "stokes"        - Verify Stokes' theorem on a tetrahedron face
      "whitney_edge"  - Symbolic construction of  w_e = w_m ∇w_n − w_n ∇w_m
                        and verification of the duality  ∫_{e'} τ·w_e = δ
      "maxwell"       - Verify  ∂_t B + rot E = 0  and  ∇·B = 0  from
                        the potentials A, ϕ (Faraday + Gauss-magnetic)
      "hodge"         - Hodge star ⋆ in R^3 and in Minkowski R^{1,3};
                        symbolic constitutive relations d = ε e, b = μ h
      "weakform_hodge" - Committed `.wls`: weak form as Hodge pairing,
                         pulled-back material/Kelvin weights, transformation
                         optics, metric-vs-curvature checks
      "hodograph"     - Committed `.wls`: Kelvin/Clebsch/Chaplygin 3-axis
                        hodograph backbone and no-fold saturation condition
      "canonical"     - Committed `.wls`: Hamiltonian/Legendre/canonical
                        transform reading of flux lines and the hodograph
      "surface_derham" - Committed `.wls`: HOIBC as surface topology
                         plus analytic exterior DtN/Steklov approximation
      "differential_geometry" - Runbook/index for the committed hodograph
                                differential-geometry `.wls` suite
      "tex"           - TeXForm output for paper writing
      "lorentz"       - Lorentz boost of the EM 2-form F  (E ↔ B mixing)
      "hex_dga"       - Codecasa-Specogna-Trevisan basis functions for
                        polyhedral / hexahedral grids (DGA energetic method)
      "kelvin"        - Kelvin transformation: derive kelvin_factor
                        R/|y| via conformal weight λ^((n-2k)/2);
                        cover k=0 (scalar pot), k=1 (vector pot),
                        k=2 (B field).  Pairs with `radia.axifem`
                        Kelvin shell + `radia_mcp.radia_ngsolve` Kelvin
                        identification workflow.
      "maxwell_stress" - Maxwell stress tensor σ_em = b⊗h̃ - {h̃·b - ρΨ}I:
                        symbolic derivation for linear isotropic, then
                        anisotropic / magnetostrictive extension.  Pairs
                        with `differential_forms_forces` (theory side),
                        gives executable Mathematica for hand-derivation-
                        free σ_em construction.  Useful for NGSolve /
                        Radia C++ codegen.
    """
    return get_mathematica_recipes_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def differential_forms_starter(level: str = "intro") -> str:
    """Suggest a study order through the tools.

    Args:
      level: "intro"        - For users new to differential forms
             "fem"          - For users with FEM background
             "research"     - For users doing PDE / numerical analysis
    """
    if level == "intro":
        return (
            "Suggested learning path through this MCP server:\n"
            "  1. differential_forms_basics('tangent_vectors')\n"
            "  2. differential_forms_basics('wedge_product')\n"
            "  3. differential_forms_basics('exterior_derivative')\n"
            "  4. differential_forms_basics('stokes')\n"
            "  5. differential_forms_maxwell('maxwell')\n"
            "  6. differential_forms_maxwell('hodge')\n"
            "  7. differential_forms_maxwell('twisted')\n"
            "Goal: understand how Maxwell's equations become four lines\n"
            "       (dF = 0, dG = J, F = ⋆G, plus material law).\n"
        )
    if level == "fem":
        return (
            "Suggested learning path for FEM-experienced reader:\n"
            "  1. differential_forms_de_rham('complex')\n"
            "  2. differential_forms_de_rham('weak')\n"
            "  3. differential_forms_de_rham('maxwell_house')\n"
            "  4. differential_forms_whitney('overview')\n"
            "  5. differential_forms_whitney('edge')\n"
            "  6. differential_forms_whitney('face')\n"
            "  7. differential_forms_whitney('properties')\n"
            "  8. differential_forms_homology('tree_cotree')\n"
            "Goal: connect classical edge / RT / BDM elements to the\n"
            "      Bossavit Whitney-complex picture; understand why\n"
            "      mixed FEM needs structure preservation.\n"
        )
    if level == "research":
        return (
            "Suggested deep-dive for PDE / numerical analyst:\n"
            "  1. differential_forms_feec('overview')\n"
            "  2. differential_forms_feec('cochain_projection')\n"
            "  3. differential_forms_feec('koszul')\n"
            "  4. differential_forms_feec('hodge_laplacian')\n"
            "  5. differential_forms_homology('chain_complex')\n"
            "  6. differential_forms_homology('euler_poincare')\n"
            "  7. differential_forms_bibliography()\n"
            "  → Then read Arnold-Falk-Winther 2006 in full (155 p).\n"
        )
    return (
        f"Unknown level '{level}'. Available: 'intro', 'fem', 'research'."
    )


@mcp.prompt()
def verify_with_mathematica(identity: str = "dsquared") -> str:
    """Workflow to symbolically verify a differential-form identity
    using the sibling `radia_mcp.mathematica` server.

    The differential_forms server provides THEORY; the mathematica
    server provides COMPUTATION (Wolfram Language subprocess bridge).
    Together they let an AI verify identities, not just describe them.

    Args:
      identity:  which identity to verify.  Values:
        "dsquared"      - d² = 0  (rot grad = 0, div rot = 0)
        "stokes"        - Stokes on a tetrahedron face
        "whitney_edge"  - the Whitney edge element formula
        "maxwell"       - dF = 0  (Faraday + ∇·B = 0)
        "hodge"         - ⋆⋆ = ±id, constitutive ⋆-laws
        "weakform_hodge" - pulled-back Hodge/material weight
        "hodograph"     - Kelvin/Clebsch/Chaplygin symbolic backbone
        "canonical"     - Hamiltonian/Legendre symbolic backbone
        "surface_derham" - HOIBC surface de Rham + analytic DtN
        "lorentz"       - Lorentz covariance of F
    """
    return (
        f"To verify '{identity}' symbolically:\n"
        "  1. THEORY: call `differential_forms_mathematica_recipes("
        f"'{identity}')` — read the Wolfram Language recipe and the\n"
        "     theoretical context.\n"
        "  2. COMPUTATION: copy the Wolfram code block from the recipe\n"
        "     and send it to `mathematica_evaluate` (from the\n"
        "     `radia_mcp.mathematica` server) via:\n"
        "       mathematica_evaluate(code=<the recipe code>, timeout=60)\n"
        "  3. INTERPRET: parse the returned 'result' field.  For\n"
        "     identities, expected values are {0,0,0} or 0 / True.\n"
        "  4. PAPER WRITING: for LaTeX output, wrap with\n"
        "     `mathematica_to_tex(<expression>)` instead.\n"
        "\n"
        "Cross-reference:\n"
        "  - `radia_mcp.differential_forms` server: theory + recipes\n"
        "  - `radia_mcp.mathematica` server: 10 tools "
        "(evaluate / status / simplify / to_tex / check_identity / "
        "vector_calc / unit_convert / solve / integrate / differentiate)\n"
        "  - Companion notes: `notes_fem_hcurl.md` inside the\n"
        "    mathematica server covers Zaglmayr H(curl) hierarchical\n"
        "    bases — directly relevant when extending Whitney elements\n"
        "    to higher polynomial order p > 1.\n"
    )


# ============================================================
# Entry point
# ============================================================



register_status_tool(
    mcp,
    server_name='mcp-server-differential-forms',
    description='Differential forms / exterior calculus for EM: de Rham complex, cohomology, EM forces theory',
    subpackage='radia_mcp.differential_forms',
    related_servers=["fem", "mathematica"],
)


def main():
    if "--selftest" in sys.argv:
        print("Differential-forms MCP server self-test:")
        tools = [
            ("basics", differential_forms_basics),
            ("maxwell", differential_forms_maxwell),
            ("whitney", differential_forms_whitney),
            ("de_rham", differential_forms_de_rham),
            ("homology", differential_forms_homology),
            ("feec", differential_forms_feec),
            ("forces", differential_forms_forces),
            ("bibliography", differential_forms_bibliography),
            ("mathematica_recipes", differential_forms_mathematica_recipes),
        ]
        for name, fn in tools:
            n_chars = len(fn("all"))
            print(f"  {name}('all'): {n_chars} chars")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
