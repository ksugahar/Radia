"""MCP Server: radia_mcp.comsol_converter -- COMSOL <-> IR <-> NGSolve.

mcp-server-comsol-converter (tool prefix ``cc_``). Hub-and-spoke model
translation for low-frequency magnetics (magnetostatic A-formulation, v1).
Reuses interop (COMSOL Java API + lab tips), fem/bem (NGSolve formulation),
team_benchmark (reference solutions); coordinates the EXTERNAL COMSOL MCP
(mph_reader / Java gen-run) and radia-ngsolve for the actual solve.

Realistic stance: these tools provide the IR schema, the COMSOL<->IR<->NGSolve
mapping knowledge, and code-generation templates. The model-specific
translation is performed by the calling LLM using that scaffolding -- full
automatic translation of arbitrary multiphysics is out of scope. v1 hardens
the magnetostatic lab core (TEAM 6/13/20) and validates via team_benchmark.

Usage:
    mcp-server-comsol-converter            # stdio
    mcp-server-comsol-converter --selftest # self-test
"""
import json
import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from . import ir as _ir

mcp = FastMCP("mcp-server-comsol-converter")


# ============================================================
# IR schema
# ============================================================

@mcp.tool()
def cc_ir_schema(example: str = "team20") -> str:
    """Neutral IR schema (COMSOL <-> IR <-> NGSolve hub) + a worked example.

    The IR carries physics intent, not solver syntax. Build/parse it on
    either side. v1 formulation = ``magnetostatic_A`` (B = curl A).

    Args:
        example: "team20" (3-D static force) or "none".

    Returns the field reference + (optionally) a concrete example IR as JSON.
    """
    lines = [
        "# comsol-converter IR  (version %s)" % _ir.IR_VERSION,
        "",
        "ModelIR fields:",
        "  name, formulation('magnetostatic_A'), length_unit",
        "  geometry: dict  (OCC/Netgen spec OR source-native ref; geometry",
        "            translation is the hardest layer -- see cc_mapping('geometry'))",
        "  materials[]: {name, domain, kind(%s), mu_r, bh_curve, conductivity,"
        % "|".join(_ir.MATERIAL_KINDS),
        "               br, magnetization_dir}",
        "  sources[]:   {name, kind(%s), domain, n_turns, current, j_density,"
        % "|".join(_ir.SOURCE_KINDS),
        "               direction}",
        "  boundary_conditions[]: {boundary, kind(%s)}" % "|".join(_ir.BC_KINDS),
        "  mesh: {element_order, max_size, symmetry_fraction}",
        "        ** symmetry_fraction is load-bearing: a 1/N model needs",
        "           force/flux scaled by 1/fraction -- see cc_mapping('symmetry') **",
        "  study: {kind(stationary|frequency), frequency, nonlinear}",
        "  outputs[]: {name, kind(%s), target, point}" % "|".join(_ir.OUTPUT_KINDS),
        "  provenance: {source_tool, source_file, ...}  (round-trip audit)",
    ]
    if example == "team20":
        lines += ["", "## Example IR (TEAM problem 20, 3-D static force):", "```json",
                  json.dumps(_ir.team20_example().to_dict(), ensure_ascii=False, indent=1),
                  "```"]
    return "\n".join(lines)


# ============================================================
# COMSOL <-> IR <-> NGSolve mapping knowledge
# ============================================================

_MAPPING = {
    "overview": (
        "magnetostatic_A correspondence (COMSOL <-> IR <-> NGSolve):\n"
        "  physics      COMSOL 'Magnetic Fields (mf)' / Ampere's Law\n"
        "               <-> IR formulation='magnetostatic_A'\n"
        "               <-> NGSolve HCurl vector potential A, B = curl(A)\n"
        "  unknown      A (magnetic vector potential), Coulomb-gauged in COMSOL;\n"
        "               NGSolve HCurl is gauge-stable with a small regularization\n"
        "               or a tree-cotree / grad-grad term.\n"
        "Sub-topics: materials, sources, boundary_conditions, symmetry, force,\n"
        "study, geometry."
    ),
    "materials": (
        "MATERIALS:\n"
        "  air/linear:  COMSOL 'Ampere's Law' Relative permeability mur\n"
        "      <-> IR Material(kind='linear'|'air', mu_r) <-> NGSolve\n"
        "      nu = 1/(mu0*mu_r) as a domain-wise CoefficientFunction.\n"
        "  nonlinear iron: COMSOL 'Ampere's Law' B-H curve (interpolation)\n"
        "      <-> IR Material(kind='nonlinear_BH', bh_curve=ref) <-> NGSolve\n"
        "      nu(|B|) CoefficientFunction from the H(B) table; Newton needs\n"
        "      the differential reluctivity dH/dB (see fem knowledge).\n"
        "  permanent magnet: COMSOL remanent flux Br + direction\n"
        "      <-> IR Material(kind='permanent_magnet', br, magnetization_dir)\n"
        "      <-> NGSolve source term curl-side (Hc) or remanence."
    ),
    "sources": (
        "SOURCES:\n"
        "  COMSOL 'Coil' (Homogenized multiturn, N turns, coil current I)\n"
        "      <-> IR Source(kind='coil', n_turns=N, current=I)\n"
        "      <-> NGSolve impressed current density J0 over the coil region\n"
        "          (J0 = N*I/area * direction), RHS = (J0, v).\n"
        "  COMSOL External current density  <-> IR Source(kind='current_density',\n"
        "      j_density=[Jx,Jy,Jz])  <-> NGSolve (J, v).\n"
        "  TRAP (interop comsol_lab_tips 'homogenized_multiturn_coil'): the\n"
        "  (3)x(4) input gotcha -- confirm N*I total, not per-turn, matches."
    ),
    "boundary_conditions": (
        "BOUNDARY CONDITIONS:\n"
        "  magnetic_insulation (n x A = 0): COMSOL 'Magnetic Insulation'\n"
        "      <-> NGSolve Dirichlet on the tangential A (dirichlet=boundary).\n"
        "  perfect_magnetic_conductor (n x H = 0): COMSOL 'PMC'\n"
        "      <-> NGSolve natural BC (do nothing).\n"
        "  symmetry: choose Magnetic Insulation vs PMC per the field symmetry\n"
        "      (interop comsol_lab_tips 'symmetry_plane_bc' has the decision).\n"
        "  infinity: COMSOL Infinite Element Domain <-> NGSolve Kelvin transform\n"
        "      (radia-ngsolve kelvin_transformation) or a large truncation box."
    ),
    "symmetry": (
        "SYMMETRY SCALING (the single biggest COMSOL trip-wire):\n"
        "  interop comsol_lab_tips 'coil_current_scaling' -- a 1/N symmetry\n"
        "  model needs the post-processed force/flux/energy multiplied by N\n"
        "  (1/8 model -> x8). COMSOL does NOT warn; the wrong answer looks\n"
        "  'almost right'. IR.mesh.symmetry_fraction records the fraction so\n"
        "  BOTH generated sides apply the SAME scaling. cc_ir_to_ngsolve and\n"
        "  cc_ir_to_comsol both emit  scale = 1/symmetry_fraction  on outputs."
    ),
    "force": (
        "FORCE (TEAM 20 quantity of interest):\n"
        "  COMSOL: Maxwell surface stress tensor integrated over a surface\n"
        "      enclosing the body (or built-in Force calculation feature).\n"
        "  IR: Output(kind='force_maxwell', target=body).\n"
        "  NGSolve: Maxwell stress tensor T = (1/mu)(B (x) B) - (|B|^2/2mu) I,\n"
        "      integrate T.n over a closed surface in air around the body;\n"
        "      cross-check with radia force methods. Validate vs\n"
        "      team_force_motion('problem_20'). Apply symmetry scaling."
    ),
    "study": (
        "STUDY:\n"
        "  stationary linear: COMSOL Stationary <-> NGSolve direct/CG solve.\n"
        "  stationary nonlinear (iron): COMSOL Stationary + Fully Coupled /\n"
        "      Newton (bump iteration limit 25->250, comsol_lab_tips) <->\n"
        "      NGSolve Newton with differential reluctivity.\n"
        "  frequency (eddy, future): COMSOL Frequency Domain <-> NGSolve\n"
        "      complex HCurl A-formulation (sigma != 0)."
    ),
    "geometry": (
        "GEOMETRY (hardest layer -- handle explicitly):\n"
        "  COMSOL geometry -> IR: prefer exporting STEP from COMSOL (or read\n"
        "      the geometry sequence via mph_reader) -> IR.geometry={'step': path}.\n"
        "  IR -> NGSolve: load STEP via netgen.occ OCCGeometry, or build with\n"
        "      build123d/Cubit (radia-interop any_step_to_cubit_hex) and import\n"
        "      the mesh. IR -> COMSOL: import the same STEP in the Java model.\n"
        "  Using a SHARED STEP/mesh on both sides removes geometry as a\n"
        "  translation variable -- recommended for v1 / TEAM problems."
    ),
}


@mcp.tool()
def cc_mapping(topic: str = "overview") -> str:
    """COMSOL <-> IR <-> NGSolve correspondence knowledge (magnetostatic_A).

    Args:
        topic: overview | materials | sources | boundary_conditions |
               symmetry | force | study | geometry | all
    """
    if topic == "all":
        return "\n\n".join(f"== {k} ==\n{v}" for k, v in _MAPPING.items())
    return _MAPPING.get(topic, f"Unknown topic '{topic}'. Options: "
                         + ", ".join(_MAPPING) + ", all.")


# ============================================================
# COMSOL -> IR
# ============================================================

@mcp.tool()
def cc_comsol_to_ir(comsol_extract: str = "") -> str:
    """Guide building an IR from a COMSOL model extraction.

    Feed the output of the EXTERNAL COMSOL MCP (mph_reader / mph_to_java):
    the physics interface, materials, coil/current features, boundary
    conditions, and study. This returns the IR skeleton to fill plus the
    mapping pointers; the calling LLM completes the IR using cc_mapping.

    Args:
        comsol_extract: text/JSON dump from mph_reader or the .java export
                        (empty -> just returns the skeleton + instructions).
    """
    skeleton = _ir.ModelIR(name="<from_comsol>",
                           provenance={"source_tool": "comsol_mph_reader"}).to_dict()
    out = [
        "Build the IR from the COMSOL extraction below by mapping each feature:",
        "  - Magnetic Fields interface  -> formulation='magnetostatic_A'",
        "  - each Material / Ampere's Law mur or B-H  -> materials[] (cc_mapping('materials'))",
        "  - each Coil / External current density      -> sources[]  (cc_mapping('sources'))",
        "  - Magnetic Insulation / PMC / symmetry      -> boundary_conditions[] (cc_mapping('boundary_conditions'))",
        "  - symmetry fraction of the model            -> mesh.symmetry_fraction (cc_mapping('symmetry'))",
        "  - Stationary (+Newton if iron)              -> study",
        "  - Force / flux / point probes               -> outputs[] (cc_mapping('force'))",
        "  - geometry: export STEP, set geometry={'step': <path>} (cc_mapping('geometry'))",
        "",
        "IR skeleton to populate (then validate with cc_ir_schema):",
        "```json", json.dumps(skeleton, ensure_ascii=False, indent=1), "```",
    ]
    if comsol_extract.strip():
        out += ["", "--- supplied COMSOL extract (map the features above) ---",
                comsol_extract[:4000]]
    return "\n".join(out)


# ============================================================
# IR -> NGSolve  (code generation)
# ============================================================

def _gen_ngsolve(ir: "_ir.ModelIR") -> str:
    diri = "|".join(b.boundary for b in ir.boundary_conditions
                    if b.kind in ("magnetic_insulation",))
    order = ir.mesh.element_order
    scale = (1.0 / ir.mesh.symmetry_fraction) if ir.mesh.symmetry_fraction else 1.0
    mats = "\n".join(
        f"#   {m.domain}: kind={m.kind} "
        + (f"mu_r={m.mu_r}" if m.kind in ("air", "linear")
           else f"bh_curve={m.bh_curve}" if m.kind == "nonlinear_BH"
           else f"Br={m.br} dir={m.magnetization_dir}")
        for m in ir.materials)
    srcs = "\n".join(
        f"#   {s.name}: {s.kind} domain={s.domain} "
        + (f"NI={s.n_turns * s.current}" if s.kind == "coil"
           else f"J={s.j_density}")
        for s in ir.sources)
    nonlinear = ir.study.nonlinear
    L = [
        "# Auto-generated by comsol-converter cc_ir_to_ngsolve",
        f"# IR: {ir.name}  formulation={ir.formulation}  unit={ir.length_unit}",
        "# NGSolve magnetostatic A-formulation (B = curl A).",
        "from ngsolve import *",
        "from netgen.occ import OCCGeometry",
        "mu0 = 4e-7*pi",
        "",
        "# --- geometry / mesh ---",
        f"# IR.geometry = {json.dumps(ir.geometry, ensure_ascii=False)}",
        "# Recommended: shared STEP on both sides (cc_mapping('geometry')).",
        "# geo = OCCGeometry(r'<team20.step>'); mesh = Mesh(geo.GenerateMesh(maxh=...))",
        "mesh = Mesh('<team20.vol>')   # <-- supply mesh; region names must match IR domains",
        "",
        "# --- materials (domain-wise reluctivity nu = 1/(mu0*mu_r)) ---",
        f"{mats}",
        "# linear example: nu = mesh.MaterialCF({'air':1/mu0, 'pole':1/(mu0*1000)}, default=1/mu0)",
        ("# nonlinear: nu = nu_of_Bnorm(Norm(curl(gfu))) from the B-H table; Newton."
         if nonlinear else "nu = CoefficientFunction(1/mu0)  # fill per-material"),
        "",
        f"fes = HCurl(mesh, order={order}, dirichlet=\"{diri}\", nograds=False)",
        "u, v = fes.TnT()",
        "gfu = GridFunction(fes)   # the vector potential A",
        "",
        "# --- impressed source current density J0 (from IR.sources) ---",
        f"{srcs}",
        "# J0 = CoefficientFunction((...))   # NI/area * winding direction over coil region",
        "",
        "# --- weak form: (nu*curl A, curl v) + eps*(A,v) = (J0, v) ---",
        "a = BilinearForm(fes)",
        "a += (nu*curl(u)*curl(v) + 1e-6*u*v)*dx",
        "f = LinearForm(fes)",
        "f += J0*v*dx",
    ]
    if nonlinear:
        L += [
            "",
            "# nonlinear iron -> Newton (nu depends on |curl A|):",
            "# a = BilinearForm(fes); a += (nu_of(Norm(curl(u)))*curl(u)*curl(v))*dx + ...",
            "# from ngsolve.solvers import Newton ; Newton(a, gfu)",
        ]
    else:
        L += [
            "",
            "a.Assemble(); f.Assemble()",
            "gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse='sparsecholesky')*f.vec",
        ]
    L += [
        "",
        "B = curl(gfu)",
        f"SYM_SCALE = {scale}   # 1/symmetry_fraction (cc_mapping('symmetry'))",
        "",
        "# --- outputs ---",
    ]
    for o in ir.outputs:
        if o.kind == "B_at_point" and o.point:
            L.append(f"print('{o.name} B@{o.point}:', B(mesh(*{o.point})))")
        elif o.kind == "force_maxwell":
            L += [
                f"# {o.name}: Maxwell stress on '{o.target}', integrate T.n over a",
                "# closed air surface around it; multiply by SYM_SCALE.",
                "# T = (1/mu0)*(OuterProduct(B,B) - 0.5*(B*B)*Id(3)); force = Integrate(T*n, surf)*SYM_SCALE",
            ]
    return "\n".join(L)


@mcp.tool()
def cc_ir_to_ngsolve(ir_json: str = "") -> str:
    """Generate an NGSolve magnetostatic A-formulation script from an IR.

    Args:
        ir_json: a ModelIR as JSON (see cc_ir_schema). Empty -> use the
                 TEAM-20 example IR.

    Returns a runnable script SKELETON: HCurl A-formulation, domain-wise
    reluctivity, coil source, Dirichlet (magnetic insulation) BCs, linear
    solve or Newton (nonlinear iron), B = curl A, and outputs (B-probe,
    Maxwell-stress force) with the symmetry scaling applied. Geometry/mesh
    is a supply-point (cc_mapping('geometry') recommends a shared STEP).
    """
    try:
        ir = _ir.team20_example() if not ir_json.strip() else _ir.ModelIR.from_dict(json.loads(ir_json))
    except Exception as e:
        return f"Bad IR JSON: {type(e).__name__}: {e}\nUse cc_ir_schema to see the shape."
    return "```python\n" + _gen_ngsolve(ir) + "\n```"


# ============================================================
# IR -> COMSOL Java  (code generation scaffold)
# ============================================================

@mcp.tool()
def cc_ir_to_comsol(ir_json: str = "") -> str:
    """Generate a COMSOL Model Builder Java skeleton from an IR.

    Pairs with interop_comsol_livelink('java_api') for the com.comsol.model
    idioms and comsol_lab_tips for the gotchas (symmetry scaling, coil
    setup, infinite-element open boundary). Run the emitted Java via the
    external COMSOL MCP. Empty ir_json -> TEAM-20 example.
    """
    try:
        ir = _ir.team20_example() if not ir_json.strip() else _ir.ModelIR.from_dict(json.loads(ir_json))
    except Exception as e:
        return f"Bad IR JSON: {type(e).__name__}: {e}"
    scale = (1.0 / ir.mesh.symmetry_fraction) if ir.mesh.symmetry_fraction else 1.0
    L = [
        "// Auto-generated COMSOL Java skeleton by cc_ir_to_comsol",
        f"// IR: {ir.name}  ({ir.formulation})",
        "// Fill idioms from interop_comsol_livelink('java_api'); gotchas from",
        "// interop_comsol_lab_tips (symmetry/coil/infinite element).",
        "import com.comsol.model.*; import com.comsol.model.util.*;",
        "public class GenModel {",
        "  public static Model run() {",
        "    Model model = ModelUtil.create(\"Model\");",
        "    model.component().create(\"comp1\", true);",
        "    model.component(\"comp1\").geom().create(\"geom1\", 3);",
        "    // geometry: import shared STEP (cc_mapping('geometry'))",
        "    // model.component(\"comp1\").geom(\"geom1\").create(\"imp1\",\"Import\"); ... .set(\"filename\", \"<team20.step>\");",
        "    model.component(\"comp1\").physics().create(\"mf\",\"MagneticFields\",\"geom1\");",
    ]
    for m in ir.materials:
        L.append(f"    // material '{m.name}' on {m.domain}: kind={m.kind} "
                 + (f"mur={m.mu_r}" if m.kind in ("air", "linear") else f"BH={m.bh_curve}"))
    for s in ir.sources:
        L.append(f"    // coil '{s.name}': N={s.n_turns} I={s.current} (NI={s.n_turns*s.current})")
    for b in ir.boundary_conditions:
        L.append(f"    // BC '{b.boundary}': {b.kind}")
    L += [
        f"    // study: {ir.study.kind} nonlinear={ir.study.nonlinear}",
        "    model.study().create(\"std1\"); model.study(\"std1\").create(\"stat\",\"Stationary\");",
        f"    // OUTPUT scaling: multiply force/flux by {scale} (1/symmetry_fraction!)",
        "    return model;",
        "  }",
        "}",
    ]
    return "```java\n" + "\n".join(L) + "\n```"


# ============================================================
# Verification (team-benchmark cross-check)
# ============================================================

@mcp.tool()
def cc_verify(team_problem: int = 20) -> str:
    """Plan a COMSOL-vs-NGSolve cross-check against a TEAM reference solution.

    Args:
        team_problem: 20 (3-D static force, default), 13 (nonlinear
                      magnetostatic), 6 (sphere, analytic), 7 (eddy current).
    """
    ref = {
        20: "team_force_motion('problem_20')  -- 3-D static force, lab core",
        13: "team_magnetostatic('problem_13') -- 3-D nonlinear C-yoke",
        6:  "team_magnetostatic('problem_06') -- sphere, analytic reference",
        7:  "team_eddy_current('problem_07')  -- asymmetrical conductor w/ hole",
    }.get(team_problem, f"team_catalog('overview') -- problem {team_problem} not pre-wired")
    return "\n".join([
        f"# cc_verify: TEAM problem {team_problem}",
        f"reference: {ref}",
        "",
        "Procedure (IR-centered cross-check):",
        " 1. Build ONE IR for the problem (cc_ir_schema; share a STEP/mesh).",
        " 2. cc_ir_to_ngsolve  -> run via radia-ngsolve  -> quantity Q_ngs.",
        " 3. cc_ir_to_comsol   -> run via external COMSOL MCP -> quantity Q_comsol.",
        " 4. Pull the published reference value from the team_benchmark tool above.",
        " 5. Assert |Q_ngs - ref| and |Q_comsol - ref| within tolerance, and",
        "    |Q_ngs - Q_comsol| small. Apply symmetry scaling on BOTH (cc_mapping('symmetry')).",
        " 6. Record as a regression case.",
    ])


register_status_tool(
    mcp,
    server_name="mcp-server-comsol-converter",
    description=("COMSOL <-> IR <-> NGSolve model translation (magnetostatic "
                 "A-formulation v1). IR hub + mapping knowledge + NGSolve/COMSOL "
                 "codegen + team-benchmark cross-check."),
    subpackage="radia_mcp.comsol_converter",
    related_servers=["radia-interop", "radia-ngsolve", "team-benchmark", "fem", "bem"],
)


def main():
    if "--selftest" in sys.argv:
        print("comsol_converter MCP server self-test:")
        print("  IR version:", _ir.IR_VERSION)
        print("  team20 example materials:", len(_ir.team20_example().materials))
        print("  mapping topics:", ", ".join(_MAPPING))
        print("  ngsolve codegen chars:", len(_gen_ngsolve(_ir.team20_example())))
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
