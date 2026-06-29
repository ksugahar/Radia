"""MCP server: Hiruma 3-term Cauer Ladder Network with Kelvin truncation.

REPOSITORY POLICY (2026-05-05): Hiruma 3-term recurrence is the SOLE
CLN extraction method exposed. The accumulator-CLN
(Sugahara-Tanimoto J-side Schmidt with running A_pot accumulator) has
been DROPPED because in finite-dimensional FE space its A modes are
not exactly (1/μ)-orthogonal: per-stage L_dot vs L_B2 disagree by
15-38%, and this disagreement is mesh-refinement-INDEPENDENT
(empirically verified v22b vs v22f digit-for-digit). Hiruma 3-term
recurrence (Lanczos-like short recurrence on (G, C) system) gives
ONE canonical L per stage and is structurally tridiagonal in any
discretization.

Reference (canonical): Hiruma, S. and Igarashi, H., "Model order
reduction for linear time-invariant system with symmetric
positive-definite matrices: Synthesis of Cauer-equivalent circuit",
IEEE Trans. Magn., vol. 56, no. 3, 2020.

Tools:
    hiruma_3term_axisym_disk
        Hiruma 3-term recurrence on axisymmetric Cu disk + Kelvin.
        Returns per-stage λ_{2n-1} (= 1/R_n), λ_{2n} (= L_n),
        τ_n = λ_{2n} × λ_{2n-1}, A-drift diagnostic.

    kameari_kelvin_sphere_benchmark   (stub)
        Kameari A-Ω_r magnetic-sphere benchmark with Kelvin. Requires
        block-coupled HCurl × H¹ formulation that no NGSolve reference
        provides; deferred.

Removed 2026-05-05:
    kameari_kelvin_cuboid_breakdown   (accumulator-CLN, v15)
        Per repository policy. v15 script retained for historical
        breakdown demonstration but no longer exposed via MCP.

Backend: NGSolve (Python, in-process via subprocess).

Formulation (cuboid tool, cuboid_521_kameari_kelvin_v15_canonical.py):
    Geometry: Sugahara two-sphere Kelvin (radia.kelvin_geometry) — inner
        physical sphere R_K centered at origin contains cuboid + air gap;
        outer Kelvin sphere at offset, radius R_K, periodic-identified at
        inner boundary, GND vertex at offset (image of infinity).
    FE space: HCurl, nograds=True, dirichlet_bbnd="GND", Periodic. Plus
        BFS spanning-tree gauge mask on tree edges.
    Material: ν' = (ρ'/R)² ν_0 in Kelvin region (Nagamine CEFC 2026).
    Source: A_s = (B_0/2)(-y, x, 0) in inner, Convention B reduced-A
        -(ρ'/R)² × A_phys(local) in Kelvin (radia.kelvin_material).
    Iteration (Tanimoto-Sugahara TMAG3304725 eq. 2-4):
        R_n = 1 / <J_n, J_n/σ>_cond
        L_n = R_n × <J_n, A_pot_accum>_cond
        Schmidt: J_{n+1} = J_n - σ A_pot_accum / L_n

Date: 2026-05-05
"""
import asyncio
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

from mcp.server import Server, NotificationOptions
from mcp.server.models import InitializationOptions
from mcp.server.stdio import stdio_server
import mcp.types as types

# Path to the validated v15 driver (cuboid CLN).
V15_DRIVER = (Path(__file__).parent /
              "cuboid_521_kameari_kelvin_v15_canonical.py")

server = Server("kameari-kelvin")


def _run_cuboid_taun(a_mm, b_mm, c_mm, sigma, B0,
                      R_K_mm, h_cond_mm, h_air_mm,
                      order, n_stages, bonus_intorder):
    """Drive the v15 NGSolve script with parameter overrides via env JSON.

    The current v15 hard-codes the cuboid 5x2x1 Cu @ R_K=25mm; this
    wrapper writes a small driver that imports v15's helpers and runs
    with caller-provided parameters. Returns the parsed JSON result.
    """
    raise NotImplementedError(
        "Parameterized cuboid driver pending — v15 currently hard-codes "
        "5x2x1 Cu / R_K=25mm. To unblock: extract the parameter constants "
        "in v15 main() into a function signature, then call it here."
    )


@server.list_tools()
async def list_tools() -> list[types.Tool]:
    return [
        # NOTE: accumulator-CLN tool (kameari_kelvin_cuboid_breakdown, v15)
        # was REMOVED 2026-05-05. In finite-dim FE space, accumulator-CLN
        # produces L_dot ≠ L_B2 by 15-38% (cross-term contamination),
        # mesh-refinement-INDEPENDENT (v22b vs v22f digit-for-digit).
        # Hiruma 3-term recurrence is the canonical replacement: structurally
        # tridiagonal in discrete space, single L per stage. See v22 → v23
        # comparison study and Hiruma & Igarashi IEEE TMag 56(3) 2020.
        types.Tool(
            name="kameari_kelvin_sphere_benchmark",
            description=(
                "[STUB] Kameari A-Omega_r magnetic sphere benchmark with "
                "Kelvin transformation. Reference: Kameari 2025/10/14 slides, "
                "0.001%% error on coarse mesh for mu_r=1000, B_0=1T, theory "
                "Bz0=2.9940 T. Requires HCurl x H1 block-coupled formulation "
                "with boundary integral source ∫ N·(H_s × n) ds at the "
                "conductor surface — different code path from the cuboid "
                "CLN tool. Implementation pending."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "mu_r": {"type": "number", "default": 1000.0},
                    "B0": {"type": "number", "default": 1.0},
                    "a_m": {"type": "number", "default": 1.0,
                            "description": "Sphere radius [m]"},
                    "R_K_m": {"type": "number", "default": 2.0,
                              "description": "Kelvin sphere radius [m]"},
                    "order": {"type": "integer", "default": 3}
                },
                "required": []
            }
        ),
        types.Tool(
            name="hiruma_3term_3d_disk",
            description=(
                "3D HCurl Hiruma 3-term Cauer extraction on Cu CYLINDER ONLY + "
                "Sugahara two-sphere Kelvin (vacuum BC). "
                "**WARNING (2026-05-09)**: gauge collapse on rectangular "
                "(square prism, cuboid 5x2x1): K^-1(Mw) lives entirely in "
                "tree-cotree near-kernel subspace, recurrence locks onto a "
                "rank-2 sub-dominant Foster mode (e.g., 65 μs on A1 vs BEM "
                "212 μs leading). Use `tanimoto_AT_3d_cuboid` instead for any "
                "rectangular conductor — Tanimoto A-T avoids the K^-1·M step "
                "by curl-based propagation (J=curl T, B=curl A) and recovers "
                "BEM leading τ within 3-12% on cylinder/A1/cuboid_521. "
                "Cylinder-only validation 2026-05-08 vs axisym v23: 5/6 stages "
                "within 5%, 3/6 within 0.5% (k=0: 0.02%, k=1: 0.03%, k=2: "
                "0.48%, k=3: 2.31%, k=4: 2.07%, k=5: 92%). "
                "Critical implementation features (do NOT modify lightly): "
                "(a) algebraic Helmholtz-Hodge projection to remove gauge "
                "contamination from K^-1 step (gauge_frac is ~80%, only ~15%% "
                "of M-energy is physical signal); (b) ORDER_PHI = ORDER (=2) "
                "in the H-H Poisson space (ORDER_PHI=1 leaks higher-order "
                "grad-zero modes through K^-1 and inflates τ_pair[5] by 5x); "
                "(c) Kahan-Parlett MGS-twice K-orth (free insurance, drift "
                "→ machine ε but does NOT affect τ — Lanczos drift was a "
                "red herring per the f3c0a10 diagnostic). "
                "Backend: disk_3d_kameari_hiruma_v5_orderphi.py "
                "(commit af9c681). ORDER=2 HCurl, Pardiso. ~1-2 min/stage."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "R_disk_mm": {"type": "number", "default": 10.0},
                    "thickness_mm": {"type": "number", "default": 2.0},
                    "sigma": {"type": "number", "default": 5.8e7},
                    "n_stages": {"type": "integer", "default": 6},
                    "order_phi": {"type": "integer", "default": 2,
                                  "description": "H-H Poisson space order. "
                                                 "Use 2 for paper-grade results "
                                                 "(default), 1 for regression."},
                    "reorth_mode": {"type": "string", "default": "mgs_twice",
                                    "enum": ["none", "single", "mgs_twice"],
                                    "description": "Re-orth mode (insurance, "
                                                   "does not affect τ)."}
                },
                "required": []
            }
        ),
        types.Tool(
            name="tanimoto_AT_3d_cuboid",
            description=(
                "3D FEM Cauer ladder extraction on arbitrary 3D conductor "
                "(cylinder, square prism, cuboid) using Tanimoto A-T "
                "alternating Poisson with Helmholtz-Hodge projected source. "
                "**Solves rectangular pathology** of v5/Hiruma (gauge collapse) "
                "by replacing M = ∫σuv mass matrix with curl-based propagation "
                "(J=curl T, B=curl A — curl annihilates HCurl gauge subspace) "
                "AND projecting source A_s onto A_s'·n=0 subspace before T-eq "
                "(Sugahara diagnostic 2026-05-09: original `f_T = ∫σA_s·W` is "
                "ill-posed on rectangular because A_s·n ≠ 0 on x/y faces "
                "violates div(σA)=0 across boundary, leaks noise into the "
                "iteration). Reference: Tanimoto修論 CLN_AT.ipynb (W:/00_CAE/"
                "NGSolve/谷本/修論/) wire skin-effect base recipe + Sugahara "
                "H-H source projection adaptation for uniform B_z eddy "
                "physics. Validation 2026-05-09 (3 geometries): "
                "cylinder R=10mm t=2mm: stage 1 τ=225.8 μs vs BEM-Foster 219 "
                "(+3.0%); square prism A1 17.72×17.72×2 (V=cyl): 219.0 μs vs "
                "BEM 212 (+3.3%); cuboid 5×2×1: 9.69 μs vs BEM 11 (-12%). "
                "Stage 0 always spurious (τ ≈ 5× leading, initialization "
                "artifact). Compare in offset-by-1 fashion vs BEM/Hiruma "
                "indexing. Stages 1-3 monotone, 4+ noise floor. Backend: "
                "tanimoto_AT_geom_sweep.py with H-H projection (ORDER+1 H1, "
                "one DOF pinned). ORDER=2 HCurl, Pardiso. ~5-10s/stage."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "geometry": {
                        "type": "string",
                        "enum": ["cylinder", "square", "cuboid_521"],
                        "default": "cuboid_521",
                        "description": "cylinder (R=10mm,t=2mm), square (A1: "
                                       "17.72x17.72x2mm), cuboid_521 (5x2x1mm)"
                    },
                    "n_stages": {"type": "integer", "default": 6},
                    "order": {"type": "integer", "default": 2,
                              "description": "HCurl order"},
                    "h_cond_mm": {"type": "number", "default": 1.0,
                                  "description": "Mesh size in conductor (mm). "
                                                 "Use 0.3 for cuboid_521."}
                },
                "required": []
            }
        ),
        types.Tool(
            name="hiruma_3term_axisym_disk",
            description=(
                "Hiruma 3-term recurrence (Lanczos-like) for CLN extraction "
                "on axisymmetric Cu disk + Kelvin truncation. "
                "Reference: Hiruma & Igarashi, IEEE Trans. Magn. 56(3), 2020 "
                "'Model order reduction for linear time-invariant system "
                "with symmetric positive-definite matrices: Synthesis of "
                "Cauer-equivalent circuit'. "
                "Algorithm: w_{2n+1} = w_{2n-1} - A·w_{2n}/λ_{2n} (3-term, "
                "subtracts from PREVIOUS-PREVIOUS, not previous). Builds "
                "Cauer ladder via short Lanczos recurrence on (G, C) system, "
                "structurally tridiagonal in any FE discretization. "
                "OUTPUT: per-stage λ_{2n-1} (= 1/R_n) and λ_{2n} (= L_n), "
                "τ_n = λ_{2n} × λ_{2n-1}. Single canonical L per stage "
                "(no L_dot vs L_B2 ambiguity). "
                "Adopted as the CLN method of choice in this server "
                "(accumulator-CLN dropped 2026-05-05 due to "
                "discrete-space orthogonality degradation). "
                "Backend: cuboid_521_kameari_kelvin_v23_hiruma_3term.py "
                "(parameterized only via script edit currently)."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "R_disk_mm": {"type": "number", "default": 10.0},
                    "thickness_mm": {"type": "number", "default": 2.0},
                    "sigma": {"type": "number", "default": 5.8e7},
                    "B0": {"type": "number", "default": 1.0},
                    "n_stages": {"type": "integer", "default": 6}
                },
                "required": []
            }
        ),
    ]


@server.call_tool()
async def call_tool(name: str, arguments: dict[str, Any]) -> list[types.TextContent]:
    # NOTE: kameari_kelvin_cuboid_breakdown (accumulator-CLN, v15) was
    # removed 2026-05-05. Server-wide policy: Hiruma 3-term only.

    if name == "kameari_kelvin_sphere_benchmark":
        return [types.TextContent(
            type="text",
            text=("kameari_kelvin_sphere_benchmark is a STUB. Implementation "
                  "requires writing the HCurl×H¹ block-coupled A-Ω_r driver "
                  "from Kameari slide-deck weak forms (no NGSolve reference "
                  "exists). See docs/kelvin/KELVIN_TRANSFORMATION.md §7.4.1.")
        )]

    if name == "hiruma_3term_3d_disk":
        defaults = {"R_disk_mm": 10.0, "thickness_mm": 2.0,
                    "sigma": 5.8e7, "n_stages": 6,
                    "order_phi": 2, "reorth_mode": "mgs_twice"}
        non_default = {k: arguments.get(k) for k in defaults
                       if k in arguments and arguments[k] != defaults[k]}
        v5_driver = SCRIPTS_DIR / "disk_3d_kameari_hiruma_v5_orderphi.py"
        if not v5_driver.exists():
            return [types.TextContent(
                type="text",
                text=("v5 driver not found at " + str(v5_driver))
            )]
        cmd = [sys.executable, str(v5_driver)]
        if "order_phi" in arguments:
            cmd.extend(["--order-phi", str(arguments["order_phi"])])
        if "reorth_mode" in arguments:
            cmd.extend(["--mode", str(arguments["reorth_mode"])])
        if non_default and any(k in non_default
                               for k in ["R_disk_mm", "thickness_mm",
                                         "sigma", "n_stages"]):
            return [types.TextContent(
                type="text",
                text=("Geometry/material parameterization not yet wired "
                      "through (R, thickness, sigma, n_stages). Currently "
                      "fixed: R=10mm, t=2mm Cu (sigma=5.8e7). To enable: "
                      "parameterize disk_3d_kameari_hiruma_v5_orderphi.py.")
            )]
        result = subprocess.run(
            cmd, capture_output=True, text=True, encoding="utf-8",
            timeout=1800,  # 30 min for 6 stages with N=2 HCurl
        )
        # v5 saves JSON with --out-tag suffix
        out_tag = f"orderphi{arguments.get('order_phi', 2)}"
        json_path = SCRIPTS_DIR / (
            f"2026-05-08-disk_3d_kameari_hiruma_v5_{out_tag}_results.json"
        )
        if not json_path.exists():
            return [types.TextContent(
                type="text",
                text=("v5 driver did not produce expected JSON. Stdout tail:\n"
                      + (result.stdout[-2000:] if result.stdout else ""))
            )]
        data = json.loads(json_path.read_text())
        # Format depends on v5 schema; keep robust by dumping the stages list
        stages = data.get("stages", data.get("tau_pair", []))
        if isinstance(stages, list) and stages:
            stages_text = "\n".join(
                f"  k={i}: τ_pair = {s if isinstance(s, (int, float)) else s.get('tau_us', s)} μs"
                for i, s in enumerate(stages)
            )
        else:
            stages_text = "  (see JSON for full results)"
        return [types.TextContent(
            type="text",
            text=("3D HCurl Hiruma 3-term + H-H projection (ORDER_PHI=2) + "
                  "MGS-twice re-orth on Cu cylinder + Kelvin:\n"
                  f"{stages_text}\n\n"
                  "Validated 2026-05-08 vs axisym v23: 5/6 stages within 5%, "
                  "3/6 within 0.5% (commit af9c681). "
                  "Reference: Hiruma & Igarashi IEEE TMag 56(3) 2020.\n"
                  f"Full results: {json_path}")
        )]

    if name == "tanimoto_AT_3d_cuboid":
        defaults = {"geometry": "cuboid_521", "n_stages": 6, "order": 2,
                    "h_cond_mm": 1.0}
        geometry = arguments.get("geometry", defaults["geometry"])
        n_stages = arguments.get("n_stages", defaults["n_stages"])
        order = arguments.get("order", defaults["order"])
        h_cond_mm = arguments.get("h_cond_mm", defaults["h_cond_mm"])
        # Recommend smaller mesh for cuboid_521 if user did not override
        if geometry == "cuboid_521" and "h_cond_mm" not in arguments:
            h_cond_mm = 0.3
        driver = SCRIPTS_DIR / "tanimoto_AT_geom_sweep.py"
        if not driver.exists():
            return [types.TextContent(
                type="text",
                text=("Tanimoto A-T driver not found at " + str(driver))
            )]
        out_tag = f"mcp_{geometry}_p{order}"
        cmd = [sys.executable, str(driver),
               "--geom", str(geometry),
               "--stages", str(n_stages),
               "--order", str(order),
               "--h-cond-mm", str(h_cond_mm),
               "--out-tag", out_tag]
        result = subprocess.run(
            cmd, capture_output=True, text=True, encoding="utf-8",
            timeout=600,
        )
        json_path = SCRIPTS_DIR / (
            f"2026-05-09-tanimoto_AT_{out_tag}_results.json"
        )
        if not json_path.exists():
            return [types.TextContent(
                type="text",
                text=("Tanimoto A-T driver did not produce expected JSON. "
                      "Stdout tail:\n" + (result.stdout[-2000:]
                                            if result.stdout else ""))
            )]
        data = json.loads(json_path.read_text())
        stages = data.get("stages", [])
        stages_text = "\n".join(
            f"  stage {s['stage']}: R_{s['R_2n_idx']} = {s['R_2n']:.4e}, "
            f"L_{s['L_2nplus1_idx']} = {s['L_2nplus1']:.4e}, "
            f"τ_pair = {s['tau_pair_us']:.4f} μs"
            for s in stages
        )
        bem_ref = {"cylinder": 219.3, "square": 212.1, "cuboid_521": 11.0}
        ref_msg = (f"BEM-Foster leading τ for {geometry}: {bem_ref[geometry]} μs"
                   if geometry in bem_ref else "")
        return [types.TextContent(
            type="text",
            text=(f"Tanimoto A-T + Helmholtz-Hodge source projection on "
                  f"{geometry} (ORDER={order}, h_cond={h_cond_mm} mm, "
                  f"{n_stages} stages, no Kelvin sphere needed):\n"
                  f"{stages_text}\n\n"
                  f"Stage 0 is initialization-spurious (always offset by -1 "
                  f"vs BEM/v5 indexing). Stage 1 = leading physical Foster "
                  f"pair. {ref_msg}\n"
                  f"Validation 2026-05-09: cylinder 3.0%%, square 3.3%%, "
                  f"cuboid 12%% vs BEM-Foster.\n"
                  f"Full results: {json_path}")
        )]

    if name == "hiruma_3term_axisym_disk":
        defaults = {"R_disk_mm": 10.0, "thickness_mm": 2.0,
                    "sigma": 5.8e7, "B0": 1.0, "n_stages": 6}
        non_default = {k: arguments.get(k) for k in defaults
                       if k in arguments and arguments[k] != defaults[k]}
        if non_default:
            return [types.TextContent(
                type="text",
                text=("Parameterization not yet wired through. Custom args: "
                      f"{non_default}. Currently fixed: R=10mm, t=2mm Cu, "
                      "Order=3, N=6 stages. To enable: parameterize "
                      "cuboid_521_kameari_kelvin_v23_hiruma_3term.py.")
            )]

        v23_driver = SCRIPTS_DIR / "cuboid_521_kameari_kelvin_v23_hiruma_3term.py"
        result = subprocess.run(
            [sys.executable, str(v23_driver)],
            capture_output=True, text=True, encoding="utf-8",
            timeout=300,
        )
        json_path = v23_driver.with_suffix(".json")
        if not json_path.exists():
            return [types.TextContent(
                type="text",
                text=("v23 driver did not produce JSON output. Stdout tail:\n"
                      + (result.stdout[-2000:] if result.stdout else ""))
            )]
        data = json.loads(json_path.read_text())
        stages_text = "\n".join(
            f"  n={d['n']}: λ_{2*d['n']-1} = {d['lam_2nplus1_R_next']:.4e} "
            f"(= 1/R_n), λ_{2*d['n']} = {d['lam_2n_L']:.4e} (= L_n), "
            f"τ_n = {d['tau_us']:.4f} μs, "
            f"A-drift = {d['a_drift']:.3e}"
            for d in data["stages"]
        )
        return [types.TextContent(
            type="text",
            text=(f"Hiruma 3-term recurrence on Cu disk + Kelvin (axisym):\n"
                  f"  λ_1 (R_0) = {data['lam_1_R0']:.4e}\n"
                  f"{stages_text}\n\n"
                  "Reference: Hiruma & Igarashi, IEEE Trans. Magn. 56(3), 2020. "
                  "Per-stage canonical L_n (no L_dot/L_B2 ambiguity). "
                  "A-drift growth indicates breakdown stage in finite FE.")
        )]

    raise ValueError(f"Unknown tool: {name}")


async def main():
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            InitializationOptions(
                server_name="kameari-kelvin",
                server_version="0.1.0",
                capabilities=server.get_capabilities(
                    notification_options=NotificationOptions(),
                    experimental_capabilities={}
                ),
            ),
        )


if __name__ == "__main__":
    asyncio.run(main())
