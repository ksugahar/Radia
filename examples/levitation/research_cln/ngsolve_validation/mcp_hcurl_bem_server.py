"""MCP server exposing HCurl BEM + closed-form Newton potential tools.

Tools provided:
    phi_box        - Banerjee-Gupta cuboid Newton potential (axis-aligned)
    phi_hex        - Arbitrary 8-node hex Newton potential (5-tet decomposition)
    hcurl_basis    - Schöberl-Zaglmayr H(curl) basis evaluation
    hex_bem_eigen  - HCurl BEM K matrix + eigenvalues for cuboid

Backend: Mathematica via wolframscript (works on Windows with Wolfram Engine).
Precision: 30 digits (Banerjee-Gupta), ~6 digits (HCurl BEM at WP=15).

Usage with Claude Code:
    Add to .claude/mcp.json or claude_desktop_config.json:
    {
      "mcpServers": {
        "hcurl-bem": {
          "command": "python",
          "args": ["w:/30_CauerLadderNetwork/2026_04_01_long-shaped-CLN/ngsolve_validation/mcp_hcurl_bem_server.py"]
        }
      }
    }

Date: 2026-05-03
"""

import asyncio
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from mcp.server import Server, NotificationOptions
from mcp.server.models import InitializationOptions
from mcp.server.stdio import stdio_server
import mcp.types as types

WOLFRAMSCRIPT = r"C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe"
SCRIPTS_DIR = Path(__file__).parent
PHI_WLS = SCRIPTS_DIR / "foster_cln_hex_phi.wls"
BASIS_WLS = SCRIPTS_DIR / "schoberl_zaglmayr_basis.wls"

server = Server("hcurl-bem")


def run_wolfram(code: str, timeout: float = 60.0) -> str:
    """Execute Mathematica code via wolframscript and return stdout."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".wls", delete=False, encoding="utf-8") as f:
        f.write(code)
        temp_path = f.name
    try:
        result = subprocess.run(
            [WOLFRAMSCRIPT, "-file", temp_path],
            capture_output=True, text=True, timeout=timeout, encoding="utf-8"
        )
        return result.stdout.strip()
    finally:
        Path(temp_path).unlink(missing_ok=True)


@server.list_tools()
async def list_tools() -> list[types.Tool]:
    return [
        types.Tool(
            name="phi_box",
            description=(
                "Banerjee-Gupta closed-form Newton potential of axis-aligned cuboid: "
                "Phi(P) = int_box 1/|r-P| dV. 30-digit precision."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "P": {
                        "type": "array", "items": {"type": "number"},
                        "minItems": 3, "maxItems": 3,
                        "description": "Observation point [Px, Py, Pz]"
                    },
                    "x0": {"type": "number"}, "x1": {"type": "number"},
                    "y0": {"type": "number"}, "y1": {"type": "number"},
                    "z0": {"type": "number"}, "z1": {"type": "number"},
                    "precision": {"type": "integer", "default": 30}
                },
                "required": ["P", "x0", "x1", "y0", "y1", "z0", "z1"]
            }
        ),
        types.Tool(
            name="phi_hex",
            description=(
                "Newton potential of arbitrary 8-node hex via 5-tet decomposition + "
                "Wilton-Rao-Glisson triangle face integrals. 30-digit precision."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "verts": {
                        "type": "array",
                        "items": {"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3},
                        "minItems": 8, "maxItems": 8,
                        "description": "8 vertices in standard hex ordering: 0..3 bottom CCW, 4..7 top CCW"
                    },
                    "P": {
                        "type": "array", "items": {"type": "number"},
                        "minItems": 3, "maxItems": 3,
                        "description": "Observation point [Px, Py, Pz]"
                    },
                    "precision": {"type": "integer", "default": 30}
                },
                "required": ["verts", "P"]
            }
        ),
        types.Tool(
            name="hcurl_basis",
            description=(
                "Schöberl-Zaglmayr H(curl) basis evaluation at a point. "
                "Returns list of {value vector, curl vector/scalar} for each basis function. "
                "Element types: TRIG (2D), QUAD (2D), TET (3D), HEX (3D), PRISM (3D), PYRAMID (3D)."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "element": {
                        "type": "string",
                        "enum": ["TRIG", "QUAD", "TET", "HEX", "PRISM", "PYRAMID"],
                        "description": "Element type"
                    },
                    "order": {"type": "integer", "minimum": 1, "default": 1, "description": "Polynomial order"},
                    "vnums": {
                        "type": "array", "items": {"type": "integer"},
                        "description": "Vertex numbers (for orientation; e.g. [0,1,2,3] for QUAD)"
                    },
                    "point": {
                        "type": "array", "items": {"type": "number"},
                        "description": "Reference-element evaluation point (2D for TRIG/QUAD, 3D otherwise)"
                    }
                },
                "required": ["element", "vnums", "point"]
            }
        ),
        types.Tool(
            name="hex_bem_eigenvalue_cuboid",
            description=(
                "Compute HCurl order-1 BEM K matrix (54 DOFs) on cuboid [0,a]x[0,b]x[0,c] and "
                "return generalized eigenvalues lambda such that K v = lambda M v. "
                "Convert to time constants tau_n = sigma * lambda. ~12 seconds. "
                "WARNING (2026-05-08): K assembly uses Mathematica NIntegrate (WP=15) which "
                "does NOT converge against 1/r self-singularity. tau_lead is artifactually "
                "inflated ~80%. For 5x2x1 mm Cu this returns ~19.43 us; true value (4-method "
                "cross-validated: radia-vim Phase F-4 Spherical Duffy + ELF time-domain + "
                "ELF Prony + ELF Foster N=3..5) is tau_lead = 10.9-11.6 us. "
                "Use radia-vim (S:/Radia/01_GitHub/src/ext/radia_vim) for quantitative work."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "a": {"type": "number", "description": "x-dimension [m]"},
                    "b": {"type": "number", "description": "y-dimension [m]"},
                    "c": {"type": "number", "description": "z-dimension [m]"},
                    "sigma": {"type": "number", "default": 5.8e7, "description": "Conductivity [S/m]"},
                    "n_modes": {"type": "integer", "default": 5, "description": "Number of leading modes to return"}
                },
                "required": ["a", "b", "c"]
            }
        ),
    ]


@server.call_tool()
async def call_tool(name: str, arguments: dict[str, Any]) -> list[types.TextContent]:
    if name == "phi_box":
        P = arguments["P"]
        x0, x1 = arguments["x0"], arguments["x1"]
        y0, y1 = arguments["y0"], arguments["y1"]
        z0, z1 = arguments["z0"], arguments["z1"]
        prec = arguments.get("precision", 30)
        code = f"""
Get["{PHI_WLS.as_posix()}"];
P = SetPrecision[{{{P[0]}, {P[1]}, {P[2]}}}, {prec}];
result = FosterClnHexPhi`PhiBoxClosed[P,
  SetPrecision[{x0}, {prec}], SetPrecision[{x1}, {prec}],
  SetPrecision[{y0}, {prec}], SetPrecision[{y1}, {prec}],
  SetPrecision[{z0}, {prec}], SetPrecision[{z1}, {prec}]];
Print[ToString[N[result, {prec - 2}], InputForm]];
"""
        result = run_wolfram(code)
        return [types.TextContent(type="text", text=f"Phi_box(P) = {result}")]

    if name == "phi_hex":
        verts = arguments["verts"]
        P = arguments["P"]
        prec = arguments.get("precision", 30)
        verts_str = ",".join(
            f"{{{v[0]},{v[1]},{v[2]}}}" for v in verts
        )
        code = f"""
Get["{PHI_WLS.as_posix()}"];
verts = SetPrecision[{{{verts_str}}}, {prec}];
P = SetPrecision[{{{P[0]}, {P[1]}, {P[2]}}}, {prec}];
result = FosterClnHexPhi`PhiHex[verts, P];
Print[ToString[N[result, {prec - 2}], InputForm]];
"""
        result = run_wolfram(code)
        return [types.TextContent(type="text", text=f"Phi_hex(P) = {result}")]

    if name == "hcurl_basis":
        elem = arguments["element"]
        order = arguments.get("order", 1)
        vnums = arguments["vnums"]
        point = arguments["point"]
        vnums_str = "{" + ",".join(map(str, vnums)) + "}"
        if len(point) == 2:
            point_subst = f"x -> {point[0]}, y -> {point[1]}"
        else:
            point_subst = f"x -> {point[0]}, y -> {point[1]}, z -> {point[2]}"
        coords = "x, y" if len(point) == 2 else "x, y, z"
        code = f"""
Get["{BASIS_WLS.as_posix()}"];
shapes = SchoberlZaglmayrBasis`Shapes{elem}[{order}, {vnums_str}, {coords}];
evaluated = Table[{{shape[[1]] /. {{{point_subst}}}, shape[[2]] /. {{{point_subst}}}}}, {{shape, shapes}}];
Print["NumBasis: ", Length[shapes]];
Print[ToString[N[evaluated, 12], InputForm]];
"""
        result = run_wolfram(code)
        return [types.TextContent(type="text", text=result)]

    if name == "hex_bem_eigenvalue_cuboid":
        a = arguments["a"]
        b = arguments["b"]
        c = arguments["c"]
        sigma_val = arguments.get("sigma", 5.8e7)
        n_modes = arguments.get("n_modes", 5)
        code = f"""
Get["{BASIS_WLS.as_posix()}"];
ClearAll[x, y, z, xp, yp, zp];
$MaxExtraPrecision = 50;
sigma = SetPrecision[{sigma_val}, 30];
mu0 = 4 Pi SetPrecision[10^-7, 30];
{{aBox, bBox, cBox}} = SetPrecision[{{{a}, {b}, {c}}}, 30];
ord = 1;
shapesRef = SchoberlZaglmayrBasis`ShapesHEX[ord, Range[0, 7], x, y, z];
nDOF = Length[shapesRef];

fastIntegrate3D[poly_, vars_] := Module[{{p = Expand[poly], rules}},
  If[p === 0, Return[0]];
  rules = CoefficientRules[p, vars];
  Total[Function[r, r[[2]] / Times @@ Map[(# + 1) &, r[[1]]]] /@ rules]
];

componentMassMatrix[shapes_, comp_] := Module[{{n = Length[shapes], M, exprs}},
  exprs = Table[Expand[shapes[[i, 1, comp]]], {{i, 1, n}}];
  M = Table[
    If[j < i, 0, fastIntegrate3D[exprs[[i]] exprs[[j]], {{x, y, z}}]],
    {{i, 1, n}}, {{j, 1, n}}
  ];
  M + Transpose[M] - DiagonalMatrix[Diagonal[M]]
];

MrefXX = componentMassMatrix[shapesRef, 1];
MrefYY = componentMassMatrix[shapesRef, 2];
MrefZZ = componentMassMatrix[shapesRef, 3];
Mphys = N[(bBox cBox / aBox) MrefXX + (aBox cBox / bBox) MrefYY + (aBox bBox / cBox) MrefZZ, 30];

Kentry[i_, j_] := Module[{{vi, vj, integrand, factor, scaled}},
  vi = shapesRef[[i, 1]];
  vj = shapesRef[[j, 1]] /. {{x -> xp, y -> yp, z -> zp}};
  integrand = (vi[[1]] vj[[1]] / aBox^2 + vi[[2]] vj[[2]] / bBox^2 + vi[[3]] vj[[3]] / cBox^2);
  factor = aBox^2 bBox^2 cBox^2;
  scaled = factor / Sqrt[aBox^2 (x - xp)^2 + bBox^2 (y - yp)^2 + cBox^2 (z - zp)^2];
  (mu0 / (4 Pi)) NIntegrate[integrand scaled,
    {{x, 0, 1}}, {{y, 0, 1}}, {{z, 0, 1}}, {{xp, 0, 1}}, {{yp, 0, 1}}, {{zp, 0, 1}},
    Method -> "GlobalAdaptive", WorkingPrecision -> 15,
    PrecisionGoal -> 4, AccuracyGoal -> 6, MaxRecursion -> 8]
];

K = Table[0, {{nDOF}}, {{nDOF}}];
Do[K[[i, j]] = Kentry[i, j]; K[[j, i]] = K[[i, j]],
   {{i, 1, nDOF}}, {{j, i, nDOF}}];
Kmat = N[K, 30];

eigvals = Eigenvalues[{{Kmat, Mphys}}];
realEigs = Re[eigvals];
sortedEigs = Reverse[Sort[Select[realEigs, # > 10^-30 &]]];
tausUs = 10^6 sigma * sortedEigs;
top = tausUs[[1 ;; Min[{n_modes}, Length[tausUs]]]];
Print[ToString[N[top, 8], InputForm]];
"""
        result = run_wolfram(code, timeout=120.0)
        return [types.TextContent(
            type="text",
            text=f"Top {n_modes} time constants tau_n [us] (cuboid {a}x{b}x{c} m):\n{result}"
        )]

    raise ValueError(f"Unknown tool: {name}")


async def main():
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            InitializationOptions(
                server_name="hcurl-bem",
                server_version="0.1.0",
                capabilities=server.get_capabilities(
                    notification_options=NotificationOptions(),
                    experimental_capabilities={}
                ),
            ),
        )


if __name__ == "__main__":
    asyncio.run(main())
