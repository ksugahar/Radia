# HCurl BEM MCP Server

Mathematica + NGSolve Schöberl-Zaglmayr H(curl) basis + closed-form Newton
potential exposed as MCP tools.

## Tools

### `phi_box(P, x0, x1, y0, y1, z0, z1, precision=30)`
Banerjee-Gupta closed-form Newton potential of axis-aligned cuboid:

    Phi(P) = int_box 1/|r-P| dV

30-digit precision. ~100 ms per call.

### `phi_hex(verts, P, precision=30)`
Newton potential of arbitrary 8-node hex via 5-tet decomposition + WRG triangle face integrals. Verified against Banerjee-Gupta to 30 digits on cuboid.

`verts` = list of 8 [x, y, z] in standard hex ordering (0-3 bottom CCW, 4-7 top).

### `hcurl_basis(element, order, vnums, point)`
Schöberl-Zaglmayr H(curl) basis values + curls at a reference point.

Element types: `TRIG`, `QUAD`, `TET`, `HEX`, `PRISM`, `PYRAMID`.

Returns: list of `{value vector, curl}` for each basis function.

### `hex_bem_eigenvalue_cuboid(a, b, c, sigma=5.8e7, n_modes=5)`
Compute HCurl order-1 BEM K matrix (54 DOFs) on cuboid `[0,a]×[0,b]×[0,c]`,
solve generalized eigenvalue, return top `n_modes` time constants tau_n [us].

> **WARNING (2026-05-08 update)**: The K matrix in this tool uses general-purpose
> Mathematica `NIntegrate` (WP=15, PG=3, AG=5, MR=6) which does NOT converge against
> the 1/r self-singularity. The leading time constant returned is **artifactually
> inflated by ~80%**. For Cu 5×2×1 mm this tool returns ~19.43 us, but the true
> value (cross-validated by 4 methods: radia-vim Phase F-4 with Spherical Duffy +
> ELF time-domain single-exp tail fit + ELF Prony + ELF Foster N=3..5) is
> **τ_lead ≈ 10.9-11.6 us**.
>
> The HDiv-div-free + Spherical Duffy implementation in `S:/Radia/01_GitHub/packages/
> radia-vim` (Phase F-4) gives the correct answer. This tool's eigenvalue output
> should be treated as illustrative, not quantitative.

~12 seconds per call.

## Setup

Backend uses `wolframscript` (Mathematica). On Windows the path is hardcoded as
`C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe`.

### Requirements
- Python 3.12+ with `mcp` Python SDK (`pip install mcp`)
- Wolfram Mathematica (or Wolfram Engine free tier)
- Project files: `foster_cln_hex_phi.wls`, `schoberl_zaglmayr_basis.wls`

### Register with Claude Desktop / Code
Add to your MCP config (see `mcp_server_config_example.json`):

```json
{
  "mcpServers": {
    "hcurl-bem": {
      "command": "python",
      "args": ["w:/30_CauerLadderNetwork/2026_04_01_long-shape-CLN/ngsolve_validation/mcp_hcurl_bem_server.py"]
    }
  }
}
```

Restart Claude. Verify tools appear via `/mcp` or by asking Claude to use the `hcurl-bem` tools.

## Verification

Manual smoke test:
```bash
python mcp_hcurl_bem_server.py
# Then send MCP messages via stdio
```

Or invoke tools through Claude Code by asking:
- "Use phi_box with P=(10mm,0,0) on a 5x2x1 mm cuboid centered at origin"
- "Compute hex_bem_eigenvalue_cuboid for a=5e-3, b=2e-3, c=1e-3"

Expected: this tool returns leading tau ≈ 19.43 us, but that value is **artifactually
inflated** (see WARNING above). True τ_lead ≈ 10.9-11.6 μs (4-method cross-validated).
