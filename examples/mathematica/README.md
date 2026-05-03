# Mathematica reference notebooks

Canonical Mathematica notation for selected Radia / radia-mcp
mathematical artifacts.  Mathematica is the *reference* layer; the
matching Python ports (CI-verified) live elsewhere in the repo.

## Contents

| File | Topic | Python port |
|------|-------|-------------|
| `RadiaBasis.m` | FEM basis functions (Tri / Tet, H1 Lagrange + RWG/HDiv RT₀ + L2) | `tests/basis/test_basis_functions.py` |

## Why Mathematica is the reference

- Symbolic differentiation / integration is exact and cheap.
- Notation matches the standard FEM textbooks (Brezzi-Marini, Bossavit,
  Monk), making the source equations directly auditable.
- The radia-mcp `basis_functions` tool returns the Mathematica string
  alongside SymPy/NumPy ports so the AI assistant always sees the
  authoritative form first.

## Running the notebook

Open `RadiaBasis.m` in Mathematica 13+.  The package executes in-place
and exposes the basis function list; copy any expression with
`InputForm[...]` to paste into a SymPy / Python translation when adding
a new test case.

## See also

- `radia-mcp basis_functions` MCP tool
- `tests/basis/test_basis_functions.py`
