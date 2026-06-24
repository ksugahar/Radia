---
name: release-qud
description: Four-machine Radia release gate. Use when the user asks for release-qud, release_quad, release_qud, post-release deploy, or the Definition Of Done for Radia releases. Coordinates the three PyPI packages across LAB, 100号機, mdx, and hibino via tools/release_qud.py.
---

# release-qud

## Canonical Entry Point

Use only `tools/release_qud.py`; release work must go through QUD.

```powershell
python tools/release_qud.py preflight
python tools/release_qud.py phase0
python tools/release_qud.py phase8 --target lab,100,hibino
python tools/release_qud.py phase8e
python tools/release_qud.py phase9
python tools/release_qud.py all
python tools/release_qud.py done
```

## Machine Policy

| Machine | Install tier | Release command path |
|---|---|---|
| LAB | NAS editable | `phase8 --target lab` |
| 100号機 | NAS editable over SSH | `phase8 --target 100` |
| hibino | PyPI wheel consumer over `ssh hibino` | `phase8 --target hibino` |
| mdx | PyPI wheel consumer, no `radia-mcp` | `phase8e` |

`phase9` is the hard gate: LAB / 100号機 / mdx / hibino must agree on
versions, compatibility constants, and tracked file hashes. mdx reports
`radia-mcp` as `N/A`; that is intentional and is excluded from drift
comparison.

## Rules

- Do not call a release done until `python tools/release_qud.py done`
  exits 0.
- Use `ssh hibino` for hibino. For multi-line remote PowerShell, follow
  the repository SSH policy: pipe a script into
  `ssh hibino 'pwsh -ExecutionPolicy Bypass -Command -'`.
- Keep `packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/`
  and `packages/radia-mcp/docs/TOOLS.md` in sync when changing release
  or deploy knowledge.
