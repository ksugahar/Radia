"""Test that all 35 radia_mcp.* subpackages import + have status_tool wired.

This is the canonical "is the package healthy" test — run on every PR
+ release. Catches:
  - import-time errors in any subpackage (broken refactor)
  - missing register_status_tool call in new servers (policy violation)
  - catalog.py drift from actual subpackage list
"""

from pathlib import Path

import pytest


def test_meta_health_all_subpackages_import():
    """All 35 subpackages must import cleanly."""
    from radia_mcp.meta.server import radia_mcp_health
    h = radia_mcp_health()
    assert h["all_healthy"], (
        f"{h['n_servers_total'] - h['n_servers_healthy']} subpackages "
        f"failed to import: "
        f"{[r['name'] + ': ' + r.get('error', '?')[:80] for r in h['results'] if not r['import_ok']]}"
    )


def test_meta_catalog_has_at_least_30_servers():
    """Sanity floor — catalog should not silently shrink."""
    from radia_mcp.meta import catalog
    assert len(catalog.CATALOG) >= 30, \
        f"catalog only has {len(catalog.CATALOG)} entries"


def test_every_cataloged_server_has_register_status_tool():
    """Policy: every server.py in catalog must wire register_status_tool.

    Exception: the meta server itself (it IS the catalog, no self-status).
    """
    from radia_mcp.meta import catalog
    repo_root = Path(__file__).resolve().parent.parent
    src_root = repo_root / "src" / "radia_mcp"

    missing = []
    for name, info in catalog.CATALOG.items():
        subpkg_dir = info["subpackage"].replace("radia_mcp.", "")
        server_py = src_root / subpkg_dir / "server.py"
        if not server_py.exists():
            missing.append(f"{name}: server.py missing at {server_py}")
            continue
        src = server_py.read_text(encoding="utf-8")
        if "register_status_tool" not in src:
            missing.append(f"{name}: register_status_tool not wired")

    assert not missing, "Status tool wiring missing:\n" + "\n".join(missing)


def test_meta_overview_returns_expected_shape():
    """Smoke-test the overview tool shape."""
    from radia_mcp.meta.server import radia_mcp_overview
    ov = radia_mcp_overview()
    assert "n_servers" in ov
    assert "servers" in ov
    assert "tags_available" in ov
    assert ov["n_servers"] == len(ov["servers"])
    # Each server entry has required keys
    for srv in ov["servers"]:
        for k in ("name", "subpackage", "entry_point",
                   "description", "primary_tools", "tags"):
            assert k in srv, f"{srv.get('name','?')} missing key {k}"


def test_meta_by_tag_optimization_finds_at_least_4():
    """Sanity: at least 4 optimization servers (optuna / bayesian-opt /
    evolutionary / mcmc / data-assimilation / topology-optimization)."""
    from radia_mcp.meta.server import radia_mcp_by_tag
    result = radia_mcp_by_tag("optimization")
    assert result["n_matches"] >= 4, \
        f"only {result['n_matches']} optimization servers"


def test_meta_related_to_mcmc_includes_optuna():
    """Sanity: mcmc → optuna cross-link."""
    from radia_mcp.meta.server import radia_mcp_related
    result = radia_mcp_related("mcmc")
    names = [r["name"] for r in result["related"]]
    assert "optuna" in names, f"mcmc related: {names}"
