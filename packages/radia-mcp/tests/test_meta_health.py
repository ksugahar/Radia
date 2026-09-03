"""Test that all cataloged radia_mcp.* subpackages import + have status_tool wired.

This is the canonical "is the package healthy" test — run on every PR
+ release. Catches:
  - import-time errors in any subpackage (broken refactor)
  - missing register_status_tool call in new servers (policy violation)
  - catalog.py drift from actual subpackage list
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest


def test_meta_health_all_subpackages_import():
    """All cataloged subpackages must import cleanly."""
    from radia_mcp.meta.server import radia_mcp_health
    h = radia_mcp_health()
    assert h["all_healthy"], (
        f"{h['n_servers_total'] - h['n_servers_healthy']} subpackages "
        f"failed to import: "
        f"{[r['name'] + ': ' + r.get('error', '?')[:80] for r in h['results'] if not r['import_ok']]}"
    )


def test_meta_health_caches_repeated_probe_in_process():
    """Repeated MCP health calls should not re-probe the whole fleet."""
    from radia_mcp.meta.server import radia_mcp_health

    fresh = radia_mcp_health(force_refresh=True)
    assert fresh["from_cache"] is False
    assert fresh["elapsed_ms"] >= 0.0
    assert fresh["cache_ttl_s"] > 0.0

    cached = radia_mcp_health()
    assert cached["from_cache"] is True
    assert cached["cache_age_s"] >= 0.0
    assert cached["n_servers_total"] == fresh["n_servers_total"]
    assert cached["n_servers_healthy"] == fresh["n_servers_healthy"]

    refreshed = radia_mcp_health(force_refresh=True)
    assert refreshed["from_cache"] is False


def test_status_optional_dep_probe_avoids_runtime_import(monkeypatch):
    """Status tools should use import metadata, not import heavy deps."""
    from radia_mcp.common import status as status_mod

    status_mod._probe_dep.cache_clear()

    def _boom(name):
        raise AssertionError(f"unexpected import of optional dependency {name}")

    monkeypatch.setattr(status_mod.importlib, "import_module", _boom)
    payload = status_mod.build_status_payload(
        server_name="mcp-server-test",
        description="test",
        subpackage="radia_mcp.test",
        optional_deps=["json"],
    )

    assert payload["optional_deps"]["json"]["installed"] is True
    assert "version" in payload["optional_deps"]["json"]


def test_meta_catalog_has_at_least_30_servers():
    """Sanity floor — catalog should not silently shrink."""
    from radia_mcp.meta import catalog
    assert len(catalog.CATALOG) >= 30, \
        f"catalog only has {len(catalog.CATALOG)} entries"


def test_meta_catalog_does_not_require_local_mcp_client_config():
    """The package catalog must not depend on an ignored workstation file."""
    from radia_mcp.meta import catalog

    assert ".mcp.json" not in (catalog.__doc__ or "")


def test_every_cataloged_server_has_register_status_tool():
    """Policy: every server.py in catalog must wire register_status_tool.

    Applies uniformly to all servers including `meta` itself — meta has its
    own status tool reporting on the catalog server's runtime state.
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


def test_meta_golden_gate_passes_catalog_contracts():
    """Meta server should expose a compact machine-readable golden gate."""
    from radia_mcp.meta.server import radia_mcp_golden_gate

    gate = radia_mcp_golden_gate()
    assert gate["level"] == "golden", gate
    assert gate["all_passed"] is True
    assert gate["n_servers"] >= 30
    assert any(c["name"] == "optuna_external_boundary" and c["ok"]
               for c in gate["checks"])
    assert "python scripts/gen_tools_doc.py --check" in gate["full_gate_commands"]


def test_generated_tools_inventory_is_current():
    """docs/TOOLS.md must stay synchronized with the live MCP tool registry."""
    pkg_root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, "scripts/gen_tools_doc.py", "--check"],
        cwd=pkg_root,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=60,
    )
    assert result.returncode == 0, (
        "docs/TOOLS.md is stale; run `python scripts/gen_tools_doc.py`.\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )


def test_readme_mcp_entrypoints_are_cataloged_or_external():
    """README discovery snippets must not advertise removed server commands."""
    from radia_mcp.meta import catalog

    pkg_root = Path(__file__).resolve().parents[1]
    readme = (pkg_root / "README.md").read_text(encoding="utf-8")
    mentioned = set(re.findall(r"\bmcp-server-[a-z0-9-]+\b", readme))
    allowed = {info["entry_point"] for info in catalog.CATALOG.values()}
    allowed.update(
        str(info.get("entry_point", ""))
        for info in catalog.EXTERNAL_PACKAGES.values()
        if str(info.get("entry_point", "")).startswith("mcp-server-")
    )
    unknown = sorted(mentioned - allowed)
    assert not unknown, f"README advertises unknown MCP entry points: {unknown}"


def test_meta_by_tag_optimization_finds_at_least_4():
    """Sanity: at least 4 optimization servers (bayesian-opt /
    evolutionary / data-assimilation / topology-optimization)."""
    from radia_mcp.meta.server import radia_mcp_by_tag
    result = radia_mcp_by_tag("optimization")
    assert result["n_matches"] >= 4, \
        f"only {result['n_matches']} optimization servers"


def test_meta_related_to_chart2d_includes_paper_writing():
    """Sanity: chart2d -> paper-writing cross-link (the canonical
    "use case A naturally suggests use case B" pattern).
    chart2d inherits the figure profile + gate stack, which as of
    2026-07-18 is served by mcp-server-paper-writing (figure merged
    in; the standalone figure server was retired).  This pins that
    the merged figure tooling stays discoverable from chart2d.
    """
    from radia_mcp.meta.server import radia_mcp_related
    result = radia_mcp_related("chart2d")
    names = [r["name"] for r in result["related"]]
    assert "paper-writing" in names, f"chart2d related: {names}"


def test_meta_related_exposes_external_optuna_mcp_without_catalog_import():
    """Official optuna-mcp is external but should still be discoverable."""
    from radia_mcp.meta.server import radia_mcp_related

    from_optuna = radia_mcp_related("optuna-mcp")
    from_optuna_names = [r["name"] for r in from_optuna["related"]]
    assert {"bayesian-opt", "evolutionary", "topology-optimization",
            "radia-streamfunction"}.issubset(from_optuna_names)

    from_bayes = radia_mcp_related("bayesian-opt")
    optuna = [r for r in from_bayes["related"] if r["name"] == "optuna-mcp"]
    assert optuna, f"bayesian-opt related: {[r['name'] for r in from_bayes['related']]}"
    assert optuna[0]["external"] is True
    assert optuna[0]["entry_point"] == "optuna-mcp"

    from_matlab = radia_mcp_related("radia-matlab")
    upstream = [
        row for row in from_matlab["related"] if row["name"] == "optuna-mcp"
    ]
    assert upstream
    assert upstream[0]["owns"].startswith("Every shared Optuna operation")
    assert upstream[0]["radia_difference_tool"].endswith(
        "matlab_optuna_mcp_route"
    )


def test_meta_related_mesh_chain_points_to_radia_ngsolve_registry():
    """CAD/mesh servers should point agents toward radia-ngsolve validation."""
    from radia_mcp.meta.server import radia_mcp_related

    for name in ("cubit", "build123d", "gmsh"):
        related = radia_mcp_related(name)
        names = [r["name"] for r in related["related"]]
        assert "radia-ngsolve" in names, f"{name} related: {names}"

    reverse = radia_mcp_related("radia-ngsolve")
    reverse_names = [r["name"] for r in reverse["related"]]
    for name in ("cubit", "build123d", "gmsh"):
        assert name in reverse_names, f"radia-ngsolve related: {reverse_names}"


def test_meta_catalog_exposes_selftest_and_heavy_audit_commands():
    """Agents should discover lightweight health checks separately from audits."""
    from radia_mcp.meta.server import radia_mcp_get, radia_mcp_overview

    overview = radia_mcp_overview()
    by_name = {entry["name"]: entry for entry in overview["servers"]}

    for name, entry in by_name.items():
        assert entry["selftest_command"] == f"{entry['entry_point']} --selftest", name

    cubit = radia_mcp_get("cubit")
    gmsh = radia_mcp_get("mcp-server-gmsh")
    assert cubit["audit_command"] == "mcp-server-cubit --selftest --audit-repo"
    assert gmsh["audit_command"] == "mcp-server-gmsh --selftest --audit-repo"
    assert "audit_command" not in radia_mcp_get("radia-ngsolve")


def test_all_related_links_are_bidirectional():
    """If A's related list contains B, then B's must contain A.

    LLMs navigating the catalog expect cross-server links to work in
    both directions. The 2026-05-24 thorough review found 33 one-way
    edges that silently broke navigation (e.g. `fem` -> `bem` but
    `bem` had no `fem` back-link). This test pins the invariant.
    """
    from radia_mcp.meta import catalog
    adj = {n: set(info.get("related", []))
           for n, info in catalog.CATALOG.items()}
    asymmetric = []
    for a, bs in adj.items():
        for b in bs:
            if b not in catalog.CATALOG:
                asymmetric.append(f"{a} -> {b} (target missing from catalog)")
            elif a not in adj.get(b, set()):
                asymmetric.append(
                    f"{a} -> {b}, but {b}'s related = "
                    f"{sorted(adj.get(b, []))} (missing reverse '{a}')"
                )
    assert not asymmetric, (
        "Found one-way related edges:\n  " + "\n  ".join(asymmetric)
    )


def test_topics_dispatcher_sync_for_dispatcher_servers():
    """For every dispatcher-style server (one tool with `topic` param +
    TOPICS dict in knowledge), verify that calling the dispatcher with
    each declared TOPICS key returns substantial content (not the
    "Unknown topic" error message).

    Catches drift where a topic gets renamed in the dispatcher's if-chain
    but the TOPICS dict in knowledge.py is not updated, or vice versa.
    The 2026-05-24 review identified this as a high-value safety net.
    """
    import asyncio
    import importlib
    import inspect
    from radia_mcp.meta import catalog

    drift = []
    skipped = []
    for name, info in catalog.CATALOG.items():
        subpkg = info["subpackage"]
        try:
            mod = importlib.import_module(f"{subpkg}.server")
        except Exception as e:
            skipped.append(f"{name}: server import failed ({type(e).__name__})")
            continue
        mcp = getattr(mod, "mcp", None)
        if mcp is None:
            continue
        # Find the dispatcher tool: a tool with a `topic` parameter,
        # excluding the auto-wired *_topics and *_status helpers.
        tools = list(mcp._tool_manager._tools.items())
        dispatchers = []
        for tname, tool in tools:
            if tname.endswith("_topics") or tname.endswith("_status"):
                continue
            fn = getattr(tool, "fn", None)
            if fn is None:
                continue
            try:
                sig = inspect.signature(fn)
            except (TypeError, ValueError):
                continue
            if "topic" in sig.parameters:
                dispatchers.append((tname, fn))
        if not dispatchers:
            continue  # multi-tool server with no `topic`-arg dispatcher
        if len(dispatchers) > 1:
            # Multiple dispatchers — accept (e.g. magnetic_materials),
            # check each.
            pass

        # Find the TOPICS dict in the same subpackage's knowledge module
        topics_dict = None
        for kn_mod_name in (subpkg + ".knowledge",
                            subpkg + "." + name.replace("-", "_") + "_knowledge",
                            subpkg + ".em_knowledge"):  # electromagnet
            try:
                kn = importlib.import_module(kn_mod_name)
            except ImportError:
                continue
            for attr in ("TOPICS",):
                t = getattr(kn, attr, None)
                if isinstance(t, dict) and t:
                    topics_dict = t
                    break
            if topics_dict:
                break
        if topics_dict is None:
            continue  # no TOPICS exposed (multi-tool / multi-knowledge-file server)

        # Call the dispatcher for each TOPICS key and verify
        for dname, dfn in dispatchers:
            for tk in topics_dict:
                try:
                    result = dfn(topic=tk)
                except Exception as e:
                    drift.append(f"{name}::{dname}(topic={tk!r}) raised "
                                  f"{type(e).__name__}: {str(e)[:80]}")
                    continue
                if not isinstance(result, str):
                    continue  # non-text response, skip
                head = result[:80].lower()
                if "unknown topic" in head or "not found" in head:
                    drift.append(
                        f"{name}::{dname}(topic={tk!r}) returned "
                        f"unknown-topic error although TOPICS lists it"
                    )

    if skipped:
        import warnings
        warnings.warn(f"Servers skipped: {skipped[:3]}")
    assert not drift, (
        "TOPICS / dispatcher chain out of sync:\n  " +
        "\n  ".join(drift)
    )


def test_all_tags_in_canonical_keep_set():
    """Tags must come from the 9-tag canonical set.

    The 2026-05-24 review found 48 distinct tags, 36 of which were
    1-server-limited (zero filter value). Consolidation drops the
    long tail and pins the keep-set so `radia_mcp_by_tag(...)` returns
    meaningful buckets.
    """
    from radia_mcp.meta import catalog
    KEEP = {"application", "optimization", "theory",
            "cad", "mesh", "fem", "solver", "ml", "meta"}
    seen = set()
    bad = []
    for n, info in catalog.CATALOG.items():
        for t in info.get("tags", []):
            seen.add(t)
            if t not in KEEP:
                bad.append(f"{n}: tag '{t}' outside keep-set")
    assert not bad, (
        f"Non-canonical tags present ({len(bad)}):\n  " +
        "\n  ".join(bad) +
        f"\nKeep-set: {sorted(KEEP)}"
    )
    # Also report tags that exist in keep-set but are NEVER used —
    # useful signal that the keep-set itself should shrink.
    unused = KEEP - seen
    if unused:
        # Not a hard failure, just a warning surfaced in pytest -v.
        import warnings
        warnings.warn(f"Keep-set tags never used: {sorted(unused)}")
