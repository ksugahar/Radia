"""Test that all cataloged radia_mcp.* subpackages import + have status_tool wired.

This is the canonical "is the package healthy" test — run on every PR
+ release. Catches:
  - import-time errors in any subpackage (broken refactor)
  - missing register_status_tool call in new servers (policy violation)
  - catalog.py drift from actual subpackage list
"""

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


def test_meta_catalog_has_at_least_30_servers():
    """Sanity floor — catalog should not silently shrink."""
    from radia_mcp.meta import catalog
    assert len(catalog.CATALOG) >= 30, \
        f"catalog only has {len(catalog.CATALOG)} entries"


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


def test_meta_by_tag_optimization_finds_at_least_4():
    """Sanity: at least 4 optimization servers (bayesian-opt /
    evolutionary / data-assimilation / topology-optimization)."""
    from radia_mcp.meta.server import radia_mcp_by_tag
    result = radia_mcp_by_tag("optimization")
    assert result["n_matches"] >= 4, \
        f"only {result['n_matches']} optimization servers"


def test_meta_related_to_chart2d_includes_figure():
    """Sanity: chart2d → figure cross-link (the canonical
    "use case A naturally suggests use case B" pattern).
    Previously used mcmc → optuna, but mcmc was removed from the
    catalog 2026-05-26; chart2d/figure is the most stable pair to
    pin here (both shipped 2026-05 and not at risk of removal).
    """
    from radia_mcp.meta.server import radia_mcp_related
    result = radia_mcp_related("chart2d")
    names = [r["name"] for r in result["related"]]
    assert "figure" in names, f"chart2d related: {names}"


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
