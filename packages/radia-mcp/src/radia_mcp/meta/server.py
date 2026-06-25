"""
MCP Server: radia_mcp.meta

Cross-server catalog + health check for the radia_mcp
ecosystem. Use this server as the FIRST stop when you don't know
which other server to call.

Usage:
    mcp-server-radia-meta              # stdio
    mcp-server-radia-meta --selftest   # self-test
"""

import importlib
import shutil
import sys

from mcp.server.fastmcp import FastMCP

from . import catalog, bug_patterns
from ..common import register_status_tool

mcp = FastMCP("mcp-server-radia-meta")


# ============================================================
# Catalog tools
# ============================================================

@mcp.tool()
def radia_mcp_overview() -> dict:
    """Authoritative catalog of all radia_mcp.* servers.

    Use this as your FIRST tool call if you don't know which radia_mcp
    server has the knowledge you need. Returns a dict with:
        - n_servers: how many subpackages exist
        - servers: list of {name, subpackage, entry_point,
                              description, primary_tools, related, tags,
                              selftest_command, optional audit_command}

    Filter by tag with `radia_mcp_by_tag`; drill into a single server
    with `radia_mcp_get`.
    """
    servers = catalog.list_all()
    external = catalog.list_external()
    return {
        "n_servers": len(servers),
        "servers": servers,
        "tags_available": sorted({
            t for s in servers for t in s.get("tags", [])
        }),
        "external_packages": external,
        "n_external_packages": len(external),
        "next_step_hint":
            "Call <server>_status() on a specific server for full "
            "introspection + dependency probe. External packages "
            "(optuna-mcp / elf / comsol / mcp-server-document) ship from their own "
            "repos — see entries in `external_packages` for install paths.",
    }


@mcp.tool()
def radia_mcp_get(name: str) -> dict:
    """Look up one server by short name (e.g. 'bayesian-opt', 'ih', 'kelvin')."""
    info = catalog.get(name)
    if info is None:
        return {
            "error": f"Unknown server '{name}'",
            "available": sorted(catalog.CATALOG.keys()),
        }
    return {"name": name, **info}


@mcp.tool()
def radia_mcp_by_tag(tag: str) -> dict:
    """Servers tagged with `tag`.

    Common tags: 'fem', 'bem', 'application', 'theory', 'optimization',
    'ml', 'cad', 'mesh', 'solver', 'materials', 'meta', 'benchmark'.
    """
    matches = catalog.find_by_tag(tag)
    return {
        "tag": tag,
        "n_matches": len(matches),
        "servers": matches,
    }


@mcp.tool()
def radia_mcp_related(name: str) -> dict:
    """Servers that pair well with `name` (e.g. radia_mcp_related('bayesian-opt')
    returns evolutionary + topology-optimization + related ML servers)."""
    related = catalog.find_related(name)
    return {
        "name": name,
        "n_related": len(related),
        "related": related,
    }


# ============================================================
# Health check
# ============================================================

@mcp.tool()
def radia_mcp_health() -> dict:
    """Probe importability of every radia_mcp.* subpackage.

    Returns per-subpackage import_ok status. Surfaces broken installs
    (missing entry points, half-applied editable install, etc.) early.
    """
    out = []
    for name, info in catalog.CATALOG.items():
        subpkg = info["subpackage"]
        result = {"name": name, "subpackage": subpkg}
        try:
            importlib.import_module(subpkg)
            result["import_ok"] = True
        except ImportError as e:
            result["import_ok"] = False
            result["error"] = str(e)
        except Exception as e:
            result["import_ok"] = False
            result["error"] = f"{type(e).__name__}: {e}"
        out.append(result)

    healthy = sum(1 for r in out if r["import_ok"])
    return {
        "n_servers_total": len(out),
        "n_servers_healthy": healthy,
        "all_healthy": healthy == len(out),
        "results": out,
        "python_version": (f"{sys.version_info.major}."
                            f"{sys.version_info.minor}."
                            f"{sys.version_info.micro}"),
    }


@mcp.tool()
def radia_mcp_golden_gate(check_path: bool = False) -> dict:
    """Machine-readable golden-quality gate for the radia-mcp server fleet.

    This is intentionally lightweight: it audits catalog, discovery, related
    links, external-MCP boundaries, and optional entry-point visibility without
    launching every server selftest. For the full gate, run the commands
    returned under ``full_gate_commands``.
    """
    servers = catalog.list_all()
    catalog_names = set(catalog.CATALOG)
    external_names = set(catalog.EXTERNAL_PACKAGES)
    checks: list[dict] = []

    def add(name: str, ok: bool, detail: str, evidence: dict | None = None) -> None:
        row = {"name": name, "ok": bool(ok), "detail": detail}
        if evidence is not None:
            row["evidence"] = evidence
        checks.append(row)

    missing_required = []
    missing_selftest = []
    bad_entrypoints = []
    for entry in servers:
        name = entry["name"]
        for key in ("subpackage", "entry_point", "description", "primary_tools", "tags"):
            if key not in entry or entry[key] in ("", [], None):
                missing_required.append(f"{name}.{key}")
        entry_point = entry.get("entry_point", "")
        if not str(entry_point).startswith("mcp-server-"):
            bad_entrypoints.append(f"{name}: {entry_point}")
        expected_selftest = f"{entry_point} --selftest"
        if entry.get("selftest_command") != expected_selftest:
            missing_selftest.append(f"{name}: {entry.get('selftest_command')}")

    add(
        "catalog_required_fields",
        not missing_required,
        "Every cataloged server declares subpackage, entry point, description, tools, and tags.",
        {"missing": missing_required[:10], "n_missing": len(missing_required)},
    )
    add(
        "entrypoint_naming",
        not bad_entrypoints,
        "Every public radia-mcp server uses the mcp-server-* console-script convention.",
        {"bad": bad_entrypoints[:10], "n_bad": len(bad_entrypoints)},
    )
    add(
        "selftest_command_coverage",
        not missing_selftest,
        "Every catalog entry exposes a lightweight --selftest command.",
        {"bad": missing_selftest[:10], "n_bad": len(missing_selftest)},
    )

    missing_related = []
    asymmetric = []
    for name, info in catalog.CATALOG.items():
        related = set(info.get("related", []))
        for other in related:
            if other not in catalog_names and other not in external_names:
                missing_related.append(f"{name} -> {other}")
            elif other in catalog_names and name not in set(catalog.CATALOG[other].get("related", [])):
                asymmetric.append(f"{name} -> {other}")
    add(
        "related_links_resolve",
        not missing_related,
        "All related links point to a cataloged server or declared external MCP package.",
        {"missing": missing_related[:10], "n_missing": len(missing_related)},
    )
    add(
        "internal_related_links_bidirectional",
        not asymmetric,
        "Internal related links are bidirectional for reliable agent navigation.",
        {"asymmetric": asymmetric[:10], "n_asymmetric": len(asymmetric)},
    )

    external_entrypoints = {
        name: info.get("entry_point", "")
        for name, info in catalog.EXTERNAL_PACKAGES.items()
    }
    add(
        "optuna_external_boundary",
        "optuna" not in catalog_names and "optuna-mcp" in external_names,
        "Optuna Study/Trial operation stays in the official external optuna-mcp server.",
        {
            "catalog_has_optuna": "optuna" in catalog_names,
            "external_has_optuna_mcp": "optuna-mcp" in external_names,
            "optuna_entry_point": external_entrypoints.get("optuna-mcp"),
        },
    )

    if check_path:
        missing_cli = [
            entry["entry_point"]
            for entry in servers
            if shutil.which(entry["entry_point"]) is None
        ]
        add(
            "entrypoints_on_path",
            not missing_cli,
            "All cataloged console scripts are visible on PATH in this environment.",
            {"missing": missing_cli[:10], "n_missing": len(missing_cli)},
        )
    else:
        add(
            "entrypoints_on_path",
            True,
            "Skipped PATH probing; pass check_path=True for local editable-install verification.",
            {"skipped": True},
        )

    all_passed = all(row["ok"] for row in checks)
    return {
        "level": "golden" if all_passed else "needs_attention",
        "all_passed": all_passed,
        "n_servers": len(servers),
        "n_external_packages": len(external_names),
        "checks": checks,
        "full_gate_commands": [
            "python tools/policy_lint.py --tracked-only",
            "python scripts/gen_tools_doc.py --check",
            "pytest tests/test_meta_health.py tests/test_each_server_selftest.py tests/test_policy_lint.py",
        ],
    }


# ============================================================
# Bug-pattern catalog -- learned anti-patterns the lab has hit
# (call BEFORE writing new code in the affected area).
# ============================================================

@mcp.tool()
def bug_patterns_lookup(topic: str = "",
                        severity: str = "",
                        recent_days: int = 0) -> dict:
    """Query the learned bug-pattern catalog.

    USE THIS BEFORE writing new panel / calc_*.py / release / Cubit
    plugin / license-handling code -- the catalog records every
    bug class the lab has hit in real incidents, with the prevention
    rule for each.  Saves repeating the same mistakes.

    Args:
        topic: substring filter against tags + id + title.
            Examples: "panel", "release", "cubit-license", "taskmanager".
            Empty = no filter.
        severity: "high" | "medium" | "low".  Empty = any.
        recent_days: only patterns observed within the last N days.
            0 = no age filter.

    Returns:
        ``{"matched": N, "patterns": [...]}`` where each entry has
        ``id`` / ``title`` / ``severity`` / ``what`` / ``root_cause``
        / ``prevention`` / ``detection`` / ``related``.

    Example use cases:
      - Before editing radia_<topic>.py: ``bug_patterns_lookup("panel")``
      - Before tagging a release: ``bug_patterns_lookup("release")``
      - Before touching Cubit licensing: ``bug_patterns_lookup("cubit-license")``
      - Weekly digest: ``bug_patterns_lookup(recent_days=14)``
    """
    matched = bug_patterns.lookup(
        topic=topic or None,
        severity=severity or None,
        recent_days=recent_days or None,
    )
    return {
        "matched": len(matched),
        "patterns": matched,
        "hint": ("Read 'prevention' first; that's the rule to follow.  "
                 "Check 'detection' to know which test/audit catches "
                 "the regression if you forget."),
    }


@mcp.tool()
def bug_patterns_stats() -> dict:
    """Counts of catalogued bug patterns by severity + topic.

    Quick health check for the bug catalog itself.  Use to see
    which areas of the codebase have the most learned anti-patterns
    (= the areas that historically bite the most)."""
    return bug_patterns.stats()


# ============================================================
# Self-introspection (uniform with other radia_mcp servers)
# ============================================================

register_status_tool(
    mcp,
    server_name="mcp-server-radia-meta",
    description="Cross-server catalog of all radia_mcp.* servers "
                "(★ recommended first call when picking a tool)",
    subpackage="radia_mcp.meta",
    related_servers=["literature-index"],
    optional_deps=[],
)


# ============================================================
# Entry point
# ============================================================

def main():
    if "--selftest" in sys.argv:
        print("radia-mcp meta server self-test:")
        ov = radia_mcp_overview()
        print(f"  servers cataloged: {ov['n_servers']}")
        print(f"  tags available:    {len(ov['tags_available'])}")
        # Quick health (lighter — don't fail selftest if some
        # subpackages are unbuildable on this machine)
        h = radia_mcp_health()
        print(f"  health: {h['n_servers_healthy']}/{h['n_servers_total']} import OK")
        if not h["all_healthy"]:
            unhealthy = [r for r in h["results"] if not r["import_ok"]]
            print(f"  unhealthy subpackages ({len(unhealthy)}):")
            for r in unhealthy[:5]:
                print(f"    - {r['name']}: {r.get('error', '?')[:60]}")
            if len(unhealthy) > 5:
                print(f"    ... and {len(unhealthy) - 5} more")
        # Tag query
        opt = radia_mcp_by_tag("optimization")
        print(f"  servers tagged 'optimization': {opt['n_matches']}")
        # Related query
        rel = radia_mcp_related("bayesian-opt")
        print(f"  related to bayesian-opt: {[r['name'] for r in rel['related']]}")
        gate = radia_mcp_golden_gate()
        print(f"  golden gate: {gate['level']} ({len(gate['checks'])} checks)")
        assert gate["all_passed"], gate
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
