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
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
