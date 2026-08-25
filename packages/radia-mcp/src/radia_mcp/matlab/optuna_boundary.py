"""Ownership and routing contract for upstream Optuna MCP and radia-optuna.

This module is deliberately declarative.  It does not import Optuna, proxy the
official MCP server, or recreate any of its Study/Trial tools.  radia-mcp owns
only the MATLAB/Simulink adaptation and differential-oracle boundary.
"""

from __future__ import annotations

from typing import Any


_TOPICS = (
    "overview",
    "shared",
    "matlab",
    "differential",
    "composition",
    "stewardship",
)


def matlab_optuna_mcp_route(topic: str = "overview") -> dict[str, Any]:
    """Return the checked owner for an Optuna/MATLAB MCP capability.

    The official server's live ``tools/list`` is authoritative for shared
    Optuna behavior.  This contract records the verified 4.9.0/0.2.0 snapshot
    but never shadows those public tools inside radia-mcp.
    """
    key = str(topic or "overview").strip().lower().replace("-", "_")
    aliases = {
        "upstream": "shared",
        "official": "shared",
        "radia": "matlab",
        "simulink": "matlab",
        "oracle": "differential",
        "routing": "composition",
        "license": "stewardship",
        "trademark": "stewardship",
        "upstream_stewardship": "stewardship",
    }
    key = aliases.get(key, key)
    if key not in _TOPICS:
        raise ValueError(f"unknown topic {topic!r}; expected {list(_TOPICS)}")

    shared = {
        "owner": "optuna/optuna-mcp",
        "server": "optuna-mcp",
        "authority": "the official server's live MCP tools/list response",
        "capabilities": [
            "Study creation and discovery",
            "ask/tell and Trial state operation",
            "Study and Trial attributes and metric names",
            "best-trial and Pareto-trial queries",
            "Optuna visualization and Dashboard launch",
            "sampler selection exposed by the official MCP tool contract",
        ],
        "verified_snapshot": {
            "optuna": "4.9.0",
            "optuna_mcp": "0.2.0",
            "transport": "stdio",
            "fixture": "tests/matlab/fixtures/optuna49_mcp_oracle.json",
            "sampler_seed_exposed": False,
        },
        "radia_policy": (
            "Do not proxy, rename, or reimplement an operation present in the "
            "official server's live tools/list response."
        ),
    }
    matlab = {
        "owner": "radia-mcp/radia-matlab",
        "server": "mcp-server-radia-matlab",
        "distribution": "radia-optuna",
        "execution_owner": "MathWorks official MATLAB MCP Server",
        "tools": [
            "matlab_optuna_mcp_route",
            "matlab_optuna_simulink_contract",
            "matlab_optimize_build",
            "matlab_optimize_resume",
        ],
        "capabilities": [
            "MATLAB namespace and installed radia-optuna distribution contract",
            "table/MAT progress persistence and resume code generation",
            "Simulink block, Scope/XY monitor, progress, and failure telemetry",
            "MATLAB-only parallel trial execution with client-owned Study tables",
            "required standalone optuna_mex lifecycle and numerical kernels",
            "Radia CAE, LTspice, and sheet-metal trial artifact adapters",
        ],
        "does_not_own": [
            "generic Python Optuna Study/Trial operation",
            "official Optuna MCP visualization or Dashboard tools",
            "a second Optuna MCP server or optuna-mcp proxy",
        ],
    }
    differential = {
        "owner": "radia-mcp test and compatibility layer",
        "behavioral_oracle": "optuna==4.9.0",
        "public_mcp_oracle": "optuna/optuna-mcp over a real MCP transport",
        "seeded_numeric_route": (
            "Run pinned upstream Optuna directly because optuna-mcp 0.2.0 "
            "set_sampler does not expose a seed."
        ),
        "matlab_difference_scope": [
            "table/MAT storage",
            "Simulink execution and monitoring",
            "standalone optuna_mex",
            "MATLAB parallel execution",
            "MATLAB/Radia CAE artifact contracts",
        ],
        "fixtures": [
            "tests/matlab/fixtures/optuna49_oracle.json",
            "tests/matlab/fixtures/optuna49_mcp_oracle.json",
            "tests/matlab/fixtures/optuna49_public_api.json",
        ],
        "disagreement_policy": (
            "Regenerate the pinned upstream fixture first, then fix MATLAB; "
            "never make MATLAB output the compatibility truth."
        ),
    }
    composition = {
        "order": [
            "Discover the official optuna-mcp live tools/list response.",
            "Route every supported shared operation to optuna-mcp.",
            "Use radia-mcp only for MATLAB/Simulink adaptation or a recorded gap.",
            "Execute generated MATLAB through the official MathWorks MATLAB MCP Server.",
            "Use pinned upstream Optuna fixtures for differential verification.",
        ],
        "no_duplicate_tool_names": True,
        "no_optuna_runtime_dependency_in_radia_mcp": True,
        "no_python_per_simulink_step": True,
    }
    stewardship = {
        "project_relationship": (
            "radia-optuna is independent and unofficial; it is not affiliated "
            "with, sponsored by, or endorsed by Preferred Networks, Inc. or "
            "the Optuna project."
        ),
        "trademark_attribution": (
            "Optuna, the Optuna logo and any related marks are trademarks of "
            "Preferred Networks, Inc."
        ),
        "visual_identity": "Do not use the Optuna logo or imply official status.",
        "upstream_licenses": [
            {
                "project": "Optuna",
                "license": "MIT",
                "copyright": "Copyright (c) 2018 Preferred Networks, Inc.",
                "url": "https://github.com/optuna/optuna/blob/master/LICENSE",
            },
            {
                "project": "optuna-mcp",
                "license": "MIT",
                "copyright": "Copyright (c) 2025 Preferred Networks, Inc.",
                "url": "https://github.com/optuna/optuna-mcp/blob/main/LICENSE",
            },
        ],
        "notice_policy": (
            "Preserve the applicable MIT copyright and permission notice in "
            "every copy or substantial portion of upstream software, and "
            "record copied or adapted source provenance."
        ),
        "distribution_notice": "radia_optuna/matlab/THIRD_PARTY_NOTICES.md",
        "upstream_runtime_bundled": False,
        "validation_operation": {
            "transport": "local stdio",
            "storage": (
                "a fresh per-run temporary SQLite database; C:/temp on LAB"
            ),
            "shared_or_production_storage": False,
            "dashboard_in_automated_tests": False,
            "live_tools_list": (
                "query deliberately for fixture regeneration, then use the "
                "checked fixture in routine tests"
            ),
            "automatic_upstream_issue_or_pr_creation": False,
        },
    }
    routes = {
        "shared": shared,
        "matlab": matlab,
        "differential": differential,
        "composition": composition,
        "stewardship": stewardship,
    }
    if key != "overview":
        return {
            "schema": "radia-mcp.matlab-optuna-mcp-route/v2",
            "status": "ready",
            "topic": key,
            "route": routes[key],
        }
    return {
        "schema": "radia-mcp.matlab-optuna-mcp-route/v2",
        "status": "ready",
        "policy": "upstream for shared behavior; radia-mcp for MATLAB differences",
        "topics": list(_TOPICS),
        "routes": routes,
    }


__all__ = ["matlab_optuna_mcp_route"]
