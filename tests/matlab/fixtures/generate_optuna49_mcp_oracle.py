"""Generate MATLAB parity fixtures through the official Optuna MCP server.

The seeded sampler oracle lives in ``generate_optuna49_oracle.py`` because
optuna-mcp 0.2.0 does not expose a sampler seed.  This fixture exercises the
public MCP Study/Trial surface over a real stdio session instead of importing
the server implementation in-process.
"""

from __future__ import annotations

import csv
import importlib.metadata
import io
import json
import shutil
from pathlib import Path
from typing import Any

import anyio
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

EXPECTED_OPTUNA_VERSION = "4.9.0"
EXPECTED_OPTUNA_MCP_VERSION = "0.2.0"


def _float_distribution(low: float, high: float) -> dict[str, object]:
    return {
        "name": "FloatDistribution",
        "attributes": {"step": None, "low": low, "high": high, "log": False},
    }


def _int_distribution(low: int, high: int) -> dict[str, object]:
    return {
        "name": "IntDistribution",
        "attributes": {"step": 1, "low": low, "high": high, "log": False},
    }


def _categorical_distribution(choices: list[str]) -> dict[str, object]:
    return {
        "name": "CategoricalDistribution",
        "attributes": {"choices": choices},
    }


def _trial_to_add(
    *,
    x: float,
    values: list[float] | None,
    state: str,
    user_attrs: dict[str, Any] | None = None,
    system_attrs: dict[str, Any] | None = None,
) -> dict[str, object]:
    return {
        "params": {"x": x},
        "distributions": {"x": _float_distribution(-1.0, 2.0)},
        "values": values,
        "state": state,
        "user_attrs": user_attrs,
        "system_attrs": system_attrs,
    }


async def _call(
    session: ClientSession, name: str, arguments: dict[str, Any] | None = None
) -> dict[str, Any] | list[Any] | str | None:
    result = await session.call_tool(name, arguments or {})
    if result.isError:
        details = "\n".join(
            str(getattr(content, "text", content)) for content in result.content
        )
        raise RuntimeError(f"optuna-mcp tool {name!r} failed: {details}")
    if result.structuredContent is not None:
        return result.structuredContent
    if not result.content:
        return None
    return str(getattr(result.content[0], "text", result.content[0]))


async def _call_expect_error(
    session: ClientSession, name: str, arguments: dict[str, Any] | None = None
) -> bool:
    result = await session.call_tool(name, arguments or {})
    return bool(result.isError)


def _stable_trial_summary(csv_payload: str) -> dict[str, object]:
    prefix = "Trials: \n"
    if not csv_payload.startswith(prefix):
        raise RuntimeError("optuna-mcp get_trials returned an unexpected payload.")
    rows = list(csv.DictReader(io.StringIO(csv_payload[len(prefix) :])))

    def column(name: str, default: str = "") -> list[str]:
        return [row.get(name, default) or default for row in rows]

    return {
        "numbers": [int(row["number"]) for row in rows],
        "states": column("state"),
        "value_present": [bool(row.get("value")) for row in rows],
        "values": [float(row["value"]) if row.get("value") else 0.0 for row in rows],
        "params_x": [float(row["params_x"]) for row in rows],
        "params_n": [
            int(float(value)) if value else 0 for value in column("params_n")
        ],
        "params_mode": column("params_mode"),
        "user_attrs_source": column("user_attrs_source"),
    }


async def _build_oracle() -> dict[str, object]:
    optuna_version = importlib.metadata.version("optuna")
    optuna_mcp_version = importlib.metadata.version("optuna-mcp")
    if optuna_version != EXPECTED_OPTUNA_VERSION:
        raise RuntimeError(
            f"Expected optuna=={EXPECTED_OPTUNA_VERSION}, found {optuna_version}."
        )
    if optuna_mcp_version != EXPECTED_OPTUNA_MCP_VERSION:
        raise RuntimeError(
            "Expected optuna-mcp=="
            f"{EXPECTED_OPTUNA_MCP_VERSION}, found {optuna_mcp_version}."
        )

    command = shutil.which("optuna-mcp")
    if command is None:
        raise RuntimeError("The official optuna-mcp entry point is not on PATH.")

    server = StdioServerParameters(command=command, args=[])
    async with stdio_client(server) as (read, write):  # noqa: SIM117
        async with ClientSession(read, write) as session:
            initialized = await session.initialize()
            tool_result = await session.list_tools()
            tools = {tool.name: tool for tool in tool_result.tools}
            sampler_schema = tools["set_sampler"].inputSchema
            sampler_properties = sorted(sampler_schema.get("properties", {}))

            single_name = "radia-matlab-mcp-single"
            create_single = await _call(
                session,
                "create_study",
                {"study_name": single_name, "directions": ["minimize"]},
            )
            sampler = await _call(session, "set_sampler", {"name": "RandomSampler"})
            single_directions = await _call(session, "get_directions")
            search_space = {
                "x": _float_distribution(1.25, 1.25),
                "n": _int_distribution(3, 3),
                "mode": _categorical_distribution(["A"]),
            }
            asked = await _call(session, "ask", {"search_space": search_space})
            await _call(
                session,
                "set_trial_user_attr",
                {"trial_number": 0, "key": "source", "value": "official-mcp"},
            )
            user_attrs = await _call(
                session, "get_trial_user_attrs", {"trial_number": 0}
            )
            told = await _call(
                session, "tell", {"trial_number": 0, "values": 2.5}
            )
            added_one = await _call(
                session,
                "add_trial",
                {
                    "trial": _trial_to_add(
                        x=-0.5,
                        values=[0.5],
                        state="COMPLETE",
                        user_attrs={"source": "archive"},
                        system_attrs={"origin": "official-mcp"},
                    )
                },
            )
            added_many = await _call(
                session,
                "add_trials",
                {
                    "trials": [
                        _trial_to_add(x=0.0, values=None, state="PRUNED"),
                        _trial_to_add(x=0.75, values=None, state="FAIL"),
                    ]
                },
            )
            single_best = await _call(session, "best_trial")
            single_trials_csv = await _call(session, "get_trials")
            if not isinstance(single_trials_csv, str):
                raise TypeError("get_trials did not return its documented CSV text.")

            multi_name = "radia-matlab-mcp-multi"
            create_multi = await _call(
                session,
                "create_study",
                {
                    "study_name": multi_name,
                    "directions": ["minimize", "maximize"],
                },
            )
            metric_names = await _call(
                session, "set_metric_names", {"metric_names": ["loss", "score"]}
            )
            metric_names_read = await _call(session, "get_metric_names")
            multi_directions = await _call(session, "get_directions")
            multi_asked = await _call(
                session,
                "ask",
                {"search_space": {"x": _float_distribution(0.0, 0.0)}},
            )
            multi_told = await _call(
                session, "tell", {"trial_number": 0, "values": [1.0, 1.0]}
            )
            multi_added = await _call(
                session,
                "add_trials",
                {
                    "trials": [
                        _trial_to_add(x=0.5, values=[0.5, 0.5], state="COMPLETE"),
                        _trial_to_add(x=2.0, values=[2.0, 2.0], state="COMPLETE"),
                        _trial_to_add(x=1.5, values=[1.5, 0.0], state="COMPLETE"),
                    ]
                },
            )
            pareto = await _call(session, "best_trials")
            best_trial_is_error = await _call_expect_error(session, "best_trial")
            all_studies = await _call(session, "get_all_study_names")

            return {
                "schema": "radia.test.optuna-upstream-mcp-oracle.v1",
                "transport": "stdio",
                "optuna_version": optuna_version,
                "optuna_mcp_version": optuna_mcp_version,
                "mcp_server_name": initialized.serverInfo.name,
                "mcp_server_reported_version": initialized.serverInfo.version,
                "tools": sorted(tools),
                "set_sampler_arguments": sampler_properties,
                "sampler_seed_supported": "seed" in sampler_properties,
                "single": {
                    "study_name": single_name,
                    "create_study": create_single,
                    "set_sampler": sampler,
                    "directions": single_directions,
                    "ask": asked,
                    "user_attrs": user_attrs,
                    "tell": told,
                    "add_trial": added_one,
                    "add_trials": added_many,
                    "best_trial": single_best,
                    "trials": _stable_trial_summary(single_trials_csv),
                },
                "multi": {
                    "study_name": multi_name,
                    "create_study": create_multi,
                    "directions": multi_directions,
                    "set_metric_names": metric_names,
                    "get_metric_names": metric_names_read,
                    "ask": multi_asked,
                    "tell": multi_told,
                    "add_trials": multi_added,
                    "best_trials": pareto,
                    "best_trial_is_error": best_trial_is_error,
                },
                "all_study_names": all_studies,
            }


def build_oracle() -> dict[str, object]:
    return anyio.run(_build_oracle)


def main() -> None:
    destination = Path(__file__).with_name("optuna49_mcp_oracle.json")
    destination.write_text(
        json.dumps(build_oracle(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(destination)


if __name__ == "__main__":
    main()
