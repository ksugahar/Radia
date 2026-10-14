# Cubit-owned implementation; intentionally maintained independently of Radia MCP.
# Derived support retains BSD-3-Clause terms: see LICENSE-BSD-3-Clause.txt.
"""Coarse MCP entry points for large families of validation tools.

Radia accumulated many one-purpose validation gates.  Registering every gate
as a top-level MCP tool makes ``tools/list`` large and forces the client model
to consider hundreds of schemas on every turn.  ``CoarseToolRegistry`` keeps
the Python callables directly testable while exposing two stable MCP tools:
a searchable catalog and one dictionary-driven runner.

The default ``core`` profile exposes only those coarse entry points.  Set
``CUBIT_MCP_TOOL_PROFILE=full`` or pass ``--tool-profile full`` to retain the
historical individual tools while migrating existing clients.
"""

from __future__ import annotations

import asyncio
import inspect
import os
import sys
from collections.abc import Callable
from dataclasses import dataclass
from functools import cached_property
from typing import Any

from mcp.server.fastmcp.utilities.func_metadata import func_metadata
from pydantic import ValidationError

from cubit_mesh_export.mcp._support.server_hardening import ANN_READONLY, ANN_WRITES, error_payload

_VALID_PROFILES = frozenset({"core", "full"})


def selected_tool_profile(argv: list[str] | None = None) -> str:
    """Return the requested MCP tool-surface profile."""

    args = list(sys.argv[1:] if argv is None else argv)
    requested = os.environ.get("CUBIT_MCP_TOOL_PROFILE", "core")
    for index, arg in enumerate(args):
        if arg.startswith("--tool-profile="):
            requested = arg.split("=", 1)[1]
        elif arg == "--tool-profile" and index + 1 < len(args):
            requested = args[index + 1]
    profile = requested.strip().lower()
    if profile not in _VALID_PROFILES:
        choices = ", ".join(sorted(_VALID_PROFILES))
        raise ValueError(f"invalid MCP tool profile {requested!r}; choose {choices}")
    return profile


def _first_line(doc: str | None) -> str:
    for line in (doc or "").splitlines():
        if line.strip():
            return line.strip()
    return "(no description)"


def _parameter_summary(function: Callable[..., Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for parameter in inspect.signature(function).parameters.values():
        annotation = parameter.annotation
        annotation_name = (
            "Any" if annotation is inspect.Parameter.empty else str(annotation)
        )
        default = parameter.default
        default_name = (
            "required" if default is inspect.Parameter.empty else repr(default)
        )
        rows.append(
            {
                "name": parameter.name,
                "type": annotation_name,
                "default": default_name,
            }
        )
    return rows


@dataclass(frozen=True)
class _Entry:
    name: str
    function: Callable[..., Any]
    description: str
    decorator_args: tuple[Any, ...]
    decorator_kwargs: dict[str, Any]

    @cached_property
    def argument_metadata(self):
        # Match SDK input semantics without converting the domain return value.
        # Build lazily so catalog discovery does not compile every gate schema.
        return func_metadata(self.function, structured_output=False)


class CoarseToolRegistry:
    """Group many related callables behind catalog and run MCP tools."""

    def __init__(
        self,
        mcp: Any,
        *,
        namespace: str,
        category: str = "validation",
        profile: str | None = None,
        min_group_size: int = 3,
    ) -> None:
        self.mcp = mcp
        self.namespace = namespace.replace("-", "_")
        self.category = category.replace("-", "_")
        self.profile = profile or selected_tool_profile()
        if self.profile not in _VALID_PROFILES:
            raise ValueError(f"unsupported tool profile: {self.profile}")
        if min_group_size < 1:
            raise ValueError("min_group_size must be positive")
        self.min_group_size = min_group_size
        self._entries: dict[str, _Entry] = {}
        self._installed = False

    def tool(self, *decorator_args: Any, **decorator_kwargs: Any):
        """Record a callable and expose it individually in ``full`` mode."""

        def decorate(function: Callable[..., Any]) -> Callable[..., Any]:
            name = str(decorator_kwargs.get("name") or function.__name__)
            if name in self._entries:
                raise ValueError(f"duplicate grouped MCP tool: {name}")
            self._entries[name] = _Entry(
                name=name,
                function=function,
                description=_first_line(function.__doc__),
                decorator_args=tuple(decorator_args),
                decorator_kwargs=dict(decorator_kwargs),
            )
            if self.profile == "full":
                self.mcp.tool(*decorator_args, **decorator_kwargs)(function)
            return function

        return decorate

    def install(self) -> None:
        """Register the stable catalog/runner pair once."""

        if self._installed:
            raise RuntimeError("coarse tool registry already installed")
        self._installed = True
        if len(self._entries) < self.min_group_size:
            if self.profile == "core":
                for entry in self._entries.values():
                    self.mcp.tool(
                        *entry.decorator_args,
                        **entry.decorator_kwargs,
                    )(entry.function)
            self._append_summary()
            return
        catalog_name = f"{self.namespace}_{self.category}_catalog"
        run_name = f"{self.namespace}_{self.category}_run"

        def catalog(query: str = "", limit: int = 50) -> dict[str, Any]:
            """Search grouped operations without loading individual schemas."""

            needle = query.casefold().strip()
            bounded_limit = max(1, min(int(limit), 200))
            matches = [
                entry
                for entry in self._entries.values()
                if not needle
                or needle in entry.name.casefold()
                or needle in entry.description.casefold()
            ]
            selected = sorted(matches, key=lambda item: item.name)[:bounded_limit]
            return {
                "namespace": self.namespace,
                "category": self.category,
                "profile": self.profile,
                "matched": len(matches),
                "returned": len(selected),
                "operations": [
                    {
                        "name": entry.name,
                        "description": entry.description,
                        "parameters": _parameter_summary(entry.function),
                    }
                    for entry in selected
                ],
                "run_tool": run_name,
            }

        async def run(name: str, arguments: dict[str, Any] | None = None) -> Any:
            """Run one grouped operation with named arguments."""

            entry = self._entries.get(name)
            if entry is None:
                available = ", ".join(sorted(self._entries)[:12])
                return error_payload(
                    run_name,
                    f"unknown {self.category} operation {name!r}; "
                    f"query {catalog_name} first (first names: {available})",
                    kind="input",
                )
            if arguments is not None and not isinstance(arguments, dict):
                return error_payload(run_name, "arguments must be a dictionary",
                                     kind="input", hint=f"Query {catalog_name} for signatures.")
            kwargs = {} if arguments is None else arguments
            # Retain rejection of unknown keywords (the SDK model ignores extras).
            # Bind only, never execute a tool inside the input-error boundary.
            try:
                inspect.signature(entry.function).bind(**kwargs)
            except TypeError as exc:
                return error_payload(run_name, str(exc), kind="input",
                                     hint=f"Query {catalog_name} for {name}'s signature.")
            metadata = entry.argument_metadata
            try:
                parsed = metadata.pre_parse_json(kwargs)
                kwargs = metadata.arg_model.model_validate(parsed).model_dump_one_level()
            except ValidationError as exc:
                return error_payload(run_name, str(exc), kind="input",
                                     hint=f"Query {catalog_name} for {name}'s signature.")
            if inspect.iscoroutinefunction(entry.function):
                result = await entry.function(**kwargs)
            else:
                result = await asyncio.to_thread(entry.function, **kwargs)
                if inspect.isawaitable(result):
                    result = await result
            return result

        catalog.__name__ = catalog_name
        catalog.__doc__ = (
            f"Search the {self.namespace} {self.category} operation catalog. "
            f"Use the returned name with {run_name}."
        )
        run.__name__ = run_name
        run.__doc__ = (
            f"Run one {self.namespace} {self.category} operation. Pass the "
            f"operation name and a dictionary of named arguments. Query "
            f"{catalog_name} first for signatures."
        )
        self.mcp.tool(
            title=f"{self.namespace} {self.category} catalog",
            annotations=ANN_READONLY,
        )(catalog)
        self.mcp.tool(
            title=f"Run {self.namespace} {self.category}",
            # The selected operation may execute a solver or write an
            # artifact. The dispatcher cannot truthfully advertise itself as
            # read-only even when many individual gates are pure checks.
            annotations=ANN_WRITES,
        )(run)

        self._append_summary()

    def _append_summary(self) -> None:
        groups = getattr(self.mcp, "_cubit_tool_groups", None)
        if groups is None:
            groups = []
            setattr(self.mcp, "_cubit_tool_groups", groups)
        groups.append(self.summary())

    def summary(self) -> dict[str, Any]:
        """Return profile metadata for status/introspection tools."""

        grouped = len(self._entries) >= self.min_group_size
        return {
            "namespace": self.namespace,
            "category": self.category,
            "profile": self.profile,
            "mode": "catalog-runner" if grouped else "direct-small-group",
            "candidate_operations": len(self._entries),
            "grouped_operations": len(self._entries) if grouped else 0,
            "individual_tools_exposed": (
                len(self._entries) if self.profile == "full" or not grouped else 0
            ),
        }
