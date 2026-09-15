"""Small lazy-call adapter for optional and validation-only subsystems."""

from __future__ import annotations

from importlib import import_module
from importlib.util import resolve_name
from typing import Any, Callable


def lazy_callable(
    module_name: str,
    attribute: str,
    package: str,
) -> Callable[..., Any]:
    """Return a callable that resolves its target only when invoked.

    The target attribute is looked up on every call.  Python still caches the
    imported module, while Radia's MCP hot-reload support can replace the
    attribute without leaving this adapter bound to the old function object.
    """

    resolved_module = resolve_name(module_name, package)

    def call(*args: Any, **kwargs: Any) -> Any:
        target = getattr(import_module(resolved_module), attribute)
        return target(*args, **kwargs)

    call.__name__ = attribute
    call.__qualname__ = attribute
    call.__doc__ = f"Lazy proxy for {resolved_module}.{attribute}."
    setattr(call, "__radia_lazy_target__", f"{resolved_module}:{attribute}")
    return call
