"""Inventory the pinned Optuna 4.9.0 public API for MATLAB closure work."""

from __future__ import annotations

import importlib
import inspect
import json
import re
from pathlib import Path
from types import ModuleType
from typing import Any

import optuna


EXPECTED_VERSION = "4.9.0"
MODULE_NAMES = [
    "optuna",
    "optuna.artifacts",
    "optuna.distributions",
    "optuna.exceptions",
    "optuna.importance",
    "optuna.integration",
    "optuna.logging",
    "optuna.pruners",
    "optuna.samplers",
    "optuna.samplers.nsgaii",
    "optuna.search_space",
    "optuna.storages",
    "optuna.storages.journal",
    "optuna.study",
    "optuna.terminator",
    "optuna.trial",
    "optuna.visualization",
    "optuna.visualization.matplotlib",
]
PUBLIC_WITHOUT_ALL = {
    "optuna.distributions": [
        "BaseDistribution",
        "CategoricalChoiceType",
        "CategoricalDistribution",
        "DISTRIBUTION_CLASSES",
        "DiscreteUniformDistribution",
        "FloatDistribution",
        "IntDistribution",
        "IntLogUniformDistribution",
        "IntUniformDistribution",
        "LogUniformDistribution",
        "UniformDistribution",
        "check_distribution_compatibility",
        "distribution_to_json",
        "json_to_distribution",
    ],
    "optuna.exceptions": [
        "CLIUsageError",
        "DuplicatedStudyError",
        "ExperimentalWarning",
        "OptunaError",
        "StorageInternalError",
        "TrialPruned",
        "UpdateFinishedTrialError",
    ],
}
PUBLIC_EXTRAS = {
    "optuna.logging": [
        "create_default_formatter",
        "disable_default_handler",
        "disable_propagation",
        "enable_default_handler",
        "enable_propagation",
        "get_logger",
        "get_verbosity",
        "set_verbosity",
    ],
}


def _stable_signature(value: Any) -> str | None:
    try:
        signature = str(inspect.signature(value))
    except (TypeError, ValueError):
        return None
    return re.sub(r"0x[0-9A-Fa-f]+", "0x...", signature)


def _kind(value: Any) -> str:
    if inspect.ismodule(value):
        return "module"
    if inspect.isclass(value):
        return "class"
    if inspect.isfunction(value) or inspect.isbuiltin(value):
        return "function"
    return "constant-or-type-alias"


def _class_members(value: type[Any]) -> list[dict[str, object]]:
    names = {
        name for base in value.__mro__ for name in vars(base) if not name.startswith("_")
    }
    names.update(name for name in getattr(value, "__annotations__", {}) if not name.startswith("_"))
    names.update(getattr(value, "__members__", {}).keys())
    members: list[dict[str, object]] = []
    for name in sorted(names):
        try:
            member = inspect.getattr_static(value, name)
            resolved = getattr(value, name)
        except Exception:
            members.append({"kind": "descriptor", "name": name, "signature": None})
            continue
        if isinstance(member, property):
            kind = "property"
        elif isinstance(member, (classmethod, staticmethod)) or callable(resolved):
            kind = "method"
        else:
            kind = "constant-or-field"
        members.append(
            {"kind": kind, "name": name, "signature": _stable_signature(resolved)}
        )
    return members


def _public_names(module: ModuleType) -> list[str]:
    configured = PUBLIC_WITHOUT_ALL.get(module.__name__)
    if configured is not None:
        names = set(configured)
    else:
        names = set(getattr(module, "__all__", []))
    names.update(PUBLIC_EXTRAS.get(module.__name__, []))
    return sorted(names)


def _module_inventory(module_name: str) -> dict[str, object]:
    module = importlib.import_module(module_name)
    symbols: list[dict[str, object]] = []
    for name in _public_names(module):
        if module_name == "optuna.integration":
            # The integration namespace is deliberately lazy. Resolving these
            # names imports optional third-party stacks and can block or make
            # the inventory depend on the host environment. Their exported
            # names are the pinned public contract; inventory them unresolved.
            symbols.append(
                {
                    "kind": "lazy-optional",
                    "members": [],
                    "name": name,
                    "resolvable": False,
                    "signature": None,
                }
            )
            continue
        try:
            value = getattr(module, name)
        except Exception:
            symbols.append(
                {
                    "kind": "lazy-optional",
                    "members": [],
                    "name": name,
                    "resolvable": False,
                    "signature": None,
                }
            )
            continue
        kind = _kind(value)
        symbols.append(
            {
                "kind": kind,
                "members": _class_members(value) if kind == "class" else [],
                "name": name,
                "resolvable": True,
                "signature": _stable_signature(value),
            }
        )
    return {"module": module_name, "symbols": symbols}


def build_inventory() -> dict[str, object]:
    if optuna.__version__ != EXPECTED_VERSION:
        raise RuntimeError(
            f"Expected optuna=={EXPECTED_VERSION}, found {optuna.__version__}."
        )
    modules = [_module_inventory(name) for name in MODULE_NAMES]
    symbol_count = sum(len(module["symbols"]) for module in modules)
    member_count = sum(
        len(symbol["members"])
        for module in modules
        for symbol in module["symbols"]
    )
    return {
        "schema": "radia.test.optuna49-public-api.v1",
        "optuna_version": optuna.__version__,
        "module_count": len(modules),
        "symbol_count": symbol_count,
        "class_member_count": member_count,
        "modules": modules,
    }


def main() -> None:
    destination = Path(__file__).with_name("optuna49_public_api.json")
    destination.write_text(
        json.dumps(build_inventory(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(destination)


if __name__ == "__main__":
    main()
