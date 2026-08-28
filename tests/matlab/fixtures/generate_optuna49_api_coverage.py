"""Compare the Optuna 4.9.0 public inventory with the MATLAB package surface."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
INVENTORY_PATH = Path(__file__).with_name("optuna49_public_api.json")
MATLAB_DIRECTORY = ROOT / "matlab" / "+radia" / "+optuna"
DESTINATION = ROOT / "matlab" / "optuna49_api_coverage.json"
FUNCTION_PATTERN = re.compile(
    r"(?m)^\s*function\s+(?:\[[^]]*\]|\w+)\s*=\s*(\w+)\s*\(|"
    r"^\s*function\s+(\w+)\s*\("
)
MODULE_EQUIVALENTS = {
    "artifacts": False,
    "distributions": True,
    "exceptions": True,
    "importance": True,
    "integration": False,
    "logging": False,
    "pruners": True,
    "samplers": True,
    "search_space": True,
    "storages": True,
    "study": True,
    "trial": True,
    "version": False,
    "visualization": False,
}
VERIFIED_FUNCTIONS = {
    "copy_study",
    "check_distribution_compatibility",
    "create_study",
    "create_trial",
    "delete_study",
    "get_all_study_names",
    "get_all_study_summaries",
    "get_param_importances",
    "intersection_search_space",
    "distribution_to_json",
    "json_to_distribution",
    "load_study",
}
VERIFIED_MEMBERS = {
    "BasePruner": {"prune"},
    "BaseTrial": {
        "datetime_start",
        "distributions",
        "number",
        "params",
        "report",
        "set_system_attr",
        "set_user_attr",
        "should_prune",
        "suggest_categorical",
        "suggest_discrete_uniform",
        "suggest_float",
        "suggest_int",
        "suggest_loguniform",
        "suggest_uniform",
        "system_attrs",
        "user_attrs",
    },
    "Study": {
        "add_trial",
        "add_trials",
        "ask",
        "best_params",
        "best_trial",
        "best_trials",
        "best_value",
        "direction",
        "directions",
        "enqueue_trial",
        "get_trials",
        "metric_names",
        "optimize",
        "set_metric_names",
        "set_system_attr",
        "set_user_attr",
        "stop",
        "system_attrs",
        "tell",
        "trials",
        "user_attrs",
    },
    "StudySummary": {"direction", "directions", "system_attrs"},
    "Trial": {
        "distributions",
        "number",
        "params",
        "report",
        "set_system_attr",
        "set_user_attr",
        "should_prune",
        "suggest_categorical",
        "suggest_discrete_uniform",
        "suggest_float",
        "suggest_int",
        "suggest_loguniform",
        "suggest_uniform",
        "system_attrs",
        "user_attrs",
    },
    "FixedTrial": {
        "distributions",
        "number",
        "params",
        "report",
        "set_system_attr",
        "set_user_attr",
        "should_prune",
        "suggest_categorical",
        "suggest_discrete_uniform",
        "suggest_float",
        "suggest_int",
        "suggest_loguniform",
        "suggest_uniform",
        "system_attrs",
        "user_attrs",
    },
    "FrozenTrial": {
        "datetime_complete",
        "datetime_start",
        "distributions",
        "duration",
        "intermediate_values",
        "last_step",
        "number",
        "params",
        "report",
        "set_system_attr",
        "set_user_attr",
        "should_prune",
        "state",
        "suggest_categorical",
        "suggest_discrete_uniform",
        "suggest_float",
        "suggest_int",
        "suggest_loguniform",
        "suggest_uniform",
        "system_attrs",
        "user_attrs",
        "value",
        "values",
    },
    "BestValueStagnationEvaluator": {"evaluate"},
    "Terminator": {"should_terminate"},
    "IntersectionSearchSpace": {"calculate"},
    "StudyDirection": {"NOT_SET", "MINIMIZE", "MAXIMIZE"},
    "TrialState": {
        "RUNNING",
        "COMPLETE",
        "PRUNED",
        "FAIL",
        "WAITING",
        "is_finished",
    },
}
PARTIALLY_VERIFIED_CLASSES = {
    "BasePruner",
    "BaseSampler",
    "BaseTrial",
    "BestValueStagnationEvaluator",
    "BruteForceSampler",
    "CategoricalDistribution",
    "CmaEsSampler",
    "FixedTrial",
    "FloatDistribution",
    "FrozenTrial",
    "GPSampler",
    "GridSampler",
    "HyperbandPruner",
    "IntDistribution",
    "UniformDistribution",
    "LogUniformDistribution",
    "DiscreteUniformDistribution",
    "IntUniformDistribution",
    "IntLogUniformDistribution",
    "IntersectionSearchSpace",
    "MaxTrialsCallback",
    "MedianPruner",
    "NopPruner",
    "NSGAIISampler",
    "NSGAIIISampler",
    "PartialFixedSampler",
    "PatientPruner",
    "PercentilePruner",
    "QMCSampler",
    "RandomSampler",
    "Study",
    "StudyDirection",
    "StudySummary",
    "SuccessiveHalvingPruner",
    "Terminator",
    "TerminatorCallback",
    "ThresholdPruner",
    "TPESampler",
    "Trial",
    "TrialState",
    "TrialPruned",
    "WilcoxonPruner",
    "BaseCrossover",
    "BLXAlphaCrossover",
    "SBXCrossover",
    "SPXCrossover",
    "UNDXCrossover",
    "UniformCrossover",
    "VSBXCrossover",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def _normalized(name: str) -> str:
    return re.sub(r"_", "", name).lower()


def _matlab_surface() -> tuple[set[str], dict[str, dict[str, set[str]]]]:
    names: set[str] = set()
    members: dict[str, dict[str, set[str]]] = {}
    superclasses: dict[str, str] = {}
    for path in sorted(MATLAB_DIRECTORY.rglob("*.m")):
        if "+internal" in path.parts:
            continue
        names.add(path.stem)
        source = path.read_text(encoding="utf-8")
        inheritance = re.search(
            r"(?m)^\s*classdef\s+(?:\([^)]*\)\s*)?(\w+)\s*<\s*([\w.]+)", source
        )
        if inheritance:
            superclasses[inheritance.group(1)] = inheritance.group(2).split(".")[-1]
        found = {
            candidate
            for groups in FUNCTION_PATTERN.findall(source)
            for candidate in groups
            if candidate
        }
        class_members = members.setdefault(path.stem, {})
        for name in found:
            class_members.setdefault(_normalized(name), set()).add(name)
        property_blocks = re.findall(
            r"(?ms)^\s*properties(?:\s*\([^)]*\))?\s*(.*?)^\s*end\s*$", source
        )
        for block in property_blocks:
            for line in block.splitlines():
                match = re.match(r"\s*([A-Za-z]\w*)\b", line)
                if match and not line.lstrip().startswith("%"):
                    name = match.group(1)
                    class_members.setdefault(_normalized(name), set()).add(name)
        enumeration_blocks = re.findall(
            r"(?ms)^\s*enumeration\s*(.*?)^\s*end\s*$", source
        )
        for block in enumeration_blocks:
            for line in block.splitlines():
                match = re.match(r"\s*([A-Za-z]\w*)\b", line)
                if match and not line.lstrip().startswith("%"):
                    name = match.group(1)
                    class_members.setdefault(_normalized(name), set()).add(name)
    _inherit_members(members, superclasses)
    return names, members


def _inherit_members(
    members: dict[str, dict[str, set[str]]], superclasses: dict[str, str]
) -> None:
    """Credit a class with what it inherits.

    MATLAB resolves BaseSampler.reseed_rng on every sampler, so a ledger that
    only looked inside each file reported the whole family as missing while
    the method was callable all along.
    """
    for name in list(members):
        seen = set()
        parent = superclasses.get(name)
        while parent and parent in members and parent not in seen:
            seen.add(parent)
            for normalized, actual in members[parent].items():
                members[name].setdefault(normalized, set()).update(actual)
            parent = superclasses.get(parent)


# Why a module is not held to the "port it" bar.  Anything absent from this
# table is `required`: it must be present AND oracle-mapped before
# full_compatibility_complete can be true.  A downgrade here is a design
# decision that has to name the thing that discharges it, so the ledger can
# never quietly excuse an unimplemented feature.
SCOPE_RULES = {
    "optuna.storages": {
        "scope": "bridged",
        "reason": (
            "MATLAB keeps its history in tables and a MAT-file rather than an "
            "optuna.storages backend. Interoperability is discharged by the "
            "explicit handoff document instead of by porting the backends."
        ),
        "discharged_by": [
            "radia.optuna.export_study",
            "radia.optuna.import_study",
            "radia_optuna.bridge",
        ],
    },
    "optuna.visualization": {
        "scope": "out-of-scope",
        "reason": (
            "Plotting is not a MATLAB Optuna concern; export the study and "
            "plot it from Python, or use MATLAB's own plotting on the tables."
        ),
        "discharged_by": ["radia_optuna.bridge"],
    },
    "optuna.integration": {
        "scope": "out-of-scope",
        "reason": (
            "Third-party Python framework callbacks (PyTorch, XGBoost, ...) "
            "have no MATLAB counterpart; run them on the bridged study."
        ),
        "discharged_by": ["radia_optuna.bridge"],
    },
    "optuna.artifacts": {
        "scope": "out-of-scope",
        "reason": (
            "Artifact stores are a Python-side file-management surface; "
            "Radia applications already own run.log / result.json artifacts."
        ),
        "discharged_by": [],
    },
    "optuna.logging": {
        "scope": "replaced",
        "reason": (
            "MATLAB reports through warning/error identifiers and disp, not "
            "through a Python logging hierarchy."
        ),
        "discharged_by": ["MATLAB warning and error identifiers"],
    },
    "optuna.exceptions": {
        "scope": "replaced",
        "reason": (
            "MATLAB signals failures with error identifiers such as "
            "radia:optuna:TrialPruned rather than exception classes."
        ),
        "discharged_by": ["radia:optuna:* error identifiers"],
    },
}


# The inventory walks __mro__, so an IntEnum drags in int's methods and an
# exception class drags in BaseException's.  Those are Python-language
# surface, not Optuna API: holding a MATLAB port to bit_count() or
# with_traceback() would be nonsense.
PYTHON_LANGUAGE_MEMBERS = frozenset(
    name
    for base in (int, BaseException)
    for name in dir(base)
    if not name.startswith("_")
)


def _scope_for(upstream: str) -> str:
    leaf = upstream.rsplit(".", 1)[-1]
    if leaf in PYTHON_LANGUAGE_MEMBERS:
        return "python-language"
    if leaf.startswith("_"):
        # A private upstream symbol is not part of the public contract.
        return "out-of-scope"
    for prefix, rule in SCOPE_RULES.items():
        if upstream == prefix or upstream.startswith(prefix + "."):
            return str(rule["scope"])
    return "required"


def _entry(
    upstream: str,
    kind: str,
    present: bool,
    matlab_name: str | None,
    oracle_status: str = "not-mapped",
) -> dict[str, object]:
    return {
        "kind": kind,
        "matlab_name": matlab_name,
        "oracle_status": oracle_status,
        "scope": _scope_for(upstream),
        "surface_status": "present" if present else "missing",
        "upstream": upstream,
    }


def build_coverage() -> dict[str, Any]:
    inventory = json.loads(INVENTORY_PATH.read_text(encoding="utf-8"))
    if inventory.get("optuna_version") != "4.9.0":
        raise RuntimeError("The public API inventory is not pinned to Optuna 4.9.0.")
    names, members = _matlab_surface()
    entries: list[dict[str, object]] = []
    for module in inventory["modules"]:
        module_name = str(module["module"])
        for symbol in module["symbols"]:
            name = str(symbol["name"])
            upstream = f"{module_name}.{name}"
            kind = str(symbol["kind"])
            if kind == "module":
                present = MODULE_EQUIVALENTS.get(name, False)
                matlab_name = f"radia.optuna ({name} namespace)" if present else None
            else:
                present = name in names
                matlab_name = f"radia.optuna.{name}" if present else None
            if present and name in VERIFIED_FUNCTIONS and kind == "function":
                oracle_status = "verified"
            elif present and name in PARTIALLY_VERIFIED_CLASSES and kind == "class":
                oracle_status = "partial"
            else:
                oracle_status = "not-mapped"
            entries.append(
                _entry(upstream, kind, present, matlab_name, oracle_status)
            )
            if kind != "class":
                continue
            for member in symbol["members"]:
                member_name = str(member["name"])
                normalized_member = _normalized(member_name)
                member_present = name in members and normalized_member in members[name]
                if member_present:
                    candidates = members[name][normalized_member]
                    actual_member = (
                        member_name
                        if member_name in candidates
                        else sorted(candidates, key=lambda value: (value.lower(), value))[0]
                    )
                else:
                    actual_member = None
                member_oracle = (
                    "verified"
                    if member_present
                    and member_name in VERIFIED_MEMBERS.get(name, set())
                    else "not-mapped"
                )
                entries.append(
                    _entry(
                        f"{upstream}.{member_name}",
                        f"class-{member['kind']}",
                        member_present,
                        f"radia.optuna.{name}.{actual_member}" if member_present else None,
                        member_oracle,
                    )
                )
    present_count = sum(entry["surface_status"] == "present" for entry in entries)
    missing_count = len(entries) - present_count
    verified_count = sum(entry["oracle_status"] == "verified" for entry in entries)
    partial_count = sum(entry["oracle_status"] == "partial" for entry in entries)
    required = [entry for entry in entries if entry["scope"] == "required"]
    required_present = [
        entry for entry in required if entry["surface_status"] == "present"
    ]
    required_mapped = [
        entry for entry in required_present if entry["oracle_status"] != "not-mapped"
    ]
    scope_counts: dict[str, int] = {}
    for entry in entries:
        scope_counts[str(entry["scope"])] = scope_counts.get(str(entry["scope"]), 0) + 1
    complete = len(required_present) == len(required) and len(required_mapped) == len(
        required
    )
    return {
        "schema": "radia.optuna49-api-coverage.v1",
        "upstream_version": "4.9.0",
        "upstream_inventory": "tests/matlab/fixtures/optuna49_public_api.json",
        "upstream_inventory_sha256": _sha256(INVENTORY_PATH),
        "closure_rule": (
            "full_compatibility_complete is true only when every entry whose scope "
            "is 'required' is present AND has an upstream differential-oracle "
            "mapping; the documented MATLAB extensions do not waive shared "
            "behavior, and a non-required scope must name what discharges it"
        ),
        "scope_rules": SCOPE_RULES,
        "scope_counts": scope_counts,
        "required_entry_count": len(required),
        "required_present_count": len(required_present),
        "required_missing_count": len(required) - len(required_present),
        "required_oracle_mapped_count": len(required_mapped),
        "required_oracle_unmapped_count": len(required) - len(required_mapped),
        "allowed_matlab_extensions": [
            "parallel execution and scheduling",
            "MATLAB table and MAT-file storage",
        ],
        "surface_entry_count": len(entries),
        "surface_present_count": present_count,
        "surface_missing_count": missing_count,
        "oracle_verified_count": verified_count,
        "oracle_partial_count": partial_count,
        "oracle_unmapped_count": len(entries)-verified_count-partial_count,
        "full_compatibility_complete": complete,
        "entries": entries,
    }


def main() -> None:
    DESTINATION.write_text(
        json.dumps(build_coverage(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(DESTINATION)


if __name__ == "__main__":
    main()
