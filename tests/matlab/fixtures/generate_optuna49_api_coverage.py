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
        "trials_dataframe",
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
    for path in sorted(MATLAB_DIRECTORY.rglob("*.m")):
        if "+internal" in path.parts:
            continue
        names.add(path.stem)
        source = path.read_text(encoding="utf-8")
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
    return names, members


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
    return {
        "schema": "radia.optuna49-api-coverage.v1",
        "upstream_version": "4.9.0",
        "upstream_inventory": "tests/matlab/fixtures/optuna49_public_api.json",
        "upstream_inventory_sha256": _sha256(INVENTORY_PATH),
        "closure_rule": (
            "full_compatibility_complete is true only when every required public "
            "symbol/member is present and has an upstream differential-oracle mapping; "
            "the two documented MATLAB extensions do not waive shared behavior"
        ),
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
        "full_compatibility_complete": False,
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
