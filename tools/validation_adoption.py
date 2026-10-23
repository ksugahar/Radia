"""Check recorded adoption arithmetic/provenance, without executing validation.

User consent, reference independence and theory applicability require an
independent reviewer. A matching string/hash does not authenticate a user.
"""
from __future__ import annotations

import ast
from datetime import date
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys
from functools import lru_cache


@lru_cache(maxsize=1)
def _bib_parser():
    """Reuse the canonical bibliography's stdlib parser without MCP imports."""
    source = (Path(__file__).resolve().parents[1] /
              "packages/radia-mcp/src/radia_mcp/bibliography/_bibparse.py")
    spec = importlib.util.spec_from_file_location("validation_record_bibparse", source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.parse_bib


def _analytic_reference(root, reference):
    """Resolve the single canonical key; URLs alone never establish identity."""
    _require("bib_key" in reference and "symbolic_verification" in reference,
             "analytic reference", sorted(reference), "canonical bib_key and hashed symbolic-verification record")
    canonical = root / "packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib"
    entries = [entry for entry in _bib_parser()(canonical.read_text(encoding="utf-8")) if entry.key]
    keys = [entry.key for entry in entries]
    _require(len(keys) == len(set(keys)), "canonical bibliography keys", keys, "unique citation keys")
    matches = [entry for entry in entries if entry.key == reference["bib_key"]]
    _require(len(matches) == 1, "analytic bibliography key", reference["bib_key"], "one existing canonical entry")
    fields = matches[0].fields
    doi = fields.get("doi", "").strip()
    isbn = re.sub(r"[-\s]", "", fields.get("isbn", ""))
    valid_isbn = (re.fullmatch(r"97[89][0-9]{10}", isbn) is not None and
                  sum(int(digit) * (1 if i % 2 == 0 else 3) for i, digit in enumerate(isbn)) % 10 == 0)
    _require(re.fullmatch(r"10\.[0-9]{4,9}/[^\s]+", doi) is not None or valid_isbn,
             "analytic bibliography identifier", {"bib_key": reference["bib_key"], "doi": doi, "isbn": isbn},
             "DOI or checksum-valid ISBN-13 (978/979); URL alone is insufficient")


def criterion_hash(criterion):
    """Hash canonical JSON; approval must identify the complete criterion."""
    return hashlib.sha256(json.dumps(criterion, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def scope_hash(record):
    """Bind consent to the named production feature and this decision record."""
    return criterion_hash({"record_id": record["id"], "feature": record["adoption"]["feature"]})


def _resolve_evidence(root, folder, record, reference, *, local_path):
    adoption = record["adoption"]
    if reference == "artifact:result":
        return {"kind": "artifact", "path": record["result"]["path"], "sha256": record["result"]["sha256"]}
    if reference == "artifact:entry":
        return {"kind": "artifact", "path": record["entry"]["script"], "sha256": record["entry"]["script_sha256"]}
    if reference == "approval:" + adoption["approval"]["approval_ref"]:
        return {"kind": "approval", "verification": "independent-review-required"}
    if reference == "review:" + adoption["criterion"]["threshold_verification"]["review_ref"]:
        return {"kind": "review", "verification": "independent-review-required"}
    if reference.startswith("record:"):
        identity = reference[len("record:"):]
        _require(identity != record["id"], "self evidence reference", identity, "a different record id")
        base = folder
        for _ in record["id"].split("/"):
            base = base.parent
        target = local_path(base, identity) / "record.json"
        raw = target.read_bytes()
        other = json.loads(raw)
        _require(other.get("id") == identity, "evidence record id", other.get("id"), identity)
        return {"kind": "record", "id": identity, "sha256": hashlib.sha256(raw).hexdigest(),
                "verification": "existence-and-identity-only; independent-review-required"}
    raise ValueError(f"adoption unresolved evidence reference: actual={reference!r}; "
                     "expected=artifact:result, artifact:entry, declared approval:/review: or existing record:<id>")


def _require(ok, label, actual, expected):
    if not ok:
        raise ValueError(f"adoption {label}: actual={actual!r}; expected={expected!r}")


def _unique(rows, label):
    result = {row["id"]: row for row in rows}
    _require(len(result) == len(rows), label, [r["id"] for r in rows], "unique ids")
    return result


def _pointer(value, pointer):
    _require(pointer.startswith("/"), "binding pointer", pointer, "JSON pointer")
    try:
        for token in pointer[1:].split("/"):
            token = token.replace("~1", "/").replace("~0", "~")
            value = value[int(token)] if isinstance(value, list) else value[token]
        return value
    except (KeyError, IndexError, ValueError, TypeError) as exc:
        raise ValueError(f"adoption invalid criterion pointer {pointer!r}") from exc


def _artifact(folder, ref, local_path, sha256):
    data = local_path(folder, ref["path"]).read_bytes()
    _require(sha256(data) == ref["sha256"], "artifact hash", sha256(data), ref["sha256"])
    return data


def validate_adoption(root, folder, record, payload, *, local_path, sha256):
    adoption = record.get("adoption")
    if adoption is None:
        _require(record["category"] != "adoption", "adoption category", record["category"], "adoption block required")
        return {"adoption_claim": False}
    _require(record["category"] == "adoption", "adoption category", record["category"], "adoption")
    rationale = adoption["rationale"]
    assumptions = _unique(rationale["assumptions"], "rationale assumption ids")
    dependencies = rationale["depends_on"]
    _require(len(set(dependencies)) == len(dependencies) and set(dependencies) == set(assumptions),
             "rationale dependencies", dependencies, sorted(assumptions))
    _require(adoption["approval"]["rationale_sha256"] == criterion_hash(rationale),
             "approved rationale hash", adoption["approval"]["rationale_sha256"], criterion_hash(rationale))
    conditions = _unique(rationale["break_conditions"], "rationale invalidation ids")
    resolved = {}
    for identity, assumption in assumptions.items():
        _require(set(assumption["invalidated_by"]) <= set(conditions), "invalidation condition ids",
                 assumption["invalidated_by"], sorted(conditions))
        resolved[identity] = {ref: _resolve_evidence(root, folder, record, ref, local_path=local_path)
                              for ref in assumption["evidence_refs"]}
    _require(adoption["approval"]["scope_sha256"] == scope_hash(record), "approved feature/record scope",
             adoption["approval"]["scope_sha256"], scope_hash(record))
    criterion = adoption["criterion"]
    approval = adoption["approval"]
    approved_on = date.fromisoformat(approval["approved_on"])
    _require(approved_on <= date.fromisoformat(record["date"]), "approval date",
             approval["approved_on"], f"not later than record date {record['date']}")
    _require(approval["criterion_sha256"] == criterion_hash(criterion), "approved criterion hash",
             approval["criterion_sha256"], criterion_hash(criterion))
    _require(record["source"]["provenance"] == "exact" and record["source"]["repo_sha"] is not None,
             "source provenance", record["source"]["provenance"], "exact, non-null source")
    result = adoption["result"]
    _require(result["source_sha"] == record["source"]["repo_sha"], "result source",
             result["source_sha"], record["source"]["repo_sha"])
    _require(result["artifact_sha256"] == record["result"]["sha256"], "result identity",
             result["artifact_sha256"], record["result"]["sha256"])
    # Artifact hash is outside the artifact: no self-referential hash field.
    expected_payload = {key: value for key, value in result.items() if key != "artifact_sha256"}
    _require(isinstance(payload, dict) and
             json.dumps(payload.get("adoption_result"), sort_keys=True, allow_nan=False) ==
             json.dumps(expected_payload, sort_keys=True, allow_nan=False),
             "artifact values", payload.get("adoption_result") if isinstance(payload, dict) else payload,
             expected_payload)
    reference = criterion["reference"]
    if reference["kind"] == "analytic":
        _analytic_reference(root, reference)
        _artifact(folder, reference["symbolic_verification"], local_path, sha256)
    elif reference["kind"] == "independent-route":
        _require("independence" in reference and "shared_components" in reference,
                 "independent reference", sorted(reference), "independence rationale and shared components")
    else:
        _require("uncertainty" in reference, "reference uncertainty", sorted(reference), "uncertainty statement")
    verification = criterion["threshold_verification"]
    if verification["mode"] == "literal-constants":
        _require("script" in verification and bool(verification.get("bindings")),
                 "literal verification", sorted(verification), "script and nonempty bindings")
        text = _artifact(folder, verification["script"], local_path, sha256).decode("utf-8-sig")
        try:
            tree = ast.parse(text)
        except SyntaxError as exc:
            raise ValueError("adoption validation source is not parseable Python") from exc
        constants = {}
        assignments = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
                assignments[node.id] = assignments.get(node.id, 0) + 1
        for node in tree.body:
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant):
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        _require(target.id not in constants, "duplicate constant", target.id, "one assignment")
                        constants[target.id] = node.value.value
        for binding in verification["bindings"]:
            _require(re.fullmatch(r"/(?:residual_threshold|series/[0-9]+/(?:minimum_rate|finest_threshold/(?:value|operator)))",
                                  binding["criterion_pointer"]) is not None,
                     "binding scope", binding["criterion_pointer"], "threshold/operator/rate pointer")
            _require(assignments.get(binding["constant"]) == 1, "literal assignment count",
                     assignments.get(binding["constant"]), "one unambiguous literal assignment")
            expected = _pointer(criterion, binding["criterion_pointer"])
            actual = constants.get(binding["constant"])
            _require(type(actual) is type(expected) and actual == expected,
                     "literal binding " + binding["constant"], actual, expected)
    else:
        _require("script" not in verification and "bindings" not in verification,
                 "reviewer verification", sorted(verification), "mode and public review_ref only")
    variables = criterion["variables"]
    _require(bool(variables), "variables", variables, "nonempty variable definitions")
    levels = _unique(criterion["levels"], "criterion levels")
    observed = _unique(result["levels"], "observed levels")
    _require(levels.keys() == observed.keys(), "level coverage", sorted(observed), sorted(levels))
    for key, level in levels.items():
        row = observed[key]
        _require(set(level["parameters"]) == set(variables), "level variables", level["parameters"], sorted(variables))
        _require(row["parameters"] == level["parameters"], "recorded parameters", row["parameters"], level["parameters"])
        _require(row["source_sha"] == result["source_sha"], "level source", row["source_sha"], result["source_sha"])
    series_ids = set()
    used = set()
    failed = any(row["relative_residual"] > criterion["residual_threshold"] for row in observed.values())
    unresolved = False
    rates = {}
    for series in criterion["series"]:
        _require(series["id"] not in series_ids, "series id", series["id"], "unique")
        series_ids.add(series["id"])
        ids = series["level_ids"]
        _require(len(set(ids)) == len(ids) and set(ids) <= levels.keys(), "series levels", ids, "distinct declared levels")
        axis = series["axis"]
        _require(axis in variables, "series axis", axis, sorted(variables))
        used.update(ids)
        rates[series["id"]] = []
        for a, b in zip(ids, ids[1:]):
            pa, pb = levels[a]["parameters"], levels[b]["parameters"]
            _require(all(pa[k] == pb[k] for k in variables if k != axis), "held-fixed variables", (pa, pb), "only axis changes")
            x, y = pa[axis], pb[axis]
            direction = series["direction"]
            _require(y > x if direction == "increasing" else x > y, "axis direction", (x, y), direction)
            ea, eb = observed[a]["error"], observed[b]["error"]
            failed |= eb >= ea
            if ea == 0 or eb == 0:
                unresolved = True
                rates[series["id"]].append(None)
                continue
            if series["rate_model"] == "log-ratio":
                _require(x > 0 and y > 0, "log-rate coordinates", (x, y), "positive")
                denominator = abs(math.log(x) - math.log(y))
            else:
                denominator = abs(y - x)
            _require(denominator > 0 and math.isfinite(denominator), "rate denominator", denominator, "positive finite")
            rate = (math.log(ea) - math.log(eb)) / denominator
            _require(math.isfinite(rate), "observed rate", rate, "finite")
            rates[series["id"]].append(rate)
            failed |= rate < series["minimum_rate"]
        final = observed[ids[-1]]["error"]
        gate = series["finest_threshold"]
        failed |= not (final < gate["value"] if gate["operator"] == "lt" else final <= gate["value"])
    _require(used == levels.keys(), "unexamined levels", sorted(levels.keys() - used), "all levels in a declared series")
    verdict = "FAIL" if failed else "UNRESOLVED" if unresolved else "PASS"
    _require(result["verdict"] == verdict, "series verdict", result["verdict"], verdict)
    expected_outcome = {"PASS": "pass", "FAIL": "fail", "UNRESOLVED": "observed"}[verdict]
    _require(record["outcome"] == expected_outcome, "record outcome", record["outcome"], expected_outcome)
    for path in adoption["protected_by"]:
        _require(path.endswith(".py"), "lock extension", path, ".py module")
        _require(path.startswith("tests/") or "/tests/" in path, "lock location", path, "tests path")
        _require(local_path(root, path).is_file(), "lock existence", path, "existing file")
        try:
            doc = ast.get_docstring(ast.parse(local_path(root, path).read_text(encoding="utf-8-sig"))) or ""
        except SyntaxError as exc:
            raise ValueError(f"adoption lock module syntax: actual={path!r}; expected=valid Python") from exc
        required = f"surrogate test of adoption record {record['id']} (convergence level 1)"
        _require(required in doc.partition("Assumes:")[2], "lock Assumes citation", path, required)

    return {"adoption_claim": verdict == "PASS" and record["status"] == "current",
            "adoption_verdict": verdict, "observed_rates": rates,
            "assumption_dependencies": {
                key: {"record_id": record["id"], "feature": adoption["feature"],
                      "evidence_refs": row["evidence_refs"], "resolved_evidence": resolved[key],
                      "invalidated_by": row["invalidated_by"],
                      "recheck": row["recheck"]} for key, row in assumptions.items()},
            "rationale_validity": "independent-review-required",
            "approval_authentication": "independent-review-required",
            "theory_applicability": "independent-review-required",
            "threshold_verification": verification["mode"], "threshold_wiring_review": "required"}


def validate_champion(root, base, folder, record, payload, *, local_path, sha256, read_json, load_record):
    """Check one-off champion provenance; never execute or rank timed work."""
    champion = record.get("champion")
    named = record["category"] in {"accuracy-champion", "speed", "scale"}
    _require(named == (champion is not None), "champion category/block", record["category"],
             "champion block exactly for accuracy-champion, speed or scale")
    if champion is None:
        return {"champion_record": False}
    _require("adoption" not in record, "champion role", "adoption" in record,
             "extension of a separate adoption decision, not a new adoption claim")
    ref = champion["adoption_record"]
    _require(ref["id"] != record["id"], "champion adoption reference", ref["id"], "different record id")
    target = local_path(base, ref["id"]) / "record.json"
    raw = target.read_bytes()
    _require(sha256(raw) == ref["sha256"], "adoption record hash", sha256(raw), ref["sha256"])
    original = read_json(raw)
    _require(original["id"] == ref["id"] and "champion" not in original and "adoption" in original,
             "adoption reference role", original.get("id"), "matching non-champion adoption record")
    review = load_record(root, target)
    _require(review["adoption_verdict"] == "PASS", "origin adoption verdict",
             review["adoption_verdict"], "PASS")
    _require(original["status"] != "stale" or record["status"] == "stale", "champion origin status",
             (original["status"], record["status"]), "stale champion when origin is stale")
    _require(record["outcome"] == "pass", "champion outcome", record["outcome"], "pass")
    criterion = original["adoption"]["criterion"]
    expected_quantity = {key: criterion["quantity"][key] for key in ("name", "unit")}
    _require(champion["quantity"] == expected_quantity, "champion quantity", champion["quantity"], expected_quantity)
    for level in champion["levels"]:
        _require(set(level["parameters"]) == set(criterion["variables"]), "champion parameter keys",
                 sorted(level["parameters"]), sorted(criterion["variables"]))
    if record["category"] == "accuracy-champion":
        errors = [level["error"] for level in champion["levels"]]
        origin_levels = {level["id"]: level for level in original["adoption"]["result"]["levels"]}
        origin_finest = min(origin_levels[series["level_ids"][-1]]["error"] for series in criterion["series"])
        _require(all(error <= origin_finest for error in errors), "champion origin error bound", errors,
                 f"each error <= finest origin error {origin_finest}")
        _require(all(b <= a for a, b in zip(errors, errors[1:])), "champion error trend", errors,
                 "nonincreasing errors in declared table order")
    _require(champion["feature"] == original["adoption"]["feature"], "champion feature",
             champion["feature"], original["adoption"]["feature"])
    _require(record["source"]["provenance"] == "exact", "champion provenance",
             record["source"]["provenance"], "exact")
    expected = {key: value for key, value in champion.items() if key != "adoption_record"}
    _require(isinstance(payload, dict) and criterion_hash(payload.get("champion_result")) == criterion_hash(expected),
             "champion artifact values", payload.get("champion_result") if isinstance(payload, dict) else None,
             expected)
    levels = _unique(champion["levels"], "champion level ids")
    _require(champion["selected_level_id"] in levels, "champion selected level",
             champion["selected_level_id"], sorted(levels))
    _require(champion["selected_level_id"] == champion["levels"][-1]["id"],
             "champion table endpoint", champion["selected_level_id"], champion["levels"][-1]["id"])
    _require(levels[champion["selected_level_id"]]["source_sha"] == record["source"]["repo_sha"],
             "champion selected source", levels[champion["selected_level_id"]]["source_sha"], record["source"]["repo_sha"])
    for level in levels.values():
        evidence = read_json(_artifact(folder, level["artifact"], local_path, sha256))
        fields = {key: value for key, value in level.items() if key != "artifact"}
        _require(isinstance(evidence, dict) and criterion_hash(evidence.get("champion_level")) == criterion_hash(fields),
                 "champion level artifact values", evidence, fields)
    return {"champion_record": True, "champion_category": record["category"],
            "adoption_origin": ref["id"], "adoption_origin_status": original["status"],
            "selected_level_id": champion["selected_level_id"], "champion_levels": len(levels),
            "champion_scope_and_best_result": "independent-review-required",
            "numerical_replay": False, "timing_used_as_gate": False}
