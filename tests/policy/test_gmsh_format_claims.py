"""Guards: tracked source and documentation state the current MSH contract.

Assumes: released changelog sections and provenance-bearing stale calculations
are historical (saved notebooks must match their recorded hash); native mesh assets have
their own structural gates. Rejection inputs are rendered from version data.
"""

from __future__ import annotations

import json
import hashlib
from datetime import datetime
from pathlib import Path
import re
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[2]
TEXT_SUFFIXES = {
    ".bat", ".c", ".cc", ".cmake", ".cmd", ".cpp", ".geo", ".h",
    ".hpp", ".html", ".ipynb", ".jou", ".js", ".json", ".m", ".md",
    ".opt", ".ps1", ".py", ".rst", ".sh", ".toml", ".txt", ".yaml",
    ".yml",
}

# Exact format tokens, not software versions or arbitrary physical values.
OLD_VERSION = r"(?:4[.]1|2[.]2)(?!\d|[.]\d)"
PUNCTUATION = r"[\s`*_()\[\]|:\"'=-]*"
FORMAT_TOKENS = re.compile(r"4[.]1|2[.]2|msh2\b", re.IGNORECASE)
FORMAT_CLAIM = re.compile(
    rf"(?:\bMSH\b|[.]msh\b){PUNCTUATION}"
    rf"(?:(?:ASCII|binary|format|version|is|only)\b{PUNCTUATION}){{0,3}}"
    rf"v?{OLD_VERSION}",
    re.IGNORECASE,
)
GMSH_VERSION = re.compile(
    rf"\bGMSH\b{PUNCTUATION}v?{OLD_VERSION}", re.IGNORECASE,
)
FORMAT_CONTEXT = re.compile(
    r"[.]msh\b|\b(?:export\w*|format\w*|mesh\w*|emit\w*|write\w*|output)\b",
    re.IGNORECASE,
)
FORMAT_SETTING = re.compile(
    rf"\b(?:Mesh[.]MshFileVersion|MshFileVersion|mesh_version|mesh_format|msh_version)"
    rf"\b[\s\"'`]*(?:=|:|,)[\s\"'`]*{OLD_VERSION}",
    re.IGNORECASE,
)
MESH_HEADER = re.compile(
    rf"\$MeshFormat(?:(?!\$EndMeshFormat).){{0,120}}?"
    rf"{OLD_VERSION}\s+[01]\s+[48](?!\d)",
    re.IGNORECASE | re.DOTALL,
)
OLD_CLI_FORMAT = re.compile(r"(?<!\w)-format\s+[\"']?msh2\b", re.IGNORECASE)
HEADING = re.compile(r"^(#{1,6})\s+(.+?)\s*$")
RELEASE = re.compile(r"^\[?v?\d+[.](?:\d+|x)(?:[.][\w.-]+)?(?:\]|\s|$)", re.I)
UNRELEASED = re.compile(r"^\[?unreleased\]?(?:\s|$)", re.IGNORECASE)


def _active_changelog(text: str) -> str:
    """Keep introduction/Unreleased text and preserve original line numbers."""
    lines = text.splitlines(keepends=True)
    unreleased_levels = {
        len(match[1])
        for line in lines
        if (match := HEADING.match(line))
        and UNRELEASED.match(match[2])
    }
    historical = False
    result = []
    for line in lines:
        heading = HEADING.match(line)
        if heading:
            level, title = len(heading[1]), heading[2]
            if UNRELEASED.match(title):
                historical = False
            elif RELEASE.match(title) and (not unreleased_levels or level in unreleased_levels):
                historical = True
        result.append("\n" if historical and line.endswith("\n") else "" if historical else line)
    return "".join(result)


def _claims(text: str) -> list[tuple[int, str]]:
    if not FORMAT_TOKENS.search(text):
        return []
    matches = list(FORMAT_CLAIM.finditer(text))
    matches += list(FORMAT_SETTING.finditer(text))
    matches += list(MESH_HEADER.finditer(text))
    matches += list(OLD_CLI_FORMAT.finditer(text))
    for match in GMSH_VERSION.finditer(text):
        start = text.rfind("\n", 0, match.start()) + 1
        end = text.find("\n", match.end())
        if FORMAT_CONTEXT.search(text[start:end if end >= 0 else len(text)]):
            matches.append(match)
    return sorted({
        (text.count("\n", 0, match.start()) + 1, match.group().replace("\n", r"\n"))
        for match in matches
    })


def _historical_calculation(text: str) -> dict | None:
    """Recognize self-described stale calculations, never current configuration."""
    try:
        record = json.loads(text)
        if not isinstance(record, dict):
            return None
        schema = record.get("schema", "")
        assert schema.startswith("radia.validation.")
        assert record.get("stale") is True and record.get("stale_reason", "").strip()
        datetime.fromisoformat(record["generated_at_utc"].replace("Z", "+00:00"))
        assert re.fullmatch(r"\d+(?:[.]\d+){1,3}", record["versions"]["radia_version"])
        assert re.fullmatch(r"\d+[.]\d+ [01] [48]", record["msh"]["mesh_format"])
        return record
    except (ValueError, KeyError, TypeError, AttributeError, AssertionError):
        return None


def _archived_notebooks(root: Path, records: list[dict]) -> set[str]:
    """Only unchanged saved notebook bytes linked by calculation provenance."""
    paths = set()
    for record in records:
        notebook = record.get("saved_notebook")
        if notebook is None:
            continue
        path = root / notebook["path"]
        assert path.suffix == ".ipynb" and path.resolve().is_relative_to(root.resolve())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == notebook["sha256"], path
        paths.add(path.relative_to(root).as_posix())
    return paths


def _document_units(path: Path, text: str):
    if path.suffix.lower() == ".ipynb":
        notebook = json.loads(text)
        for index, cell in enumerate(notebook["cells"], 1):
            source = cell.get("source", [])
            yield f"cell {index} source", "".join(source) if isinstance(source, list) else source
            for number, output in enumerate(cell.get("outputs", []), 1):
                content = dict(output.get("data", {}))
                if "text" in output:
                    content["text/plain"] = output["text"]
                for kind in ("text/plain", "text/markdown"):
                    value = content.get(kind)
                    if value is not None:
                        yield f"cell {index} output {number} {kind}", "".join(value) if isinstance(value, list) else value
    else:
        yield "", _active_changelog(text) if path.stem.casefold() == "changelog" else text


@pytest.mark.parametrize("major,minor", [(4, 1), (2, 2)])
def test_claim_detector_distinguishes_format_contracts_from_values(major, minor):
    old = f"{major}.{minor}"
    bad = (
        f"Export MSH v{old} fields",
        f"| Gmsh (.msh) | v{old}, order 1-3 |",
        f'gmsh.option.setNumber("Mesh.MshFileVersion", {old})',
        f'{{"mesh_format": "{old} 0 8"}}',
        f'$MeshFormat\\n{old} 0 8\\n$EndMeshFormat',
        f'assert "GMSH v{old} .msh output" in response',
    )
    assert all(_claims(value) for value in bad), bad
    good = (
        "Current ASCII MSH 4.2 output",
        "Gmsh 4.15.2 exports current MSH 4.2",
        "Gmsh 4.1.0 software provenance",
        "physical coefficients: 2.2 and 4.1",
        "field.msh has amplitude 2.2",
        "rejected_versions = [(4, 1), (2, 2)]",
    )
    assert all(not _claims(value) for value in good), good
    history = f"# Changelog\n## [Unreleased]\nMSH 4.2\n## [1.4.29]\nMSH v{old}\n"
    assert not _claims(_active_changelog(history))
    assert _claims(_active_changelog(history.replace("MSH 4.2", f"MSH v{old}", 1)))
    if major == 2:
        assert _claims(f"gmsh result.geo -format msh{minor}")
    image_only = json.dumps({"cells": [{"source": [], "outputs": [{"data": {"image/png": f"MSH v{old}"}}]}]})
    assert not any(_claims(unit) for _, unit in _document_units(Path("result.ipynb"), image_only))
    for cell in ({"source": [f"MSH v{old}"]}, {"outputs": [{"text": [f"MSH v{old}"]}]}):
        notebook_text = json.dumps({"cells": [cell]})
        assert any(_claims(unit) for _, unit in _document_units(Path("result.ipynb"), notebook_text))


def test_historical_calculation_requires_provenance_and_unchanged_notebook(tmp_path):
    record = {
        "schema": "radia.validation.saved_calculation.v1",
        "generated_at_utc": "2026-08-05T05:10:48Z",
        "versions": {"radia_version": "4.95.40"},
        "msh": {"mesh_format": ".".join(map(str, (4, 1))) + " 0 8"},
        "stale": True, "stale_reason": "historical calculation",
    }
    assert _historical_calculation(json.dumps(record)) == record
    for key in ("schema", "generated_at_utc", "versions", "msh", "stale", "stale_reason"):
        assert _historical_calculation(json.dumps({k: v for k, v in record.items() if k != key})) is None, key
    notebook = tmp_path / "saved.ipynb"
    notebook.write_bytes(b'{"cells": []}')
    record["saved_notebook"] = {"path": notebook.name, "sha256": hashlib.sha256(notebook.read_bytes()).hexdigest()}
    assert _archived_notebooks(tmp_path, [record]) == {notebook.name}
    notebook.write_bytes(b'{"cells": [{"source": ["changed claim"]}]}')
    with pytest.raises(AssertionError, match="saved.ipynb"):
        _archived_notebooks(tmp_path, [record])


def test_tracked_files_make_no_retired_msh_format_claims():
    inventory = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", "ls-files", "-z"],
        cwd=ROOT, capture_output=True, check=True,
    )
    tracked = [name.decode("utf-8", errors="surrogateescape") for name in inventory.stdout.split(b"\0") if name]
    assert tracked, "required tracked-file inventory is empty"
    historical = []
    archived_records = set()
    for relative in tracked:
        if Path(relative).suffix.lower() == ".json":
            record = _historical_calculation((ROOT / relative).read_text(encoding="utf-8-sig"))
            if record is not None:
                historical.append(record)
                archived_records.add(relative)
    archived_notebooks = _archived_notebooks(ROOT, historical)
    failures = []
    resolved_root = ROOT.resolve()
    for relative in tracked:
        if relative in archived_records or relative in archived_notebooks:
            continue
        path = ROOT / relative
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            assert path.resolve(strict=True).is_relative_to(resolved_root), "tracked text resolves outside repository"
            raw = path.read_bytes()
            if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
                text = raw.decode("utf-16")
            elif b"\0" in raw:
                continue  # Binary payload, not a source/prose claim.
            else:
                text = raw.decode("utf-8-sig", errors="surrogateescape")
            for location, unit in _document_units(path, text):
                for line, claim in _claims(unit):
                    failures.append(f"{relative}:{location}:{line}: {claim!r}")
        except (OSError, ValueError, KeyError, TypeError, AssertionError, RuntimeError) as exc:
            failures.append(f"{relative}: cannot lint required tracked text: {exc!r}")
    assert not failures, "Retired MSH format claims in active tracked content:\n" + "\n".join(failures)
