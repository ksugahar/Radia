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
FORMAT_TOKENS = re.compile(r"4[.]1|2[.]2|msh(?:2|22|41)\b|MshFileVersion|_v(?:41|22)(?=[_.-])", re.IGNORECASE)
FORMAT_CLAIM = re.compile(
    rf"(?:\bMSH\b|[.]msh\b){PUNCTUATION}"
    rf"(?:(?:ASCII|binary|file|format|version|is|only)\b{PUNCTUATION}){{0,3}}"
    rf"v?{OLD_VERSION}",
    re.IGNORECASE,
)
GMSH_VERSION = re.compile(
    rf"\bGMSH\b{PUNCTUATION}v?{OLD_VERSION}", re.IGNORECASE,
)
FORMAT_SETTING = re.compile(
    rf"\b(?:Mesh[.]MshFileVersion|MshFileVersion|mesh_version|mesh_format|msh_version)"
    rf"\b[\s\"'`]*(?:=|:|,)[\s\"'`]*{OLD_VERSION}",
    re.IGNORECASE,
)
MESH_HEADER = re.compile(
    rf"\$MeshFormat(?:(?!\$EndMeshFormat).){{0,120}}?"
    rf"{OLD_VERSION}\s+[01](?!\d)",
    re.IGNORECASE | re.DOTALL,
)
OLD_CLI_FORMAT = re.compile(r"(?<!\w)-format\s+[\"']?msh(?:2|22|41)\b", re.IGNORECASE)
INTEGER_SETTING = re.compile(r'\b(?:Mesh[.])?MshFileVersion[\s"\'`]*(?:=|:|,)[\s"\'`]*2(?![.\d])', re.I)
OLD_FILENAME = re.compile(r"[\w./\\-]*_v(?:4[.]?1|2[.]?2)(?=[_.-])[\w.-]*[.](?:msh|geo)\b", re.I)
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
    matches += list(INTEGER_SETTING.finditer(text))
    matches += list(OLD_FILENAME.finditer(text))
    matches += list(GMSH_VERSION.finditer(text))
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
        notebook = record["saved_notebook"]
        assert isinstance(notebook["path"], str) and notebook["path"].endswith(".ipynb")
        assert re.fullmatch(r"[0-9a-f]{64}", notebook["sha256"])
        return record
    except (ValueError, KeyError, TypeError, AttributeError, AssertionError):
        return None


def _archived_notebooks(root: Path, records: list[dict]) -> set[str]:
    """Only unchanged saved notebook bytes linked by calculation provenance."""
    paths = set()
    for record in records:
        notebook = record["saved_notebook"]
        path = root / notebook["path"]
        assert path.suffix == ".ipynb" and path.resolve().is_relative_to(root.resolve())
        payload = path.read_bytes()
        if notebook.get("text_eol") == "lf":
            payload = payload.replace(b"\r\n", b"\n")  # Git text checkout conversion only.
        assert hashlib.sha256(payload).hexdigest() == notebook["sha256"], path
        paths.add(path.relative_to(root).as_posix())
    return paths


def _active_record_text(record: dict) -> str:
    """Exclude only the archived format value, never arbitrary new prose."""
    active = dict(record, msh=dict(record["msh"]))
    active["msh"]["mesh_format"] = "historical format (original preserved in record)"
    return json.dumps(active)


def _candidate_paths(root: Path) -> set[str]:
    """Native Git prefilter; still repository-wide, including UTF-16 text."""
    # Prefilter ordinary text by format-related tokens. The binary-marked
    # suffix fallback below preserves UTF-16; image payload files are excluded.
    pattern = r"4.?[.].?1|2.?[.].?2|m.?s.?h.?(2|4|f)|v.?(4.?1|2.?2)[_.-]|s.?a.?v.?e.?d.?_.?n.?o.?t.?e.?b.?o.?o.?k"
    found = subprocess.run(
        ["git", "-c", f"safe.directory={root.as_posix()}", "grep", "--text", "-i", "-l", "-z", "-E", pattern,
         "--", *[f":(icase)*{suffix}" for suffix in sorted(TEXT_SUFFIXES)]],
        cwd=root, capture_output=True,
    )
    assert found.returncode in (0, 1), found.stderr
    paths = {name.decode("utf-8", errors="surrogateescape") for name in found.stdout.split(b"\0") if name}
    # POSIX Git regex stops at NUL despite --text. Native Git classifies
    # binary working-tree bytes; retain binary-marked text suffixes so the
    # normal decoder below still inspects UTF-16 (without PCRE dependency).
    endings = subprocess.run(["git", "-c", f"safe.directory={root.as_posix()}",
                              "ls-files", "--eol", "-z"], cwd=root,
                             capture_output=True, check=True)
    for row in endings.stdout.split(b"\0"):
        metadata, separator, raw_name = row.partition(b"\t")
        if separator and b"w/-text" in metadata:
            name = raw_name.decode("utf-8", errors="surrogateescape")
            if Path(name).suffix.lower() in TEXT_SUFFIXES:
                paths.add(name)
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
        f"## GMSH v{old}",
        f"MSH file format {old}",
        f"case_v{old}_dim2.msh",
        f"case_v{major}{minor}_dim2.geo",
        f'"$MeshFormat\\n" << "{old} 0 "',
        f"gmsh case.geo -format msh{major}{minor}",
    )
    assert all(_claims(value) for value in bad), bad
    good = (
        "Current ASCII MSH 4.2 output",
        "Gmsh 4.15.2 exports current MSH 4.2",
        "Gmsh 4.1.0 software provenance",
        "physical coefficients: 2.2 and 4.1",
        "field.msh has amplitude 2.2",
        "rejected_versions = [(4, 1), (2, 2)]",
        "test_generalization_v41.py",
        "test_gate_v22_revision",
    )
    assert all(not _claims(value) for value in good), good
    history = f"# Changelog\n## [Unreleased]\nMSH 4.2\n## [1.4.29]\nMSH v{old}\n"
    assert not _claims(_active_changelog(history))
    assert _claims(_active_changelog(history.replace("MSH 4.2", f"MSH v{old}", 1)))
    if major == 2:
        assert _claims(f"gmsh result.geo -format msh{minor}")
        assert _claims('setNumber("Mesh.MshFileVersion", ' + str(minor) + ')')
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
    assert _historical_calculation(json.dumps(record)) is None, "unlinked stale record must not be exempt"
    notebook = tmp_path / "saved.ipynb"
    notebook.write_bytes(b'{"cells": []}')
    record["saved_notebook"] = {"path": notebook.name, "sha256": hashlib.sha256(notebook.read_bytes()).hexdigest()}
    assert _historical_calculation(json.dumps(record)) == record
    for key in ("schema", "generated_at_utc", "versions", "msh", "stale", "stale_reason", "saved_notebook"):
        assert _historical_calculation(json.dumps({k: v for k, v in record.items() if k != key})) is None, key
    assert _archived_notebooks(tmp_path, [record]) == {notebook.name}
    assert not _claims(_active_record_text(record))
    record["new_prose"] = "current exporter writes MSH v" + ".".join(map(str, (4, 1)))
    assert _claims(_active_record_text(record)), "valid archival metadata must not hide active claims"
    notebook.write_bytes(b'{"cells": [{"source": ["changed claim"]}]}')
    with pytest.raises(AssertionError, match="saved.ipynb"):
        _archived_notebooks(tmp_path, [record])


def test_git_prefilter_preserves_utf16_and_cli_claims(tmp_path):
    root = tmp_path / "tracked"
    root.mkdir()
    for command in (["git", "init", "-q", str(root)],):
        subprocess.run(command, capture_output=True, check=True)
    old = ".".join(map(str, (4, 1)))
    inputs = {"legacy.py": (f'Mesh.MshFileVersion = {old}', "utf-16"),
              "command.sh": ("gmsh file.geo -format msh" + "41", "utf-8"),
              "reference.py": ('path="case_v' + "41" + '.geo"', "utf-8"),
              "current.py": ('Mesh.MshFileVersion = 4.2', "utf-8")}
    for name, (text, encoding) in inputs.items():
        (root / name).write_text(text, encoding=encoding)
    subprocess.run(["git", "-C", str(root), "add", "."], capture_output=True, check=True)
    selected = _candidate_paths(root)
    assert {"legacy.py", "command.sh", "reference.py"}.issubset(selected), selected
    assert _claims((root / "legacy.py").read_text(encoding="utf-16"))


def test_tracked_files_make_no_retired_msh_format_claims(record_testsuite_property):
    inventory = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", "ls-files", "-z"],
        cwd=ROOT, capture_output=True, check=True,
    )
    tracked = [name.decode("utf-8", errors="surrogateescape") for name in inventory.stdout.split(b"\0") if name]
    assert tracked, "required tracked-file inventory is empty"
    candidates = _candidate_paths(ROOT)
    assert candidates.issubset(set(tracked)), "grep returned nontracked input"
    record_testsuite_property("tracked_paths", len(tracked))
    record_testsuite_property("content_candidates", len(candidates))
    historical = []
    records = {}
    for relative in candidates:
        if Path(relative).suffix.lower() == ".json":
            text = (ROOT / relative).read_text(encoding="utf-8-sig")
            record = _historical_calculation(text) if "saved_notebook" in text else None
            if record is not None:
                historical.append(record)
                records[relative] = record
    archived_notebooks = _archived_notebooks(ROOT, historical)
    failures = [f"{relative}: retired format filename" for relative in tracked if OLD_FILENAME.search(relative)]
    for relative in sorted(candidates - archived_notebooks):
        path = ROOT / relative
        try:
            assert path.resolve(strict=True).is_relative_to(ROOT.resolve()), "tracked text resolves outside repository"
            raw = path.read_bytes()
            if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
                text = raw.decode("utf-16")
            elif b"\0" in raw:
                continue  # Binary payload, not a source/prose claim.
            else:
                text = raw.decode("utf-8-sig", errors="surrogateescape")
            if relative in records:
                text = _active_record_text(records[relative])
            for location, unit in _document_units(path, text):
                for line, claim in _claims(unit):
                    failures.append(f"{relative}:{location}:{line}: {claim!r}")
        except (OSError, ValueError, KeyError, TypeError, AssertionError, RuntimeError) as exc:
            failures.append(f"{relative}: cannot lint required tracked text: {exc!r}")
    assert not failures, "Retired MSH format claims in active tracked content:\n" + "\n".join(failures)
