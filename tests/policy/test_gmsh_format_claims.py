"""Guards: tracked source and documentation state the current MSH contract.

Assumes: released changelog sections are historical; native mesh assets have
their own structural gates. Rejection inputs are rendered from version data.
"""

from __future__ import annotations

import json
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


def test_tracked_files_make_no_retired_msh_format_claims():
    inventory = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", "ls-files", "-z"],
        cwd=ROOT, capture_output=True, check=True,
    )
    tracked = [name.decode("utf-8", errors="surrogateescape") for name in inventory.stdout.split(b"\0") if name]
    assert tracked, "required tracked-file inventory is empty"
    failures = []
    resolved_root = ROOT.resolve()
    for relative in tracked:
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
