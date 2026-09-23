"""Require published formulation citations and canonical bibliography entries.

The total/reduced scalar-potential method is attributed where it is defined.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BIB = (ROOT / "packages/radia-mcp/src/radia_mcp/bibliography/data"
       / "references.bib")

# The total/reduced scalar potential, and its nonlinear 3-D treatment.
KEYS = {
    "simkin1979use": "10.1002/nme.1620140308",
    "simkin1980three": "10.1049/ip-b.1980.0052",
}

# Where the formulation is defined, and must say whose it is.
DEFINITION_SITES = (
    "src/radia/static_electromagnet.py",
    "src/radia/kelvin_solver.py",
    "packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/kelvin.py",
)


def _entry(key):
    text = BIB.read_text(encoding="utf-8")
    match = re.search(r"@(\w+)\{" + re.escape(key) + r",(.*?)\n\}\n",
                      text, re.DOTALL)
    return match


@pytest.mark.parametrize("key,doi", sorted(KEYS.items()))
def test_the_cited_entry_is_in_the_canonical_bibliography(key, doi):
    """One bibliography is the parent; the citation has to resolve in it."""
    match = _entry(key)
    assert match is not None, f"{key} is not in {BIB.name}"
    kind, body = match.group(1), match.group(2)
    assert kind == "article", (
        f"{key} is @{kind}; both papers appeared in journals, and typing a "
        f"journal paper as a proceedings misreports where it was published")
    assert doi in body, f"{key} does not carry its DOI"
    for field in ("journal", "volume", "pages", "year"):
        assert re.search(rf"\n\s*{field}\s*=", body), f"{key} has no {field}"


@pytest.mark.parametrize("relative", DEFINITION_SITES)
def test_the_formulation_says_whose_method_it_is(relative):
    """Every place that defines the split names the paper it comes from."""
    text = (ROOT / relative).read_text(encoding="utf-8")
    assert "mixed total/reduced Omega" in text, (
        f"{relative} no longer defines the formulation; update this list")
    assert "Simkin" in text and "Trowbridge" in text, (
        f"{relative} defines the formulation without attributing it")
    assert any(key in text for key in KEYS), (
        f"{relative} names the authors but gives no resolvable citation key")
