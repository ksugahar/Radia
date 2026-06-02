"""Journal full name <-> IEEE / IEICE / ISO 4 abbreviation table.

The list is short and curated for the electromagnetics / power
electronics fields where this lab publishes. Extend as needed.
"""
from __future__ import annotations


# Full name → {style: abbreviation}
# Keys are normalized to lower + alphanumeric for matching.
_TABLE = {
    "ieee transactions on magnetics": {
        "ieee": "IEEE Trans. Magn.",
        "ieice": "IEEE Trans. Magn.",
        "iso4": "IEEE Trans. Magn.",
    },
    "ieee transactions on industry applications": {
        "ieee": "IEEE Trans. Ind. Appl.",
        "ieice": "IEEE Trans. Ind. Appl.",
        "iso4": "IEEE Trans. Ind. Appl.",
    },
    "ieee transactions on power electronics": {
        "ieee": "IEEE Trans. Power Electron.",
        "ieice": "IEEE Trans. Power Electron.",
        "iso4": "IEEE Trans. Power Electron.",
    },
    "ieee transactions on industrial electronics": {
        "ieee": "IEEE Trans. Ind. Electron.",
        "ieice": "IEEE Trans. Ind. Electron.",
        "iso4": "IEEE Trans. Ind. Electron.",
    },
    "ieee transactions on antennas and propagation": {
        "ieee": "IEEE Trans. Antennas Propag.",
        "ieice": "IEEE Trans. Antennas Propag.",
        "iso4": "IEEE Trans. Antennas Propag.",
    },
    "ieee transactions on microwave theory and techniques": {
        "ieee": "IEEE Trans. Microw. Theory Tech.",
        "ieice": "IEEE Trans. Microw. Theory Tech.",
        "iso4": "IEEE Trans. Microw. Theory Tech.",
    },
    "ieee transactions on electromagnetic compatibility": {
        "ieee": "IEEE Trans. Electromagn. Compat.",
        "iso4": "IEEE Trans. Electromagn. Compat.",
    },
    "ieice transactions on electronics": {
        "ieee": "IEICE Trans. Electron.",
        "ieice": "IEICE Trans. Electron.",
        "iso4": "IEICE Trans. Electron.",
    },
    "ieee magnetics letters": {
        "ieee": "IEEE Magn. Lett.",
        "iso4": "IEEE Magn. Lett.",
    },
    "journal of magnetism and magnetic materials": {
        "ieee": "J. Magn. Magn. Mater.",
        "ieice": "J. Magn. Magn. Mater.",
        "iso4": "J. Magn. Magn. Mater.",
    },
    "physica b condensed matter": {
        "ieee": "Physica B Condens. Matter",
        "iso4": "Physica B Condens. Matter",
    },
    "compel the international journal for computation and mathematics in electrical and electronic engineering": {
        "ieee": "COMPEL",
        "ieice": "COMPEL",
        "iso4": "COMPEL",
    },
    "journal of computational and applied mathematics": {
        "ieee": "J. Comput. Appl. Math.",
        "iso4": "J. Comput. Appl. Math.",
    },
    "applied physics letters": {
        "ieee": "Appl. Phys. Lett.",
        "iso4": "Appl. Phys. Lett.",
    },
    "journal of applied physics": {
        "ieee": "J. Appl. Phys.",
        "iso4": "J. Appl. Phys.",
    },
    "physical review letters": {
        "ieee": "Phys. Rev. Lett.",
        "iso4": "Phys. Rev. Lett.",
    },
    "physical review b": {
        "ieee": "Phys. Rev. B",
        "iso4": "Phys. Rev. B",
    },
    "nature": {
        "ieee": "Nature",
        "iso4": "Nature",
    },
    "science": {
        "ieee": "Science",
        "iso4": "Science",
    },
    "nature communications": {
        "ieee": "Nat. Commun.",
        "iso4": "Nat. Commun.",
    },
    "nature physics": {
        "ieee": "Nat. Phys.",
        "iso4": "Nat. Phys.",
    },
}


import re


def _normalize(name: str) -> str:
    s = (name or "").lower()
    s = re.sub(r"[^a-z0-9]+", " ", s).strip()
    return s


def lookup(full_name: str, style: str = "ieee") -> str | None:
    """Return the abbreviation for ``full_name`` under ``style``.

    ``style`` is one of ``ieee``, ``ieice``, ``iso4``. Returns None if
    the journal is not in the table.
    """
    norm = _normalize(full_name)
    entry = _TABLE.get(norm)
    if not entry:
        return None
    return entry.get(style) or entry.get("ieee")


def reverse_lookup(abbrev: str) -> str | None:
    """Find the full journal name from an abbreviation (best-effort)."""
    target = _normalize(abbrev)
    for full, styles in _TABLE.items():
        for v in styles.values():
            if _normalize(v) == target:
                return full
    return None
