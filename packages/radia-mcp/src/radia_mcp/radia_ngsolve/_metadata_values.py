"""Private metadata readers shared by domain gates, independent of MCP.

These preserve existing input conventions; parsing is not physical validation.
Keep distinct list conventions separate instead of silently broadening inputs.
"""

import re


def _norm(value):
    return str(value or "").strip().lower().replace("-", "_").replace(" ", "_")


def _first(row, names):
    for name in names:
        if name in row and row[name] is not None:
            return row[name]
    return None


def _string_list(value):
    if value in (None, ""):
        return []
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return []
        for sep in (";", "\n"):
            text = text.replace(sep, ",")
        return [item.strip() for item in text.split(",") if item.strip()]
    if isinstance(value, dict):
        # Keep the existing mapping contract, including dict subclasses.
        return [str(item).strip() for item in value.keys() if str(item).strip()]  # noqa: SIM118
    return [str(item).strip() for item in value if str(item).strip()]


def _unit_mapping(value):
    if value in (None, ""):
        return {}
    if isinstance(value, dict):
        return {
            str(key).strip(): str(unit).strip()
            for key, unit in value.items()
            if str(key).strip()
        }
    pairs = {}
    for item in _string_list(value):
        if ":" in item:
            key, unit = item.split(":", 1)
        elif "=" in item:
            key, unit = item.split("=", 1)
        else:
            continue
        key = key.strip()
        if key:
            pairs[key] = unit.strip()
    return pairs


def _coordinate_tuple(value):
    if value in (None, ""):
        return None
    if isinstance(value, dict):
        lower_keys = ("x", "y", "z")
        upper_keys = ("X", "Y", "Z")
        if all(key in value for key in lower_keys[:2]):
            coords = [value[key] for key in lower_keys if key in value]
        elif all(key in value for key in upper_keys[:2]):
            coords = [value[key] for key in upper_keys if key in value]
        else:
            raise ValueError("coordinate dictionaries must include x/y or X/Y")
    elif isinstance(value, str):
        coords = [item for item in re.split(r"[,;\s]+", value.strip()) if item]
    else:
        coords = list(value)
    if len(coords) not in (2, 3):
        raise ValueError("coordinates must contain two or three values")
    return tuple(float(coord) for coord in coords)


def _phase_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [part.strip() for part in value.replace(";", ",").split(",") if part.strip()]
    try:
        return [str(part).strip() for part in value if str(part).strip()]
    except TypeError:
        return [str(value).strip()]
