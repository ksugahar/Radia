"""Read the physics contract from a legacy JMAG ``.jproj`` archive.

JMAG project files are ZIP containers.  This module deliberately extracts only
portable validation inputs (model/study identity, coil current conditions, and
magnetic-material B-H tables).  It does not attempt to run JMAG or translate
its mesh.  Geometry remains owned by the accompanying SAT model and must be
converted through the supported Cubit -> Netgen ``.vol`` boundary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


def _text(element: ET.Element | None, child_name: str) -> str:
    if element is None:
        return ""
    child = element.find(child_name)
    return "" if child is None or child.text is None else child.text.strip()


def _definition_index(root: ET.Element) -> dict[tuple[str, str], ET.Element]:
    result: dict[tuple[str, str], ET.Element] = {}
    for element in root.iter():
        key = element.attrib.get("key")
        if key and element.tag.startswith("T"):
            if (element.tag, key) in result:
                raise ValueError(f"duplicate JMAG definition {element.tag}:{key}")
            result[(element.tag, key)] = element
    return result


def _resolve(
    definitions: dict[tuple[str, str], ET.Element], key: str, *tags: str
) -> ET.Element | None:
    for tag in tags:
        value = definitions.get((tag, key))
        if value is not None:
            return value
    return None


def _property_value(element: ET.Element) -> tuple[object, str]:
    quantity = element.find("quantity")
    if quantity is not None:
        raw = _text(quantity, "value")
        unit = _text(quantity.find("unit"), "key")
        try:
            return float(raw), unit
        except ValueError:
            return raw, unit
    raw = _text(element, "value")
    try:
        return int(raw), ""
    except ValueError:
        try:
            return float(raw), ""
        except ValueError:
            return raw, ""


def _entity_ids(
    definitions: dict[tuple[str, str], ET.Element], reference: ET.Element | None
) -> list[int]:
    if reference is None:
        return []
    key = reference.attrib.get("key", "")
    entity_set = _resolve(
        definitions, key, "TCompositeEntitySet", "TSingleEntitySet")
    if entity_set is None:
        return []
    table = entity_set.find("./local/entity_list_table")
    if table is None:
        table = entity_set.find("entity_list_table")
    if table is None or not table.text:
        return []
    return [int(value) for value in table.text.split()]


def _current_conditions(
    root: ET.Element, definitions: dict[tuple[str, str], ET.Element]
) -> list[dict]:
    wanted = {
        "turn", "current", "constant_current", "flowin_type",
        "flg_uniform_current_density",
    }
    result = []
    property_tags = (
        "TRealProperty", "TFlagProperty", "TStringProperty", "TVectorProperty",
    )
    for condition in root.iter("TConditionData"):
        condition_type = condition.find("type")
        if condition_type is None or condition_type.attrib.get("name") != (
            "condition/CurrentCtrl"
        ):
            continue
        properties = {}
        for reference in condition.findall("property"):
            definition = _resolve(
                definitions, reference.attrib.get("key", ""), *property_tags)
            if definition is None:
                continue
            name = _text(definition, "name")
            if name in wanted:
                value, unit = _property_value(definition)
                properties[name] = {"value": value, "unit": unit}
        for name, unit in (("current", "A"), ("turn", "turn")):
            if name not in properties:
                raise ValueError(f"current condition is missing {name}")
            if properties[name]["unit"] != unit:
                raise ValueError(f"unsupported {name} unit: {properties[name]['unit']!r}")
        current = float(properties["current"]["value"])
        turns = float(properties["turn"]["value"])
        if not math.isfinite(current) or not math.isfinite(turns) or turns <= 0:
            raise ValueError("current and positive turns must be finite")
        ampere_turns = current * turns
        if not math.isfinite(ampere_turns):
            raise ValueError("ampere-turn product must be finite")
        result.append({
            "title": _text(condition, "title"),
            "entity_ids": _entity_ids(definitions, condition.find("parts")),
            "properties": properties,
            "ampere_turns": ampere_turns,
        })
    return result


def _persistent_table(
    definitions: dict[tuple[str, str], ET.Element], key: str
) -> tuple[list[list[float]], list[str]]:
    table = _resolve(definitions, key, "TPersistentTable")
    if table is None:
        return [], []
    rows = int(_text(table, "rows") or 0)
    cols = int(_text(table, "cols") or 0)
    values = [float(item.text) for item in table.findall("data") if item.text]
    if rows <= 0 or cols <= 0 or len(values) != rows * cols:
        raise ValueError(
            f"invalid JMAG persistent table {key}: rows={rows}, cols={cols}, "
            f"values={len(values)}")
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"non-finite JMAG persistent table {key}")
    matrix = [values[index:index + cols] for index in range(0, len(values), cols)]
    units = [item.text or "" for item in table.findall("col_units")]
    return matrix, units


def _bh_table(
    material_data: ET.Element,
    definitions: dict[tuple[str, str], ET.Element],
) -> tuple[list[list[float]], list[str]]:
    for reference in material_data.findall("property"):
        prop = _resolve(
            definitions, reference.attrib.get("key", ""), "TPointArrayProperty")
        if prop is None or _text(prop, "name") != "bhtable":
            continue
        point_reference = prop.find("point_array")
        if point_reference is None:
            return [], []
        point_array = _resolve(
            definitions, point_reference.attrib.get("key", ""), "TPointArray")
        if point_array is None:
            return [], []
        data_reference = point_array.find("data")
        if data_reference is None:
            return [], []
        return _persistent_table(
            definitions, data_reference.attrib.get("key", ""))
    return [], []


def _materials(
    root: ET.Element, definitions: dict[tuple[str, str], ET.Element]
) -> list[dict]:
    result = []
    for material in root.iter("TMaterial"):
        reference = material.find("property_list")
        if reference is None:
            continue
        data = _resolve(
            definitions, reference.attrib.get("key", ""), "TMaterialData")
        if data is None:
            continue
        name = ""
        for prop_reference in data.findall("property"):
            prop = _resolve(
                definitions, prop_reference.attrib.get("key", ""),
                "TStringProperty")
            if prop is not None and _text(prop, "name") == "text_name":
                name = _text(prop, "value")
                break
        bh_values, bh_units = _bh_table(data, definitions)
        result.append({
            "name": name,
            "id": _text(material, "id"),
            "read_only": bool(int(_text(material, "read_only") or 0)),
            "bh_units": bh_units,
            "bh": bh_values,
        })
    return result


def inspect_jmag_project(path: str | Path, *, include_sha256: bool = True) -> dict:
    """Return the portable physics manifest stored in one JMAG project."""
    source = Path(path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    with zipfile.ZipFile(source) as archive:
        names = set(archive.namelist())
        required = {"hierarchy.xml", "project.xml"}
        missing = required - names
        if missing:
            raise ValueError(f"JMAG archive is missing: {sorted(missing)}")
        hierarchy = ET.fromstring(archive.read("hierarchy.xml"))
        project = ET.fromstring(archive.read("project.xml"))
        entry_count = len(archive.infolist())
    definitions = _definition_index(project)
    models = hierarchy.findall("model")
    if len(models) != 1 or len(models[0].findall("study")) != 1:
        raise ValueError("legacy reader requires exactly one model and one study")
    model = models[0]
    study = model.find("study")
    document = {
        "schema": "radia.validation.legacy-jmag-ffag.v1",
        "source": str(source),
        "source_size_bytes": source.stat().st_size,
        "archive_entry_count": entry_count,
        "model": {
            "title": _text(model, "title"),
            "uuid": _text(model, "uuid"),
        },
        "study": {
            "title": _text(study, "title"),
            "uuid": _text(study, "uuid"),
            "cases": int(_text(study, "cases") or 0),
        },
        "current_conditions": _current_conditions(project, definitions),
        "materials": _materials(project, definitions),
    }
    if include_sha256:
        digest = hashlib.sha256()
        with source.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        document["source_sha256"] = digest.hexdigest()
    return document


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("jproj", type=Path)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--no-sha256", action="store_true")
    options = parser.parse_args()
    result = inspect_jmag_project(
        options.jproj, include_sha256=not options.no_sha256)
    payload = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if options.output_json:
        options.output_json.parent.mkdir(parents=True, exist_ok=True)
        options.output_json.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
