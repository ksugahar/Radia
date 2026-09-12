"""Regression tests for the read-only legacy JMAG FFAG manifest reader."""

from __future__ import annotations

import zipfile

import pytest

from validation_test.ffag_topopt.legacy_jmag_ffag import inspect_jmag_project


HIERARCHY = """<?xml version="1.0" encoding="UTF-8"?>
<jproj>
  <model><title>FFAG F-D model</title><uuid>model-uuid</uuid>
    <study><title>static</title><uuid>study-uuid</uuid><cases>1</cases></study>
  </model>
</jproj>
"""


PROJECT = """<?xml version="1.0" encoding="UTF-8"?>
<designer>
  <TConditionData key="1">
    <type type="TStaticCurrentConditionType" key="2" name="condition/CurrentCtrl" />
    <parts type="TCompositeEntitySet" key="10" />
    <title>D coil</title>
    <property type="TRealProperty" key="20" />
    <property type="TRealProperty" key="21" />
    <property type="TFlagProperty" key="22" />
  </TConditionData>
  <TCompositeEntitySet key="10"><local><entity_list_table>12 13</entity_list_table></local></TCompositeEntitySet>
  <TRealProperty key="20"><name>turn</name><quantity><value>980</value><unit><key>turn</key></unit></quantity></TRealProperty>
  <TRealProperty key="21"><name>current</name><quantity><value>5</value><unit><key>A</key></unit></quantity></TRealProperty>
  <TFlagProperty key="22"><name>flowin_type</name><value>1</value></TFlagProperty>
  <TMaterial key="100"><property_list type="TMaterialData" key="101" /><id>mat-id</id><read_only>0</read_only></TMaterial>
  <TMaterialData key="101">
    <property type="TStringProperty" key="102" />
    <property type="TPointArrayProperty" key="103" />
  </TMaterialData>
  <TStringProperty key="102"><name>text_name</name><value>nkj1</value></TStringProperty>
  <TPointArrayProperty key="103"><name>bhtable</name><point_array type="TPointArray" key="104" /></TPointArrayProperty>
  <TPointArray key="104"><data type="TPersistentTable" key="105" /></TPointArray>
  <TPersistentTable key="105">
    <rows>3</rows><cols>2</cols>
    <data>0</data><data>0</data><data>40</data><data>0.13</data><data>80</data><data>0.83</data>
    <col_units>A/m</col_units><col_units>tesla</col_units>
  </TPersistentTable>
</designer>
"""


def _write_project(path, *, include_project: bool = True):
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("hierarchy.xml", HIERARCHY)
        if include_project:
            archive.writestr("project.xml", PROJECT)


def test_extracts_coil_ampere_turns_and_material_bh_table(tmp_path):
    path = tmp_path / "ffag.jproj"
    _write_project(path)

    result = inspect_jmag_project(path, include_sha256=False)

    assert result["model"]["title"] == "FFAG F-D model"
    assert result["study"] == {
        "title": "static", "uuid": "study-uuid", "cases": 1}
    assert result["current_conditions"] == [{
        "title": "D coil",
        "entity_ids": [12, 13],
        "properties": {
            "turn": {"value": 980.0, "unit": "turn"},
            "current": {"value": 5.0, "unit": "A"},
            "flowin_type": {"value": 1, "unit": ""},
        },
        "ampere_turns": 4900.0,
    }]
    material = result["materials"][0]
    assert material["name"] == "nkj1"
    assert material["bh_units"] == ["A/m", "tesla"]
    assert material["bh"] == [[0.0, 0.0], [40.0, 0.13], [80.0, 0.83]]


def test_rejects_an_archive_without_project_xml(tmp_path):
    path = tmp_path / "incomplete.jproj"
    _write_project(path, include_project=False)

    with pytest.raises(ValueError, match="project.xml"):
        inspect_jmag_project(path, include_sha256=False)


@pytest.mark.parametrize("old,new,match", [
    ('<name>current</name>', '<name>other</name>', 'missing current'),
    ('<key>A</key>', '<key>mA</key>', 'unsupported current unit'),
    ('<value>5</value>', '<value>nan</value>', 'must be finite'),
    ('<value>980</value>', '<value>0</value>', 'positive turns'),
    ('<data>0.83</data>', '<data>inf</data>', 'non-finite'),
])
def test_invalid_physics_does_not_become_a_zero_current_manifest(
        tmp_path, old, new, match):
    path = tmp_path / 'invalid.jproj'
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('hierarchy.xml', HIERARCHY)
        archive.writestr('project.xml', PROJECT.replace(old, new))
    with pytest.raises(ValueError, match=match):
        inspect_jmag_project(path)


def test_multiple_models_are_not_reported_as_one_physics_contract(tmp_path):
    path = tmp_path / 'multiple.jproj'
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('hierarchy.xml', HIERARCHY.replace('</jproj>', '<model /></jproj>'))
        archive.writestr('project.xml', PROJECT)
    with pytest.raises(ValueError, match='exactly one model'):
        inspect_jmag_project(path)


def test_source_hash_identifies_the_archive(tmp_path):
    import hashlib
    path = tmp_path / 'ffag.jproj'
    _write_project(path)
    assert inspect_jmag_project(path)['source_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
