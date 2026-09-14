"""Boundary contracts for shared, transport-independent metadata readers."""

import pytest
from radia_mcp.radia_ngsolve import _metadata_values as values


@pytest.mark.parametrize("reader,argument,expected", [
    (values._norm, None, ""),
    (values._norm, 0, ""),
    (values._norm, " AC-Field Map ", "ac_field_map"),
    (values._string_list, None, []),
    (values._string_list, " a;b\nc,d ", ["a", "b", "c", "d"]),
    (values._string_list, {" a ": 1, "": 2}, ["a"]),
    (values._string_list, [0, " b ", ""], ["0", "b"]),
    (values._unit_mapping, None, {}),
    (values._unit_mapping, {" x ": " m ", "": "s"}, {"x": "m"}),
    (values._unit_mapping, "x:m;y=N;ignored;x:mm", {"x": "mm", "y": "N"}),
    (values._coordinate_tuple, None, None),
    (values._coordinate_tuple, "", None),
    (values._coordinate_tuple, {"x": 1, "y": 2}, (1.0, 2.0)),
    (values._coordinate_tuple, {"X": 1, "Y": 2, "Z": 3}, (1.0, 2.0, 3.0)),
    (values._coordinate_tuple, "1;2 3", (1.0, 2.0, 3.0)),
    (values._coordinate_tuple, [1, 2], (1.0, 2.0)),
    (values._phase_list, None, []),
    (values._phase_list, " U;V,W ", ["U", "V", "W"]),
    (values._phase_list, 0, ["0"]),
    (values._phase_list, [" U ", "", 2], ["U", "2"]),
    # Phase names and generic lists deliberately retain different newline rules.
    (values._phase_list, "U\nV", ["U\nV"]),
])
def test_reader_contract(reader, argument, expected):
    assert reader(argument) == expected


@pytest.mark.parametrize("argument,message", [
    ({"x": 1}, "coordinate dictionaries must include x/y or X/Y"),
    ({"x": 1, "Y": 2}, "coordinate dictionaries must include x/y or X/Y"),
    ([1], "coordinates must contain two or three values"),
    ([1, 2, 3, 4], "coordinates must contain two or three values"),
    ("   ", "coordinates must contain two or three values"),
])
def test_coordinate_shape_errors(argument, message):
    with pytest.raises(ValueError, match=message):
        values._coordinate_tuple(argument)


@pytest.mark.parametrize("value", [0, False, "", [], {}])
def test_first_preserves_falsey_values(value):
    assert values._first({"skip": None, "chosen": value, "later": 9},
                         ["missing", "skip", "chosen", "later"]) == value


def test_first_returns_none_when_no_value_exists():
    assert values._first({"empty": None}, ["missing", "empty"]) is None
