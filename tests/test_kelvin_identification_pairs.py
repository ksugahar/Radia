"""Contract for reading Netgen point identifications across Netgen 6.2.2606/2607.

6.2.2607 appends the identification number to each pair.  The reader must accept
both shapes, and a changed shape must raise instead of reading as "no Kelvin
identification" (an earlier catch-all hid exactly that change).
"""
from types import SimpleNamespace

import pytest

from radia.kelvin_identify_ngsolve import _identified_point_pairs, has_kelvin_identification


def _pid(nr):
    return SimpleNamespace(nr=nr)


def test_pairs_and_triples_read_the_same_point_numbers():
    pairs = [(_pid(3), _pid(7)), (_pid(4), _pid(8))]
    triples = [(_pid(3), _pid(7), 1), (_pid(4), _pid(8), 1)]
    assert _identified_point_pairs(pairs) == [(3, 7), (4, 8)]
    assert _identified_point_pairs(triples) == [(3, 7), (4, 8)]


@pytest.mark.parametrize("entry", [(_pid(1),), (_pid(1), _pid(2), 1, "x")])
def test_other_entry_shapes_raise(entry):
    with pytest.raises(ValueError, match="identification entry"):
        _identified_point_pairs([entry])


def test_mesh_errors_propagate_instead_of_reading_as_unidentified():
    class BrokenMesh:
        def GetBoundaries(self):
            raise RuntimeError("mesh access failed")

    with pytest.raises(RuntimeError, match="mesh access failed"):
        has_kelvin_identification(BrokenMesh())
