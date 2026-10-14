"""Contract: the 2D-parity loss recipes keep the driven-conductor field term."""

import re

from radia_mcp.radia_ngsolve.knowledge.magnetostatic_2d_parity import FEMM_MAGNETICS


def _loss_argument(call: str) -> str:
    match = re.search(rf"{call}\((\w+),", FEMM_MAGNETICS)
    assert match, f"{call} example missing"
    return match.group(1)


def test_planar_loss_uses_the_full_axial_field():
    name = _loss_argument("ohmic_loss_2d")
    assert re.search(rf"{name}\s*=\s*-1j\*omega\*Az \+ Vc\b", FEMM_MAGNETICS)


def test_axisymmetric_loss_keeps_the_winding_field():
    name = _loss_argument("ohmic_loss_axi")
    assert re.search(rf"{name}\s*=\s*-1j\*omega\*Az \+ Vc/x\b", FEMM_MAGNETICS)
    assert "ohmic_loss_axi(-1j*omega*Az" not in FEMM_MAGNETICS
