r"""2D->3D machine-global scaling (AxialLength x SymmetryFactor, ONELAB #1/#2/#4).

Pure-arithmetic, fast, no FE solve.  Locks that :class:`MachineScaling` reproduces
the hand-rolled ``MM_TORQUE_TO_NM = 1e-6 * (DEPTH_MM*1e-3)`` convention used in
test_motor_cogging_torque.py, that a 1/N sector scales x N, and that
``winding_resistance_onelab`` matches ONELAB's ``Rb[]`` form verbatim.
"""
import dataclasses
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.machine_scaling import (MachineScaling,
                                                     torque_scaling_summary,
                                                     winding_resistance_onelab)


def test_torque_matches_mm_convention():
    # test_motor_cogging_torque.py: MM_TORQUE_TO_NM = 1e-6 * (DEPTH_MM*1e-3),
    # i.e. mm mesh (length_unit_m=1e-3), 1 mm stack, full machine.
    sc = MachineScaling(stack_length=1e-3, symmetry_factor=1.0, length_unit_m=1e-3)
    assert math.isclose(sc.torque(1.0), 1e-6 * 1e-3, rel_tol=1e-15)
    # back out the raw 2D value that yields the stored 7.5919e-3 N*m reference
    tau_2d = 7.5919e-3 / (1e-6 * 1e-3)
    assert math.isclose(sc.torque(tau_2d), 7.5919e-3, rel_tol=1e-12)


def test_symmetry_factor_scales_globals():
    full = MachineScaling(stack_length=0.05, symmetry_factor=1.0)
    sector8 = MachineScaling(stack_length=0.05, symmetry_factor=8.0)
    assert math.isclose(sector8.global_scale, 8.0 * full.global_scale)
    # a 1/8 anti-periodic sector torque must be x8 to recover the whole machine
    assert math.isclose(sector8.torque(2.0), 8.0 * full.torque(2.0))


def test_global_scale_is_axial_times_symmetry():
    sc = MachineScaling(stack_length=0.12, symmetry_factor=4.0)
    assert math.isclose(sc.global_scale, 0.12 * 4.0)
    assert math.isclose(sc.flux_linkage(3.0), 3.0 * 0.12 * 4.0)
    # force carries one mesh-length power, torque two
    scmm = MachineScaling(stack_length=0.02, length_unit_m=1e-3)
    assert math.isclose(scmm.force(5.0), 5.0 * 1e-3 * 0.02)
    assert math.isclose(scmm.torque(5.0), 5.0 * 1e-6 * 0.02)


def test_torque_scaling_summary_spells_out_sector_to_whole_machine():
    row = torque_scaling_summary(
        2.5e6,
        stack_length_m=0.08,
        symmetry_factor=6.0,
        length_unit_m=1.0e-3,
        label="sixth_sector",
    )
    assert row["label"] == "sixth_sector"
    assert math.isclose(row["mesh_unit_torque_scale"], 1.0e-6)
    assert math.isclose(row["torque_one_model_sector_Nm"], 0.2)
    assert math.isclose(row["whole_machine_torque_Nm"], 1.2)
    assert math.isclose(row["whole_over_sector"], 6.0)
    assert math.isclose(row["total_scaling_factor"], 4.8e-7)

    braking = torque_scaling_summary(-2.5e6, 0.08, symmetry_factor=6.0, length_unit_m=1.0e-3)
    assert math.isclose(braking["whole_machine_torque_Nm"], -1.2)

    for kwargs in (
        {"stack_length_m": 0.0},
        {"stack_length_m": 0.08, "symmetry_factor": 0.0},
        {"stack_length_m": 0.08, "length_unit_m": 0.0},
    ):
        try:
            torque_scaling_summary(1.0, **kwargs)
            assert False, f"expected ValueError for {kwargs}"
        except ValueError:
            pass


def test_winding_resistance_replicates_onelab_Rb():
    N, surf, sigma, lstk = 20, 1.5e-4, 5.9e7, 0.05
    f3d, fill = 1.15, 0.96
    expect = f3d * lstk * fill * N * N / (surf * sigma)
    assert math.isclose(winding_resistance_onelab(N, surf, sigma, lstk, f3d, fill),
                        expect, rel_tol=1e-15)
    # factor_3d_effects >= 1 inflates R (end-turn copper the 2D model can't see)
    base = winding_resistance_onelab(N, surf, sigma, lstk, 1.0, 1.0)
    assert winding_resistance_onelab(N, surf, sigma, lstk, 1.15, 1.0) > base


def test_defaults_and_frozen():
    sc = MachineScaling(stack_length=0.05)
    assert sc.symmetry_factor == 1.0 and sc.length_unit_m == 1.0
    try:
        sc.stack_length = 1.0
        assert False, "MachineScaling should be frozen"
    except dataclasses.FrozenInstanceError:
        pass


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("[OK] machine_scaling (ONELAB #1/#2/#4 scale corrections) validated")


if __name__ == "__main__":
    main()
