"""Geometry-role normalization on IHDesignSpec.

The .vol / .step / .sol path fields are the inputs users re-point most
often; a crossed pair (coil .step in wp_vol, workpiece .vol in
peec_step) must be repaired deterministically by extension at the spec
boundary, recorded, and warned -- and anything without a unique repair
must fail immediately with the expected extensions spelled out.
"""

import warnings

import pytest

from radia.ih_design import (
    HEAT_SRC_SPATIAL,
    IHDesignSpec,
    METHOD_BEMA_BEM,
    METHOD_PEEC_BEM,
    METHOD_THERMAL_AXISYM,
    METHOD_THERMAL_3D_STATIC,
)


def test_swapped_wp_vol_and_peec_step_are_repaired():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        spec = IHDesignSpec(
            method=METHOD_PEEC_BEM,
            wp_vol="coil.step",
            peec_step="wp.vol",
        )
    assert spec.wp_vol == "wp.vol"
    assert spec.peec_step == "coil.step"
    assert len(spec.geometry_role_notes) == 2
    assert any("wp_vol" in note for note in spec.geometry_role_notes)
    assert any(
        issubclass(w.category, UserWarning)
        and "reassigned by extension" in str(w.message)
        for w in caught
    )


def test_stp_extension_counts_as_step():
    with pytest.warns(UserWarning, match="reassigned by extension"):
        spec = IHDesignSpec(
            method=METHOD_PEEC_BEM,
            wp_vol="coil.STP",
            peec_step="wp.VOL",
        )
    assert spec.wp_vol == "wp.VOL"
    assert spec.peec_step == "coil.STP"


def test_matching_inputs_untouched_without_notes():
    spec = IHDesignSpec(
        method=METHOD_BEMA_BEM,
        wp_vol="wp.vol",
        coil_vol="coil.vol",
    )
    assert spec.wp_vol == "wp.vol"
    assert spec.coil_vol == "coil.vol"
    assert spec.geometry_role_notes == ()


def test_vol_gz_accepted_for_mesh_slots():
    spec = IHDesignSpec(method=METHOD_PEEC_BEM,
                        wp_vol="wp.vol.gz", peec_step="coil.step")
    assert spec.geometry_role_notes == ()


def test_step_in_wp_vol_without_step_slot_raises():
    with pytest.raises(ValueError) as excinfo:
        IHDesignSpec(
            method=METHOD_BEMA_BEM,
            wp_vol="coil.step",
            coil_vol="wp.vol",
        )
    message = str(excinfo.value)
    assert "wp_vol" in message
    assert ".vol" in message
    assert "peec_step" in message  # the hint names the right slot


def test_unknown_extension_raises_with_expected_extensions():
    with pytest.raises(ValueError) as excinfo:
        IHDesignSpec(method=METHOD_PEEC_BEM, wp_vol="wp.msh",
                     peec_step="coil.step")
    assert ".vol" in str(excinfo.value)
    assert "wp.msh" in str(excinfo.value)


def test_thermal_qsurf_and_em_vol_swap_repaired():
    with pytest.warns(UserWarning, match="reassigned by extension"):
        spec = IHDesignSpec(
            method=METHOD_THERMAL_3D_STATIC,
            heat_source=HEAT_SRC_SPATIAL,
            qsurf_sol="model_em.vol",
            em_vol="model_q.sol",
        )
    assert spec.qsurf_sol == "model_q.sol"
    assert spec.em_vol == "model_em.vol"
    assert spec.geometry_role_notes


def test_build_command_normalizes_post_construction_mutation():
    spec = IHDesignSpec(method=METHOD_PEEC_BEM,
                        wp_vol="wp.vol", peec_step="coil.step")
    spec.wp_vol, spec.peec_step = spec.peec_step, spec.wp_vol
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        command = spec.build_command()
    assert spec.wp_vol == "wp.vol"
    assert "coil.step" in command
    index = command.index("--coil-step")
    assert command[index + 1] == "coil.step"


def test_normalization_is_idempotent():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec = IHDesignSpec(method=METHOD_PEEC_BEM,
                            wp_vol="coil.step", peec_step="wp.vol")
    notes_after_repair = spec.geometry_role_notes
    assert spec.normalize_geometry_roles() == ()
    assert spec.geometry_role_notes == notes_after_repair


def test_geometry_role_notes_roundtrip_accepts_list():
    spec = IHDesignSpec(geometry_role_notes=["earlier note"])
    assert spec.geometry_role_notes == ("earlier note",)


@pytest.mark.parametrize(
    ("method", "expected_order"),
    [
        (METHOD_THERMAL_AXISYM, "2"),
        (METHOD_THERMAL_3D_STATIC, "1"),
    ],
)
def test_thermal_fes_order_uses_method_specific_default(method, expected_order):
    command = IHDesignSpec(method=method, wp_vol="workpiece.vol").build_command(
        python="python", panels_dir="panels"
    )
    assert command[command.index("--fes-order") + 1] == expected_order


@pytest.mark.parametrize("method", [METHOD_THERMAL_AXISYM, METHOD_THERMAL_3D_STATIC])
def test_thermal_fes_order_explicit_override_is_preserved(method):
    command = IHDesignSpec(
        method=method,
        wp_vol="workpiece.vol",
        thermal_fes_order=3,
    ).build_command(python="python", panels_dir="panels")
    assert command[command.index("--fes-order") + 1] == "3"


@pytest.mark.parametrize("invalid", [True, 0, -1, 1.5, "two"])
def test_thermal_fes_order_rejects_invalid_values(invalid):
    spec = IHDesignSpec(
        method=METHOD_THERMAL_AXISYM,
        wp_vol="workpiece.vol",
        thermal_fes_order=invalid,
    )
    with pytest.raises(ValueError, match="positive integer"):
        spec.build_command(python="python", panels_dir="panels")
