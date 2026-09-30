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
    METHOD_PEEC_FEM_KELVIN,
    METHOD_FEM_FULL,
    METHOD_THERMAL_AXISYM,
    METHOD_THERMAL_3D_STATIC,
)


@pytest.mark.parametrize("order", [1, 2, 3])
def test_hcurl_design_preserves_automatic_ams_selection(order):
    command = IHDesignSpec(method=METHOD_PEEC_FEM_KELVIN,
                           wp_vol="workpiece.vol", peec_step="coil.step",
                           fes_order=order).build_command(python="python", panels_dir="panels")
    assert command[command.index("--solver") + 1] == "auto"


def test_compound_av_design_uses_direct_and_rejects_unsupported_ams():
    spec = IHDesignSpec(method=METHOD_FEM_FULL)
    assert spec._fem_solver() == "sparsecholesky"
    spec.solver = "AMS (iterative, p=1)"
    with pytest.raises(ValueError, match="A-V"):
        spec._fem_solver()


def test_hcurl_design_honours_direct_choice_and_rejects_unknown_solver():
    spec = IHDesignSpec(method=METHOD_PEEC_FEM_KELVIN,
                        solver="sparsecholesky (direct)")
    assert spec._fem_solver() == "sparsecholesky"
    spec.solver = "misspelled"
    with pytest.raises(ValueError, match="Unknown FEM solver"):
        spec._fem_solver()


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
    command = IHDesignSpec(
        method=method,
        wp_vol="workpiece.vol",
        heat_flux_boundaries="heated",
        convection_boundaries="exposed",
    ).build_command(
        python="python", panels_dir="panels"
    )
    assert command[command.index("--fes-order") + 1] == expected_order


@pytest.mark.parametrize("method", [METHOD_THERMAL_AXISYM, METHOD_THERMAL_3D_STATIC])
def test_thermal_fes_order_explicit_override_is_preserved(method):
    command = IHDesignSpec(
        method=method,
        wp_vol="workpiece.vol",
        heat_flux_boundaries="heated",
        convection_boundaries="exposed",
        thermal_fes_order=3,
    ).build_command(python="python", panels_dir="panels")
    assert command[command.index("--fes-order") + 1] == "3"


def test_axisym_spatial_heat_uses_converged_phi_default_and_no_3d_axis_flag():
    command = IHDesignSpec(
        method=METHOD_THERMAL_AXISYM,
        wp_vol="workpiece_axisym.vol",
        heat_source=HEAT_SRC_SPATIAL,
        qsurf_sol="qsurf.sol",
        em_vol="em.vol",
        heat_flux_boundaries="heated",
        convection_boundaries="exposed",
    ).build_command(python="python", panels_dir="panels")

    assert command[command.index("--n-phi-samples") + 1] == "128"
    assert "--rotation-axis" not in command


@pytest.mark.parametrize("invalid", [True, 0, -1, 1.5, "two"])
def test_axisym_spatial_heat_rejects_invalid_phi_sample_count(invalid):
    spec = IHDesignSpec(
        method=METHOD_THERMAL_AXISYM,
        wp_vol="workpiece_axisym.vol",
        heat_source=HEAT_SRC_SPATIAL,
        qsurf_sol="qsurf.sol",
        em_vol="em.vol",
        n_phi_samples=invalid,
        heat_flux_boundaries="heated",
        convection_boundaries="exposed",
    )
    with pytest.raises(ValueError, match="positive integer"):
        spec.build_command(python="python", panels_dir="panels")


@pytest.mark.parametrize("method", [METHOD_THERMAL_AXISYM, METHOD_THERMAL_3D_STATIC])
def test_thermal_spatial_output_is_gmsh_only(method):
    command = IHDesignSpec(
        method=method,
        wp_vol="workpiece.vol",
        heat_flux_boundaries="heated",
        convection_boundaries="exposed",
    ).build_command(python="python", panels_dir="panels")
    assert command[command.index("--msh-output") + 1].endswith("_heat.msh")
    assert "--vtu-prefix" not in command
    assert not hasattr(IHDesignSpec(), "vtu_prefix")


@pytest.mark.parametrize("invalid", [True, 0, -1, 1.5, "two"])
def test_thermal_fes_order_rejects_invalid_values(invalid):
    spec = IHDesignSpec(
        method=METHOD_THERMAL_AXISYM,
        wp_vol="workpiece.vol",
        thermal_fes_order=invalid,
    )
    with pytest.raises(ValueError, match="positive integer"):
        spec.build_command(python="python", panels_dir="panels")


def test_thermal_boundary_roles_are_independent_cli_settings():
    spec = IHDesignSpec(
        method=METHOD_THERMAL_AXISYM,
        wp_vol="workpiece.vol",
        heat_flux_boundaries="heated_outer|heated_end",
        convection_boundaries="outer|top|bottom",
        radiation_boundaries="outer",
        emissivity="0.8",
    )
    command = spec.build_command(python="python", panels_dir="panels")

    assert command[command.index("--heat-flux-boundaries") + 1] == (
        "heated_outer|heated_end"
    )
    assert command[command.index("--convection-boundaries") + 1] == (
        "outer|top|bottom"
    )
    assert command[command.index("--radiation-boundaries") + 1] == "outer"
    assert "--surface-label" not in command


def test_thermal_boundary_roles_fail_fast_when_active_role_is_missing():
    spec = IHDesignSpec(
        method=METHOD_THERMAL_AXISYM,
        wp_vol="workpiece.vol",
        heat_flux_boundaries="heated",
    )
    assert spec.missing_required_inputs() == ["Convection boundary selector"]
    with pytest.raises(ValueError, match="convection_boundaries"):
        spec.build_command(python="python", panels_dir="panels")


def test_legacy_surface_label_is_rejected_with_migration_guidance():
    with pytest.raises(ValueError, match="heat_flux_boundaries"):
        IHDesignSpec(surface_label="outer")


def test_spatial_qsurf_high_order_is_rejected_before_command_execution():
    spec = IHDesignSpec(
        method=METHOD_THERMAL_AXISYM,
        wp_vol="workpiece.vol",
        heat_source=HEAT_SRC_SPATIAL,
        qsurf_sol="q.sol",
        em_vol="em.vol",
        qsurf_order=2,
        heat_flux_boundaries="heated",
        convection_boundaries="exposed",
    )
    with pytest.raises(ValueError, match="qsurf_order=1"):
        spec.build_command(python="python", panels_dir="panels")


def test_overheating_constraint_and_material_table_reach_the_heat_solver():
    command = IHDesignSpec(
        method=METHOD_THERMAL_AXISYM,
        wp_vol="workpiece_axisym.vol",
        heat_source=HEAT_SRC_SPATIAL,
        qsurf_sol="qsurf.sol",
        em_vol="em.vol",
        em_heat_boundaries="sibc",
        heat_flux_boundaries="heated",
        convection_boundaries="exposed",
        exposure_thresholds="850,1400",
        temperature_limit="1450",
        thermal_material_table="steel_kcp.csv",
    ).build_command(python="python", panels_dir="panels")
    assert command[command.index("--temperature-limit") + 1] == "1450"
    assert command[command.index("--exposure-thresholds") + 1] == "850,1400"
    assert command[command.index("--material-table") + 1] == "steel_kcp.csv"
    assert command[command.index("--em-heat-boundaries") + 1] == "sibc"


def test_thermal_constraint_fields_are_visible_and_halvings_are_explicit():
    from dataclasses import replace
    spec = IHDesignSpec(
        method=METHOD_THERMAL_AXISYM, wp_vol="workpiece_axisym.vol",
        heat_source=HEAT_SRC_SPATIAL, qsurf_sol="qsurf.sol", em_vol="em.vol",
        heat_flux_boundaries="heated", convection_boundaries="exposed")
    visible = spec.visible_fields()
    assert {"exposure_thresholds", "temperature_limit",
            "thermal_material_table", "em_heat_boundaries"} <= visible
    assert "max_halvings" not in visible          # linear solve: no Newton
    command = spec.build_command(python="python", panels_dir="panels")
    assert "--max-halvings" not in command
    with pytest.raises(ValueError, match="material-table"):
        replace(spec, max_halvings=2).build_command(
            python="python", panels_dir="panels")
    nonlinear = replace(spec, max_halvings=2,
                        thermal_material_table="kcp.csv")
    assert "max_halvings" in nonlinear.visible_fields()
    command = nonlinear.build_command(python="python", panels_dir="panels")
    assert command[command.index("--max-halvings") + 1] == "2"
    with pytest.raises(ValueError, match="non-negative integer"):
        replace(nonlinear, max_halvings=1.5) \
            .build_command(python="python", panels_dir="panels")


def test_rotor_states_reach_the_rotating_heat_solver():
    from dataclasses import replace
    from radia.ih_design import METHOD_THERMAL_3D_ROTATING
    spec = IHDesignSpec(
        method=METHOD_THERMAL_3D_ROTATING, wp_vol="workpiece.vol",
        heat_source=HEAT_SRC_SPATIAL, rotor_states="rotor.json",
        rotation_rpm=60.0, heat_flux_boundaries="heated",
        convection_boundaries="exposed")
    assert "rotor_states" in spec.visible_fields()
    command = spec.build_command(python="python", panels_dir="panels")
    assert command[command.index("--rotor-states") + 1] == "rotor.json"
    assert "--qsurf-sol" not in command
    for bad, match in (({"qsurf_sol": "q.sol"}, "leave qsurf_sol"),
                       ({"rotation_rpm": 0.0}, "rotation_rpm > 0"),
                       ({"method": METHOD_THERMAL_AXISYM,
                         "wp_vol": "workpiece_axisym.vol"}, "rotating 3D")):
        with pytest.raises(ValueError, match=match):
            replace(spec, **bad).build_command(python="python",
                                               panels_dir="panels")


def test_a_non_numeric_temperature_limit_fails_early():
    with pytest.raises(ValueError):
        IHDesignSpec(
            method=METHOD_THERMAL_3D_STATIC,
            wp_vol="workpiece.vol",
            heat_flux_boundaries="heated",
            convection_boundaries="exposed",
            temperature_limit="hot",
        ).build_command(python="python", panels_dir="panels")


@pytest.mark.parametrize("method", [METHOD_THERMAL_AXISYM, METHOD_THERMAL_3D_STATIC])
def test_thermal_convection_map_path_is_preserved_in_command(method):
    cooling_path = "case with spaces/cooling.json"
    spec = IHDesignSpec(method=method, wp_vol="workpiece.vol",
                        heat_flux_boundaries="heated", convection_map=cooling_path)
    command = spec.build_command(python="python", panels_dir="panels")
    assert command[command.index("--convection-map") + 1] == cooling_path
    assert "--convection-boundaries" not in command
    assert "Convection boundary selector" not in spec.missing_required_inputs()


def test_thermal_convection_map_rejects_uniform_selector_mix():
    spec = IHDesignSpec(method=METHOD_THERMAL_AXISYM, wp_vol="workpiece.vol",
                        heat_flux_boundaries="heated", convection_map="cooling.json",
                        convection_boundaries="air_cooling")
    with pytest.raises(ValueError, match="mutually exclusive"):
        spec.build_command(python="python", panels_dir="panels")
