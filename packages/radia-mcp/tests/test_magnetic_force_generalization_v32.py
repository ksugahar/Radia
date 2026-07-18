from __future__ import annotations

from test_magnetic_force_generalization_v31 import _gate, _identity_v31


_PROMOTED_CASE_IDS = (
    "v32_public_axisymmetric_force_weighted_stress_coenergy_contour_radius_mesh_mismatch",
    "v32_public_laminated_diffusion_skin_depth_complex_power_frequency_mesh_mismatch",
)


def _identity_v32():
    identity = _identity_v31()
    generation = "axisymmetric-force-closure-191"
    identity[
        "axisymmetric_weighted_stress_coenergy_contour_displacement_radius_weight_material_mesh_owner_result_identity"
    ] = {
        "force_generation": generation,
        **{
            key: generation
            for key in (
                "stress_generation",
                "coenergy_generation",
                "contour_generation",
                "displacement_generation",
                "weight_generation",
                "material_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "stress_method": "weighted_stress_tensor",
        "result_stress_method": "weighted_stress_tensor",
        "coenergy_method": "symmetric_virtual_displacement",
        "result_coenergy_method": "symmetric_virtual_displacement",
        "contour_radius_m": 0.025,
        "result_contour_radius_m": 0.025,
        "virtual_displacement_m": 1.0e-5,
        "result_virtual_displacement_m": 1.0e-5,
        "axisymmetric_weight": "2*pi*r",
        "result_axisymmetric_weight": "2*pi*r",
        "material_side": "air_gap",
        "result_material_side": "air_gap",
        "weighted_stress_force_n": 12.5,
        "coenergy_force_n": 12.5,
        "force_relative_tolerance": 1.0e-8,
        "force_mesh_sha256": "1" * 64,
        "result_force_mesh_sha256": "1" * 64,
        "force_owner": "axisymmetric/plunger/group1",
        "result_force_owner": "axisymmetric/plunger/group1",
        "force_result_sha256": "2" * 64,
        "accepted_force_result_sha256": "2" * 64,
    }
    generation = "laminated-diffusion-closure-191"
    conductivity = [[2.0e6, 0.0], [0.0, 1.0e3]]
    identity[
        "laminated_diffusion_conductivity_skin_depth_frequency_phasor_power_volume_mesh_loss_result_identity"
    ] = {
        "diffusion_generation": generation,
        **{
            key: generation
            for key in (
                "conductivity_generation",
                "skin_depth_generation",
                "frequency_generation",
                "phasor_generation",
                "power_generation",
                "volume_generation",
                "mesh_generation",
                "loss_generation",
                "result_generation",
            )
        },
        "conductivity_tensor_s_per_m": conductivity,
        "result_conductivity_tensor_s_per_m": [row[:] for row in conductivity],
        "skin_depth_m": 5.0e-4,
        "result_skin_depth_m": 5.0e-4,
        "frequency_hz": 400.0,
        "result_frequency_hz": 400.0,
        "phasor_convention": "exp(+jwt)_rms",
        "result_phasor_convention": "exp(+jwt)_rms",
        "complex_power_va_ri": [4.7, 1.2],
        "result_complex_power_va_ri": [4.7, 1.2],
        "active_volume_m3": 9.5e-4,
        "result_active_volume_m3": 9.5e-4,
        "laminated_mesh_sha256": "3" * 64,
        "result_laminated_mesh_sha256": "3" * 64,
        "laminated_loss_w": 4.7,
        "result_laminated_loss_w": 4.7,
        "laminated_result_sha256": "4" * 64,
        "accepted_laminated_result_sha256": "4" * 64,
    }
    return identity


def test_v32_public_positive_axisymmetric_force_and_laminated_diffusion_closure():
    assert _gate(_identity_v32())["status"] == "ok"


def test_v32_public_axisymmetric_force_weighted_stress_coenergy_contour_radius_mesh_mismatch():
    identity = _identity_v32()
    record = identity[
        "axisymmetric_weighted_stress_coenergy_contour_displacement_radius_weight_material_mesh_owner_result_identity"
    ]
    record.update(
        {
            "stress_generation": "axisymmetric-force-190",
            "mesh_generation": "axisymmetric-force-189",
            "result_generation": "axisymmetric-force-188",
            "result_stress_method": "contour_maxwell_stress",
            "result_coenergy_method": "forward_displacement",
            "result_contour_radius_m": 0.03,
            "result_virtual_displacement_m": 1.0e-3,
            "result_axisymmetric_weight": "1",
            "result_material_side": "steel",
            "coenergy_force_n": -10.0,
            "result_force_mesh_sha256": "9" * 64,
            "result_force_owner": "planar/group2",
            "accepted_force_result_sha256": "a" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_force_closes_weighted_stress_coenergy_contour_displacement_weight_material_mesh_owner_and_result"
    ]


def test_v32_public_laminated_diffusion_skin_depth_complex_power_frequency_mesh_mismatch():
    identity = _identity_v32()
    record = identity[
        "laminated_diffusion_conductivity_skin_depth_frequency_phasor_power_volume_mesh_loss_result_identity"
    ]
    record.update(
        {
            "conductivity_generation": "laminated-diffusion-190",
            "frequency_generation": "laminated-diffusion-189",
            "result_generation": "laminated-diffusion-188",
            "result_conductivity_tensor_s_per_m": [[1.0e3, 0.0], [0.0, 2.0e6]],
            "result_skin_depth_m": 5.0e-3,
            "result_frequency_hz": 50.0,
            "result_phasor_convention": "exp(-jwt)_peak",
            "result_complex_power_va_ri": [1.2, -4.7],
            "result_active_volume_m3": 1.0e-3,
            "result_laminated_mesh_sha256": "b" * 64,
            "result_laminated_loss_w": 9.0,
            "accepted_laminated_result_sha256": "c" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "laminated_diffusion_uses_current_conductivity_skin_depth_frequency_phasor_power_volume_mesh_loss_and_result"
    ]
