"""Map converged 2D AGE sweep rows to existing motor comparison contracts."""

from math import degrees, isfinite


def _finite_real(value):
    if isinstance(value, (bool, str)):
        raise ValueError("AGE lanes require finite real magnetostatic samples")
    try:
        number = complex(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("AGE lanes require finite real magnetostatic samples") from exc
    if number.imag != 0 or not isfinite(number.real):
        raise ValueError("AGE lanes require finite real magnetostatic samples")
    return number.real


def _sample_angles(rows, required_field, provenance_field):
    """Read converged sweep angles in radians and convert them to degrees."""
    if not isinstance(rows, list) or not rows:
        raise ValueError("provide a nonempty list of AGE sweep rows")
    angles = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"AGE row {index} must be an object")
        for key in ("converged", "theta_r", required_field, provenance_field):
            if key not in row:
                raise ValueError(f"AGE row {index} is missing {key}")
        if row["converged"] is not True:
            raise ValueError(f"AGE row {index} did not converge")
        angle = _finite_real(degrees(_finite_real(row["theta_r"])))
        if angles and angle <= angles[-1]:
            raise ValueError(f"AGE row {index} angle must be strictly increasing")
        angles.append(angle)
    return angles


def age_motor_lanes(*, back_emf_rows=None, cogging_rows=None, stack_length_m,
                    pole_pairs, phase_current_A, speed_rpm=None,
                    electrical_offset_deg=0.0, linkage_polarity=None, stack_basis=None,
                    slots=None, slot_pole_symmetry=None, torque_basis=None,
                    model_fraction=None, drop_duplicate_endpoint=False):
    """Emit comparison lanes from converged sweep rows with theta_r in radians."""
    if type(phase_current_A) not in (float, int) or phase_current_A != 0:
        raise ValueError("declare phase_current_A = 0 for the open-circuit sweep")
    if type(pole_pairs) is not int or pole_pairs < 1:
        raise ValueError("pole_pairs must be a positive integer")
    length = _finite_real(stack_length_m)
    if length <= 0:
        raise ValueError("stack_length_m must be positive")
    if type(drop_duplicate_endpoint) is not bool:
        raise ValueError("drop_duplicate_endpoint must be boolean")
    if back_emf_rows is None and cogging_rows is None:
        raise ValueError("provide back_emf_rows or cogging_rows")
    common = {"model_dimension": "2d", "stack_length_m": length,
              "source_angle_unit": "rad",
              "drop_duplicate_endpoint": drop_duplicate_endpoint}
    result = {}
    if back_emf_rows is not None:
        speed = _finite_real(speed_rpm)
        if speed <= 0 or stack_basis not in ("per_metre", "total"):
            raise ValueError("back-EMF requires positive speed_rpm and a declared stack_basis")
        if linkage_polarity not in ("positive_winding", "negative_winding"):
            raise ValueError("declare linkage_polarity for the signed abc winding")
        angles = _sample_angles(back_emf_rows, "phase_flux_wb", "linkage_stack_basis")
        flux = []
        for index, row in enumerate(back_emf_rows):
            if row["linkage_stack_basis"] != stack_basis:
                raise ValueError(f"AGE row {index} linkage_stack_basis does not match stack_basis")
            values = row["phase_flux_wb"]
            if not isinstance(values, (list, tuple)) or len(values) != 3:
                raise ValueError(f"AGE row {index} phase_flux_wb requires three signed abc linkages")
            flux.append([_finite_real(value) for value in values])
        result["back_emf_lane"] = {
            **common,
            "mechanical_angles_deg": angles, "phase_flux_wb": flux,
            "pole_pairs": pole_pairs, "speed_rpm": speed,
            "electrical_offset_deg": _finite_real(electrical_offset_deg),
            "linkage_polarity": linkage_polarity, "stack_basis": stack_basis,
            "emf_method": "central_difference_of_linkage",
        }

    if cogging_rows is not None:
        if type(slots) is not int or slots < 0 or type(slot_pole_symmetry) is not bool:
            raise ValueError("declare a nonnegative integer slots and boolean slot_pole_symmetry")
        if torque_basis not in ("per_metre", "total"):
            raise ValueError("declare torque_basis as per_metre or total")
        fraction = _finite_real(model_fraction)
        if not 0 < fraction <= 1:
            raise ValueError("model_fraction must be in (0, 1]")
        angles = _sample_angles(cogging_rows, "torque", "axial_length_m")
        expected_length = 1.0 if torque_basis == "per_metre" else length
        torque = []
        for index, row in enumerate(cogging_rows):
            axial_length = _finite_real(row["axial_length_m"])
            if axial_length != expected_length:
                raise ValueError(f"AGE row {index} axial_length_m does not match torque_basis")
            torque.append(_finite_real(row["torque"]))
        if slots == 0:
            result["cogging_unavailable"] = {
                "source_angle_unit": "rad",
                "reason": "slotless: the cogging gate requires a positive slot count",
                "mechanical_angles_deg": angles, "torque_nm": torque,
                "torque_basis": torque_basis, "stack_length_m": length,
                "model_fraction": fraction,
            }
        else:
            result["cogging_lane"] = {
                **common,
                "mechanical_angles_deg": angles, "torque_nm": torque,
                "slots": slots, "pole_count": 2 * pole_pairs,
                "slot_pole_symmetry": slot_pole_symmetry,
                "extraction_method": "maxwell_stress", "skew": "none",
                "torque_basis": torque_basis, "model_fraction": fraction,
                "excitation": {"phase_current_A": 0},
            }
    return result
