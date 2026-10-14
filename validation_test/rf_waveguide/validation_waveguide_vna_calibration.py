"""Validation-class rectangular-waveguide VNA calibration example.

The example uses a WR-90-like TE10 guide and an offset short:

* synthesize S11(f) for a shorted guide section
* recover group delay from unwrapped S11 phase
* convert group delay back to the physical short offset
* record TE/TM wave-impedance duality and dielectric-slab reflection metrics

It is intentionally an example/validation run, not a pytest test.

Run:

    python validation_test/rf_waveguide/validation_waveguide_vna_calibration.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from result_metadata import add_result_metadata  # noqa: E402

from radia_mcp.radia_ngsolve.waveguide import (  # noqa: E402
    C0,
    rectangular_waveguide_cutoff,
    reflection_metrics,
    sparameter_group_delay,
    waveguide_dielectric_slab_sparams,
    waveguide_dispersion,
    waveguide_offset_short_length_from_group_delay,
    waveguide_offset_short_s11,
    waveguide_wave_impedance,
)


OUT_JSON = HERE / "validation_waveguide_vna_calibration_summary.json"


def _sample_frequencies(start_hz=8.0e9, stop_hz=12.0e9, step_hz=10.0e6):
    n = int(round((stop_hz - start_hz) / step_hz))
    return [start_hz + i * step_hz for i in range(n + 1)]


def _complex_row(z: complex) -> dict:
    return {"re": z.real, "im": z.imag, "mag": abs(z), "phase_rad": math.atan2(z.imag, z.real)}


def main() -> int:
    width_a = 0.02286
    height_b = 0.01016
    offset_length = 0.040
    center_frequency = 10.0e9
    fc = rectangular_waveguide_cutoff(width_a, height_b, 1, 0)

    dispersion_rows = []
    eta0 = 4.0e-7 * math.pi * C0
    for f in (8.0e9, 10.0e9, 12.0e9):
        disp = waveguide_dispersion(f, fc)
        zte = waveguide_wave_impedance(f, fc, "TE")["Z"]
        ztm = waveguide_wave_impedance(f, fc, "TM")["Z"]
        dispersion_rows.append({
            "frequency": f,
            "beta": disp["beta"],
            "lambda_g": disp["lambda_g"],
            "v_group_over_c": disp["v_group"] / C0,
            "Z_TE": zte,
            "Z_TM": ztm,
            "Z_product_rel_error": abs(zte * ztm - eta0 * eta0) / (eta0 * eta0),
        })

    freqs = _sample_frequencies()
    trace = [waveguide_offset_short_s11(f, width_a, offset_length)["S11"] for f in freqs]
    delays = sparameter_group_delay(freqs, trace)
    recovered = [
        waveguide_offset_short_length_from_group_delay(f, width_a, tau)["offset_length"]
        for f, tau in zip(freqs, delays)
    ]
    interior_errors = [
        abs(length - offset_length) / offset_length
        for length in recovered[1:-1]
    ]
    sample_indices = [0, len(freqs) // 4, len(freqs) // 2, 3 * len(freqs) // 4, len(freqs) - 1]
    offset_rows = [{
        "frequency": freqs[i],
        "group_delay": delays[i],
        "recovered_offset_length": recovered[i],
        "relative_length_error": abs(recovered[i] - offset_length) / offset_length,
        "S11": _complex_row(trace[i]),
    } for i in sample_indices]

    center = waveguide_offset_short_s11(center_frequency, width_a, offset_length)
    center_delay = sparameter_group_delay(
        [center_frequency - 1.0e6, center_frequency, center_frequency + 1.0e6],
        [waveguide_offset_short_s11(f, width_a, offset_length)["S11"]
         for f in (center_frequency - 1.0e6, center_frequency, center_frequency + 1.0e6)],
    )[1]
    center_recovered = waveguide_offset_short_length_from_group_delay(
        center_frequency, width_a, center_delay)

    slab = waveguide_dielectric_slab_sparams(center_frequency, width_a, 2.2, 0.010)
    slab_metrics = reflection_metrics(slab["S11"])

    summary = {
        "kind": "waveguide_vna_calibration_validation",
        "validation_class": True,
        "guide": {
            "name": "WR-90-like",
            "width_a": width_a,
            "height_b": height_b,
            "te10_cutoff": fc,
        },
        "offset_short": {
            "true_length": offset_length,
            "n_frequency_samples": len(freqs),
            "frequency_start": freqs[0],
            "frequency_stop": freqs[-1],
            "frequency_step": freqs[1] - freqs[0],
            "max_interior_length_rel_error": max(interior_errors),
            "center_frequency": center_frequency,
            "center_analytic_group_delay": center["group_delay"],
            "center_sampled_group_delay": center_delay,
            "center_recovered_length": center_recovered["offset_length"],
            "sample_rows": offset_rows,
        },
        "dispersion_rows": dispersion_rows,
        "dielectric_slab_10ghz": {
            "eps_r": 2.2,
            "length": 0.010,
            "S11": _complex_row(slab["S11"]),
            "S21": _complex_row(slab["S21"]),
            "unitarity": slab["unitarity"],
            "return_loss_db": slab_metrics["return_loss_db"],
            "delivered_power_fraction": slab_metrics["delivered_power_fraction"],
            "vswr": slab_metrics["vswr"],
        },
    }
    summary = add_result_metadata(summary, __file__)
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"[guide] TE10 cutoff = {fc / 1e9:.6f} GHz")
    print(
        f"[offset short] length={offset_length:.6f} m, "
        f"max interior recovery error={max(interior_errors):.3e}"
    )
    print(
        f"[center] analytic tau={center['group_delay']:.12e} s, "
        f"sampled tau={center_delay:.12e} s, "
        f"recovered d={center_recovered['offset_length']:.12f} m"
    )
    print(
        f"[slab] |S11|={abs(slab['S11']):.6f}, |S21|={abs(slab['S21']):.6f}, "
        f"unitarity={slab['unitarity']:.12f}, return_loss={slab_metrics['return_loss_db']:.6f} dB"
    )
    print(f"[OK] wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
