"""Panel-smoke scenarios catalogue (consumed by
.claude/skills/panel-smoke/smoke.py).

Each scenario drives a Layer-3 panel with specific widget values,
runs the subprocess command the panel builds, and verifies the JSON
output against a golden band.

Add new scenarios here as panels gain features.  The reference
scenarios MUST stay green (they are the paper-citable / release-
critical cases — see SKILL.md "Reference scenarios" table).
"""
from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
SAMPLES = REPO / "src" / "radia" / "panels" / "samples"


SCENARIOS = [
    # ------------------------------------------------------------------
    # IH panel — PEEC + BEM weak coupling, scalar ESIM
    # Reference: sweep_v2 scalar at I=100 A, f=50 kHz.
    # ------------------------------------------------------------------
    {
        "id": "ih_peec_bem_scalar_I100",
        "panel": "radia_ih",
        "window_class": "IHWindow",
        "method_const": "METHOD_PEEC_BEM",
        "vol": str(SAMPLES / "ih_bem_sample_p1.vol"),
        "set": {
            "impedance_model": "Nonlinear ESIM (experimental, WIP)",
            "peec_step": str(SAMPLES / "ih_fem_kelvin_demo_coil.step"),
            "bh_file":   str(SAMPLES / "em_sample_bh.txt"),
            "freq":            "50000",
            "current":         "100.0",
            "wp_sigma":        "2e6",
            "mu_r":            "100",
            "half_thickness":  "0.005",
            "esim_per_panel":  False,
            "esim_anderson_m": 5,
            "esim_max_iter":   30,
            "esim_relax":      "0.5",
            "fes_order":       1,
        },
        "builder": "_build_peec_bem_command",
        "golden": {
            "P_wp_W":          {"value": 30.51, "tol_rel": 0.03},
            "L_coil_nH":       {"value": 100.67, "tol_rel": 0.01},
            "esim_iterations": {"value": 6, "tol_abs": 3},
            "esim_converged":  {"value": True},
            "esim_per_panel":  {"value": False},
        },
        "timeout_s": 600,
    },
    # ------------------------------------------------------------------
    # IH panel — PEEC + BEM weak coupling, PER-ELEMENT ESIM
    # Reference: sweep_v2 per-panel at I=100 A, f=50 kHz.
    # This IS the IGTE 2026 paper Fig. 1 representative cell.
    # ------------------------------------------------------------------
    {
        "id": "ih_peec_bem_per_element_I100",
        "panel": "radia_ih",
        "window_class": "IHWindow",
        "method_const": "METHOD_PEEC_BEM",
        "vol": str(SAMPLES / "ih_bem_sample_p1.vol"),
        "set": {
            "impedance_model": "Nonlinear ESIM (experimental, WIP)",
            "peec_step": str(SAMPLES / "ih_fem_kelvin_demo_coil.step"),
            "bh_file":   str(SAMPLES / "em_sample_bh.txt"),
            "freq":            "50000",
            "current":         "100.0",
            "wp_sigma":        "2e6",
            "mu_r":            "100",
            "half_thickness":  "0.005",
            "esim_per_panel":  True,
            "esim_anderson_m": 5,
            "esim_max_iter":   30,
            "esim_relax":      "0.5",
            # fes_order is auto-clamped to 1 by _on_impedance_changed when
            # per_panel is on — leave unset to verify the clamp logic works.
        },
        "builder": "_build_peec_bem_command",
        "golden": {
            "P_wp_W":          {"value": 18.75, "tol_rel": 0.02},
            "L_coil_nH":       {"value": 100.67, "tol_rel": 0.01},
            "esim_iterations": {"value": 7, "tol_abs": 3},
            "esim_converged":  {"value": True},
            "esim_per_panel":  {"value": True},
            "esim_anderson_m": {"value": 5},
        },
        "timeout_s": 600,
    },
]
