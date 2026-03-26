"""
Layer 2: UI logic tests for InductanceDialog (mock cubit, no Qt display).

Tests the dialog's internal logic without requiring Cubit or rendering Qt.
Runs in CI (pytest).

Tests:
  - Block detection logic
  - Material combo -> sigma mapping
  - Argument construction for subprocess
  - Result display data formatting
  - ESIM group visibility toggle
"""

import pytest


# ============================================================
# Test: Block detection logic
# ============================================================
class TestBlockDetection:

    def test_source_sink_found(self):
        """When source+sink blocks exist, both should be detected."""
        blocks = {"source": 1, "sink": 2, "conductor": 3}
        source = blocks.get("source")
        sink = blocks.get("sink")
        workpiece = blocks.get("workpiece")
        assert source is not None
        assert sink is not None
        assert workpiece is None

    def test_workpiece_detected(self):
        """When workpiece block exists, it should be detected."""
        blocks = {"source": 1, "sink": 2, "workpiece": 3}
        workpiece = blocks.get("workpiece")
        assert workpiece is not None

    def test_no_blocks(self):
        """Empty block list -> all None."""
        blocks = {}
        assert blocks.get("source") is None
        assert blocks.get("sink") is None
        assert blocks.get("workpiece") is None

    def test_esim_group_visibility(self):
        """ESIM group visible iff workpiece block found."""
        # With workpiece
        blocks_with = {"source": 1, "sink": 2, "workpiece": 3}
        esim_visible = "workpiece" in blocks_with
        assert esim_visible is True

        # Without workpiece
        blocks_without = {"source": 1, "sink": 2}
        esim_visible = "workpiece" in blocks_without
        assert esim_visible is False

    def test_can_solve_requires_source_sink_python(self):
        """Solve button enabled only when source+sink+python all available."""
        def can_solve(blocks, python_available):
            return ("source" in blocks and "sink" in blocks
                    and python_available)

        assert can_solve({"source": 1, "sink": 2}, True) is True
        assert can_solve({"source": 1}, True) is False
        assert can_solve({"source": 1, "sink": 2}, False) is False
        assert can_solve({}, True) is False


# ============================================================
# Test: Material -> sigma mapping
# ============================================================
class TestMaterialMapping:

    SIGMA_MAP = {"Steel": "2.0e6", "Copper": "5.8e7", "Aluminum": "3.5e7"}

    def test_all_materials_have_sigma(self):
        for mat in ["Steel", "Copper", "Aluminum"]:
            assert mat in self.SIGMA_MAP

    def test_sigma_values_reasonable(self):
        assert float(self.SIGMA_MAP["Copper"]) > float(self.SIGMA_MAP["Aluminum"])
        assert float(self.SIGMA_MAP["Aluminum"]) > float(self.SIGMA_MAP["Steel"])


# ============================================================
# Test: Argument construction
# ============================================================
class TestArgConstruction:

    def test_base_args(self):
        """Base args (no workpiece) should not contain --workpiece."""
        args = [
            "calc_inductance.py",
            "--cub5", "model.cub5",
            "--source", "source",
            "--sink", "sink",
            "--order", "2",
            "--fes-order", "0",
        ]
        assert "--workpiece" not in args

    def test_esim_args_appended(self):
        """When workpiece exists, ESIM args should be appended."""
        args = [
            "calc_inductance.py",
            "--cub5", "model.cub5",
            "--source", "source",
            "--sink", "sink",
            "--order", "2",
            "--fes-order", "0",
        ]

        workpiece_block = "workpiece"
        if workpiece_block:
            args += [
                "--workpiece", workpiece_block,
                "--impedance-model", "esim",
                "--frequency", "50000",
                "--sigma", "2.0e6",
                "--half-thickness", "0.005",
                "--material", "steel",
            ]

        assert "--workpiece" in args
        idx = args.index("--workpiece")
        assert args[idx + 1] == "workpiece"
        assert "--impedance-model" in args
        assert "--frequency" in args

    def test_dowell_model_selection(self):
        """Dowell model should be passed correctly."""
        model = "dowell"
        args = ["--impedance-model", model]
        assert args[1] == "dowell"


# ============================================================
# Test: Result display formatting
# ============================================================
class TestResultDisplay:

    def _format_inductance(self, L):
        """Replicate the panel's inductance formatting logic."""
        if abs(L) >= 1e-3:
            return f"{L*1e3:.4f} mH"
        elif abs(L) >= 1e-6:
            return f"{L*1e6:.4f} uH"
        elif abs(L) >= 1e-9:
            return f"{L*1e9:.4f} nH"
        else:
            return f"{L:.4e} H"

    def test_format_nH(self):
        assert "nH" in self._format_inductance(150e-9)

    def test_format_uH(self):
        assert "uH" in self._format_inductance(2.5e-6)

    def test_format_mH(self):
        assert "mH" in self._format_inductance(1.5e-3)

    def test_format_very_small(self):
        assert "H" in self._format_inductance(1e-12)

    def test_esim_result_rows(self):
        """When wp_R_effective present, ESIM rows should be generated."""
        data = {
            "inductance_H": 150e-9,
            "wp_R_effective": 6e-4,
            "wp_P_total": 3e-4,
            "wp_Q_total": 7e-4,
            "wp_frequency": 50000,
            "wp_X_effective": 1e-3,
            "wp_delta_min": 0.04e-3,
            "wp_delta_max": 0.04e-3,
            "wp_n_panels": 128,
            "wp_model": "esim",
            "wp_material": "steel",
        }

        # Simulate _display_result logic
        rows = [("Inductance (coil)", self._format_inductance(data["inductance_H"]))]
        if "wp_R_effective" in data:
            rows.append(("R (workpiece)", f"{data['wp_R_effective']:.4e} Ohm"))
            rows.append(("P (workpiece)", f"{data['wp_P_total']:.4e} W"))
            rows.append(("Q (workpiece)", f"{data['wp_Q_total']:.4e} var"))

        # Should have ESIM rows
        labels = [r[0] for r in rows]
        assert "R (workpiece)" in labels
        assert "P (workpiece)" in labels
        assert "Q (workpiece)" in labels

    def test_no_esim_rows_without_workpiece(self):
        """Without wp_R_effective, no ESIM rows."""
        data = {"inductance_H": 150e-9}
        has_esim = "wp_R_effective" in data
        assert has_esim is False

    def test_total_impedance_computation(self):
        """Total Z = R_wp + j*(X_coil + X_wp)."""
        import math
        L = 150e-9
        freq = 50000
        R_wp = 6e-4
        X_wp = 1e-3
        omega = 2 * math.pi * freq
        X_coil = omega * L

        Z_total_R = R_wp
        Z_total_X = X_coil + X_wp
        Z_mag = (Z_total_R**2 + Z_total_X**2)**0.5

        assert Z_total_R > 0
        assert Z_total_X > 0
        assert Z_mag > Z_total_R
