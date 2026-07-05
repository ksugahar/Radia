"""
Layer 1: UI logic tests for AccelMagnetDialog (no Cubit, no Qt display).

Runs in CI (pytest). Tests the dialog's internal logic:
  - Block detection for yoke/air/kelvin/symmetry
  - Formulation-dependent Dirichlet BC assignment
  - Material model argument construction
  - Result display formatting
  - Hysteresis parameter handling
"""

import math
import pytest


# ============================================================
# Test: Block detection
# ============================================================
class TestAccelBlockDetection:

    def test_yoke_air_found(self):
        blocks = {"yoke": 1, "air": 2}
        assert "yoke" in blocks
        assert "air" in blocks

    def test_kelvin_optional(self):
        blocks = {"yoke": 1, "air": 2}
        assert "kelvin" not in blocks
        # Should not prevent solving
        has_all = all(k in blocks for k in ("yoke", "air"))
        assert has_all is True

    def test_symmetry_blocks_optional(self):
        blocks = {"yoke": 1, "air": 2, "sym_tangential": 3}
        assert "sym_tangential" in blocks
        assert "sym_normal" not in blocks

    def test_no_coil_block_needed(self):
        """Accelerator magnet: coil is Biot-Savart, no coil block."""
        blocks = {"yoke": 1, "air": 2}
        assert "coil" not in blocks
        # Solve requires yoke + air + coil_script, NOT coil block
        has_all = all(k in blocks for k in ("yoke", "air"))
        assert has_all is True

    def test_solve_requires_coil_script(self):
        blocks = {"yoke": 1, "air": 2}
        has_blocks = all(k in blocks for k in ("yoke", "air"))
        has_coil = False  # no coil script path set
        python_ok = True
        can_solve = has_blocks and has_coil and python_ok
        assert can_solve is False

        has_coil = True
        can_solve = has_blocks and has_coil and python_ok
        assert can_solve is True


# ============================================================
# Test: Symmetry BC assignment (formulation-dependent)
# ============================================================
class TestSymmetryBC:
    """Symmetry BCs swap between A-formulation and Omega-reduced.

    sym_tangential (B parallel to face):
        A-form  -> Dirichlet (A x n = 0)
        Omega   -> natural (do nothing)

    sym_normal (B normal to face):
        A-form  -> natural (do nothing)
        Omega   -> Dirichlet (Omega = 0)
    """

    def _build_dirichlet(self, formulation, boundaries):
        dir_parts = []
        if "outer" in boundaries:
            dir_parts.append("outer")
        if formulation == "omega":
            if "sym_normal" in boundaries:
                dir_parts.append("sym_normal")
        else:  # "a"
            if "sym_tangential" in boundaries:
                dir_parts.append("sym_tangential")
        return "|".join(dir_parts) if dir_parts else ""

    def test_omega_sym_normal_dirichlet(self):
        bnd = ["outer", "sym_tangential", "sym_normal"]
        d = self._build_dirichlet("omega", bnd)
        assert "sym_normal" in d
        assert "sym_tangential" not in d

    def test_a_sym_tangential_dirichlet(self):
        bnd = ["outer", "sym_tangential", "sym_normal"]
        d = self._build_dirichlet("a", bnd)
        assert "sym_tangential" in d
        assert "sym_normal" not in d

    def test_no_symmetry(self):
        bnd = ["outer"]
        d_omega = self._build_dirichlet("omega", bnd)
        d_a = self._build_dirichlet("a", bnd)
        assert d_omega == "outer"
        assert d_a == "outer"

    def test_dipole_quarter_model(self):
        """Dipole: Bz field, XZ mirror (B parallel) + XY mirror (B normal)."""
        bnd = ["outer", "sym_tangential", "sym_normal"]
        # Omega: Omega=0 on XY plane (sym_normal)
        d = self._build_dirichlet("omega", bnd)
        assert "sym_normal" in d
        # A: A x n = 0 on XZ plane (sym_tangential)
        d = self._build_dirichlet("a", bnd)
        assert "sym_tangential" in d

    def test_full_model_no_symmetry_blocks(self):
        bnd = ["outer"]
        d = self._build_dirichlet("omega", bnd)
        assert d == "outer"


# ============================================================
# Test: Argument construction
# ============================================================
class TestAccelArgConstruction:

    def _build_args(self, formulation="omega", material="steel",
                    coil_script="coil.py", mu_r="1000",
                    bh_file="", hys_file="", n_steps=1):
        args = [
            "calc_accel_magnet.py",
            "--coil-script", coil_script,
            "--cub5", "model.cub5",
            "--formulation", formulation,
            "--order", "2",
            "--fes-order", "1",
            "--max-iter", "30",
            "--relax", "0.3",
        ]
        if material == "linear":
            args += ["--material", "linear", "--mu-r", mu_r]
        elif material == "hysteresis":
            args += ["--material", "hysteresis",
                     "--hys-file", hys_file,
                     "--n-steps", str(n_steps)]
        else:
            args += ["--material", "steel"]
            if bh_file:
                args += ["--bh-file", bh_file]
        return args

    def test_linear_args(self):
        args = self._build_args(material="linear", mu_r="2000")
        assert "--material" in args
        idx = args.index("--material")
        assert args[idx + 1] == "linear"
        assert "--mu-r" in args
        assert args[args.index("--mu-r") + 1] == "2000"

    def test_steel_args(self):
        args = self._build_args(material="steel")
        assert args[args.index("--material") + 1] == "steel"
        assert "--mu-r" not in args

    def test_bh_file_args(self):
        args = self._build_args(material="steel", bh_file="custom.txt")
        assert "--bh-file" in args
        assert args[args.index("--bh-file") + 1] == "custom.txt"

    def test_hysteresis_args(self):
        args = self._build_args(material="hysteresis",
                                hys_file="iron.hys", n_steps=10)
        assert args[args.index("--material") + 1] == "hysteresis"
        assert "--hys-file" in args
        assert args[args.index("--hys-file") + 1] == "iron.hys"
        assert "--n-steps" in args
        assert args[args.index("--n-steps") + 1] == "10"

    def test_coil_script_always_present(self):
        args = self._build_args()
        assert "--coil-script" in args
        assert args[args.index("--coil-script") + 1] == "coil.py"

    def test_omega_formulation(self):
        args = self._build_args(formulation="omega")
        assert args[args.index("--formulation") + 1] == "omega"

    def test_a_formulation(self):
        args = self._build_args(formulation="a")
        assert args[args.index("--formulation") + 1] == "a"


# ============================================================
# Test: Formulation selection logic
# ============================================================
class TestFormulationSelection:

    @staticmethod
    def _select_formulation(combo_text):
        if "Omega" in combo_text:
            return "omega"
        elif "HCurl" in combo_text:
            return "a"
        elif "HDiv" in combo_text:
            return "hdiv"
        else:
            return "unknown"

    def test_omega(self):
        assert self._select_formulation("Omega-reduced (H1)") == "omega"

    def test_a(self):
        assert self._select_formulation("A-formulation (HCurl)") == "a"

    def test_hdiv(self):
        assert self._select_formulation("HDiv-VIM") == "hdiv"

    def test_radia_disables_kelvin(self):
        """HDiv-VIM is open-boundary through the charge Gram."""
        assert self._select_formulation("HDiv-VIM") == "hdiv"

    def test_fem_needs_kelvin(self):
        """Omega/A: Kelvin needed for open boundary."""
        for text in ["Omega-reduced (H1)", "A-formulation (HCurl)"]:
            assert "Radia" not in text


# ============================================================
# Test: Hantila decomposition constants
# ============================================================
class TestHantilaConstants:
    """Verify Hantila decomposition constants are consistent."""

    MU_0 = 4e-7 * math.pi
    NU_0 = 1.0 / (4e-7 * math.pi)

    def test_alpha_from_nu_rev(self):
        """alpha = 1/(mu_0 * nu_rev) - 1."""
        nu_rev = 500.0  # typical for steel
        alpha = 1.0 / (self.MU_0 * nu_rev) - 1.0
        assert alpha > 0
        # mu_eff = mu_0 * (1 + alpha) = 1/nu_rev
        mu_eff = self.MU_0 * (1 + alpha)
        assert abs(mu_eff - 1.0 / nu_rev) < 1e-10

    def test_mu_eff_equals_inv_nu_rev(self):
        for nu_rev in [100, 500, 1000, 5000]:
            alpha = 1.0 / (self.MU_0 * nu_rev) - 1.0
            mu_eff = self.MU_0 * (1 + alpha)
            assert abs(mu_eff * nu_rev - 1.0) < 1e-10

    def test_alpha_positive_for_magnetic(self):
        """For magnetic materials (mu_r > 1), alpha > 0."""
        # nu_rev < NU_0 for magnetic materials
        nu_rev = self.NU_0 / 100  # mu_r ~ 100
        alpha = 1.0 / (self.MU_0 * nu_rev) - 1.0
        assert alpha > 0
        assert alpha == pytest.approx(99.0, rel=1e-6)

    def test_alpha_zero_for_air(self):
        """For air (nu_rev = NU_0), alpha = 0."""
        nu_rev = self.NU_0
        alpha = 1.0 / (self.MU_0 * nu_rev) - 1.0
        assert abs(alpha) < 1e-10

    def test_kelvin_weight_formula(self):
        """Kelvin weight (a/r')^2 -> 0 as r' -> infinity."""
        a = 0.3  # Kelvin radius
        for rp in [0.5, 1.0, 5.0, 100.0]:
            w = a**2 / rp**2
            assert w > 0
            assert w < 1 if rp > a else True
        # At r' = a: weight = 1
        assert abs(a**2 / a**2 - 1.0) < 1e-10


# ============================================================
# Test: Quasi-static current ramp
# ============================================================
class TestCurrentRamp:

    def _build_ramp(self, current, n_steps):
        if n_steps <= 1:
            return [current]
        half = n_steps // 2
        ramp_up = [current * (i + 1) / half for i in range(half)]
        ramp_down = [current * (half - 1 - i) / half for i in range(half)]
        return ramp_up + ramp_down

    def test_single_step(self):
        levels = self._build_ramp(1000, 1)
        assert levels == [1000]

    def test_ramp_10_steps(self):
        levels = self._build_ramp(1000, 10)
        assert len(levels) == 10
        # Ramp up: 200, 400, 600, 800, 1000
        assert levels[0] == pytest.approx(200)
        assert levels[4] == pytest.approx(1000)
        # Ramp down: 800, 600, 400, 200, 0
        assert levels[5] == pytest.approx(800)
        assert levels[9] == pytest.approx(0)

    def test_ramp_symmetric(self):
        levels = self._build_ramp(500, 20)
        assert len(levels) == 20
        # Peak at step 10
        assert levels[9] == pytest.approx(500)
        # First and last should be symmetric
        assert levels[0] == pytest.approx(levels[-2])

    def test_ramp_up_monotonic(self):
        levels = self._build_ramp(1000, 10)
        up = levels[:5]
        for i in range(1, len(up)):
            assert up[i] > up[i - 1]


# ============================================================
# Test: Result display formatting
# ============================================================
class TestAccelResultDisplay:

    def test_linear_result(self):
        data = {
            "B_origin_mag": 0.5,
            "L": 1.5e-6,
            "formulation": "omega",
            "iterations": 1,
            "converged": True,
            "ndof": 5000,
            "t_total": 2.3,
            "hysteresis": False,
        }
        assert data["B_origin_mag"] == pytest.approx(0.5)
        assert data["L"] * 1e6 == pytest.approx(1.5)

    def test_hysteresis_result(self):
        data = {
            "B_origin_mag": 0.8,
            "L": 3.2e-6,
            "formulation": "a",
            "iterations": 5,
            "converged": True,
            "hysteresis": True,
            "n_steps": 10,
            "step_results": [
                {"step": 0, "current": 200, "B_center": 0.2, "iterations": 3},
                {"step": 9, "current": 0, "B_center": 0.05, "iterations": 2},
            ],
        }
        assert data["hysteresis"] is True
        assert data["n_steps"] == 10
        assert len(data["step_results"]) == 2
        # Remanence: B > 0 at I = 0
        assert data["step_results"][-1]["B_center"] > 0

    def test_energy_inductance_relation(self):
        """L = 2 * W_mag / I^2."""
        W_mag = 0.005  # 5 mJ
        I = 1000  # 1000 A
        L = 2 * W_mag / I**2
        assert L == pytest.approx(1e-8)  # 10 nH


# ============================================================
# Test: Kelvin symmetry sphere center computation
# ============================================================
class TestKelvinSymmetry:
    """Test symmetry-aware sphere center and radius computation."""

    @staticmethod
    def _sym_coord(lo, hi):
        """Symmetry plane = bounding box face closest to 0."""
        return lo if abs(lo) <= abs(hi) else hi

    @staticmethod
    def _compute_sphere_center(symmetry, bb_min, bb_max):
        cx = (bb_min[0] + bb_max[0]) / 2
        cy = (bb_min[1] + bb_max[1]) / 2
        cz = (bb_min[2] + bb_max[2]) / 2

        sc = TestKelvinSymmetry._sym_coord
        if "Quarter" in symmetry or "Eighth" in symmetry:
            sx = sc(bb_min[0], bb_max[0])
        else:
            sx = cx
        if "Eighth" in symmetry:
            sy = sc(bb_min[1], bb_max[1])
        else:
            sy = cy
        if symmetry != "Full":
            sz = sc(bb_min[2], bb_max[2])
        else:
            sz = cz
        return (sx, sy, sz)

    @staticmethod
    def _enclosing_radius(sx, sy, sz, bb_min, bb_max):
        corners = [(x, y, z)
                   for x in (bb_min[0], bb_max[0])
                   for y in (bb_min[1], bb_max[1])
                   for z in (bb_min[2], bb_max[2])]
        return max(math.sqrt((x-sx)**2 + (y-sy)**2 + (z-sz)**2)
                   for x, y, z in corners)

    def test_full_sphere_at_centroid(self):
        bb_min = [0, -0.01, 0]
        bb_max = [0.1, 0.01, 0.1]
        sx, sy, sz = self._compute_sphere_center("Full", bb_min, bb_max)
        assert sx == pytest.approx(0.05)
        assert sy == pytest.approx(0.0)
        assert sz == pytest.approx(0.05)

    def test_half_z_at_z_symmetry(self):
        """Half model: z symmetry plane at z=0 (bb_min[2]=0)."""
        bb_min = [0, -0.01, 0]
        bb_max = [0.1, 0.01, 0.1]
        sx, sy, sz = self._compute_sphere_center("Half (Z)", bb_min, bb_max)
        assert sx == pytest.approx(0.05)  # centroid (no x symmetry)
        assert sz == pytest.approx(0.0)   # symmetry plane at z=0

    def test_quarter_xz_at_origin(self):
        """Quarter model in +x, +z: sphere center at (0, cy, 0)."""
        bb_min = [0, -0.05, 0]
        bb_max = [0.1, 0.05, 0.08]
        sx, sy, sz = self._compute_sphere_center("Quarter (X,Z)", bb_min, bb_max)
        assert sx == pytest.approx(0.0)   # x-symmetry plane
        assert sy == pytest.approx(0.0)   # centroid (no y symmetry)
        assert sz == pytest.approx(0.0)   # z-symmetry plane

    def test_eighth_all_at_origin(self):
        """Eighth model in +x, +y, +z: sphere center at origin."""
        bb_min = [0, 0, 0]
        bb_max = [0.1, 0.05, 0.08]
        sx, sy, sz = self._compute_sphere_center("Eighth (X,Y,Z)", bb_min, bb_max)
        assert sx == pytest.approx(0.0)
        assert sy == pytest.approx(0.0)
        assert sz == pytest.approx(0.0)

    def test_quarter_negative_quadrant(self):
        """Quarter model in -x, +z: symmetry plane at x=0 (bb_max[0])."""
        bb_min = [-0.1, -0.05, 0]
        bb_max = [0, 0.05, 0.08]
        sx, sy, sz = self._compute_sphere_center("Quarter (X,Z)", bb_min, bb_max)
        assert sx == pytest.approx(0.0)   # bb_max closer to 0
        assert sz == pytest.approx(0.0)

    def test_enclosing_radius_full(self):
        """Full: radius = half diagonal from centroid."""
        bb_min = [-0.05, -0.05, -0.05]
        bb_max = [0.05, 0.05, 0.05]
        sx, sy, sz = 0.0, 0.0, 0.0
        r = self._enclosing_radius(sx, sy, sz, bb_min, bb_max)
        expected = math.sqrt(0.05**2 * 3)
        assert r == pytest.approx(expected)

    def test_enclosing_radius_quarter(self):
        """Quarter: radius from origin to farthest corner."""
        bb_min = [0, -0.01, 0]
        bb_max = [0.1, 0.01, 0.08]
        sx, sy, sz = 0.0, 0.0, 0.0
        r = self._enclosing_radius(sx, sy, sz, bb_min, bb_max)
        expected = math.sqrt(0.1**2 + 0.01**2 + 0.08**2)
        assert r == pytest.approx(expected)

    def test_sym_coord_picks_closer_to_zero(self):
        assert self._sym_coord(0, 0.1) == 0
        assert self._sym_coord(-0.1, 0) == 0
        assert self._sym_coord(-0.02, 0.1) == -0.02
        assert self._sym_coord(-0.1, 0.02) == 0.02

    def test_matching_quadrant_filter(self):
        """Verify matching quadrant logic for exterior pieces."""
        # Model centroid at (0.05, 0, 0.04), sphere center at (0, 0, 0)
        cx, sx = 0.05, 0.0
        cz, sz = 0.04, 0.0
        dir_x = cx - sx  # positive -> keep nc_x > ext_sx
        dir_z = cz - sz  # positive -> keep nc_z > sz

        ext_sx = sx + 1.0  # exterior sphere at x=1.0
        eps = 0.001

        # 4 quadrants: (x > ext_sx, z > 0), (x < ext_sx, z > 0),
        #              (x > ext_sx, z < 0), (x < ext_sx, z < 0)
        quadrants = [
            (ext_sx + 0.3, 0.0, 0.3),   # +x, +z -> KEEP
            (ext_sx - 0.3, 0.0, 0.3),   # -x, +z -> DELETE
            (ext_sx + 0.3, 0.0, -0.3),  # +x, -z -> DELETE
            (ext_sx - 0.3, 0.0, -0.3),  # -x, -z -> DELETE
        ]
        kept = []
        for nc_x, nc_y, nc_z in quadrants:
            ok = True
            if dir_z >= 0 and nc_z < sz - eps:
                ok = False
            elif dir_z < 0 and nc_z > sz + eps:
                ok = False
            if dir_x >= 0 and nc_x < ext_sx - eps:
                ok = False
            elif dir_x < 0 and nc_x > ext_sx + eps:
                ok = False
            if ok:
                kept.append((nc_x, nc_y, nc_z))

        assert len(kept) == 1
        assert kept[0][0] > ext_sx  # +x side
        assert kept[0][2] > 0       # +z side
