"""Integration tests for IH (BEM) and IH (FEM) solver pipelines.

Tests the calc scripts directly (no Cubit, using OCC fallback meshes).
Requires: NGSolve + ngsolve.bem.

Run: python -m pytest tests/panels/test_ih_solvers.py -v
"""

import json
import math
import os
import subprocess
import sys
import tempfile

import pytest

ngsolve = pytest.importorskip("ngsolve")

_repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_panels_dir = os.path.join(_repo, "src", "radia", "panels")
_radia_src = os.path.join(_repo, "src", "radia")

# Add to path for direct imports
for p in [_panels_dir, _radia_src]:
    if p not in sys.path:
        sys.path.insert(0, p)

try:
    from ngsolve.bem import BiotSavartCF  # noqa: F401
    _HAS_BEM = True
except ImportError:
    _HAS_BEM = False


def _run_script(script_name, args, timeout=120):
    """Run a calc script as subprocess and return parsed JSON."""
    script = os.path.join(_panels_dir, script_name)
    cmd = [sys.executable, script] + args

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)

    # Parse JSON from stdout (last JSON line)
    result = None
    for line in proc.stdout.strip().splitlines():
        line = line.strip()
        if line.startswith("{"):
            try:
                result = json.loads(line)
            except json.JSONDecodeError:
                continue

    if result is None and proc.returncode != 0:
        pytest.fail(f"{script_name} failed (exit {proc.returncode}):\n"
                    f"STDERR: {proc.stderr[-500:]}")
    return result, proc.stderr


# ============================================================
# IH (BEM): calc_inductance.py with OCC fallback
# ============================================================
class TestIHBEM:
    """Test BEM inductance + workpiece impedance pipeline."""

    @pytest.fixture(scope="class")
    def occ_surface_mesh(self):
        """Create a torus surface mesh via OCC (no Cubit).

        Returns path to .vol file + metadata.
        """
        from ngsolve import Mesh, TaskManager
        from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Axis,
                                 OCCGeometry, Glue)

        R = 0.03
        a = 0.003
        gap_deg = 10

        wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
        circle = wp.Circle(a).Face()
        torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)),
                               360 - gap_deg)

        # Label faces
        expected_area = math.pi * a**2
        for f in torus.faces:
            area = f.mass
            if abs(area - expected_area) / expected_area < 0.3:
                angle = math.atan2(f.center.y, f.center.x)
                if abs(angle) < math.radians(gap_deg):
                    f.name = "source"
                else:
                    f.name = "sink"
            else:
                f.name = "default"
            f.maxh = a * 2

        geo = OCCGeometry(Glue(torus.faces))
        with TaskManager():
            ngmesh = geo.GenerateMesh(maxh=a * 3)
        mesh = Mesh(ngmesh)
        mesh.Curve(2)

        tmpdir = tempfile.mkdtemp(prefix="radia_test_bem_")
        vol_path = os.path.join(tmpdir, "coil.vol")
        mesh.ngmesh.Save(vol_path)

        return {
            "vol_path": vol_path,
            "tmpdir": tmpdir,
            "R": R,
            "a": a,
            "n_elements": mesh.GetNE(ngsolve.BND),
        }

    @pytest.mark.skipif(not _HAS_BEM, reason="ngsolve.bem not available")
    def test_inductance_direct_call(self, occ_surface_mesh):
        """Test BEM solver called directly (not as subprocess)."""
        from bem_inductance import compute_inductance_source_sink
        from netgen.meshing import Mesh as NetgenMesh

        vol_path = occ_surface_mesh["vol_path"]
        ngmesh = NetgenMesh()
        ngmesh.Load(vol_path)
        mesh = ngsolve.Mesh(ngmesh)

        boundaries = list(set(mesh.GetBoundaries()))
        if "source" not in boundaries or "sink" not in boundaries:
            pytest.skip("OCC mesh missing source/sink labels")

        sol = compute_inductance_source_sink(mesh, "source", "sink", 0)

        assert "error" not in sol, f"Solver error: {sol.get('error')}"
        L = sol["L"]
        assert L > 0, f"Inductance must be positive, got {L}"
        # Torus analytical: L ~ mu0*R*(ln(8R/a)-2)
        mu0 = 4e-7 * math.pi
        L_analytical = mu0 * occ_surface_mesh["R"] * (
            math.log(8 * occ_surface_mesh["R"] / occ_surface_mesh["a"]) - 2)
        # Allow large tolerance (BEM on coarse mesh)
        assert abs(L - L_analytical) / L_analytical < 0.5, \
            f"L={L*1e9:.1f} nH vs analytical {L_analytical*1e9:.1f} nH"


# ============================================================
# IH (FEM): calc_fem_kelvin.py with OCC fallback
# ============================================================
class TestIHFEM:
    """Test FEM-SIBC solver pipeline."""

    @pytest.fixture(scope="class")
    def occ_fem_mesh(self):
        """Create 3D IH mesh via OCC (no Cubit).

        Uses build_occ_ih_mesh_3d from calc_common.
        """
        from calc_common import build_occ_ih_mesh_3d

        mesh, info = build_occ_ih_mesh_3d(
            R_coil=0.030, a_coil=0.003, gap_deg=10,
            R_wp=0.010, H_wp=0.020,
            a_kelvin=0.060, kelvin_offset=0.300,
            maxh_coil=0.008, maxh_air=0.020,
            order=1)

        tmpdir = tempfile.mkdtemp(prefix="radia_test_fem_")
        vol_path = os.path.join(tmpdir, "ih_model.vol")
        mesh.ngmesh.Save(vol_path)

        return {
            "vol_path": vol_path,
            "tmpdir": tmpdir,
            "mesh": mesh,
            "info": info,
        }

    def test_occ_mesh_has_materials(self, occ_fem_mesh):
        """OCC mesh should have coil + air (hole approach: no workpiece)."""
        mesh = occ_fem_mesh["mesh"]
        mats = set(mesh.GetMaterials())
        assert "coil" in mats, f"Missing 'coil', got {mats}"
        assert "air" in mats or any("default" in m for m in mats), \
            f"Missing 'air', got {mats}"
        # Hole approach: workpiece NOT in materials
        assert "workpiece" not in mats, \
            f"Hole approach: 'workpiece' should NOT be in materials, got {mats}"

    def test_occ_mesh_has_boundaries(self, occ_fem_mesh):
        """OCC mesh should have source + sink + sibc boundaries."""
        mesh = occ_fem_mesh["mesh"]
        bnds = set(mesh.GetBoundaries())
        assert "source" in bnds, f"Missing 'source', got {bnds}"
        assert "sink" in bnds, f"Missing 'sink', got {bnds}"
        assert "sibc" in bnds, f"Missing 'sibc' (hole approach), got {bnds}"

    def test_t0_source_computation(self, occ_fem_mesh):
        """T0 source/sink technique should produce nonzero J."""
        from calc_common import compute_T0_source
        import numpy as np

        mesh = occ_fem_mesh["mesh"]
        bnds = set(mesh.GetBoundaries())
        if "source" not in bnds or "sink" not in bnds:
            pytest.skip("OCC mesh missing source/sink")

        J_source, gf_T0 = compute_T0_source(mesh, order=1, I_total=1.0)
        # T0 should range from 0 to 1
        vals = gf_T0.vec.FV().NumPy()
        assert vals.max() > 0.5, f"T0 max should be ~1.0, got {vals.max()}"

    def test_fem_solve_kelvin_hole(self, occ_fem_mesh):
        """End-to-end FEM-SIBC solve with hole approach + Kelvin.

        Validates: L within 10% of analytical torus formula.
        Analytical: L = mu0 * R * (ln(8R/a) - 2) * (1 - gap/360)
        """
        from calc_fem_kelvin import solve_fem
        from em_material import EMMaterial

        vol_path = occ_fem_mesh["vol_path"]
        mat = EMMaterial.from_name("copper")

        result = solve_fem(
            vol_file=vol_path,
            fes_order=1,
            frequency=50000,
            mat=mat,
            impedance_model="sibc",
            formulation="total",
            I_total=1.0,
            half_thickness=0.010,
            solver="pardiso",
            reg=1e-6,
            nthreads=0)

        assert "error" not in result, f"Solver error: {result.get('error')}"
        L = result["L"]
        assert L > 0, f"Inductance must be positive, got {L}"

        # Analytical: mu0 * R * (ln(8R/a) - 2) * (1 - gap/360)
        mu0 = 4e-7 * math.pi
        R, a = 0.030, 0.003
        gap_deg = 10
        L_analytical = mu0 * R * (math.log(8 * R / a) - 2) * (1 - gap_deg / 360)
        # Coarse mesh: allow 10% tolerance
        err = abs(L - L_analytical) / L_analytical
        assert err < 0.10, \
            f"L={L*1e9:.1f} nH vs analytical {L_analytical*1e9:.1f} nH ({err*100:.1f}%)"

        # Hole approach: H_t should be available
        assert result.get("H_t_rms") is not None, "H_t_rms should be available"
        assert result["H_t_rms"] > 0, f"H_t should be positive, got {result['H_t_rms']}"

        # Kelvin should be detected
        assert result.get("has_kelvin") is True, "Kelvin not detected"


class TestCalcCommonProtocol:
    """Test calc_common protocol used by both IH solvers."""

    def test_progress_format(self, capsys):
        """progress() writes TAG:msg to stderr."""
        from calc_common import progress
        progress("BEM-SIBC", "wp R=10.0mm, sigma=2e6")
        captured = capsys.readouterr()
        assert "BEM-SIBC:wp R=10.0mm, sigma=2e6" in captured.err

    def test_sigma_for_material(self):
        """sigma_for_material returns reasonable values."""
        from calc_common import sigma_for_material
        assert sigma_for_material("steel") == 2e6
        assert sigma_for_material("copper") == 5.8e7
        assert sigma_for_material("aluminum") == 3.5e7
        assert sigma_for_material("unknown") == 2e6  # default

    def test_get_bh_curve_steel(self):
        """Built-in steel BH curve should be monotonic."""
        from calc_common import get_bh_curve
        bh, mu_r = get_bh_curve("steel")
        assert bh is not None
        assert mu_r is None
        # Check monotonicity
        for i in range(1, len(bh)):
            assert bh[i][0] >= bh[i-1][0], "H not monotonic"
            assert bh[i][1] >= bh[i-1][1], "B not monotonic"

    def test_get_bh_curve_copper(self):
        """Copper has no BH curve, just mu_r=1."""
        from calc_common import get_bh_curve
        bh, mu_r = get_bh_curve("copper")
        assert bh is None
        assert mu_r == 1.0


class TestEMMaterial:
    """Test EMMaterial class and its integration."""

    def test_from_name_steel(self):
        from em_material import EMMaterial
        mat = EMMaterial.from_name("steel")
        assert mat.sigma == 2e6
        assert mat.mu_r == 100.0
        assert mat.bh_curve is not None
        assert mat.is_magnetic
        assert mat.is_nonlinear
        assert mat.is_conducting

    def test_from_name_copper(self):
        from em_material import EMMaterial
        mat = EMMaterial.from_name("copper")
        assert mat.sigma == 5.8e7
        assert mat.mu_r == 1.0
        assert mat.bh_curve is None
        assert not mat.is_magnetic
        assert not mat.is_nonlinear

    def test_from_name_aluminum(self):
        from em_material import EMMaterial
        mat = EMMaterial.from_name("aluminum")
        assert mat.sigma == 3.5e7
        assert mat.mu_r == 1.0

    def test_from_name_unknown(self):
        from em_material import EMMaterial
        with pytest.raises(ValueError, match="Unknown material"):
            EMMaterial.from_name("unobtanium")

    def test_custom(self):
        from em_material import EMMaterial
        mat = EMMaterial.custom(sigma=1e7, mu_r=50)
        assert mat.name == "custom"
        assert mat.sigma == 1e7
        assert mat.mu_r == 50

    def test_skin_depth(self):
        from em_material import EMMaterial
        mat = EMMaterial.from_name("copper")
        delta = mat.skin_depth(50000)
        # Copper at 50 kHz: delta ~ 0.3 mm
        assert 0.0002 < delta < 0.0005

    def test_skin_depth_dc(self):
        from em_material import EMMaterial
        mat = EMMaterial.from_name("copper")
        assert mat.skin_depth(0) == float("inf")

    def test_dowell_Zs(self):
        from em_material import EMMaterial
        mat = EMMaterial.from_name("copper")
        Z = mat.dowell_Zs(50000, 0.01)
        assert Z.real > 0
        assert Z.imag > 0

    def test_dowell_Zs_dc(self):
        from em_material import EMMaterial
        mat = EMMaterial.from_name("copper")
        Z = mat.dowell_Zs(0, 0.01)
        assert Z == 0

    def test_rho(self):
        from em_material import EMMaterial
        mat = EMMaterial.from_name("copper")
        assert abs(mat.rho - 1.0 / 5.8e7) < 1e-15

    def test_get_bh_curve_method(self):
        from em_material import EMMaterial
        mat = EMMaterial.from_name("steel")
        bh, mu = mat.get_bh_curve()
        assert bh is not None
        assert mu is None

        mat2 = EMMaterial.from_name("copper")
        bh2, mu2 = mat2.get_bh_curve()
        assert bh2 is None
        assert mu2 == 1.0

    def test_from_args_preset(self):
        """--material copper should auto-set sigma=5.8e7."""
        import argparse
        from em_material import EMMaterial, add_material_args
        parser = argparse.ArgumentParser()
        add_material_args(parser)
        args = parser.parse_args(["--material", "copper"])
        mat = EMMaterial.from_args(args)
        assert mat.sigma == 5.8e7
        assert mat.mu_r == 1.0
        assert mat.name == "copper"

    def test_from_args_preset_override_sigma(self):
        """--material steel --sigma 5e6 should override sigma only."""
        import argparse
        from em_material import EMMaterial, add_material_args
        parser = argparse.ArgumentParser()
        add_material_args(parser)
        args = parser.parse_args(["--material", "steel", "--sigma", "5e6"])
        mat = EMMaterial.from_args(args)
        assert mat.sigma == 5e6
        assert mat.bh_curve is not None  # steel BH still present

    def test_from_args_custom(self):
        """--material custom --sigma 1e7 --mu-r 50."""
        import argparse
        from em_material import EMMaterial, add_material_args
        parser = argparse.ArgumentParser()
        add_material_args(parser, include_custom=True)
        args = parser.parse_args(["--material", "custom",
                                  "--sigma", "1e7", "--mu-r", "50"])
        mat = EMMaterial.from_args(args)
        assert mat.sigma == 1e7
        assert mat.mu_r == 50

    def test_from_args_default_steel(self):
        """Default (no args) should be steel."""
        import argparse
        from em_material import EMMaterial, add_material_args
        parser = argparse.ArgumentParser()
        add_material_args(parser)
        args = parser.parse_args([])
        mat = EMMaterial.from_args(args)
        assert mat.name == "steel"
        assert mat.sigma == 2e6

    def test_linear_mode(self):
        """--material linear for accelerator yoke."""
        from em_material import EMMaterial
        mat = EMMaterial.from_name("linear")
        assert mat.mu_r == 1000.0
        assert mat.bh_curve is None
        assert not mat.is_nonlinear
