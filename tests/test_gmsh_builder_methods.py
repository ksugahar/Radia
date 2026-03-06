"""
Tests for GmshBuilder methods.

Requires: pip install gmsh
"""

import os
import sys
import tempfile
import pytest

# Skip entire module if gmsh is not installed
gmsh = pytest.importorskip('gmsh')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'radia'))
from radia.gmsh_builder import GmshBuilder


class TestContextManager:
    """Test context manager lifecycle."""

    def test_enter_exit(self):
        with GmshBuilder(verbose=False) as hb:
            assert hb._initialized is True
        assert hb._initialized is False

    def test_not_initialized_raises(self):
        hb = GmshBuilder(verbose=False)
        with pytest.raises(RuntimeError, match="context manager"):
            hb.add_box([0, 0, 0], [0.1, 0.1, 0.1])


class TestAddBox:
    """Test box creation."""

    def test_single_box(self):
        with GmshBuilder(verbose=False) as hb:
            vid = hb.add_box([0, 0, 0], [0.1, 0.02, 0.03])
            assert vid == 1
            assert vid in hb.get_volumes()
            assert len(hb.get_gmsh_tags(vid)) == 1

    def test_two_boxes(self):
        with GmshBuilder(verbose=False) as hb:
            v1 = hb.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = hb.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            assert v1 != v2
            assert len(hb.get_volumes()) == 2


class TestAddCylinder:
    """Test cylinder creation."""

    def test_cylinder_z(self):
        with GmshBuilder(verbose=False) as hb:
            vid = hb.add_cylinder([0, 0, 0], 'z', 0.01, 0.05)
            assert vid in hb.get_volumes()

    def test_cylinder_invalid_axis(self):
        with GmshBuilder(verbose=False) as hb:
            with pytest.raises(ValueError, match="axis"):
                hb.add_cylinder([0, 0, 0], 'w', 0.01, 0.05)


class TestAddSphere:
    """Test sphere creation."""

    def test_sphere(self):
        with GmshBuilder(verbose=False) as hb:
            vid = hb.add_sphere([0, 0, 0], 0.01)
            assert vid in hb.get_volumes()


class TestBooleanOps:
    """Test boolean operations."""

    def test_fuse(self):
        with GmshBuilder(verbose=False) as hb:
            v1 = hb.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = hb.add_box([0.05, 0, 0], [0.1, 0.1, 0.1])
            v_fused = hb.fuse([v1, v2])
            assert v_fused in hb.get_volumes()
            assert v1 not in hb.get_volumes()
            assert v2 not in hb.get_volumes()

    def test_fuse_requires_two(self):
        with GmshBuilder(verbose=False) as hb:
            v1 = hb.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            with pytest.raises(ValueError, match="at least 2"):
                hb.fuse([v1])

    def test_cut(self):
        with GmshBuilder(verbose=False) as hb:
            v1 = hb.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = hb.add_box([0.03, 0, 0], [0.04, 0.04, 0.04])
            v_result = hb.cut(v1, v2)
            assert v_result in hb.get_volumes()
            assert v1 not in hb.get_volumes()
            assert v2 not in hb.get_volumes()


class TestWebcut:
    """Test webcut plane operation."""

    def test_webcut_x(self):
        with GmshBuilder(verbose=False) as hb:
            v = hb.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            new_ids = hb.webcut_plane(v, 'x', 0.0)
            assert len(new_ids) == 2
            assert v not in hb.get_volumes()
            for nid in new_ids:
                assert nid in hb.get_volumes()

    def test_webcut_y(self):
        with GmshBuilder(verbose=False) as hb:
            v = hb.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            new_ids = hb.webcut_plane(v, 'y', 0.0)
            assert len(new_ids) == 2

    def test_webcut_z(self):
        with GmshBuilder(verbose=False) as hb:
            v = hb.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            new_ids = hb.webcut_plane(v, 'z', 0.0)
            assert len(new_ids) == 2

    def test_double_webcut(self):
        """Cut a box into 4 pieces with two perpendicular cuts."""
        with GmshBuilder(verbose=False) as hb:
            v = hb.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            parts_x = hb.webcut_plane(v, 'x', 0.0)
            assert len(parts_x) == 2

            # Cut first piece along y
            parts_y = hb.webcut_plane(parts_x[0], 'y', 0.0)
            assert len(parts_y) == 2

            # Total: 3 volumes (1 uncut + 2 from second cut)
            assert len(hb.get_volumes()) == 3

    def test_webcut_invalid_axis(self):
        with GmshBuilder(verbose=False) as hb:
            v = hb.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            with pytest.raises(ValueError, match="axis"):
                hb.webcut_plane(v, 'w', 0.0)

    def test_webcut_nonexistent_volume(self):
        with GmshBuilder(verbose=False) as hb:
            with pytest.raises(KeyError):
                hb.webcut_plane(999, 'x', 0.0)


class TestMeshGeneration:
    """Test mesh generation."""

    def test_generate_hex_single_box(self):
        with GmshBuilder(verbose=False) as hb:
            hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            hb.set_divisions_all(3, 3, 3)
            hb.generate()
            info = hb.get_mesh_info()
            assert info['meshed'] is True
            assert info['n_hex'] == 27  # 3*3*3
            assert info['n_tet'] == 0

    def test_generate_hex_with_mesh_size(self):
        with GmshBuilder(verbose=False) as hb:
            hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            hb.set_mesh_size(0.01)
            hb.generate()
            info = hb.get_mesh_info()
            assert info['meshed'] is True
            assert info['n_hex'] == 27  # 0.03/0.01 = 3 per axis

    def test_generate_hex_per_volume_divisions(self):
        with GmshBuilder(verbose=False) as hb:
            v = hb.add_box([0, 0, 0], [0.04, 0.02, 0.01])
            hb.set_divisions(v, 4, 2, 1)
            hb.generate()
            info = hb.get_mesh_info()
            assert info['n_hex'] == 8  # 4*2*1

    def test_generate_tet(self):
        with GmshBuilder(verbose=False) as hb:
            hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            hb.set_mesh_size(0.01)
            hb.generate(element_type='tet')
            info = hb.get_mesh_info()
            assert info['meshed'] is True
            assert info['n_tet'] > 0

    def test_generate_hex_after_webcut(self):
        """Hex mesh on webcut volumes."""
        with GmshBuilder(verbose=False) as hb:
            v = hb.add_box([0, 0, 0], [0.06, 0.03, 0.03])
            new_ids = hb.webcut_plane(v, 'x', 0.0)

            hb.set_divisions_all(3, 3, 3)
            hb.generate()

            info = hb.get_mesh_info()
            assert info['meshed'] is True
            assert info['n_hex'] == 54  # 2 * (3*3*3) = 54


class TestExport:
    """Test mesh export."""

    def test_export_msh(self):
        with GmshBuilder(verbose=False) as hb:
            hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            hb.set_divisions_all(2, 2, 2)

            tmp = os.path.join(tempfile.gettempdir(), '_test_hex_export.msh')
            try:
                hb.export(tmp)
                assert os.path.exists(tmp)
                # Check file has content
                size = os.path.getsize(tmp)
                assert size > 100
            finally:
                if os.path.exists(tmp):
                    os.unlink(tmp)


class TestToRadia:
    """Test Radia integration."""

    def test_to_radia(self):
        try:
            import radia as rad
        except ImportError:
            pytest.skip("radia not installed")

        rad.UtiDelAll()

        with GmshBuilder(verbose=False) as hb:
            hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            hb.set_divisions_all(2, 2, 2)
            container = hb.to_radia(mu_r=1000)

        assert container > 0


class TestGetMeshInfo:
    """Test mesh info query."""

    def test_info_before_mesh(self):
        with GmshBuilder(verbose=False) as hb:
            hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            info = hb.get_mesh_info()
            assert info['meshed'] is False

    def test_info_after_mesh(self):
        with GmshBuilder(verbose=False) as hb:
            hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            hb.set_divisions_all(2, 2, 2)
            hb.generate()
            info = hb.get_mesh_info()
            assert info['meshed'] is True
            assert info['n_nodes'] > 0
            assert info['n_hex'] == 8
            assert info['n_volumes'] == 1


class TestSurfaceMesh:
    """Test surface mesh generation."""

    def test_generate_surface(self):
        with GmshBuilder(verbose=False) as hb:
            hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            hb.generate_surface(mesh_size=0.01)
            assert hb._meshed is True

    def test_surface_cylinder(self):
        with GmshBuilder(verbose=False) as hb:
            hb.add_cylinder([0, 0, 0], 'z', 0.01, 0.05)
            hb.generate_surface(mesh_size=0.005)
            assert hb._meshed is True


class TestPhysicalGroups:
    """Test physical group assignment."""

    def test_set_physical_group(self):
        with GmshBuilder(verbose=False) as hb:
            v1 = hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            v2 = hb.add_box([0.05, 0, 0], [0.03, 0.03, 0.03])
            hb.set_physical_group(v1, 'iron_core')
            hb.set_physical_group(v2, 'coil')
            # Verify no error raised

    def test_physical_group_nonexistent(self):
        with GmshBuilder(verbose=False) as hb:
            with pytest.raises(KeyError):
                hb.set_physical_group(999, 'test')


class TestTransformations:
    """Test translate, rotate, copy."""

    def test_translate(self):
        with GmshBuilder(verbose=False) as hb:
            v = hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            hb.translate(v, 0.1, 0, 0)
            # Verify volume still exists
            assert v in hb.get_volumes()

    def test_copy(self):
        with GmshBuilder(verbose=False) as hb:
            v1 = hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            v2 = hb.copy(v1)
            assert v2 != v1
            assert len(hb.get_volumes()) == 2
            hb.translate(v2, 0.05, 0, 0)

    def test_rotate(self):
        import math
        with GmshBuilder(verbose=False) as hb:
            v = hb.add_box([0.05, 0, 0], [0.03, 0.03, 0.03])
            hb.rotate(v, [0, 0, 0], [0, 0, 1], math.pi / 4)
            assert v in hb.get_volumes()


class TestSetOrder:
    """Test mesh element order."""

    def test_set_order_2(self):
        with GmshBuilder(verbose=False) as hb:
            hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            hb.set_divisions_all(2, 2, 2)
            hb.generate()
            hb.set_order(2)
            # After set_order(2), nodes should increase (Hex20 vs Hex8)
            info = hb.get_mesh_info()
            assert info['n_hex'] == 8


class TestMeshSizeAt:
    """Test per-entity mesh size."""

    def test_set_mesh_size_at(self):
        with GmshBuilder(verbose=False) as hb:
            v = hb.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            hb.set_mesh_size_at(v, 0.005)
            # Should not raise
