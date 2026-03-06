"""
Tests for GmshBuilder - GMSH OCC wrapper.

Tests all modules: geometry, boolean, transforms, mesh control,
mesh generation, mesh quality, mesh access, physical groups, query, export.

Requires: pip install gmsh
"""

import os
import sys
import math
import tempfile
import pytest

gmsh = pytest.importorskip('gmsh')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'radia'))
from radia.gmsh_builder import GmshBuilder


# ================================================================
# Core
# ================================================================

class TestCore:
    """Test context manager and state management."""

    def test_enter_exit(self):
        with GmshBuilder(verbose=False) as cm:
            assert cm._initialized is True
        assert cm._initialized is False

    def test_not_initialized_raises(self):
        cm = GmshBuilder(verbose=False)
        with pytest.raises(RuntimeError, match="context manager"):
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])

    def test_reset(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            assert len(cm.get_volumes()) == 1
            cm.reset()
            assert len(cm.get_volumes()) == 0

    def test_compress_ids(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            v3 = cm.add_box([0.4, 0, 0], [0.1, 0.1, 0.1])
            cm.remove_volume(v2)
            id_map = cm.compress_ids()
            vols = cm.get_volumes()
            assert vols == [1, 2]
            assert len(id_map) == 2


# ================================================================
# Geometry
# ================================================================

class TestGeometry:
    """Test geometry creation methods."""

    def test_add_box(self):
        with GmshBuilder(verbose=False) as cm:
            vid = cm.add_box([0, 0, 0], [0.1, 0.02, 0.03])
            assert vid == 1
            assert vid in cm.get_volumes()
            assert len(cm.get_gmsh_tags(vid)) == 1

    def test_add_cylinder(self):
        with GmshBuilder(verbose=False) as cm:
            vid = cm.add_cylinder([0, 0, 0], 'z', 0.01, 0.05)
            assert vid in cm.get_volumes()

    def test_add_cylinder_invalid_axis(self):
        with GmshBuilder(verbose=False) as cm:
            with pytest.raises(ValueError, match="axis"):
                cm.add_cylinder([0, 0, 0], 'w', 0.01, 0.05)

    def test_add_sphere(self):
        with GmshBuilder(verbose=False) as cm:
            vid = cm.add_sphere([0, 0, 0], 0.01)
            assert vid in cm.get_volumes()

    def test_add_cone(self):
        with GmshBuilder(verbose=False) as cm:
            vid = cm.add_cone([0, 0, 0], 'z', 0.02, 0.01, 0.05)
            assert vid in cm.get_volumes()

    def test_add_cone_pointed(self):
        with GmshBuilder(verbose=False) as cm:
            vid = cm.add_cone([0, 0, 0], 'z', 0.02, 0.0, 0.05)
            assert vid in cm.get_volumes()

    def test_add_torus(self):
        with GmshBuilder(verbose=False) as cm:
            vid = cm.add_torus([0, 0, 0], 'z', 0.05, 0.01)
            assert vid in cm.get_volumes()

    def test_add_torus_x_axis(self):
        with GmshBuilder(verbose=False) as cm:
            vid = cm.add_torus([0, 0, 0], 'x', 0.05, 0.01)
            assert vid in cm.get_volumes()

    def test_add_wedge(self):
        with GmshBuilder(verbose=False) as cm:
            vid = cm.add_wedge([0, 0, 0], [0.1, 0.05, 0.03])
            assert vid in cm.get_volumes()

    def test_add_two_boxes(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            assert v1 != v2
            assert len(cm.get_volumes()) == 2

    def test_low_level_point_line(self):
        with GmshBuilder(verbose=False) as cm:
            p1 = cm.add_point(0, 0, 0)
            p2 = cm.add_point(0.1, 0, 0)
            line = cm.add_line(p1, p2)
            assert p1 > 0
            assert p2 > 0
            assert line > 0

    def test_add_rectangle(self):
        with GmshBuilder(verbose=False) as cm:
            surf = cm.add_rectangle([0, 0, 0], 0.1, 0.05)
            assert surf > 0

    def test_add_disk(self):
        with GmshBuilder(verbose=False) as cm:
            surf = cm.add_disk([0, 0, 0], 0.05)
            assert surf > 0

    def test_heal_geometry(self):
        with GmshBuilder(verbose=False) as cm:
            vid = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.heal_geometry(vid)  # Should not raise


# ================================================================
# Boolean operations
# ================================================================

class TestBooleanOps:
    """Test boolean operations."""

    def test_fuse(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.05, 0, 0], [0.1, 0.1, 0.1])
            v_fused = cm.fuse([v1, v2])
            assert v_fused in cm.get_volumes()
            assert v1 not in cm.get_volumes()
            assert v2 not in cm.get_volumes()

    def test_fuse_requires_two(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            with pytest.raises(ValueError, match="at least 2"):
                cm.fuse([v1])

    def test_cut(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.03, 0, 0], [0.04, 0.04, 0.04])
            v_result = cm.cut(v1, v2)
            assert v_result in cm.get_volumes()

    def test_intersect(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.03, 0, 0], [0.1, 0.1, 0.1])
            v_int = cm.intersect(v1, v2)
            assert v_int in cm.get_volumes()

    def test_webcut_x(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            new_ids = cm.webcut_plane(v, 'x', 0.0)
            assert len(new_ids) == 2
            assert v not in cm.get_volumes()
            for nid in new_ids:
                assert nid in cm.get_volumes()

    def test_webcut_y(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            new_ids = cm.webcut_plane(v, 'y', 0.0)
            assert len(new_ids) == 2

    def test_webcut_z(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            new_ids = cm.webcut_plane(v, 'z', 0.0)
            assert len(new_ids) == 2

    def test_double_webcut(self):
        """Cut a box into 4 pieces."""
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            parts_x = cm.webcut_plane(v, 'x', 0.0)
            assert len(parts_x) == 2
            parts_y = cm.webcut_plane(parts_x[0], 'y', 0.0)
            assert len(parts_y) == 2
            assert len(cm.get_volumes()) == 3

    def test_webcut_invalid_axis(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            with pytest.raises(ValueError, match="axis"):
                cm.webcut_plane(v, 'w', 0.0)

    def test_webcut_nonexistent(self):
        with GmshBuilder(verbose=False) as cm:
            with pytest.raises(KeyError):
                cm.webcut_plane(999, 'x', 0.0)

    def test_webcut_cylinder(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            new_ids = cm.webcut_cylinder(v, [0, 0, 0], 'z', 0.03)
            assert len(new_ids) >= 2

    def test_webcut_sphere(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            new_ids = cm.webcut_sphere(v, [0, 0, 0], 0.03)
            assert len(new_ids) >= 2

    def test_webcut_box(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            new_ids = cm.webcut_box(v, [0, 0, 0], [0.04, 0.04, 0.04])
            assert len(new_ids) >= 2

    def test_remove_volume(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            cm.remove_volume(v1)
            assert v1 not in cm.get_volumes()
            assert v2 in cm.get_volumes()

    def test_imprint_and_merge(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.add_box([0.1, 0, 0], [0.1, 0.1, 0.1])
            cm.imprint_and_merge()
            assert len(cm.get_volumes()) >= 2


# ================================================================
# Transforms
# ================================================================

class TestTransforms:
    """Test transform operations."""

    def test_translate(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.translate(v, 0.1, 0, 0)
            assert v in cm.get_volumes()

    def test_copy(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            v2 = cm.copy(v1)
            assert v2 != v1
            assert len(cm.get_volumes()) == 2
            cm.translate(v2, 0.05, 0, 0)

    def test_rotate(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0.05, 0, 0], [0.03, 0.03, 0.03])
            cm.rotate(v, [0, 0, 0], [0, 0, 1], math.pi / 4)
            assert v in cm.get_volumes()

    def test_mirror(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0.05, 0, 0], [0.03, 0.03, 0.03])
            cm.mirror(v, 'x')
            assert v in cm.get_volumes()

    def test_scale(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.scale(v, 2.0)
            bb = cm.get_bounding_box(v)
            assert abs(bb['dx'] - 0.06) < 1e-6

    def test_copy_translate(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            v2 = cm.copy_translate(v1, 0.1, 0, 0)
            assert v2 != v1
            assert len(cm.get_volumes()) == 2

    def test_copy_rotate(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0.05, 0, 0], [0.03, 0.03, 0.03])
            v2 = cm.copy_rotate(v1, [0, 0, 0], [0, 0, 1], math.pi / 2)
            assert v2 != v1
            assert len(cm.get_volumes()) == 2

    def test_copy_mirror(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0.05, 0, 0], [0.03, 0.03, 0.03])
            v2 = cm.copy_mirror(v1, 'x')
            assert v2 != v1
            assert len(cm.get_volumes()) == 2

    def test_array_linear(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            new_ids = cm.array_linear(v, [1, 0, 0], 3, 0.05)
            assert len(new_ids) == 3
            assert len(cm.get_volumes()) == 4

    def test_array_circular(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0.05, 0, 0], [0.02, 0.02, 0.02])
            new_ids = cm.array_circular(v, [0, 0, 0], [0, 0, 1], 3)
            assert len(new_ids) == 3
            assert len(cm.get_volumes()) == 4


# ================================================================
# Mesh control
# ================================================================

class TestMeshControl:
    """Test mesh control methods."""

    def test_set_mesh_size(self):
        with GmshBuilder(verbose=False) as cm:
            cm.set_mesh_size(0.01)
            assert cm._mesh_size == 0.01

    def test_set_divisions(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions(v, 3, 3, 3)

    def test_set_divisions_all(self):
        with GmshBuilder(verbose=False) as cm:
            cm.set_divisions_all(3, 3, 3)
            assert cm._global_divisions == (3, 3, 3)

    def test_set_algorithm(self):
        with GmshBuilder(verbose=False) as cm:
            cm.set_algorithm(6, 1)  # Should not raise

    def test_set_min_max_size(self):
        with GmshBuilder(verbose=False) as cm:
            cm.set_min_max_size(0.001, 0.01)  # Should not raise

    def test_set_recombine_all(self):
        with GmshBuilder(verbose=False) as cm:
            cm.set_recombine_all()  # Should not raise

    def test_add_field_box(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            field = cm.add_field_box(-0.02, 0.02, -0.02, 0.02,
                                      -0.02, 0.02, 0.005, 0.02)
            assert field > 0

    def test_add_field_math(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            field = cm.add_field_math("0.01")
            assert field > 0

    def test_set_curve_divisions(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            curves = cm.get_curves(cm.get_surfaces(v)[0])
            cm.set_curve_divisions(curves[0], 5)

    def test_set_bias(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            curves = cm.get_curves(cm.get_surfaces(v)[0])
            cm.set_bias(curves[0], 5, 1.5)

    def test_set_double_bias(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            curves = cm.get_curves(cm.get_surfaces(v)[0])
            cm.set_double_bias(curves[0], 5, 0.5)


# ================================================================
# Mesh generation
# ================================================================

class TestMeshGeneration:
    """Test mesh generation."""

    def test_generate_hex_single_box(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(3, 3, 3)
            cm.generate()
            info = cm.get_mesh_info()
            assert info['meshed'] is True
            assert info['n_hex'] == 27
            assert info['n_tet'] == 0

    def test_generate_hex_with_mesh_size(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_mesh_size(0.01)
            cm.generate()
            info = cm.get_mesh_info()
            assert info['meshed'] is True
            assert info['n_hex'] == 27

    def test_generate_hex_per_volume(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.04, 0.02, 0.01])
            cm.set_divisions(v, 4, 2, 1)
            cm.generate()
            info = cm.get_mesh_info()
            assert info['n_hex'] == 8

    def test_generate_tet(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_mesh_size(0.01)
            cm.generate(element_type='tet')
            info = cm.get_mesh_info()
            assert info['meshed'] is True
            assert info['n_tet'] > 0

    def test_generate_hex_after_webcut(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.06, 0.03, 0.03])
            cm.webcut_plane(v, 'x', 0.0)
            cm.set_divisions_all(3, 3, 3)
            cm.generate()
            info = cm.get_mesh_info()
            assert info['meshed'] is True
            assert info['n_hex'] == 54

    def test_generate_surface(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.generate_surface(mesh_size=0.01)
            assert cm._meshed is True

    def test_surface_cylinder(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_cylinder([0, 0, 0], 'z', 0.01, 0.05)
            cm.generate_surface(mesh_size=0.005)
            assert cm._meshed is True

    def test_delete_mesh(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            assert cm._meshed is True
            cm.delete_all_mesh()
            assert cm._meshed is False

    def test_set_order_2(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            cm.set_order(2)
            info = cm.get_mesh_info()
            assert info['n_hex'] == 8

    def test_refine(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_mesh_size(0.015)
            cm.generate(element_type='tet')
            info1 = cm.get_mesh_info()
            cm.refine()
            info2 = cm.get_mesh_info()
            # After refinement, should have more elements
            assert info2['n_tet'] > info1['n_tet']


# ================================================================
# Mesh quality
# ================================================================

class TestMeshQuality:
    """Test mesh quality evaluation."""

    def test_get_quality(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            q = cm.get_quality()
            assert len(q) == 8
            assert all(q > 0)

    def test_get_quality_stats(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            stats = cm.get_quality_stats()
            assert stats['n_elements'] == 8
            assert stats['min'] > 0
            assert stats['max'] > 0
            assert stats['n_negative'] == 0

    def test_get_quality_histogram(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_mesh_size(0.01)
            cm.generate(element_type='tet')
            hist = cm.get_quality_histogram(n_bins=5)
            assert len(hist['counts']) == 5
            assert sum(hist['counts']) > 0

    def test_get_scaled_jacobian(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            sj = cm.get_scaled_jacobian()
            assert len(sj) == 8
            # Perfect cubes should have SJ = 1.0
            assert all(sj > 0.9)

    def test_check_inverted(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            check = cm.check_inverted_elements()
            assert check['n_inverted'] == 0
            assert check['n_total'] == 8

    def test_get_worst_elements(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            worst = cm.get_worst_elements(n=3)
            assert len(worst) == 3
            assert all('quality' in w for w in worst)


# ================================================================
# Mesh access
# ================================================================

class TestMeshAccess:
    """Test mesh data access."""

    def test_get_node_count(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            n = cm.get_node_count()
            assert n == 27  # 3*3*3

    def test_get_all_node_coordinates(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            coords = cm.get_all_node_coordinates()
            assert coords.shape == (27, 3)

    def test_get_hex_count(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            assert cm.get_hex_count() == 8

    def test_get_elements(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            elems = cm.get_elements(3)
            assert len(elems) == 8
            assert all(len(e['nodes']) == 8 for e in elems)

    def test_get_closest_node(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            result = cm.get_closest_node(0, 0, 0)
            assert result['distance'] < 0.02

    def test_get_elements_as_numpy(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            data = cm.get_elements_as_numpy()
            assert 'node_coords' in data
            assert 'hex' in data
            assert data['hex'].shape == (8, 8)

    def test_get_element_centroid(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            elems = cm.get_elements(3)
            c = cm.get_element_centroid(elems[0]['tag'])
            assert len(c) == 3

    def test_get_volume_elements(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions(v, 2, 2, 2)
            cm.generate()
            elems = cm.get_volume_elements(v)
            assert len(elems) == 8


# ================================================================
# Physical groups
# ================================================================

class TestPhysicalGroups:
    """Test physical group management."""

    def test_set_physical_group(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            v2 = cm.add_box([0.05, 0, 0], [0.03, 0.03, 0.03])
            cm.set_physical_group(v1, 'iron_core')
            cm.set_physical_group(v2, 'coil')

    def test_physical_group_nonexistent(self):
        with GmshBuilder(verbose=False) as cm:
            with pytest.raises(KeyError):
                cm.set_physical_group(999, 'test')

    def test_add_block(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            v2 = cm.add_box([0.05, 0, 0], [0.03, 0.03, 0.03])
            pg = cm.add_block([v1, v2], name='iron')
            assert pg > 0

    def test_add_sideset(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            surfs = cm.get_surfaces(v)
            pg = cm.add_sideset(surfs[:2], name='inlet')
            assert pg > 0

    def test_get_physical_groups(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_physical_group(v, 'test_group')
            groups = cm.get_physical_groups(dim=3)
            assert len(groups) >= 1
            assert groups[0]['name'] == 'test_group'

    def test_set_material(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_material(v, 'iron', mu_r=1000)
            mats = cm.get_materials()
            assert v in mats
            assert mats[v]['mu_r'] == 1000

    def test_auto_assign_blocks(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.add_box([0.05, 0, 0], [0.03, 0.03, 0.03])
            result = cm.auto_assign_blocks()
            assert len(result) == 2


# ================================================================
# Query
# ================================================================

class TestQuery:
    """Test geometry query methods."""

    def test_get_bounding_box(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.06, 0.04])
            bb = cm.get_bounding_box(v)
            assert abs(bb['dx'] - 0.1) < 1e-6
            assert abs(bb['dy'] - 0.06) < 1e-6
            assert abs(bb['dz'] - 0.04) < 1e-6

    def test_get_center(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0.1, 0.2, 0.3], [0.04, 0.04, 0.04])
            c = cm.get_center(v)
            assert abs(c[0] - 0.1) < 1e-6
            assert abs(c[1] - 0.2) < 1e-6
            assert abs(c[2] - 0.3) < 1e-6

    def test_get_volume_measure(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            vol = cm.get_volume_measure(v)
            assert abs(vol - 1e-3) < 1e-6

    def test_get_surfaces(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            surfs = cm.get_surfaces(v)
            assert len(surfs) == 6

    def test_get_curves(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            surfs = cm.get_surfaces(v)
            curves = cm.get_curves(surfs[0])
            assert len(curves) == 4

    def test_get_surface_area(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            surfs = cm.get_surfaces(v)
            area = cm.get_surface_area(surfs[0])
            assert area > 0

    def test_get_curve_length(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            surfs = cm.get_surfaces(v)
            curves = cm.get_curves(surfs[0])
            length = cm.get_curve_length(curves[0])
            assert abs(length - 0.1) < 1e-6

    def test_get_entity_count(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            assert cm.get_entity_count(3) == 1

    def test_list_entities(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            entities = cm.list_entities()
            assert 3 in entities
            assert entities[3]['count'] == 1

    def test_set_get_name(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_name(v, 'my_box')
            assert cm.get_name(v) == 'my_box'

    def test_find_by_name(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            cm.set_name(v1, 'iron_core')
            cm.set_name(v2, 'copper_coil')
            found = cm.find_volumes_by_name('iron')
            assert v1 in found
            assert v2 not in found

    def test_get_adjacent_volumes(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.1, 0, 0], [0.1, 0.1, 0.1])
            cm.imprint_and_merge()
            vols = cm.get_volumes()
            if len(vols) >= 2:
                adj = cm.get_adjacent_volumes(vols[0])
                assert len(adj) >= 1

    def test_get_point_on_curve(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            surfs = cm.get_surfaces(v)
            curves = cm.get_curves(surfs[0])
            point = cm.get_point_on_curve(curves[0], 0.5)
            assert len(point) == 3

    def test_is_meshed(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            assert not cm.is_meshed()
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            assert cm.is_meshed()


# ================================================================
# Export
# ================================================================

class TestExport:
    """Test export methods."""

    def test_export_msh(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)

            tmp = os.path.join(tempfile.gettempdir(), '_test_cm_export.msh')
            try:
                cm.export(tmp)
                assert os.path.exists(tmp)
                assert os.path.getsize(tmp) > 100
            finally:
                if os.path.exists(tmp):
                    os.unlink(tmp)

    def test_to_radia(self):
        try:
            import radia as rad
        except ImportError:
            pytest.skip("radia not installed")

        rad.UtiDelAll()

        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)
            container = cm.to_radia(mu_r=1000)

        assert container > 0

    def test_export_vtk(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(2, 2, 2)

            tmp = os.path.join(tempfile.gettempdir(), '_test_cm_export.vtk')
            try:
                cm.export_vtk(tmp)
                assert os.path.exists(tmp)
                assert os.path.getsize(tmp) > 100
            finally:
                if os.path.exists(tmp):
                    os.unlink(tmp)

    def test_export_stl(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])

            tmp = os.path.join(tempfile.gettempdir(), '_test_cm_export.stl')
            try:
                cm.export_stl(tmp)
                assert os.path.exists(tmp)
                assert os.path.getsize(tmp) > 100
            finally:
                if os.path.exists(tmp):
                    os.unlink(tmp)


# ================================================================
# Integration
# ================================================================

class TestIntegration:
    """Integration tests combining multiple features."""

    def test_c_magnet_workflow(self):
        """C-type magnet: 3 boxes + divisions + hex mesh."""
        with GmshBuilder(verbose=False) as cm:
            yoke = cm.add_box([0, 0, 0], [0.1, 0.03, 0.02])
            left = cm.add_box([-0.04, 0, 0.02], [0.02, 0.03, 0.02])
            right = cm.add_box([0.04, 0, 0.02], [0.02, 0.03, 0.02])

            cm.set_divisions(yoke, 10, 3, 2)
            cm.set_divisions(left, 2, 3, 2)
            cm.set_divisions(right, 2, 3, 2)

            cm.generate()
            info = cm.get_mesh_info()
            assert info['n_hex'] == (10 * 3 * 2) + 2 * (2 * 3 * 2)
            assert info['n_tet'] == 0

    def test_webcut_mesh_workflow(self):
        """Box -> webcut -> mesh all pieces."""
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.06, 0.06, 0.06])
            parts = cm.webcut_plane(v, 'x', 0.0)
            assert len(parts) == 2
            parts2 = cm.webcut_plane(parts[0], 'z', 0.0)
            assert len(cm.get_volumes()) == 3

            cm.set_divisions_all(3, 3, 3)
            cm.generate()
            info = cm.get_mesh_info()
            assert info['n_hex'] == 3 * 27

    def test_geometry_query_workflow(self):
        """Create geometry, query properties."""
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_name(v, 'test_box')

            bb = cm.get_bounding_box(v)
            assert abs(bb['dx'] - 0.1) < 1e-6

            center = cm.get_center(v)
            assert abs(center[0]) < 1e-6

            vol = cm.get_volume_measure(v)
            assert abs(vol - 1e-3) < 1e-6

            surfs = cm.get_surfaces(v)
            assert len(surfs) == 6

            for s in surfs:
                curves = cm.get_curves(s)
                assert len(curves) == 4

    def test_material_assignment_workflow(self):
        """Assign materials and verify."""
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.15, 0, 0], [0.1, 0.1, 0.1])

            cm.set_material(v1, 'iron', mu_r=1000)
            cm.set_material(v2, 'copper', mu_r=1.0)

            mats = cm.get_materials()
            assert mats[v1]['name'] == 'iron'
            assert mats[v2]['name'] == 'copper'

    def test_quality_workflow(self):
        """Generate mesh and check quality."""
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.03, 0.03, 0.03])
            cm.set_divisions_all(3, 3, 3)
            cm.generate()

            stats = cm.get_quality_stats()
            assert stats['n_elements'] == 27
            assert stats['min'] > 0.9  # Perfect cubes
            assert stats['n_negative'] == 0

            check = cm.check_inverted_elements()
            assert not check['has_inverted']


# ================================================================
# New methods tests (Phase 2: 176 -> 353)
# ================================================================

class TestCoreNew:
    """Test new core methods."""

    def test_set_get_verbose(self):
        with GmshBuilder(verbose=True) as cm:
            assert cm.get_verbose() is True
            cm.set_verbose(False)
            assert cm.get_verbose() is False

    def test_volume_count(self):
        with GmshBuilder(verbose=False) as cm:
            assert cm.get_volume_count() == 0
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            assert cm.get_volume_count() == 1

    def test_all_volume_ids(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            ids = cm.get_all_volume_ids()
            assert ids == [v1, v2]

    def test_has_volume(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            assert cm.has_volume(v)
            assert not cm.has_volume(999)

    def test_get_gmsh_volume_tags(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            tags = cm.get_gmsh_volume_tags(v)
            assert len(tags) >= 1

    def test_model_name(self):
        with GmshBuilder(model_name='test_model', verbose=False) as cm:
            assert cm.get_model_name() == 'test_model'

    def test_is_initialized(self):
        cm = GmshBuilder(verbose=False)
        assert not cm.is_initialized()
        with cm:
            assert cm.is_initialized()
        assert not cm.is_initialized()

    def test_snapshot(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            snap = cm.snapshot()
            assert snap['volume_count'] == 1
            assert snap['initialized'] is True


class TestGeometryNew:
    """Test new geometry creation methods."""

    def test_add_box_at(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box_at(0, 0, 0, 0.1, 0.1, 0.1)
            assert v >= 1

    def test_add_cylinder_vector(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_cylinder_vector([0, 0, 0], [0, 0, 1], 0.05, 0.1)
            assert v >= 1

    def test_add_cone_vector(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_cone_vector([0, 0, 0], [0, 0, 1], 0.05, 0.02, 0.1)
            assert v >= 1

    def test_add_ellipsoid(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_ellipsoid([0, 0, 0], 0.1, 0.05, 0.03)
            assert v >= 1

    def test_add_hollow_cylinder(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_hollow_cylinder([0, 0, 0], 'z', 0.03, 0.05, 0.1)
            assert v >= 1

    def test_add_hollow_sphere(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_hollow_sphere([0, 0, 0], 0.04, 0.05)
            assert v >= 1

    def test_add_cylinder_between(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_cylinder_between([0, 0, 0], [0.1, 0.1, 0], 0.01)
            assert v >= 1

    def test_add_polygon(self):
        with GmshBuilder(verbose=False) as cm:
            pts = [[0, 0], [0.1, 0], [0.1, 0.1], [0, 0.1]]
            surf = cm.add_polygon(pts)
            assert surf >= 1

    def test_add_prism(self):
        with GmshBuilder(verbose=False) as cm:
            pts = [[0, 0], [0.1, 0], [0.05, 0.1]]
            v = cm.add_prism(pts, 0.05)
            assert v >= 1

    def test_extrude_surface(self):
        with GmshBuilder(verbose=False) as cm:
            rect = cm.add_rectangle([0, 0, 0], 0.1, 0.1)
            vids = cm.extrude_surface(rect, 0, 0, 0.05)
            assert len(vids) >= 1

    def test_add_circle(self):
        with GmshBuilder(verbose=False) as cm:
            tag = cm.add_circle([0, 0, 0], 'z', 0.05)
            assert tag >= 1

    def test_copy_geometry(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.copy_geometry(v)
            assert v2 != v

    def test_import_mesh_nonexistent(self):
        with GmshBuilder(verbose=False) as cm:
            with pytest.raises(FileNotFoundError):
                cm.import_mesh('nonexistent.msh')


class TestBooleanOpsNew:
    """Test new boolean operation methods."""

    def test_fragment(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.05, 0, 0], [0.1, 0.1, 0.1])
            new_ids = cm.fragment([v1, v2])
            assert len(new_ids) >= 2

    def test_cut_keep_tool(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_sphere([0, 0, 0], 0.04)
            result = cm.cut_keep_tool(v1, v2)
            assert result >= 1
            assert cm.has_volume(v2)  # Tool preserved

    def test_subtract_list(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_sphere([0.03, 0, 0], 0.02)
            v3 = cm.add_sphere([-0.03, 0, 0], 0.02)
            result = cm.subtract_list(v1, [v2, v3])
            assert result >= 1

    def test_fuse_all(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.add_box([0.05, 0, 0], [0.1, 0.1, 0.1])
            result = cm.fuse_all()
            assert result >= 1
            assert cm.get_volume_count() == 1

    def test_remove_volumes(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            v3 = cm.add_box([0.4, 0, 0], [0.1, 0.1, 0.1])
            cm.remove_volumes([v1, v3])
            assert cm.get_volume_count() == 1
            assert cm.has_volume(v2)

    def test_keep_only(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            v3 = cm.add_box([0.4, 0, 0], [0.1, 0.1, 0.1])
            cm.keep_only([v2])
            assert cm.get_volume_count() == 1
            assert cm.has_volume(v2)

    def test_separate_volumes(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            ids = cm.separate_volumes(v)
            assert len(ids) >= 1


class TestTransformsNew:
    """Test new transform methods."""

    def test_translate_to(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.translate_to(v, [0.5, 0.5, 0.5])
            c = cm.get_center(v)
            assert abs(c[0] - 0.5) < 1e-6
            assert abs(c[1] - 0.5) < 1e-6
            assert abs(c[2] - 0.5) < 1e-6

    def test_rotate_degrees(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0.1, 0, 0], [0.02, 0.02, 0.02])
            cm.rotate_degrees(v, [0, 0, 0], [0, 0, 1], 90)
            c = cm.get_center(v)
            assert abs(c[0]) < 1e-6
            assert abs(c[1] - 0.1) < 1e-6

    def test_copy_scale(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.copy_scale(v, 2.0)
            assert v2 != v

    def test_array_grid(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.01, 0.01, 0.01])
            ids = cm.array_grid(v, 1, 1, 0, 0.02, 0.02, 0.02)
            assert len(ids) == 3  # 2*2*1 - 1 = 3

    def test_center_at_origin(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([1, 2, 3], [0.1, 0.1, 0.1])
            cm.center_at_origin(v)
            c = cm.get_center(v)
            assert abs(c[0]) < 1e-6
            assert abs(c[1]) < 1e-6
            assert abs(c[2]) < 1e-6

    def test_translate_all(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            cm.translate_all(1, 0, 0)
            for vid in cm.get_all_volume_ids():
                c = cm.get_center(vid)
                assert c[0] > 0.5

    def test_symmetrize(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0.05, 0, 0], [0.1, 0.1, 0.1])
            result = cm.symmetrize(v, 'x')
            assert result >= 1

    def test_swap_volumes(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.5, 0, 0], [0.2, 0.2, 0.2])
            tags1 = cm.get_gmsh_volume_tags(v1)
            tags2 = cm.get_gmsh_volume_tags(v2)
            cm.swap_volumes(v1, v2)
            assert cm.get_gmsh_volume_tags(v1) == tags2
            assert cm.get_gmsh_volume_tags(v2) == tags1


class TestMeshControlNew:
    """Test new mesh control methods."""

    def test_add_field_ball(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            f = cm.add_field_ball(0, 0, 0, 0.05, 0.001, 0.01)
            assert f >= 1

    def test_add_field_max(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            f1 = cm.add_field_math("0.01")
            f2 = cm.add_field_math("0.02")
            f3 = cm.add_field_max([f1, f2])
            assert f3 >= 1

    def test_get_field_count(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            assert cm.get_field_count() == 0
            cm.add_field_math("0.01")
            assert cm.get_field_count() == 1

    def test_remove_field(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            f = cm.add_field_math("0.01")
            cm.remove_field(f)
            assert cm.get_field_count() == 0

    def test_set_transfinite_surface(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            surfs = cm.get_all_surfaces()
            # Should not raise
            cm.set_transfinite_surface(surfs[0])

    def test_get_mesh_options(self):
        with GmshBuilder(verbose=False) as cm:
            opts = cm.get_mesh_options()
            assert 'Algorithm' in opts
            assert 'Algorithm3D' in opts

    def test_set_algorithm_2d_3d(self):
        with GmshBuilder(verbose=False) as cm:
            cm.set_algorithm_2d(5)
            cm.set_algorithm_3d(4)
            opts = cm.get_mesh_options()
            assert opts['Algorithm'] == 5
            assert opts['Algorithm3D'] == 4

    def test_set_mesh_size_factor(self):
        with GmshBuilder(verbose=False) as cm:
            cm.set_mesh_size_factor(2.0)
            # Should not raise


class TestMeshGenerateNew:
    """Test new mesh generation methods."""

    def test_generate_tet(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.generate_tet(mesh_size=0.05)
            assert cm.get_tet_count() > 0

    def test_generate_hex_convenience(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate_hex()
            assert cm.get_hex_count() == 8

    def test_generate_and_optimize(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate_and_optimize()
            assert cm.get_hex_count() == 8

    def test_generate_surface_quad(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.generate_surface_quad(mesh_size=0.05)
            assert cm.get_quad_count() > 0

    def test_convert_first_order(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            cm.convert_to_second_order()
            n2 = cm.get_node_count()
            cm.convert_to_first_order()
            n1 = cm.get_node_count()
            assert n1 < n2

    def test_renumber_nodes(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            cm.renumber_nodes()
            cm.renumber_elements()
            # Should not raise

    def test_get_mesh_statistics(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            stats = cm.get_mesh_statistics()
            assert stats['n_nodes'] > 0
            assert stats['n_elements'] > 0

    def test_embed_point(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            p = cm.add_point(0.05, 0.05, 0.05)
            vols = cm.get_all_surfaces()
            # embed point in first volume's gmsh tag
            tags = cm.get_gmsh_volume_tags(1)
            cm.embed_point_in_volume(p, tags[0])


class TestMeshQualityNew:
    """Test new mesh quality methods."""

    def test_get_quality_range(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            qmin, qmax = cm.get_quality_range()
            assert qmin > 0
            assert qmax <= 1.0 + 1e-6

    def test_count_below_threshold(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            count = cm.count_elements_below_threshold(threshold=0.1)
            assert count == 0  # Perfect cubes

    def test_get_quality_summary(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            summary = cm.get_quality_summary()
            assert 'minSICN' in summary

    def test_get_volume_quality(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            stats = cm.get_volume_quality(v)
            assert stats['n_elements'] == 8

    def test_get_edge_length_range(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            emin, emax = cm.get_edge_length_range()
            assert abs(emin - 0.05) < 1e-6
            assert abs(emax - 0.05) < 1e-6

    def test_smooth_mesh(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            cm.smooth_mesh(n_iterations=1)
            # Should not raise


class TestMeshAccessNew:
    """Test new mesh access methods."""

    def test_get_node_coordinates_batch(self):
        import numpy as np
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            nodes = cm.get_nodes()
            coords = cm.get_node_coordinates_batch(nodes[:5])
            assert coords.shape == (5, 3)

    def test_get_pyramid_count(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            assert cm.get_pyramid_count() == 0
            assert cm.get_prism_count() == 0

    def test_get_element_face_count(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            elems = cm.get_elements(dim=3)
            faces = cm.get_element_face_count(elems[0]['tag'])
            assert faces == 6  # Hexahedron

    def test_get_nodes_on_surface(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            surfs = cm.get_all_surfaces()
            nodes = cm.get_nodes_on_surface(surfs[0])
            assert len(nodes) > 0

    def test_get_internal_nodes(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(3, 3, 3)
            cm.generate()
            internal = cm.get_internal_nodes()
            assert len(internal) > 0  # 3x3x3 has internal nodes

    def test_get_node_to_element_map(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            n2e = cm.get_node_to_element_map()
            assert len(n2e) > 0

    def test_get_edge_nodes(self):
        import numpy as np
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            edges = cm.get_edge_nodes()
            assert edges.shape[1] == 2
            assert edges.shape[0] > 0

    def test_get_elements_in_sphere(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(3, 3, 3)
            cm.generate()
            elems = cm.get_elements_in_sphere([0, 0, 0], 0.05)
            assert len(elems) > 0
            assert len(elems) < 27  # Not all elements


class TestPhysicalGroupsNew:
    """Test new physical group methods."""

    def test_add_block_by_name(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.add_block_by_name([v], 'iron_core')
            assert cm.get_block_count() >= 1

    def test_get_all_blocks(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.add_block_by_name([v], 'block1')
            blocks = cm.get_all_blocks()
            assert len(blocks) >= 1
            assert blocks[0]['name'] == 'block1'

    def test_material_management(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_material(v, 'iron', mu_r=1000)
            assert cm.has_material(v)
            mat = cm.get_material(v)
            assert mat['name'] == 'iron'
            cm.clear_materials()
            assert not cm.has_material(v)

    def test_set_material_by_name(self):
        with GmshBuilder(verbose=False) as cm:
            v1 = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            v2 = cm.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            cm.set_material_by_name([v1, v2], 'copper',
                                     {'mu_r': 1.0, 'sigma': 5.8e7})
            assert cm.has_material(v1)
            assert cm.has_material(v2)

    def test_auto_assign_all(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.add_box([0.2, 0, 0], [0.1, 0.1, 0.1])
            cm.auto_assign_all()
            assert cm.get_block_count() >= 2

    def test_rename_physical_group(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.add_block_by_name([v], 'old_name')
            blocks = cm.get_all_blocks()
            cm.rename_physical_group(blocks[0]['tag'], 3, 'new_name')
            blocks2 = cm.get_all_blocks()
            assert blocks2[0]['name'] == 'new_name'


class TestQueryNew:
    """Test new query methods."""

    def test_get_all_surfaces(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            surfs = cm.get_all_surfaces()
            assert len(surfs) == 6

    def test_get_all_curves(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            curves = cm.get_all_curves()
            assert len(curves) == 12

    def test_get_all_points(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            pts = cm.get_all_points()
            assert len(pts) == 8

    def test_get_bounding_box_all(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.add_box([1, 0, 0], [0.1, 0.1, 0.1])
            bb = cm.get_bounding_box_all()
            assert bb['xmin'] < -0.04
            assert bb['xmax'] > 1.04

    def test_get_center_of_mass(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            com = cm.get_center_of_mass(v)
            assert abs(com[0]) < 1e-6
            assert abs(com[1]) < 1e-6
            assert abs(com[2]) < 1e-6

    def test_get_mass(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            # Volume of 0.1^3 = 0.001 m^3
            mass = cm.get_mass(v, density=1000)
            assert abs(mass - 1.0) < 0.01

    def test_get_surface_count(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            assert cm.get_surface_count(v) == 6

    def test_get_model_summary(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            summary = cm.get_model_summary()
            assert summary['volumes'] == 1
            assert summary['surfaces'] == 6
            assert summary['curves'] == 12
            assert summary['points'] == 8

    def test_get_entities(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            vols = cm.get_entities(3)
            assert len(vols) == 1

    def test_get_volume_topology(self):
        with GmshBuilder(verbose=False) as cm:
            v = cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            topo = cm.get_volume_topology(v)
            assert 'surfaces' in topo
            assert len(topo['surfaces']) == 6


class TestExportNew:
    """Test new export methods."""

    def test_export_msh_v4(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            with tempfile.NamedTemporaryFile(suffix='.msh', delete=False) as f:
                path = f.name
            try:
                cm.export_msh(path, version=4.1)
                assert os.path.exists(path)
                assert os.path.getsize(path) > 0
            finally:
                os.unlink(path)

    def test_export_geo(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            with tempfile.NamedTemporaryFile(suffix='.geo_unrolled',
                                              delete=False) as f:
                path = f.name
            try:
                cm.export_geo(path)
                assert os.path.exists(path)
            finally:
                os.unlink(path)

    def test_to_dict(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            d = cm.to_dict()
            assert 'nodes' in d
            assert 'elements' in d
            assert len(d['nodes']) > 0

    def test_get_export_formats(self):
        with GmshBuilder(verbose=False) as cm:
            fmts = cm.get_export_formats()
            assert 'msh' in fmts
            assert 'vtk' in fmts
            assert 'step' in fmts

    def test_export_auto(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            cm.set_divisions_all(2, 2, 2)
            cm.generate()
            with tempfile.NamedTemporaryFile(suffix='.vtk',
                                              delete=False) as f:
                path = f.name
            try:
                cm.export_auto(path)
                assert os.path.exists(path)
            finally:
                os.unlink(path)

    def test_export_step(self):
        with GmshBuilder(verbose=False) as cm:
            cm.add_box([0, 0, 0], [0.1, 0.1, 0.1])
            with tempfile.NamedTemporaryFile(suffix='.step',
                                              delete=False) as f:
                path = f.name
            try:
                cm.export_step(path)
                assert os.path.exists(path)
            finally:
                os.unlink(path)
