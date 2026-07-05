"""
Test ObjTetrahedron and ObjHexahedron APIs.

These tests verify that the tetrahedron and hexahedron creation APIs
work correctly for creating magnetic elements.
"""
import sys
import os
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/radia'))
import radia as rad


class TestObjTetrahedron:
    """Tests for ObjTetrahedron API."""

    def setup_method(self):
        rad.UtiDelAll()

    def test_basic_creation(self):
        """Test basic tetrahedron creation."""
        vertices = [
            [0, 0, 0], [0.1, 0, 0], [0.05, 0.0866, 0], [0.05, 0.0289, 0.0816]
        ]
        tetra = rad.ObjTetrahedron(vertices, [0, 0, 954930])
        assert tetra > 0

    def test_creation_without_magnetization(self):
        """Test tetrahedron creation without magnetization (defaults to zero)."""
        vertices = [
            [0, 0, 0], [0.1, 0, 0], [0.05, 0.0866, 0], [0.05, 0.0289, 0.0816]
        ]
        tetra = rad.ObjTetrahedron(vertices, [0, 0, 0])
        assert tetra > 0

    def test_field_computation(self):
        """Test field computation from tetrahedron."""
        vertices = [
            [0, 0, 0], [0.1, 0, 0], [0.05, 0.0866, 0], [0.05, 0.0289, 0.0816]
        ]
        tetra = rad.ObjTetrahedron(vertices, [0, 0, 954930])
        B = rad.Fld(tetra, 'b', [0.05, 0.03, 0.2])
        assert len(B) == 3
        assert any(abs(b) > 0 for b in B)

    def test_wrong_vertex_count_3(self):
        """Test that 3 vertices raises error."""
        vertices = [[0, 0, 0], [1, 0, 0], [0.5, 0.866, 0]]
        with pytest.raises(Exception) as excinfo:
            rad.ObjTetrahedron(vertices, [0, 0, 954930])
        assert "4 vertices" in str(excinfo.value)

    def test_wrong_vertex_count_5(self):
        """Test that 5 vertices raises error."""
        vertices = [
            [0, 0, 0], [1, 0, 0], [0.5, 0.866, 0], [0.5, 0.289, 0.816], [0, 0, 1]
        ]
        with pytest.raises(Exception) as excinfo:
            rad.ObjTetrahedron(vertices, [0, 0, 954930])
        assert "4 vertices" in str(excinfo.value)

    def test_field_consistency(self):
        """Test that field computation is consistent across multiple calls."""
        vertices = [
            [0, 0, 0], [0.1, 0, 0], [0.05, 0.0866, 0], [0.05, 0.0289, 0.0816]
        ]

        rad.UtiDelAll()
        tetra = rad.ObjTetrahedron(vertices, [0, 0, 954930])
        B1 = rad.Fld(tetra, 'b', [0.05, 0.03, 0.2])
        B2 = rad.Fld(tetra, 'b', [0.05, 0.03, 0.2])

        np.testing.assert_allclose(B1, B2, rtol=1e-10)


class TestObjHexahedron:
    """Tests for ObjHexahedron API."""

    def setup_method(self):
        rad.UtiDelAll()

    def test_basic_creation(self):
        """Test basic hexahedron creation."""
        s = 0.05
        vertices = [
            [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
            [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s]
        ]
        hex_obj = rad.ObjHexahedron(vertices, [0, 0, 954930])
        assert hex_obj > 0

    def test_creation_without_magnetization(self):
        """Test hexahedron creation without magnetization (defaults to zero)."""
        s = 0.05
        vertices = [
            [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
            [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s]
        ]
        hex_obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        assert hex_obj > 0

    def test_field_computation(self):
        """Test field computation from hexahedron."""
        s = 0.05
        vertices = [
            [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
            [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s]
        ]
        hex_obj = rad.ObjHexahedron(vertices, [0, 0, 954930])
        B = rad.Fld(hex_obj, 'b', [0, 0, 0.15])
        assert len(B) == 3
        assert abs(B[2]) > 0  # Z-component should be non-zero

    def test_wrong_vertex_count_4(self):
        """Test that 4 vertices raises error."""
        vertices = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
        with pytest.raises(Exception) as excinfo:
            rad.ObjHexahedron(vertices, [0, 0, 954930])
        assert "8 vertices" in str(excinfo.value)

    def test_wrong_vertex_count_10(self):
        """Test that 10 vertices raises error."""
        s = 0.05
        vertices = [
            [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
            [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s],
            [0, 0, 0], [0, 0, 1]  # Extra vertices
        ]
        with pytest.raises(Exception) as excinfo:
            rad.ObjHexahedron(vertices, [0, 0, 954930])
        assert "8 vertices" in str(excinfo.value)

    def test_field_consistency(self):
        """Test that field computation is consistent across multiple calls."""
        s = 0.05
        vertices = [
            [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
            [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s]
        ]

        rad.UtiDelAll()
        hex_obj = rad.ObjHexahedron(vertices, [0, 0, 954930])
        B1 = rad.Fld(hex_obj, 'b', [0, 0, 0.15])
        B2 = rad.Fld(hex_obj, 'b', [0, 0, 0.15])

        np.testing.assert_allclose(B1, B2, rtol=1e-10)


class TestMaterialApplication:
    """Test material application to tetrahedra and hexahedra."""

    def setup_method(self):
        rad.UtiDelAll()

    def test_tetrahedron_with_linear_material_requires_mesh_backed_hdiv(self):
        """Mesh-less soft iron must be modeled through the mesh-backed HDiv route."""
        vertices = [
            [0, 0, 0], [0.1, 0, 0], [0.05, 0.0866, 0], [0.05, 0.0289, 0.0816]
        ]
        tetra = rad.ObjTetrahedron(vertices, [0, 0, 0])

        # Apply linear material (mu_r = 1000)
        mat = rad.MatLin(1000)
        rad.MatApl(tetra, mat)

        with pytest.raises(RuntimeError, match="[Mm]esh-less soft iron"):
            rad.Solve(tetra, 0.0001, 100, 0)

    def test_hexahedron_soft_iron_requires_mesh_backed_hdiv(self):
        """ObjHexahedron + MatLin soft iron is not a supported mesh-less solve path."""
        rad.set_demag_backend("auto")
        s = 0.05
        vertices = [
            [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
            [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s]
        ]
        hex_obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        rad.MatApl(hex_obj, rad.MatLin(1000))
        with pytest.raises(RuntimeError, match="[Mm]esh-less soft iron"):
            rad.Solve(hex_obj, 0.0001, 100, 0)


class TestContainer:
    """Test adding tetrahedra and hexahedra to containers."""

    def setup_method(self):
        rad.UtiDelAll()

    def test_tetrahedra_in_container(self):
        """Test creating container with multiple tetrahedra."""
        vertices1 = [[0, 0, 0], [0.1, 0, 0], [0.05, 0.0866, 0], [0.05, 0.0289, 0.0816]]
        vertices2 = [[0.2, 0, 0], [0.3, 0, 0], [0.25, 0.0866, 0], [0.25, 0.0289, 0.0816]]

        t1 = rad.ObjTetrahedron(vertices1, [0, 0, 954930])
        t2 = rad.ObjTetrahedron(vertices2, [0, 0, 954930])

        container = rad.ObjCnt([t1, t2])
        assert container > 0

        # Field from container should be sum of individual fields
        B = rad.Fld(container, 'b', [0.15, 0, 0.2])
        assert len(B) == 3

    def test_mixed_container(self):
        """Test container with both tetrahedra and hexahedra."""
        tetra_verts = [[0, 0, 0], [0.1, 0, 0], [0.05, 0.0866, 0], [0.05, 0.0289, 0.0816]]
        hex_verts = [
            [0.2, -0.05, -0.05], [0.3, -0.05, -0.05], [0.3, 0.05, -0.05], [0.2, 0.05, -0.05],
            [0.2, -0.05, 0.05], [0.3, -0.05, 0.05], [0.3, 0.05, 0.05], [0.2, 0.05, 0.05]
        ]

        tetra = rad.ObjTetrahedron(tetra_verts, [0, 0, 954930])
        hex_obj = rad.ObjHexahedron(hex_verts, [0, 0, 954930])

        container = rad.ObjCnt([tetra, hex_obj])
        assert container > 0

        B = rad.Fld(container, 'b', [0.15, 0, 0.2])
        assert len(B) == 3


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
