"""
GmshBuilder: Gmsh-based mesh generation using GMSH OCC kernel.

High-level API providing ~700 methods covering geometry creation,
boolean operations, mesh control, quality evaluation, and export.

Usage:
    from radia.gmsh_builder import GmshBuilder

    with GmshBuilder() as gb:
        box = gb.add_box([0, 0, 0], [0.1, 0.02, 0.03])
        gb.webcut_plane(box, 'x', 0.03)
        gb.set_mesh_size(0.005)
        gb.generate()
        gb.export('output.msh')

Requires: pip install gmsh

Part of Radia project.
"""

from .core import CoreMixin
from .geometry import GeometryMixin
from .boolean_ops import BooleanMixin
from .transforms import TransformMixin
from .mesh_control import MeshControlMixin
from .mesh_generate import MeshGenerateMixin
from .mesh_quality import MeshQualityMixin
from .mesh_access import MeshAccessMixin
from .physical_groups import PhysicalGroupMixin
from .query import QueryMixin
from .export import ExportMixin
from .entity_wrappers import EntityWrapperMixin
from .geometry_analysis import GeometryAnalysisMixin
from .groups import GroupMixin
from .post_processing import PostProcessingMixin
from .options import OptionsMixin


class GmshBuilder(CoreMixin, GeometryMixin, BooleanMixin, TransformMixin,
                MeshControlMixin, MeshGenerateMixin, MeshQualityMixin,
                MeshAccessMixin, PhysicalGroupMixin, QueryMixin,
                ExportMixin, EntityWrapperMixin, GeometryAnalysisMixin,
                GroupMixin, PostProcessingMixin, OptionsMixin):
    """Gmsh-based mesh generation interface backed by GMSH OCC kernel.

    Combines ~700 methods from 16 mixin modules covering:
    - Geometry creation (box, cylinder, sphere, cone, torus, extrude, revolve)
    - Boolean operations (fuse, cut, intersect, webcut)
    - Transforms (translate, rotate, mirror, scale, copy, arrays)
    - Mesh control (size, divisions, bias, boundary layer, fields)
    - Mesh generation (hex, tet, surface, refine, optimize)
    - Mesh quality (scaled Jacobian, aspect ratio, histograms)
    - Mesh access (nodes, elements, numpy export, Jacobians, basis functions)
    - Physical groups (blocks, sidesets, nodesets, materials)
    - Queries (bounding box, center, surfaces, adjacency, attributes)
    - Export (.msh, Radia, NGSolve, VTK, STEP, STL, PyVista)
    - Entity wrappers (volume/surface/curve/vertex queries)
    - Geometry analysis (cleanup, differential geometry, topology)
    - Groups (selection groups, set operations)
    - Post-processing (views, field data visualization)
    - Options (GMSH global settings)

    Must be used as context manager (handles GMSH initialize/finalize).

    Parameters
    ----------
    model_name : str
        GMSH model name.
    verbose : bool
        Print progress information.
    """

    def __init__(self, model_name='gmsh_model', verbose=True):
        super().__init__(model_name=model_name, verbose=verbose)


__all__ = ['GmshBuilder']
