#!/usr/bin/env python
"""
radia_ngsolve_utils.py - Utilities for Radia-NGSolve integration

This module provides high-level utilities for working with Radia and NGSolve together,
including mesh import and VTK export functionality.

Author: Radia Development Team
Created: 2025-11-22
Version: 1.0.0

Key Features
------------
- Import Netgen/OCC tetrahedral meshes to Radia
- Export Radia geometry to VTK for visualization
- Export NGSolve solutions to VTK
- Unified interface for both solvers

Functions
---------
create_radia_from_mesh : Import Netgen mesh to Radia with VTK export
export_ngsolve_vtk : Export NGSolve solution to VTK
export_radia_vtk : Export Radia geometry to VTK

Example
-------
>>> import radia as rad
>>> from ngsolve import Mesh
>>> from netgen.occ import Box, OCCGeometry
>>> from radia_ngsolve_utils import create_radia_from_mesh, export_radia_vtk
>>>
>>> # Setup
>>> rad.FldUnits('m')
>>>
>>> # Create mesh
>>> geo = OCCGeometry(Box((-0.5, -0.5, -0.5), (0.5, 0.5, 0.5)))
>>> mesh = Mesh(geo.GenerateMesh(maxh=0.2))
>>>
>>> # Import to Radia and export VTK
>>> mag = create_radia_from_mesh(
...     mesh,
...     material={'magnetization': [0, 0, 1.2]},
...     vtk_filename='magnet'
... )
>>>
>>> # Apply material and solve
>>> rad.MatApl(mag, rad.MatStd('NdFeB', 1.2))
>>> rad.Solve(mag, 0.0001, 10000)
>>>
>>> # Export final geometry
>>> export_radia_vtk(mag, 'magnet_solved')
"""

import os
import radia as rad
from .netgen_mesh_import import netgen_mesh_to_radia
from .radia_vtk_export import exportGeometryToVTK


def create_radia_from_mesh(mesh, material=None, units='m', combine=True,
                           verbose=False, material_filter=None, vtk_filename=None):
    """
    Import Netgen mesh to Radia and optionally export to VTK.

    This is a convenience wrapper around netgen_mesh_to_radia that adds
    automatic VTK export functionality.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        NGSolve mesh object to import
    material : dict, optional
        Material properties. Should contain 'magnetization' key with [Mx, My, Mz] in Tesla.
        Default: {'magnetization': [0, 0, 0]}
    units : str, optional
        Length units: 'm' (meters) or 'mm' (millimeters)
        Default: 'm' (recommended for NGSolve integration)
    combine : bool, optional
        If True, combine all tetrahedra into single container
        If False, return list of individual tetrahedra
        Default: True
    verbose : bool, optional
        Print detailed import information
        Default: False
    material_filter : str, optional
        Filter mesh by material name (e.g., 'magnetic', 'air')
        If None, import all materials
        Default: None
    vtk_filename : str, optional
        If provided, export geometry to VTK file with this name (without .vtk extension)
        Default: None (no VTK export)

    Returns
    -------
    int or list
        Radia object ID (if combine=True) or list of object IDs (if combine=False)

    Example
    -------
    >>> from netgen.occ import Box, OCCGeometry
    >>> from ngsolve import Mesh
    >>> import radia as rad
    >>> from radia_ngsolve_utils import create_radia_from_mesh
    >>>
    >>> rad.FldUnits('m')
    >>> geo = OCCGeometry(Box((0, 0, 0), (0.1, 0.1, 0.1)))
    >>> mesh = Mesh(geo.GenerateMesh(maxh=0.03))
    >>>
    >>> magnet = create_radia_from_mesh(
    ...     mesh,
    ...     material={'magnetization': [0, 0, 1.2]},
    ...     vtk_filename='my_magnet'
    ... )
    >>> # Creates my_magnet.vtk in current directory
    """
    # Import mesh to Radia
    radia_obj = netgen_mesh_to_radia(
        mesh=mesh,
        material=material,
        units=units,
        combine=combine,
        verbose=verbose,
        material_filter=material_filter
    )

    # Export to VTK if requested
    if vtk_filename is not None:
        export_radia_vtk(radia_obj, vtk_filename)

    return radia_obj


def export_radia_vtk(radia_obj, filename):
    """
    Export Radia geometry to VTK file.

    Parameters
    ----------
    radia_obj : int
        Radia object ID
    filename : str
        Output filename without .vtk extension

    Returns
    -------
    str
        Full path to created VTK file

    Example
    -------
    >>> import radia as rad
    >>> from radia_ngsolve_utils import export_radia_vtk
    >>>
    >>> magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
    >>> vtk_path = export_radia_vtk(magnet, 'my_magnet')
    >>> print(f"VTK exported to: {vtk_path}")
    """
    # Ensure .vtk extension
    if not filename.endswith('.vtk'):
        base_filename = filename
        filename = filename + '.vtk'
    else:
        base_filename = filename[:-4]

    # Export using radia_vtk_export
    exportGeometryToVTK(radia_obj, base_filename)

    # Get absolute path
    vtk_path = os.path.abspath(filename)

    return vtk_path


def export_ngsolve_vtk(mesh, gridfunction, filename, field_name='solution'):
    """
    Export NGSolve solution to VTK file.

    Uses NGSolve's built-in VTK export functionality.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        NGSolve mesh
    gridfunction : ngsolve.GridFunction
        Solution to export (e.g., H-field, B-field, potential)
    filename : str
        Output filename without .vtk extension
    field_name : str, optional
        Name for the field in VTK file
        Default: 'solution'

    Returns
    -------
    str
        Full path to created VTK file

    Example
    -------
    >>> from ngsolve import *
    >>> from netgen.occ import Box, OCCGeometry
    >>> from radia_ngsolve_utils import export_ngsolve_vtk
    >>>
    >>> geo = OCCGeometry(Box((0, 0, 0), (1, 1, 1)))
    >>> mesh = Mesh(geo.GenerateMesh(maxh=0.2))
    >>> fes = H1(mesh, order=2)
    >>> gfu = GridFunction(fes)
    >>> # ... solve for gfu ...
    >>>
    >>> vtk_path = export_ngsolve_vtk(mesh, gfu, 'solution', 'phi')
    >>> print(f"NGSolve VTK exported to: {vtk_path}")
    """
    from ngsolve import VTKOutput

    # Ensure .vtk extension
    if filename.endswith('.vtk'):
        filename = filename[:-4]

    # Create VTK output
    vtk = VTKOutput(
        ma=mesh,
        coefs=[gridfunction],
        names=[field_name],
        filename=filename,
        subdivision=2
    )
    vtk.Do()

    # Get absolute path
    vtk_path = os.path.abspath(filename + '.vtk')

    return vtk_path


def export_comparison_vtk(mesh_ngsolve, gfu_ngsolve, radia_obj,
                         prefix='comparison', ngsolve_field_name='H_ngsolve'):
    """
    Export both NGSolve and Radia results to VTK for comparison.

    Creates two VTK files:
    - {prefix}_ngsolve.vtk : NGSolve solution
    - {prefix}_radia.vtk : Radia geometry

    Parameters
    ----------
    mesh_ngsolve : ngsolve.Mesh
        NGSolve mesh
    gfu_ngsolve : ngsolve.GridFunction
        NGSolve solution
    radia_obj : int
        Radia object ID
    prefix : str, optional
        Prefix for output filenames
        Default: 'comparison'
    ngsolve_field_name : str, optional
        Name for NGSolve field in VTK
        Default: 'H_ngsolve'

    Returns
    -------
    tuple of str
        (ngsolve_vtk_path, radia_vtk_path)

    Example
    -------
    >>> from radia_ngsolve_utils import export_comparison_vtk
    >>>
    >>> # After solving both NGSolve and Radia...
    >>> ngsolve_vtk, radia_vtk = export_comparison_vtk(
    ...     mesh, gfu, radia_system,
    ...     prefix='cube_benchmark',
    ...     ngsolve_field_name='H_field'
    ... )
    >>> print(f"NGSolve: {ngsolve_vtk}")
    >>> print(f"Radia:   {radia_vtk}")
    """
    ngsolve_filename = f"{prefix}_ngsolve"
    radia_filename = f"{prefix}_radia"

    ngsolve_vtk = export_ngsolve_vtk(mesh_ngsolve, gfu_ngsolve,
                                     ngsolve_filename, ngsolve_field_name)
    radia_vtk = export_radia_vtk(radia_obj, radia_filename)

    return ngsolve_vtk, radia_vtk


# Module-level constants
__version__ = '1.0.0'
__all__ = [
    'create_radia_from_mesh',
    'export_radia_vtk',
    'export_ngsolve_vtk',
    'export_comparison_vtk'
]
