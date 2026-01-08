# Radia Python package
# This module re-exports all symbols from the C++ extension module (radia.pyd)
# so that 'import radia' works correctly when installed via pip

__version__ = "1.4.4"

# Add package directory to DLL search path (Windows)
# This is needed for finding Intel MKL DLL (mkl_rt.2.dll)
import os
import sys

_package_dir = os.path.dirname(os.path.abspath(__file__))

# Add DLL directory for Windows (Python 3.8+)
if sys.platform == 'win32':
    if hasattr(os, 'add_dll_directory'):
        os.add_dll_directory(_package_dir)
        # Add Intel MKL path if available
        mkl_bin = r'C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin'
        if os.path.isdir(mkl_bin):
            os.add_dll_directory(mkl_bin)
        # Add Intel Compiler runtime path if available
        icx_bin = r'C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin'
        if os.path.isdir(icx_bin):
            os.add_dll_directory(icx_bin)
    # Also add to PATH as fallback for older methods
    if _package_dir not in os.environ.get('PATH', ''):
        os.environ['PATH'] = _package_dir + os.pathsep + os.environ.get('PATH', '')

# Import all symbols from the C++ extension module
try:
    from radia.radia import *
except ImportError:
    # Fallback for development: try importing from the same directory
    try:
        from .radia import *
    except ImportError as e:
        raise ImportError(
            "Failed to import radia C++ extension module (radia.pyd). "
            "Ensure the package was built correctly with Build.ps1 before installation. "
            f"Package directory: {_package_dir}"
        ) from e

# ESIM (Effective Surface Impedance Method) for induction heating analysis
# Import convenience functions for ESIM workpiece creation
try:
    from .esim_cell_problem import (
        ESIMCellProblemSolver,
        BHCurveInterpolator,
        ComplexPermeabilityInterpolator,
        ESITable,
        generate_esi_table_from_bh_curve,
    )
    from .esim_workpiece import (
        ESIMWorkpiece,
        SurfacePanel,
        create_esim_block,
        create_esim_cylinder,
    )
    from .esim_coupled_solver import (
        InductionHeatingCoil,
        ESIMCoupledSolver,
        solve_induction_heating,
    )
    from .esim_vtk_export import (
        ESIMVTKOutput,
        export_esim_workpiece_vtk,
        export_esim_coil_field_vtk,
        export_esim_combined_vtk,
    )
    ESIM_AVAILABLE = True
except ImportError:
    # ESIM requires scipy, which may not be installed
    ESIM_AVAILABLE = False
