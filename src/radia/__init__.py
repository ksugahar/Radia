# Radia Python package
# This module re-exports all symbols from the C++ extension module (_radia_pybind.pyd)
# so that 'import radia' works correctly when installed via pip
#
# Module naming: The C++ extension is named '_radia_pybind' (with underscore) so that
# 'import radia' uses this Python package with __init__.py, which then imports
# from _radia_pybind. This follows NGSolve's pattern.
#
# pybind11 Migration Complete (2026-01):
# All bindings now use pybind11 exclusively.

__version__ = "1.8.1"

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

# Import all symbols from the pybind11 C++ extension module (_radia_pybind.pyd)
try:
    from ._radia_pybind import *
except ImportError as e:
    raise ImportError(
        "Failed to import radia pybind11 module (_radia_pybind.pyd). "
        "Ensure the package was built correctly with BuildMSVC.ps1 before installation. "
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
        # WPT (Wireless Power Transfer) analysis
        WPTCoupledSolver,
        compute_mutual_inductance,
        compute_coupling_coefficient,
        analyze_coil_coupling,
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

# RWG-EFIE solver for 3D surface element analysis
# The Python implementation has been migrated to C++ with OpenMP parallelization.
# Access via rad.RwgMeshRect(), rad.RwgMeshDisk(), rad.RwgMeshCylinder(),
# rad.RwgMeshSpiral(), rad.RwgMeshLoop(), rad.RwgSolverCreate(), etc.
# See docs/API_REFERENCE.md for usage.

# VTK Export: Use rad.FldVTS() (C++ implementation)
# See docs/API_REFERENCE.md for usage

# Beam Tracking: NumPy-based particle trajectory computation
# Import beam_tracking submodule for accelerator physics applications
try:
    from . import beam_tracking
    BEAM_TRACKING_AVAILABLE = True
except ImportError:
    # beam_tracking requires numpy
    BEAM_TRACKING_AVAILABLE = False

# Analysis Framework: Static, Frequency Response, Transient (CLN)
# Unified interface for electromagnetic analysis
try:
    from .analysis import (
        # Analysis types
        AnalysisType,
        SolverType,
        # Result classes
        AnalysisResult,
        StaticResult,
        FrequencyResult,
        TransientResult,
        # PEEC solver classes
        PEECAnalysisSolver,
        UnifiedAnalysis,
        # MMM result classes
        MMMStaticResult,
        MMMFrequencyResult,
        # MMM solver classes
        MMMAnalysisSolver,
        UnifiedMMMAnalysis,
        # MMM utility functions
        build_magnetic_circuit_from_mmm,
        # Convenience waveform generators
        step_voltage,
        pulse_voltage,
        sinusoidal_voltage,
        ramp_voltage,
    )
    ANALYSIS_AVAILABLE = True
except ImportError:
    # Analysis requires numpy
    ANALYSIS_AVAILABLE = False
