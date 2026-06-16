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

__version__ = "4.90.2"

# Compatibility window with the cubit-mesh-export package. The .ccm/.ccl
# plugin binaries built there must match this radia minor series; the
# 2026-04-14 incident (stale .ccl on 100号機) is the canonical reason
# this matters. cubit-plugin-install enforces this at deploy time.
COMPAT_CUBIT_MESH_EXPORT_MIN = "0.5.0"
COMPAT_CUBIT_MESH_EXPORT_MAX = "0.999.999"  # bumped on next radia minor

# DLL loader for Windows
# MKL DLLs are installed via pip dependency (mkl>=2024.2.0)
# at {sys.prefix}/Library/bin/ (following NGSolve pattern)
import os
import sys

_package_dir = os.path.dirname(os.path.abspath(__file__))

if sys.platform == 'win32':
    _dirs_to_add = []

    # 1. Package directory (for _radia_pybind.pyd and other .pyd files)
    _dirs_to_add.append(_package_dir)

    # 2. NGSolve DLLs (ngcore.dll, libngsolve.dll, etc.)
    #    _radia_pybind links against ngstd, so NGSolve DLLs must be loadable.
    #    Detect via ngsolve Python package location → install_root/bin/
    try:
        import ngsolve as _ngsolve_mod
        _ngsolve_pkg_dir = os.path.dirname(_ngsolve_mod.__file__)
        # pip-installed: site-packages/ngsolve/ → bin/ is at ../../Library/bin/
        # source-installed: .../Lib/site-packages/ngsolve/ → bin/ is at ../../../bin/
        for _up in [
            os.path.join(_ngsolve_pkg_dir, '..', '..', '..', 'bin'),     # source install
            os.path.join(_ngsolve_pkg_dir, '..', '..', 'Library', 'bin'),# pip install
        ]:
            _ng_bin = os.path.normpath(_up)
            if os.path.isdir(_ng_bin) and os.path.isfile(os.path.join(_ng_bin, 'ngcore.dll')):
                _dirs_to_add.append(_ng_bin)
                break
        del _ngsolve_mod, _ngsolve_pkg_dir
    except ImportError:
        pass

    # 3. MKL DLLs from pip install mkl (mkl_rt.2.dll, etc.)
    _mkl_bin = os.path.join(sys.prefix, "Library", "bin")
    if os.path.isdir(_mkl_bin):
        _dirs_to_add.append(_mkl_bin)

    # 4. Intel oneAPI (fallback for development builds)
    for _intel_path in [
        os.path.join(os.environ.get("MKLROOT", ""), "bin"),
        r"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin",
        r"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin",
    ]:
        if os.path.isdir(_intel_path):
            _dirs_to_add.append(_intel_path)

    # Register directories with OS DLL loader
    for _d in _dirs_to_add:
        if hasattr(os, 'add_dll_directory'):
            os.add_dll_directory(_d)
        if _d not in os.environ.get('PATH', ''):
            os.environ['PATH'] = _d + os.pathsep + os.environ.get('PATH', '')

    del _dirs_to_add, _mkl_bin

# High-order mesh curving is handled by the Cubit C++ plugin (ACIS kernel).
# Netgen fork (SetGeomInfo) is no longer required.

# Import all symbols from the pybind11 C++ extension module (_radia_pybind.pyd)
try:
    from ._radia_pybind import *
except ImportError as e:
    raise ImportError(
        "Failed to import radia pybind11 module (_radia_pybind.pyd). "
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
        # WPT (Wireless Power Transfer) analysis
        WPTCoupledSolver,
        compute_mutual_inductance,
        compute_coupling_coefficient,
        analyze_coil_coupling,
    )
    # VTK export is not maintained by Radia — use NGSolve VTKOutput instead.
    # esim_vtk_export.py is kept for backwards compatibility but not re-exported.
    pass
    ESIM_AVAILABLE = True
except ImportError:
    # ESIM requires scipy, which may not be installed
    ESIM_AVAILABLE = False

# NOTE: Old conductor API (CndLoop, CndRecBlock, CplMag*, Rwg*) removed (2026-02-13).
# Use PEEC topology solver (peec_topology.py) and coupled solver (peec_coupled.py).

# NOTE: FldVTS() and beam_tracking removed (2026-03-22).
# Use NGSolve + GmshPostExport for visualization, CERN Xsuite for tracking.

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

# Post-hoc Kelvin Periodic identification on an NGSolve mesh (no Cubit /
# OCC needed).  See `radia.kelvin_identify_ngsolve` for full docs.
try:
    from .kelvin_identify_ngsolve import (
        add_kelvin_identification,
        detect_kelvin_offset,
        has_kelvin_identification,
    )
    KELVIN_IDENTIFY_AVAILABLE = True
except ImportError:
    # Requires ngsolve + scipy (greedy fallback) -- both should be
    # present in any Radia install, so this is just defensive.
    KELVIN_IDENTIFY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Demag backend selection (runtime, API-OPTIONAL) -- the yano-type MSC backend is DEPRECATED.
#
# The Yano MSC demag backend (rad.Solve on ObjHexahedron / ObjWedge soft iron) is being
# superseded by the FEEC HDiv-VIM (radia.hdiv_vim).  It is exposed here as an OPTIONAL, selectable
# backend so it can be removed cleanly later: callers that pin "yano" keep today's behaviour (with a
# one-time DeprecationWarning), and when the HDiv-VIM is wired into rad.Solve the default flips to
# "hdiv" -- deleting yano then changes only this default, not the public API surface.  Permanent-magnet
# and MMM (tetrahedron) solves do not use the yano-type MSC path and are unaffected.
# ---------------------------------------------------------------------------
import warnings as _warnings

_DEMAG_BACKEND = "yano"            # legacy default during migration; flips to "hdiv" when wired in
_yano_deprecation_emitted = False


def set_demag_backend(name):
    """Select the demag backend for rad.Solve: "yano" (the legacy Yano MSC for hex/wedge soft
    iron -- DEPRECATED, the only backend rad.Solve currently dispatches) or "hdiv" (the FEEC HDiv-VIM in
    radia.hdiv_vim, not yet wired into rad.Solve).  Returns the previous value."""
    global _DEMAG_BACKEND
    if name not in ("yano", "hdiv"):
        raise ValueError("demag_backend must be 'yano' or 'hdiv', got %r" % (name,))
    prev, _DEMAG_BACKEND = _DEMAG_BACKEND, name
    return prev


def get_demag_backend():
    """Return the current demag backend ("yano" legacy / "hdiv" = radia.hdiv_vim)."""
    return _DEMAG_BACKEND


if "Solve" in globals():
    _cpp_Solve = globals()["Solve"]

    def Solve(*args, **kwargs):   # noqa: F811  (intentional thin wrapper over the C++ Solve)
        """Radia relaxation solve.  Thin wrapper over the C++ solver adding the demag-backend gate
        (see set_demag_backend): the yano-type MSC backend is DEPRECATED -- migrate hex/wedge soft-iron
        demag to the FEEC HDiv-VIM (radia.hdiv_vim).  Permanent-magnet / MMM-tet solves are unaffected."""
        global _yano_deprecation_emitted
        if _DEMAG_BACKEND == "hdiv":
            # Dispatch the FEEC HDiv-VIM backend.  Requires the iron body to have been built via
            # radia.vim.soft_iron_from_mesh (which registers the NGSolve mesh + material); the bridge
            # solves with hdiv_demag_solve and writes per-element M back so rad.Fld/rad.ObjM reflect it.
            from radia.vim import _radsolve
            return _radsolve.dispatch(*args, **kwargs)
        if not _yano_deprecation_emitted:
            _yano_deprecation_emitted = True
            _warnings.warn(
                "The yano-type MSC demag backend (rad.Solve on ObjHexahedron / ObjWedge soft iron) is "
                "DEPRECATED and will be removed; it is superseded by the FEEC HDiv-VIM (radia.hdiv_vim).  "
                "Permanent-magnet and MMM (tetrahedron) solves are unaffected.  This notice fires once per "
                "session; call radia.set_demag_backend('hdiv') to forbid the yano path.",
                DeprecationWarning, stacklevel=2)
        return _cpp_Solve(*args, **kwargs)


if "SolverConfig" in globals():
    _cpp_SolverConfig = globals()["SolverConfig"]

    def SolverConfig(**kwargs):   # noqa: F811  (adds demag_backend on top of the C++ SolverConfig)
        """Unified solver config.  Adds the demag_backend selector ("yano" | "hdiv", see
        set_demag_backend); all other kwargs (hacapk_eps, bicgstab_tol, relax_param, newton_method, ...)
        pass through to the C++ SolverConfig."""
        if "demag_backend" in kwargs:
            set_demag_backend(kwargs.pop("demag_backend"))
        if kwargs:
            _cpp_SolverConfig(**kwargs)


if "UtiDelAll" in globals():
    _cpp_UtiDelAll = globals()["UtiDelAll"]

    def UtiDelAll(*args, **kwargs):   # noqa: F811  (clears the HDiv mesh<->container registry too)
        """Delete all Radia objects.  Also clears the HDiv-VIM mesh<->container registry
        (radia.vim.soft_iron_from_mesh), whose container handles are invalidated here."""
        try:
            from radia.vim import _radsolve
            _radsolve.clear_registry()
        except Exception:
            pass
        return _cpp_UtiDelAll(*args, **kwargs)
