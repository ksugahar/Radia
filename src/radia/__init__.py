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

__version__ = "4.92.0"

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
# Demag backend: BOTH the collocation surface-charge (yano-type MSC) AND the FEEC HDiv-VIM
# (radia.vim) are kept (decision 2026-06-19).  They are complementary -- yano-MSC (6 sigma DOF/hex,
# higher per-element order) reaches accuracy with fewer hex elements; HDiv-VIM is loop-free
# (mesh/mu_r-independent iteration count).
#
# DEFAULT = "auto" (API-split): the API you use selects the method.
#   - mesh-LESS soft iron (ObjHexahedron/ObjWedge + MatLin/MatSatIsoTab + rad.Solve) -> yano-MSC (C++).
#   - mesh-BACKED soft iron (radia.vim.soft_iron_from_mesh(mesh, mu_r=/bh_table=) + rad.Solve) -> HDiv-VIM.
#   - tetrahedron (MMM) and permanent-magnet solves -> C++ solver (unchanged).
# set_demag_backend("yano"|"hdiv") OVERRIDES the auto split (a soft_iron_from_mesh container carries both
# representations, so either backend can solve it); set_demag_backend("auto"/None) restores the split.
# ---------------------------------------------------------------------------

_demag_backend = None   # None/"auto" = API-split default; "yano" or "hdiv" = forced override


def set_demag_backend(name):
    """Select the soft-iron demag backend.  "yano" = collocation surface-charge MSC; "hdiv" = FEEC
    HDiv-VIM; "auto"/None = API-split default (mesh-less -> yano, soft_iron_from_mesh -> HDiv).
    The choice is consulted by rad.Solve.  Returns the effective backend string."""
    global _demag_backend
    if name in (None, "auto"):
        _demag_backend = None
    elif name in ("yano", "hdiv"):
        _demag_backend = name
    else:
        raise ValueError("demag_backend must be 'yano', 'hdiv', or 'auto'/None (got %r)" % (name,))
    return _demag_backend or "auto"


def get_demag_backend():
    """The selected soft-iron demag backend: "yano", "hdiv", or "auto" (the API-split default)."""
    return _demag_backend or "auto"


if "ObjCnt" in globals():
    _cpp_ObjCnt = globals()["ObjCnt"]

    def ObjCnt(*args, **kwargs):   # noqa: F811  (record Python-built containers for safe HDiv lookup)
        """Create a Radia container and record its direct members for HDiv-VIM dispatch.

        The underlying C++ ObjCntStuf helper is not safe to probe on arbitrary non-container handles, so
        the Solve wrapper uses this Python-side record when deciding whether a container includes a
        soft_iron_from_mesh body.
        """
        h = _cpp_ObjCnt(*args, **kwargs)
        members = args[0] if args else kwargs.get("objs", None)
        if members is not None:
            try:
                from radia.vim import _radsolve
                _radsolve.register_container(h, list(members))
            except Exception:
                pass
        return h


if "Solve" in globals():
    _cpp_Solve = globals()["Solve"]

    def Solve(*args, **kwargs):   # noqa: F811  (thin wrapper: pick the soft-iron demag backend)
        """Radia relaxation solve with the API-split demag backend (see set_demag_backend):
          - mesh-BACKED soft iron (radia.vim.soft_iron_from_mesh) -> FEEC HDiv-VIM (default), or
            yano-MSC if demag_backend='yano';
          - mesh-LESS hex/wedge soft iron (ObjHexahedron/ObjWedge + MatLin) -> yano-type MSC (C++);
          - tetrahedron (MMM) and permanent magnets -> C++ solver.
        A per-call demag_backend=('yano'|'hdiv'|'auto') overrides the global set_demag_backend choice."""
        backend = kwargs.pop("demag_backend", None) or _demag_backend   # per-call > global > auto
        top = args[0] if args else None
        registered = False
        if top is not None:
            try:
                from radia.vim import _radsolve
                registered = _radsolve.is_registered(top)
            except Exception:
                registered = False
        if registered:
            if backend == "yano":
                return _cpp_Solve(*args, **kwargs)          # yano-MSC on the mesh-built elements
            from radia.vim import _radsolve                 # auto/hdiv -> FEEC HDiv-VIM
            return _radsolve.dispatch(*args, **kwargs)
        if backend == "hdiv":
            raise ValueError(
                "demag_backend='hdiv' needs a mesh-backed soft iron built via "
                "radia.vim.soft_iron_from_mesh(mesh, mu_r=/bh_table=); this body is mesh-less.  "
                "Build it via soft_iron_from_mesh, or use the yano backend (the default for mesh-less iron).")
        return _cpp_Solve(*args, **kwargs)                  # mesh-less -> yano-MSC (or MMM/PM)


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


# Intent-based user-facing API (2-layer API; see CLAUDE.md "Reduce Proprietary API Surface").
# radia.SoftIron("yoke.vol", mu_r=) unifies the yano-MSC and HDiv-VIM soft-iron paths behind one
# .vol-driven object; the ObjHexahedron/... primitives become an internal representation detail.
from .soft_iron import SoftIron  # noqa: E402,F401
