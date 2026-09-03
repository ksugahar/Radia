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

__version__ = "4.95.77"

# Compatibility window with the cubit-mesh-export package. The .ccm/.pyd
# plugin binaries built there must match this radia minor series. The
# 2026-04-14 stale Cubit-plugin incident is the canonical reason this
# matters. cubit-plugin-install enforces this at deploy time.
COMPAT_CUBIT_MESH_EXPORT_MIN = "0.5.0"
COMPAT_CUBIT_MESH_EXPORT_MAX = "0.999.999"  # bumped on next radia minor

# DLL loader for Windows
# MKL DLLs are installed via pip dependency (mkl>=2026,<2027)
# at {sys.prefix}/Library/bin/ (following NGSolve pattern)
import os
import sys

_package_dir = os.path.dirname(os.path.abspath(__file__))

if sys.platform == 'win32':
    # oneMKL otherwise selects Intel OpenMP, while PyTorch bundles a separate
    # libiomp5md.dll. Loading both aborts the process with OMP Error #15. TBB
    # is a supported threaded oneMKL backend and is installed by the mkl wheel.
    # Preserve an explicit application-level choice when one is already set.
    os.environ.setdefault('MKL_THREADING_LAYER', 'TBB')

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

    # 3. MKL 2026 DLLs from pip install mkl (mkl_rt.3.dll, etc.)
    _mkl_bin = os.path.join(sys.prefix, "Library", "bin")
    if os.path.isdir(_mkl_bin):
        _dirs_to_add.append(_mkl_bin)

    # 4. Explicit external MKL override. There is deliberately no machine-wide
    #    oneAPI fallback: the selected Python environment is the runtime owner.
    _mklroot_bin = os.path.join(os.environ.get("MKLROOT", ""), "bin")
    if os.path.isdir(_mklroot_bin):
        _dirs_to_add.append(_mklroot_bin)

    # Register directories with OS DLL loader
    for _d in _dirs_to_add:
        if hasattr(os, 'add_dll_directory'):
            os.add_dll_directory(_d)
        if _d not in os.environ.get('PATH', ''):
            os.environ['PATH'] = _d + os.pathsep + os.environ.get('PATH', '')

    del _dirs_to_add, _mkl_bin, _mklroot_bin

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
    ESIM_AVAILABLE = True
except ImportError:
    # ESIM requires scipy, which may not be installed
    ESIM_AVAILABLE = False

# NOTE: Old conductor API (CndLoop, CndRecBlock, CplMag*, Rwg*) removed (2026-02-13).
# Use the PEEC topology solver for conductor-only PEEC. Magnetic-material
# coupling is handled through the HDiv-VIM / reduced-FEM route.

# NOTE: FldVTS() and the old in-tree beam_tracking engine were removed
# (2026-03-22). Use NGSolve + GmshPostExport for visualization and the thin
# radia.xsuite_bridge adapter to CERN Xsuite for magnetic-field tracking.

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
# Demag backend: soft iron is HDiv-VIM only in Radia.
#
# DEFAULT = "auto":
#   - mesh-BACKED pure TET / HEX / WEDGE soft iron
#     (radia.vim.MeshSoftIron(mesh, mu_r=/bh_table=, order=1|2) + rad.Solve) -> HDiv-VIM (BDM1/BDM2).
#   - mesh-BACKED unsupported element mixes fail loud until HDiv coverage is added.
#   - mesh-LESS surface-charge soft iron is retired; build a mesh and use MeshSoftIron/SoftIron.
#   - permanent-magnet field objects and legacy non-soft-iron C++ operations are unchanged.
# set_demag_backend("hdiv") is accepted for explicitness.
# set_demag_backend("auto"/None) restores the split.
# ---------------------------------------------------------------------------

_demag_backend = None   # None/"auto" = default; "hdiv" = explicit HDiv-VIM


_BACKEND_MISSING = object()


def set_demag_backend(name):
    """Select the soft-iron demag backend.

    Accepted values are ``"hdiv"`` and ``"auto"``/``None``.
    """
    global _demag_backend
    if name in (None, "auto"):
        _demag_backend = None
    elif name == "hdiv":
        _demag_backend = name
    else:
        raise ValueError("demag_backend must be 'hdiv' or 'auto'/None (got %r)" % (name,))
    return _demag_backend or "auto"


def get_demag_backend():
    """The selected soft-iron demag backend: "hdiv" or "auto"."""
    return _demag_backend or "auto"


def _normalize_demag_backend(name):
    """Normalize a per-call demag backend without mutating the global default."""
    if name in (None, "auto"):
        return None
    if name == "hdiv":
        return name
    raise ValueError("demag_backend must be 'hdiv' or 'auto'/None (got %r)" % (name,))


if "ObjCnt" in globals():
    _cpp_ObjCnt = globals()["ObjCnt"]

    def ObjCnt(*args, **kwargs):   # noqa: F811  (record Python-built containers for safe HDiv lookup)
        """Create a Radia container and record its direct members for HDiv-VIM dispatch.

        The underlying C++ ObjCntStuf helper is not safe to probe on arbitrary non-container handles, so
        the Solve wrapper uses this Python-side record when deciding whether a container includes a
        MeshSoftIron body.
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
        """Radia relaxation solve with the HDiv-VIM soft-iron backend (see set_demag_backend):
          - mesh-BACKED pure TET / HEX / WEDGE soft iron (radia.vim.MeshSoftIron)
            -> FEEC HDiv-VIM (BDM1, default);
          - mesh-BACKED unsupported element mixes fail loud until HDiv support lands;
          - mesh-LESS surface-charge soft iron is retired and rejected by the C++ relaxation layer;
          - permanent magnets and legacy non-soft-iron operations stay on the C++ path.
        A per-call demag_backend=('hdiv'|'auto') overrides the global set_demag_backend choice."""
        backend_arg = kwargs.pop("demag_backend", _BACKEND_MISSING)
        backend = _demag_backend if backend_arg is _BACKEND_MISSING else _normalize_demag_backend(backend_arg)
        top = args[0] if args else None
        registered = False
        if top is not None:
            try:
                from radia.vim import _radsolve
                registered = _radsolve.is_registered(top)
            except Exception:
                registered = False
        if registered:
            from radia.vim import _radsolve
            if not _radsolve.is_hdiv_eligible(top):
                raise ValueError(
                    "rad.Solve(auto): every mesh-backed soft iron must be HDiv-VIM eligible "
                    "(pure tet/hex/wedge only).")
            return _radsolve.dispatch(*args, **kwargs)      # auto tet/hex/wedge / explicit hdiv -> FEEC HDiv-VIM
        if backend == "hdiv":
            raise ValueError(
                "demag_backend='hdiv' needs a mesh-backed soft iron built via "
                "radia.vim.MeshSoftIron(mesh, mu_r=/bh_table=); this body is mesh-less.  "
                "Build it via radia.vim.MeshSoftIron or radia.SoftIron.")
        if "method" not in kwargs and len(args) < 4:
            # The legacy C++ relaxation path retained only its symmetric LU
            # method.  Supply that default here as well as in the native ABI so
            # an older extension binary cannot reintroduce the retired method=1.
            kwargs["method"] = 0
        return _cpp_Solve(*args, **kwargs)                  # PM / legacy non-soft-iron C++ path


if "Fld" in globals():
    _cpp_Fld = globals()["Fld"]

    def Fld(obj, *args, **kwargs):   # noqa: F811  (HDiv BDM1/BDM2 field dispatch)
        """Evaluate Radia fields.

        Solved mesh-backed HDiv-VIM objects are evaluated from their full
        BDM1/BDM2 fields by persistent C++ charge-field kernels.  A multi-body
        container sums every registered body plus its ordinary Radia sources.
        Other Radia objects call the ordinary C++ ``Fld`` unchanged.

        The caller owns ``with ngsolve.TaskManager():`` when evaluating a
        solved mesh-backed HDiv-VIM object.
        """
        try:
            from radia.vim import _radsolve
            record = _radsolve.field_solution_for(obj)
        except Exception:
            record = None
        if record is not None:
            if kwargs or len(args) != 2:
                raise TypeError("rad.Fld(HDiv): expected Fld(obj, field_type, points)")
            import numpy as _np
            from radia.vim._field_batch import field_from_solution, magnetization_from_solution
            field_type, points = args
            field_type = str(field_type).lower()
            pts = _np.asarray(points, dtype=float)
            single = pts.ndim == 1
            pts2 = pts.reshape(-1, 3)
            results = record.get("results")
            if results is None:
                results = (record["result"],)
            h_iron = field_from_solution(results[0], pts2)
            for result in results[1:]:
                h_iron = h_iron + field_from_solution(result, pts2)
            source_obj = record.get("source_object")
            if field_type in ("h", "hx", "hy", "hz"):
                value = h_iron
                if source_obj is not None:
                    value = value + _np.asarray(_cpp_Fld(source_obj, "h", pts2), float).reshape(-1, 3)
            elif field_type in ("b", "bx", "by", "bz"):
                magnetization = magnetization_from_solution(results[0], pts2)
                for result in results[1:]:
                    magnetization = magnetization + magnetization_from_solution(result, pts2)
                value = (4.0e-7 * _np.pi) * (h_iron + magnetization)
                if source_obj is not None:
                    value = value + _np.asarray(_cpp_Fld(source_obj, "b", pts2), float).reshape(-1, 3)
            elif field_type == "m":
                value = magnetization_from_solution(results[0], pts2)
                for result in results[1:]:
                    value = value + magnetization_from_solution(result, pts2)
            else:
                raise NotImplementedError(
                    "rad.Fld on an HDiv-VIM solution supports b/h/m and Cartesian components; "
                    f"{field_type!r} has no BDM1/BDM2 field contract"
                )
            if len(field_type) == 2 and field_type[1] in "xyz":
                value = value[:, "xyz".index(field_type[1])]
            if single:
                return float(value[0]) if value.ndim == 1 else value[0]
            return value
        return _cpp_Fld(obj, *args, **kwargs)


if "SolverConfig" in globals():
    _cpp_SolverConfig = globals()["SolverConfig"]

    def SolverConfig(**kwargs):   # noqa: F811  (adds demag_backend on top of the C++ SolverConfig)
        """Unified solver config.  Adds the demag_backend selector ("hdiv" | "auto", see
        set_demag_backend); nonlinear LU-state options such as relax_param,
        newton_method, and keep_magnetization pass through to the C++ SolverConfig.
        HDiv, PEEC, and BEM compression settings belong to their solver APIs."""
        if "demag_backend" in kwargs:
            set_demag_backend(kwargs.pop("demag_backend"))
        if kwargs:
            _cpp_SolverConfig(**kwargs)


if "UtiDelAll" in globals():
    _cpp_UtiDelAll = globals()["UtiDelAll"]

    def UtiDelAll(*args, **kwargs):   # noqa: F811  (clears the HDiv mesh<->container registry too)
        """Delete all Radia objects.  Also clears the HDiv-VIM mesh<->container registry
        (radia.vim.MeshSoftIron), whose container handles are invalidated here."""
        try:
            from radia.vim import _radsolve
            _radsolve.clear_registry()
        except Exception:
            pass
        return _cpp_UtiDelAll(*args, **kwargs)


# Intent-based user-facing API (2-layer API; see CLAUDE.md "Reduce Proprietary API Surface").
# radia.SoftIron("yoke.vol", mu_r=) exposes the HDiv-VIM soft-iron path behind one
# .vol-driven object; the ObjHexahedron/... primitives become an internal representation detail.
from .soft_iron import SoftIron  # noqa: E402,F401
# radia.magnet_box(center, dimensions, magnetization) -- the ObjRecMag substitute on a fixed-M
# surface-charge ObjHexahedron; a permanent magnet, no Solve.
from .magnet import magnet_box  # noqa: E402,F401
# Script-side ObjRecMag: the SOLE definition now that the C++ surface-current ObjRecMag
# constructor is retired (un-exposed from the extension). It forwards to the fixed-M
# ObjHexahedron magnet_box so production / other solvers / examples / tests
# keep working unchanged.
from .magnet import ObjRecMag  # noqa: E402,F401
from .lamination import laminated_mu_eff  # noqa: E402,F401
