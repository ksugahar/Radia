"""radia-mmm -- the LEGACY collocation MMM/MSC demag solver namespace.

Radia's demag core is migrating to the FEEC **HDiv-VIM** (``radia.hdiv_vim``), which is the DEFAULT
path for new work.  The collocation **Magnetic Moment Method** (MMM, tetrahedra) and **Magnetic Surface
Charge** (MSC, hexahedra/wedges -- including the Yano-Sugahara distortion elements) are the LEGACY demag
path, homed in this package so they can be maintained / deprecated / removed independently of the new
core.

The MMM/MSC C++ kernels remain in the ``radia`` wheel; this package re-exports the MMM/MSC-relevant API
under a clear legacy namespace, and is the home of the MMM/MSC **examples** (``packages/radia-mmm/
examples``) and **tests**.

For NEW work prefer ``radia.hdiv_vim`` (``build_demag`` / ``DemagOperator`` / ``solve_demag_newton``);
the Yano-Sugahara MSC backend is deprecated (see ``radia.set_demag_backend``).

    import radia_mmm as rm
    pm = rm.ObjHexahedron(vertices, [0, 0, 954930])   # the same C++ element as radia.ObjHexahedron
    B  = rm.Fld(pm, "b", [0, 0, 0.1])
"""
import radia as _radia

# Re-export the MMM/MSC-relevant Radia API under the legacy namespace (the kernels live in `radia`).
# getattr-guarded so a renamed/absent symbol never breaks `import radia_mmm`.
_REEXPORT = [
    "ObjRecMag", "ObjHexahedron", "ObjTetrahedron", "ObjWedge", "ObjCnt", "ObjBckg",
    "MatLin", "MatSatIsoTab", "MatStd", "MatApl",
    "Solve", "Fld", "UtiDelAll", "SolverConfig", "GetSolverConfig",
    "set_demag_backend", "get_demag_backend",
]
__all__ = []
for _n in _REEXPORT:
    _obj = getattr(_radia, _n, None)
    if _obj is not None:
        globals()[_n] = _obj
        __all__.append(_n)

# The Hantila polarization MMM solve (constant LHS, factored once) is the C++ SolverConfig path
# rm.SolverConfig(b_input_hantila=True, hantila_alpha=..., hantila_relax=...) + rm.Solve (see
# examples/hantila_solver/) -- it is reachable through the re-exported SolverConfig/Solve above, not a
# separate Python entry point.

__version__ = "0.1.0"
