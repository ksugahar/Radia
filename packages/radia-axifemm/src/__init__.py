"""radia_axifemm — FEMM/Henrotte axisymmetric FE for NGSolve.

Adds custom finite-element classes whose shape functions are linear in
{1, r^2, z, r^2 z, ...} on physical (r, z) coordinates, reproducing the
FEMM/Henrotte 1993 axisymmetric formulation while staying inside the
NGSolve API. Bring your own mesh, geometry, materials, and solvers.

Quick usage:

    from ngsolve import *
    from radia_axifemm import H1Henrotte

    mesh = Mesh(...)
    fes = H1Henrotte(mesh, order=2)   # order=2 means {1, r^2, z}
    u, v = fes.TnT()
    a = BilinearForm(...)             # standard NGSolve weak forms
    ...

Reference:
    F. Henrotte et al., "A new method for axisymmetric linear and nonlinear
    problems," IEEE Transactions on Magnetics 9(2):1352-1355, March 1993.
    D. Meeker, FEMM 4.2 axisymmetric formulation notes.
"""

# IMPORTANT: import ngsolve before the C++ module so that NGSolve's
# shared libraries are loaded into the process before ours.
import ngsolve  # noqa: F401

from .radia_axifemm import version, hello  # noqa: F401

# Phase 2-B/C/E exports.
try:
    from .radia_axifemm import (  # noqa: F401
        H1Henrotte,
        AxiHenrotteFESpace,
        AxiHenrotteStiffnessBFI,
        AxiHenrotteSigmaMassBFI,
        AxiHenrotteFE_Q1_AxisAligned,
        AxiHenrotteFE_Q2_AxisAligned,
        AxiHenrotteFE_Q3_AxisAligned,
        AxiHenrotteFE_P1_Triangle,
    )
except ImportError:
    # Phase 2-A bootstrap — FESpace/Integrators not yet built.
    pass

__version__ = "0.1.0"
