"""Time-averaged Maxwell-stress surface force for a COMPLEX time-harmonic field
(force.maxwell_surface_force_harmonic), validated by REDUCTION to the
COMSOL-cross-validated static force.maxwell_surface_force.

The time-average of the quadratic Maxwell stress over a closed surface obeys two
exact identities (n real, B = a + i b phasor):

  (1) real field:   <F> = 0.5 * static(a)                 [the cos^2 time-avg 1/2]
  (2) complex field:<F> = 0.5 * (static(a) + static(b))   [cross-terms vanish]

Both are checked here on a unit-ball surface with analytic divergence-free fields
-- no solve, pure formula consistency against the validated static integrator.
This is the force handle for the SIBC workpiece in calc_fem_kelvin (a HOLE, so
the only force is the Maxwell stress on its "sibc" boundary).
"""
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, x, y, z, TaskManager
from netgen.occ import Box, Pnt, OCCGeometry, X
from radia_mcp.radia_ngsolve.force import (maxwell_surface_force,
                                           maxwell_surface_force_harmonic)


def build_box():
    # A single OPEN face (the x=1 plane).  A closed surface around a source-free
    # div/curl-free field would give zero NET force (Gauss); one face makes the
    # net integral non-trivial so the reduction check is meaningful.  The
    # harmonic=0.5*static identity is pointwise, so it holds on any surface.
    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    box.faces.Max(X).name = "surf"
    return Mesh(OCCGeometry(box).GenerateMesh(maxh=0.12))


def main():
    mesh = build_box()

    # Two distinct divergence-free analytic fields (valid magnetostatic B):
    #   a = (x, y, -2z)   div = 1+1-2 = 0
    #   b = (2x, -y, -z)  div = 2-1-1 = 0
    a = CoefficientFunction((x, y, -2 * z))
    b = CoefficientFunction((2 * x, -y, -z))

    with TaskManager():
        Fa = maxwell_surface_force(a, mesh, "surf")
        Fb = maxwell_surface_force(b, mesh, "surf")

        # (1) purely-real phasor a -> 0.5 * static(a)
        Fh_real = maxwell_surface_force_harmonic(a, mesh, "surf")

        # (2) complex phasor B = a + i b -> 0.5 * (static(a) + static(b))
        B = CoefficientFunction((a[0] + 1j * b[0],
                                 a[1] + 1j * b[1],
                                 a[2] + 1j * b[2]))
        Fh_cplx = maxwell_surface_force_harmonic(B, mesh, "surf")

    scale = max(abs(Fa[k]) for k in range(3))
    print(f"  static(a)        = ({Fa[0]:+.4e}, {Fa[1]:+.4e}, {Fa[2]:+.4e}) N")
    print(f"  static(b)        = ({Fb[0]:+.4e}, {Fb[1]:+.4e}, {Fb[2]:+.4e}) N")
    print(f"  harmonic(a)      = ({Fh_real[0]:+.4e}, {Fh_real[1]:+.4e}, {Fh_real[2]:+.4e}) N")
    print(f"  harmonic(a+i b)  = ({Fh_cplx[0]:+.4e}, {Fh_cplx[1]:+.4e}, {Fh_cplx[2]:+.4e}) N")
    assert scale > 1e-6, "test field gave ~zero net force; check geometry/naming"

    # (1) reduction to static (the time-average 1/2)
    for k in range(3):
        err = abs(Fh_real[k] - 0.5 * Fa[k]) / scale
        assert err < 1e-6, \
            f"real-phasor reduction failed k={k}: {Fh_real[k]:.6e} vs {0.5*Fa[k]:.6e} (rel {err:.2e})"

    # (2) complex decomposition: cross-terms vanish in the time-average
    for k in range(3):
        expect = 0.5 * (Fa[k] + Fb[k])
        err = abs(Fh_cplx[k] - expect) / scale
        assert err < 1e-6, \
            f"complex decomposition failed k={k}: {Fh_cplx[k]:.6e} vs {expect:.6e} (rel {err:.2e})"

    print("\n[OK] maxwell_surface_force_harmonic validated by reduction to the "
          "COMSOL-cross-validated static force:\n"
          "     real phasor   -> 0.5*static(a)              (cos^2 time-average)\n"
          "     complex a+i b -> 0.5*(static(a)+static(b))  (oscillating cross-terms drop).\n"
          "     This is the SIBC-workpiece eddy-force handle for calc_fem_kelvin's "
          "'sibc' boundary.")


if __name__ == "__main__":
    main()
