r"""Circular waveguide cutoff (Bessel-zero 2D Helmholtz eigenmodes) -- regression test (#61).

The circular sibling of the rectangular cutoff (#53): TE = Neumann (k_c=j'_mn/a), TM = Dirichlet
(k_c=j_mn/a) on the disk; f_c=c k_c/2pi. Dominant TE11 (j'_11=1.8412). Closed-form helper
(tool-independent) + the SAME `helmholtz_cutoff_wavenumbers_2d` eigensolver on a disk mesh."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import circular_waveguide_cutoff, cutoff_frequency, C0

A = 0.0127


def validate_circular_cutoff_closed_form():
    # dominant mode TE11, j'_11 = 1.8412 -> f_c = c*1.8412/(2 pi a)
    f_te11 = circular_waveguide_cutoff(A, "TE", 1, 1)
    assert math.isclose(f_te11, C0 * 1.8412 / (2 * math.pi * A), rel_tol=1e-3)
    # mode ordering: TE11 < TM01 < TE21 < TE01=TM11
    f_tm01 = circular_waveguide_cutoff(A, "TM", 0, 1)   # j_01 = 2.4048
    f_te21 = circular_waveguide_cutoff(A, "TE", 2, 1)   # j'_21 = 3.0542
    assert f_te11 < f_tm01 < f_te21
    # the famous degeneracy TE0n / TM1n (j'_0n = j_1n)
    assert math.isclose(circular_waveguide_cutoff(A, "TE", 0, 1),
                        circular_waveguide_cutoff(A, "TM", 1, 1), rel_tol=1e-6)
    # cutoff scales as 1/radius
    assert math.isclose(circular_waveguide_cutoff(2 * A, "TE", 1, 1), f_te11 / 2, rel_tol=1e-12)
    with pytest.raises(ValueError):
        circular_waveguide_cutoff(A, "XX", 1, 1)


def validate_circular_waveguide_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane
    from radia_mcp.radia_ngsolve.waveguide import helmholtz_cutoff_wavenumbers_2d

    disk = WorkPlane().Circle(0, 0, A).Face(); disk.edges.name = "wall"
    mesh = Mesh(OCCGeometry(disk, dim=2).GenerateMesh(maxh=A / 18))

    with TaskManager():
        te = helmholtz_cutoff_wavenumbers_2d(mesh, 3, bc="neumann")     # TE11(x2), TE21
        tm = helmholtz_cutoff_wavenumbers_2d(mesh, 1, bc="dirichlet")   # TM01

    # dominant TE11 matches the Bessel-zero cutoff
    assert math.isclose(cutoff_frequency(te[0]), circular_waveguide_cutoff(A, "TE", 1, 1), rel_tol=3e-3)
    # TE21 (3rd TE wavenumber, after the TE11 degenerate pair)
    assert math.isclose(cutoff_frequency(te[2]), circular_waveguide_cutoff(A, "TE", 2, 1), rel_tol=3e-3)
    # TM01 (Dirichlet) matches
    assert math.isclose(cutoff_frequency(tm[0]), circular_waveguide_cutoff(A, "TM", 0, 1), rel_tol=3e-3)
    # TE11 appears as a degenerate pair (two near-equal lowest Neumann modes)
    assert math.isclose(te[0], te[1], rel_tol=5e-3)


if __name__ == "__main__":
    validate_circular_cutoff_closed_form()
    validate_circular_waveguide_fe()
    print("[OK] circular waveguide: disk Helmholtz eigenvalues == Bessel-zero cutoffs (TE11 dominant).")
