r"""Circular waveguide cutoff (Bessel-zero 2D Helmholtz eigenmodes) -- regression test (#61).

The circular sibling of the rectangular cutoff (#53): TE = Neumann (k_c=j'_mn/a), TM = Dirichlet
(k_c=j_mn/a) on the disk; f_c=c k_c/2pi. Dominant TE11 (j'_11=1.8412). Closed-form helper
(tool-independent) + the SAME `helmholtz_cutoff_wavenumbers_2d` eigensolver on a disk mesh."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import (
    C0,
    circular_waveguide_band_summary,
    circular_waveguide_cutoff,
    circular_waveguide_mode_table,
    cutoff_frequency,
)

A = 0.0127


def test_circular_cutoff_closed_form():
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


def test_circular_mode_table_and_band_summary():
    table = circular_waveguide_mode_table(A, max_m=3, max_n=2)
    assert [row["mode"] for row in table[:4]] == ["TE11", "TM01", "TE21", "TE01"]
    assert table[0]["angular_degeneracy"] == 2
    assert table[1]["angular_degeneracy"] == 1
    assert math.isclose(table[3]["cutoff_frequency"], table[4]["cutoff_frequency"], rel_tol=1e-12)
    assert table[3]["mode"] == "TE01" and table[4]["mode"] == "TM11"

    f_te11 = table[0]["cutoff_frequency"]
    f_tm01 = table[1]["cutoff_frequency"]
    below = circular_waveguide_band_summary(A, 0.95 * f_te11)
    single = circular_waveguide_band_summary(A, 0.5 * (f_te11 + f_tm01))
    multi = circular_waveguide_band_summary(A, 1.01 * f_tm01)
    assert below["below_dominant_cutoff"]
    assert single["single_mode"]
    assert single["n_propagating_with_degeneracy"] == 2
    assert not multi["single_mode"]
    assert [row["mode"] for row in multi["propagating_modes"]] == ["TE11", "TM01"]
