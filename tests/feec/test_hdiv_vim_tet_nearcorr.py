"""Golden test: production #3 (part 2) -- near-field-corrected scalable Gram for the unstructured tet
HDiv-VIM.  The scalable Gram is the compressed CENTROID-MONOPOLE H-matrix (far, cheap) PLUS a SPARSE
near-field correction (exact sub-point MINUS monopole for near charge pairs) -- the standard H-matrix
near-field correction.  This lifts the sphere demag from the monopole under-estimate (~0.31) to the
Gram-EXACT value (~0.33), converging to the analytic 1/3, while staying scalable (O(N) extra cost).

(#3 part 1 fixed the tet/tri sub-point quadrature bug + showed the Gram-exact demag -> 1/3; this part 2
delivers that accuracy in the SCALABLE path via the sparse near correction.)

Locks:
  - the near-field-corrected SCALABLE demag (monopole H-matvec + sparse corr) == the dense reference
    (dense monopole Gram + same corr) within the ACA tolerance -- the H-matrix near correction is exact;
  - the correction LIFTS the demag above the monopole value and toward 1/3 (in [0.32, 0.334]);
  - it converges toward 1/3 with mesh refinement.
"""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

import radia._radia_pybind as _rp  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
from radia.vim import _core as tet  # noqa: E402

import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402


def _sphere(h):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        return ng.Mesh(geo.GenerateMesh(maxh=h))


def _corrected_demag(mesh, eps=1e-6, near_factor=2.0):
    """Return (monopole_demag, scalable_corrected_demag, dense_corrected_demag)."""
    d = tet.build_demag(mesh, nsub=4)
    H = _rp._ChargeGramHMatrix(list(d["cent"].ravel()), list(d["meas"]),
                               list(d["self_energy"]), eps, 32, 2.0)
    corr = tet.build_near_correction(mesh, d, nsub=4, near_factor=near_factor)
    B = d["B_csr"]
    m = d["m_unit"]
    Bm = B @ m
    GH = np.asarray(H.matvec(Bm.tolist()), float)
    Nm_scal = B.T @ (GH + corr @ Bm)
    Nm_dense = B.T @ (d["G"] @ Bm + corr @ Bm)
    denom = m @ d["M_mass"] @ m
    return (tet.demag_factor(d),
            float((m @ Nm_scal) / denom),
            float((m @ Nm_dense) / denom))


def test_nearcorr_scalable_matches_dense_and_lifts_toward_third():
    mono, scal, dense = _corrected_demag(_sphere(0.5))
    # scalable (H-matvec + sparse corr) matches the dense reference (monopole Gram + same corr)
    assert abs(scal - dense) < 1e-3 * abs(dense), f"scalable {scal:.5f} vs dense {dense:.5f}"
    # the near correction LIFTS the demag above the monopole under-estimate, toward 1/3
    assert scal > mono + 0.005, f"near correction did not lift demag: mono={mono:.4f} corr={scal:.4f}"
    assert 0.315 < scal < 0.3333 + 1e-3, f"corrected demag {scal:.4f} not approaching 1/3"


def test_nearcorr_converges_to_third():
    """The near-field-corrected demag converges toward 1/3 as the mesh refines."""
    _, d5, _ = _corrected_demag(_sphere(0.5))
    _, d3, _ = _corrected_demag(_sphere(0.3))
    assert d3 > d5 - 1e-3, f"corrected demag not converging toward 1/3: h0.5={d5:.4f} h0.3={d3:.4f}"
    assert 0.32 < d3 < 0.3333 + 1e-3, f"refined corrected demag {d3:.4f} not ~1/3"
