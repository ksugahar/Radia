"""Cyclic (N-fold) image reduction: ONE pole reproduces the whole RING.

The physical claim.  A rotationally symmetric machine (motor pole array, magnet ring) is solved on ONE
sector; the charge Gram folds in the other N-1 poles as ROTATED images about +z.  The reduced solve must
land on the FULL RING magnetization -- not on the lone-sector value, which over-magnetizes because it
misses the neighbours' demagnetizing influence.  Closing that gap is the entire job of the reduction.

Why ROTATIONAL and not translational: a finite N-fold rotational array is an exact FINITE image sum,
while an infinite translational array is a conditionally convergent dipole lattice sum whose value
depends on the summation shape -- that shape dependence IS the demagnetizing-factor phenomenon.  So the
cyclic reduction is unconditionally well posed where a translational one would not be.

The kernel-level anchor (a pi rotation reproduces the golden mask-3 mirror to round-off) lives in the
fast lane, tests/feec/test_hdiv_vim_cyclic_image.py.  This lane is the physical end-to-end check.

Gate: the reduced model must close at least 98 % of the lone-vs-full gap.  It is deliberately a gap
FRACTION rather than an absolute tolerance because the standalone pole and the ring's pole 0 are meshed
independently (netgen is not rotation-equivariant), so a small discretization difference is expected and
is not what this lane is testing.

MEASURED 2026-08-11 (LAB, N=4, maxh 7 mm, ring 522 / sector 112 elements): <Mz> full ring 3.2157743e+05,
lone sector 3.6013978e+05 (+11.99 %), cyclic reduced 3.2157792e+05 (+0.0002 %) -- 100.00 % of the gap
closed.  A second, independent confirmation on Sculpt hex meshes (6-fold ring) is recorded in
memory/cyclic_image_reduction_design.md.

MEASURED 2026-08-26 (LAB, N=4 alternating signs, same 522 / 112 element meshes): <Mz> full ring pole 0
3.9718599e+05, lone sector 3.6013978e+05 (-9.3272 %), cyclic alternating reduced 3.9718954e+05
(+0.000895 %) -- 99.9904 % of the gap closed.
"""
import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import radia.vim as vim  # noqa: E402
import ngsolve as ng  # noqa: E402

pytestmark = pytest.mark.compute_host

N_FOLD = 4
A = 0.010                 # pole half-size
R = 0.026                 # pole centre radius: > 2A, so the rotated poles stay DISJOINT
                          # (at R = 0.016 adjacent boxes overlap by 4 % of the volume and the
                          # "ring" is no longer 4 rotated copies of the sector -- a geometry
                          # artifact that showed up as a 4 % residual, not a solver error)
MU_R = 1000.0
H0 = 1.0e5
MAXH = 0.007
CLOSE_FRACTION = 0.98     # MEASURED 2026-08-11: the reduction closes 100.00 % of the gap
                          # (<Mz> full 3.2157743e5, lone 3.6013978e5 = +11.99 %, cyclic
                          # 3.2157792e5 = +0.0002 %).  Gate at 98 % to leave meshing headroom.


def _pole_mesh(sector_indices):
    """A mesh of the listed poles, each a cube of side 2A centred at radius R, rotated k*2pi/N."""
    from netgen.occ import Box, OCCGeometry, Pnt, gp_Ax1, gp_Dir, Glue
    boxes = [Box(Pnt(R - A, -A, -A), Pnt(R + A, A, A))
             .Rotate(gp_Ax1(Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 360.0 * k / N_FOLD)
             for k in sector_indices]
    # Glue, not Fuse (`+`): the poles are disjoint bodies, and the fuse drops the material name.
    shape = Glue(boxes)
    shape.mat("pole")
    return ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=MAXH))


def _alternating_pole_mesh(sector_indices):
    """A ring whose pole materials retain their indices for an alternating drive."""
    from netgen.occ import Box, OCCGeometry, Pnt, gp_Ax1, gp_Dir, Glue
    boxes = []
    for k in sector_indices:
        box = (Box(Pnt(R - A, -A, -A), Pnt(R + A, A, A))
               .Rotate(gp_Ax1(Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 360.0 * k / N_FOLD))
        box.mat(f"pole{k}")
        boxes.append(box)
    return ng.Mesh(OCCGeometry(Glue(boxes)).GenerateMesh(maxh=MAXH))


def _mean_mz(mesh, *, H_ext=None, material="pole", **kw):
    rad.UtiDelAll()
    if H_ext is None:
        H_ext = ng.CF((0.0, 0.0, H0))
    result = vim.Solve(mesh, mu_r=MU_R, H_ext=H_ext, order=1, tol=1e-9, **kw)
    region = mesh.Materials(material)
    volume = float(ng.Integrate(ng.CF(1.0), mesh, definedon=region))
    return float(ng.Integrate(result["gfM"][2], mesh, definedon=region)) / volume


def test_cyclic_sector_reproduces_the_full_ring():
    with ng.TaskManager():
        ring = _pole_mesh(range(N_FOLD))
        sector = _pole_mesh([0])
        mz_full = _mean_mz(ring)                       # the truth: all N poles solved explicitly
        mz_lone = _mean_mz(sector)                     # no images: over-magnetized
        mz_cyclic = _mean_mz(sector, image_cyclic=N_FOLD)   # the reduced model

    gap = abs(mz_lone - mz_full)
    residual = abs(mz_cyclic - mz_full)
    assert gap / abs(mz_full) > 1e-3, (
        "the poles barely interact at this radius, so the lane cannot discriminate "
        f"(<Mz> full {mz_full:.6e}, lone {mz_lone:.6e})")
    assert residual < (1.0 - CLOSE_FRACTION) * gap, (
        f"cyclic reduction did not close the gap: <Mz> full {mz_full:.6e}, lone {mz_lone:.6e} "
        f"(gap {100.0*gap/abs(mz_full):+.3f} %), cyclic {mz_cyclic:.6e} "
        f"(residual {100.0*residual/abs(mz_full):+.3f} %)")


def test_cyclic_alternating_sector_reproduces_the_full_ring():
    """The (-1)^k cyclic lane must reproduce an explicit alternating N/S ring.

    Example 7 uses this exact symmetry class.  The kernel tests pin image signs and rotations
    separately; this end-to-end solve proves that their composition has the same physical solution as
    four explicitly meshed poles driven with the material pattern (+H0, -H0, +H0, -H0).
    """
    with ng.TaskManager():
        ring = _alternating_pole_mesh(range(N_FOLD))
        sector = _alternating_pole_mesh([0])
        alternating_drive = ring.MaterialCF({
            f"pole{k}": ng.CF((0.0, 0.0, H0 if k % 2 == 0 else -H0))
            for k in range(N_FOLD)
        })
        mz_full = _mean_mz(
            ring, H_ext=alternating_drive, material="pole0")
        mz_lone = _mean_mz(sector, material="pole0")
        mz_cyclic = _mean_mz(
            sector, material="pole0",
            image_cyclic=N_FOLD, image_cyclic_alternating=True)

    gap = abs(mz_lone - mz_full)
    residual = abs(mz_cyclic - mz_full)
    assert gap / abs(mz_full) > 1e-3, (
        "the alternating poles barely interact, so the lane cannot discriminate "
        f"(<Mz> full {mz_full:.6e}, lone {mz_lone:.6e})")
    assert residual < (1.0 - CLOSE_FRACTION) * gap, (
        f"alternating cyclic reduction did not close the gap: <Mz> full {mz_full:.6e}, "
        f"lone {mz_lone:.6e} (gap {100.0*gap/abs(mz_full):+.3f} %), "
        f"cyclic {mz_cyclic:.6e} (residual {100.0*residual/abs(mz_full):+.3f} %)")


def test_cyclic_costs_one_sector():
    """The payoff: the reduced model carries ~1/N of the elements."""
    with ng.TaskManager():
        ring = _pole_mesh(range(N_FOLD))
        sector = _pole_mesh([0])
    assert sector.ne * (N_FOLD - 1) < ring.ne, (
        f"sector {sector.ne} elements vs ring {ring.ne}: the reduction is not saving work")
