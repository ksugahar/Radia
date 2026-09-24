import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from radia.meshed_current import solve_closed_coil_current


@pytest.fixture(scope="module")
def ring():
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    parts = [Box(Pnt(1,-1,0),Pnt(2,1,1)), Box(Pnt(-2,-1,0),Pnt(-1,1,1)),
             Box(Pnt(-2,1,0),Pnt(2,2,1)), Box(Pnt(-2,-2,0),Pnt(2,-1,1))]
    for i, part in enumerate(parts):
        part.mat("drive" if i == 0 else "return")
    return ng.Mesh(OCCGeometry(Glue(parts)).GenerateMesh(maxh=.7))


def solve(mesh, current):
    drive = mesh.MaterialCF({"drive":ng.CF((0,1,0))}, default=ng.CF((0,0,0)))
    with ng.TaskManager():
        return solve_closed_coil_current(mesh, drive, current_A=current,
                    drive_length_m=2, materials=("drive","return"), inverse="pardiso")


def test_current_conservation_sign_and_zero(ring):
    positive, negative, zero = (solve(ring,i) for i in (2,-2,0))
    assert positive['stats']['relative_divergence'] < 1e-10
    assert positive['stats']['section_current_A'] == pytest.approx(2,abs=1e-12)
    np.testing.assert_allclose(negative['current'].vec.FV().NumPy(),
                              -positive['current'].vec.FV().NumPy(),atol=1e-12)
    np.testing.assert_array_equal(zero['current'].vec.FV().NumPy(),0)


def test_gradient_drive_has_no_closed_current(ring):
    with ng.TaskManager(), pytest.raises(ValueError,match="closed current path"):
        solve_closed_coil_current(ring, ng.CF((1,0,0)), current_A=1,
                                 drive_length_m=1, materials=("drive","return"))


@pytest.mark.parametrize('current,length', [(float('nan'),1),(1,0),(1,float('inf'))])
def test_invalid_normalization(ring,current,length):
    with pytest.raises(ValueError,match="finite"):
        solve_closed_coil_current(ring,ng.CF((0,1,0)),current_A=current,
                                 drive_length_m=length,materials=("drive","return"))


def test_disconnected_coils_require_separate_normalization():
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    first=Box(Pnt(0,0,0),Pnt(1,1,1))
    second=Box(Pnt(2,0,0),Pnt(3,1,1))
    first.mat('coil'); second.mat('coil')
    mesh=ng.Mesh(OCCGeometry(Glue([first,second])).GenerateMesh(maxh=1))
    with pytest.raises(ValueError,match='face-connected'):
        solve_closed_coil_current(mesh,ng.CF((0,1,0)),current_A=1,drive_length_m=1)


# --- A-phi (scalar potential with one thick cut) --------------------------------

from radia.meshed_current import solve_closed_coil_current_phi  # noqa: E402


def ring_mesh(maxh):
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    parts = [Box(Pnt(1,-1,0),Pnt(2,1,1)), Box(Pnt(-2,-1,0),Pnt(-1,1,1)),
             Box(Pnt(-2,1,0),Pnt(2,2,1)), Box(Pnt(-2,-2,0),Pnt(2,-1,1))]
    for i, part in enumerate(parts):
        part.mat("drive" if i == 0 else "return")
    return ng.Mesh(OCCGeometry(Glue(parts)).GenerateMesh(maxh=maxh))


def solve_phi(mesh, current, radius=1.0, origin=(1.5, 0.0, 0.5)):
    with ng.TaskManager():
        return solve_closed_coil_current_phi(
            mesh, current_A=current, cut_origin=origin, cut_normal=(0, 1, 0),
            cut_radius_m=radius, materials=("drive", "return"))


def test_phi_current_sign_value_and_weak_divergence(ring):
    positive, negative, zero = (solve_phi(ring, i) for i in (2, -2, 0))
    stats = positive["stats"]
    assert stats["relative_weak_divergence"] < 1e-10
    assert stats["cut_face_flux_A"] == pytest.approx(2, rel=0.05)
    drive = ring.MaterialCF({"drive": ng.CF((0, 1, 0))}, default=ng.CF((0, 0, 0)))
    section = ng.Integrate(ng.InnerProduct(positive["current"], drive), ring,
                           definedon=ring.Materials("drive")) / 2
    assert section == pytest.approx(2, rel=0.05)
    for p, n in zip(positive["components"], negative["components"]):
        np.testing.assert_allclose(n.vec.FV().NumPy(), -p.vec.FV().NumPy(), atol=1e-12)
    for z in zero["components"]:
        np.testing.assert_array_equal(z.vec.FV().NumPy(), 0)


def test_phi_current_is_orthogonal_to_all_mesh_gradients(ring):
    current = solve_phi(ring, 2)["current"]
    space = ng.H1(ring, order=1)
    w = space.TestFunction()
    region = ring.Materials("drive|return")
    load = ng.LinearForm(space)
    load += ng.InnerProduct(current, ng.grad(w)) * ng.dx(definedon=region)
    scale = ng.LinearForm(space)
    scale += ng.Norm(current) * ng.Norm(ng.grad(w)) * ng.dx(definedon=region)
    with ng.TaskManager():
        load.Assemble()
        scale.Assemble()
    assert np.linalg.norm(load.vec.FV().NumPy()) < 1e-11 * np.linalg.norm(scale.vec.FV().NumPy())


def test_phi_and_mixed_currents_converge_to_one_distribution():
    differences = []
    for maxh in (.7, .35):
        mesh = ring_mesh(maxh)
        drive = mesh.MaterialCF({"drive": ng.CF((0, 1, 0))}, default=ng.CF((0, 0, 0)))
        with ng.TaskManager():
            mixed = solve_closed_coil_current(mesh, drive, current_A=2, drive_length_m=2,
                                              materials=("drive", "return"))["current"]
        primal = solve_phi(mesh, 2)["current"]
        region = mesh.Materials("drive|return")
        delta = ng.Integrate(ng.InnerProduct(primal - mixed, primal - mixed), mesh, definedon=region)
        norm = ng.Integrate(ng.InnerProduct(mixed, mixed), mesh, definedon=region)
        differences.append(float(np.sqrt(delta / norm)))
    assert differences[1] < differences[0] < 0.5


def test_phi_cut_must_cross_and_span_the_conductor(ring):
    # A radius smaller than the unit section's half diagonal cuts only part of it.
    with pytest.raises(ValueError, match="does not span"):
        solve_phi(ring_mesh(.25), 1, radius=0.3)
    with pytest.raises(ValueError, match="does not cross"):
        solve_phi(ring, 1, origin=(1.5, 0.0, 5.0), radius=0.5)
    with pytest.raises(ValueError, match="positive cut radius"):
        solve_phi(ring, 1, radius=0.0)
