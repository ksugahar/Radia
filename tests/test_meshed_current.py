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
