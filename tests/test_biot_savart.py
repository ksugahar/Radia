import numpy as np


def test_finite_filament_coefficient_function_matches_point_formula_and_has_gradient():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.biot_savart import h_filament, h_filament_cf

    start=(-0.2,0.0,0.0); end=(0.3,0.1,-0.05)
    point=np.array([0.04,0.2,0.13]); current=17.0
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1)
    field=h_filament_cf(start,end,current=current)
    value=np.asarray(field(mesh(*point)),dtype=float)
    np.testing.assert_allclose(
        value,h_filament(start,end,point,current=current),rtol=3e-15,atol=3e-15)

    # NGSolve symbolically differentiates the same analytic expression.  The
    # nonzero lock prevents a generic, value-only CoefficientFunction from
    # silently dropping coil-kernel motion in GetTrafo shape derivatives.
    gradient=np.column_stack([
        np.asarray(field.Diff(coordinate)(mesh(*point)),dtype=float)
        for coordinate in (ng.x,ng.y,ng.z)])
    assert np.linalg.norm(gradient)>1.0
    np.testing.assert_allclose(np.trace(gradient),0.0,atol=3e-13)
