import pytest

@pytest.mark.parametrize('complex_system,iterative', [(False,False),(True,False),(False,True)])
def test_scalar_linear_solve_has_true_residual(complex_system, iterative):
    import ngsolve as ng
    from netgen.csg import unit_cube
    from radia.scalar_potential_solver import ScalarPotentialSolver
    mesh=ng.Mesh(unit_cube.GenerateMesh(maxh=.4))
    fes=ng.H1(mesh,order=2,dirichlet='.*',complex=complex_system)
    u,v=fes.TnT()
    a=ng.BilinearForm(fes,symmetric=True)
    a+=(ng.grad(u)*ng.grad(v)+(1+.3j if complex_system else 1)*u*v)*ng.dx
    f=ng.LinearForm(fes)
    f+=(1+2j if complex_system else 1)*v*ng.dx
    pre=ng.Preconditioner(a,'bddc',inverse='sparsecholesky') if iterative else None
    gf=ng.GridFunction(fes)
    solver=ScalarPotentialSolver(mesh)
    with ng.TaskManager():
        a.Assemble(); f.Assemble()
        solver._solve_system(a,f,fes,gf,pre,iterative)
        free=ng.Projector(fes.FreeDofs(),True)
        r=f.vec.CreateVector(); r.data=free*(f.vec-a.mat*gf.vec)
        assert ng.Norm(r)/ng.Norm(f.vec)<1e-8
        assert ng.Norm(gf.vec)>0
