"""Owned normalized HACApK transfer, without dense Gram reconstruction or CUDA."""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh
from radia.vim import ChargeGram


def test_owned_normalized_snapshot_preserves_the_cpu_operator():
    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=2)
        B, gram, _ = ChargeGram(fes, eps=1e-12, leafsize=16, eta=2.)
        for budget in (0, 513, True, False):
            with pytest.raises(ValueError, match="budget"):
                gram._resident_leaf_snapshot(budget)
        snapshot = gram._resident_leaf_snapshot(1)
        assert snapshot["schema"] == "radia.normalized-hacapk-leaves.v1"
        assert snapshot["normalization_active"]
        assert snapshot["owned_cpp_payload_upper_bound_bytes"] <= 2**20
        perm, scale = snapshot["permutation"], snapshot["scale"]
        assert perm.dtype == np.uint64 and scale.dtype == np.float64
        assert np.array_equal(np.sort(perm), np.arange(gram.ndof()))
        x = np.random.default_rng(20261020).standard_normal(gram.ndof())
        expected = np.asarray(gram.matvec_sym(x))
        xp = x[perm]*scale[perm]
        yp = np.zeros_like(xp)
        for leaf in snapshot["leaves"]:
            r, c = leaf["row"], leaf["column"]
            nr, nc, rank = leaf["rows"], leaf["columns"], leaf["rank"]
            rows, cols = slice(r, r+nr), slice(c, c+nc)
            if rank:
                right = leaf["a1"].reshape((nc, rank), order="F")
                left = leaf["a2"].reshape((nr, rank), order="F")
                yp[rows] += left @ (right.T @ xp[cols])
                if r != c:
                    yp[cols] += right @ (left.T @ xp[rows])
            else:
                block = leaf["a1"].reshape((nr, nc))
                yp[rows] += block @ xp[cols]
                if r != c:
                    yp[cols] += block.T @ xp[rows]
        actual = np.empty_like(xp);actual[perm] = yp*scale[perm]
        assert np.linalg.norm(actual-expected) <= 1e-12*np.linalg.norm(expected)
        # Mutating the returned transfer data cannot corrupt the native owner.
        snapshot["leaves"][0]["a1"].fill(0)
        snapshot["scale"].fill(0)
        np.testing.assert_array_equal(gram.matvec_sym(x), expected)
        fresh = gram._resident_leaf_snapshot(1)
        assert np.all(fresh["scale"] > 0)
        assert np.any(fresh["leaves"][0]["a1"] != 0)

        # Explicit accuracy oracle retains owner state and rejects malformed calls.
        face_x=np.random.default_rng(10).standard_normal(fes.ndof)
        oracle=gram._apply_configured_demag_compensated(face_x,1)
        assert np.all(np.isfinite(oracle)) and np.linalg.norm(oracle)>0
        np.testing.assert_array_equal(oracle,gram._apply_configured_demag_compensated(face_x,1))
        for budget in (0,513):
            with pytest.raises(ValueError,match="budget"):
                gram._apply_configured_demag_compensated(face_x,budget)
        with pytest.raises(ValueError,match="size"):
            gram._apply_configured_demag_compensated(face_x[:-1],1)
        np.testing.assert_array_equal(gram.matvec_sym(x),expected)


def test_explicit_compensated_pcg_preserves_original_operator_and_rejects_budget():
    from scipy.sparse import coo_matrix
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1)
    with ng.TaskManager():
        fes=ng.HDiv(mesh,order=2);_,gram,_=ChargeGram(fes,eps=1e-12,leafsize=16,eta=2.)
        u,v=fes.TnT();form=ng.BilinearForm(fes);form+=u*v*ng.dx;form.Assemble()
        ii,jj,vv=form.mat.COO();mass=coo_matrix((vv,(ii,jj)),shape=(fes.ndof,fes.ndof)).tocsr()
        gram.configure_mass_matrix_ngsolve(form.mat)
        target=np.random.default_rng(31).standard_normal(fes.ndof)
        rhs=np.ascontiguousarray(mass@target+gram._apply_configured_demag_compensated(target,16))
        for budget in (True, False):
            with pytest.raises(ValueError,match="bool"):
                gram._apply_configured_demag_compensated(target,budget)
            with pytest.raises(ValueError,match="bool"):
                gram.solve_configured_linear_material_mass_riesz(1.,rhs,5e-14,10000,compensated_memory_mb=budget)
        with pytest.raises(ValueError,match="budget"):
            gram.solve_configured_linear_material_mass_riesz(1.,rhs,5e-14,10000,compensated_memory_mb=513)
        out=gram.solve_configured_linear_material_mass_riesz(1.,rhs,5e-14,10000,compensated_memory_mb=16)
        x=np.asarray(out['m'])
        assert out['iters']<10000
        assert np.linalg.norm(x-target)<=1e-10*np.linalg.norm(target)
        assert np.linalg.norm(mass@x+gram._apply_configured_demag_compensated(x,16)-rhs)<=1e-11*np.linalg.norm(rhs)


def test_compensated_composition_matches_decimal_input_factors():
    from decimal import Decimal,localcontext
    # Four elongated cells create separated source supports. A unit cube alone
    # has only dense leaves, so it cannot certify the paired low-rank branch.
    mesh=MakeStructured3DMesh(hexes=True,nx=4,ny=1,nz=1,
                             mapping=lambda x,y,z:(60*x,y,z))
    with ng.TaskManager(),localcontext() as context:
        context.prec=60
        fes=ng.HDiv(mesh,order=2);B,gram,_=ChargeGram(fes,eps=1e-8,leafsize=16,eta=2.)
        tree=gram._resident_leaf_snapshot(16);n=tree['dimension']
        assert n<=512 and B.shape[0]%n==0
        coverage={name:0 for name in ('diagonal_dense','offdiagonal_dense','offdiagonal_lowrank')}
        for leaf in tree['leaves']:
            name=('diagonal' if leaf['row']==leaf['column'] else 'offdiagonal')+('_lowrank' if leaf['rank'] else '_dense')
            assert name in coverage, "Unexpected leaf family in composition oracle"
            coverage[name]+=1
        assert all(coverage.values()), coverage
        D=lambda value:Decimal.from_float(float(value))
        state=np.random.default_rng(31).standard_normal(fes.ndof)
        result=[Decimal(0) for _ in range(fes.ndof)]
        for component in range(B.shape[0]//n):
            physical=[]
            for a in range(n):
                row=component*n+a
                physical.append(sum((D(B.data[k])*D(state[B.indices[k]]) for k in range(B.indptr[row],B.indptr[row+1])),Decimal(0))*D(tree['scale'][a]))
            xp=[physical[int(j)] for j in tree['permutation']];yp=[Decimal(0) for _ in range(n)]
            for leaf in tree['leaves']:
                row,col,nr,nc,rank=[int(leaf[k]) for k in ('row','column','rows','columns','rank')]
                a1=[D(v) for v in leaf['a1']];a2=[D(v) for v in leaf['a2']]
                if rank:
                    for h in range(rank):
                        t=sum((a1[j+nc*h]*xp[col+j] for j in range(nc)),Decimal(0))
                        for i in range(nr):yp[row+i]+=a2[i+nr*h]*t
                        if row!=col:
                            t=sum((a2[i+nr*h]*xp[row+i] for i in range(nr)),Decimal(0))
                            for j in range(nc):yp[col+j]+=a1[j+nc*h]*t
                else:
                    for i in range(nr):yp[row+i]+=sum((a1[i*nc+j]*xp[col+j] for j in range(nc)),Decimal(0))
                    if row!=col:
                        for j in range(nc):yp[col+j]+=sum((a1[i*nc+j]*xp[row+i] for i in range(nr)),Decimal(0))
            for a,j in enumerate(tree['permutation']):physical[int(j)]=yp[a]*D(tree['scale'][int(j)])
            for a in range(n):
                row=component*n+a
                for k in range(B.indptr[row],B.indptr[row+1]):result[int(B.indices[k])]+=D(B.data[k])*physical[a]
        expected=np.array([float(v) for v in result]);actual=gram._apply_configured_demag_compensated(state,16)
        assert np.linalg.norm(actual-expected)<=2e-15*np.linalg.norm(expected)
