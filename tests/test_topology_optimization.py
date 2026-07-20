import numpy as np

from radia.topology_optimization import (
    affine_cell_self_energy_shape_derivative,
    assemble_ngsolve_hdiv_shape_tangents,
    linearize_laplace_charge_gram,
    production_hex_volume_self_block_derivatives,
    production_hex_face_self_block_derivatives,
    production_wedge_volume_self_block_derivatives,
    production_wedge_face_self_block_derivatives,
    production_wedge_charge_gram_derivatives,
    production_tet_volume_self_block_derivatives,
    production_tet_face_self_block_derivatives,
    linearize_laplace_pair_gram, linearize_vim_operator,
    VIMLinearization,
    linearize_vim_system,
    optimize_vim_lp,
    solve_lp_update,
    write_cubit_density_journal,
    linearize_production_vim_from_ngsolve,
    linearize_production_vim_matrix_free_from_ngsolve,
    sample_production_gettrafo_displacements,
)


def test_ngsolve_gettrafo_production_hex_closes_full_vim_scaling_tangent():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
        mapping=lambda x,y,z:(1.1*x+.03*y*z,.9*y+.02*x*z,1.2*z+.04*x*y))
    fes=ng.HDiv(mesh,order=1)
    with ng.TaskManager():
        cb=_charge_basis_hex(fes,cob_quad=3)
        B,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
        mode=ng.GridFunction(ng.VectorH1(mesh,order=1))
        mode.Set(ng.CF((ng.x,ng.y,ng.z)))
        result=linearize_production_vim_from_ngsolve(fes=fes,
            deformation_modes=[mode],charge_basis=cb,charge_gram=gram,
            charge_map=B,applied_coefficients=np.ones(fes.ndof),inv_chi=.2,
            family="hex")
        _,mf_charge,mf=linearize_production_vim_matrix_free_from_ngsolve(
            fes=fes,deformation_modes=[mode],charge_basis=cb,charge_gram=gram,
            charge_map=B,inv_chi=.2,family="hex",eps=1e-12,leaf=256,eta=2.0)
    np.testing.assert_allclose(result.charge_gram.jacobian[0],
        -result.charge_gram.matrix,rtol=4e-10,atol=4e-13)
    np.testing.assert_allclose(result.operator.matrix_jacobian[0],
        -result.operator.matrix,rtol=2e-9,atol=2e-11)
    np.testing.assert_allclose(result.operator.rhs_jacobian[0],
        -result.operator.rhs,rtol=2e-9,atol=2e-11)
    probe=np.linspace(-.3,.6,fes.ndof)
    np.testing.assert_allclose(mf.matvec(probe),result.operator.matrix@probe,
        rtol=2e-11,atol=2e-12)
    np.testing.assert_allclose(mf.directional_matvec(0,probe),
        result.operator.matrix_jacobian[0]@probe,rtol=2e-10,atol=2e-11)
    assert not hasattr(mf_charge,"matrix") or mf_charge.matrix is gram


def test_gettrafo_sampler_routes_tet_and_mixed_wedge_face_lattices():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis,_charge_basis_wedge
    for family,mesh,builder in (
        ("tet",MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1),
         lambda f:_charge_basis(f,4)),
        ("wedge",MakeStructured3DMesh(prism=True,nx=1,ny=1,nz=1),
         _charge_basis_wedge)):
        fes=ng.HDiv(mesh,order=1);vf=ng.VectorH1(mesh,order=1)
        mode=ng.GridFunction(vf);mode.Set(ng.CF((.01*ng.x,.02*ng.y,.03*ng.z)))
        with ng.TaskManager(): cb=builder(fes)
        sampled=sample_production_gettrafo_displacements(
            fes,[mode],cb,family=family)
        assert sampled.family==family and all(x.shape[0]==1 for x in sampled.cell)
        if family=="tet":
            assert all(x.shape[1:]==(4,3) for x in sampled.cell)
            assert all(x.shape[1:]==(3,3) for x in sampled.face)
        else:
            assert all(x.shape[1:]==(18,3) for x in sampled.cell)
            assert {x.shape[1] for x in sampled.face}=={6,9}


def test_production_wedge_self_block_python_boundaries_preserve_host_mode_order():
    class Gram:
        def wedge_volume_self_block_directional_derivative(self,host,velocity):
            return np.full((2,2),host+np.asarray(velocity)[0,0])
        def wedge_face_self_block_directional_derivative(self,host,velocity):
            return np.full((3,3),host+np.asarray(velocity)[0,0])
    volume=production_wedge_volume_self_block_derivatives(
        Gram(),[np.stack([np.zeros((18,3)),np.ones((18,3))])])
    faces=production_wedge_face_self_block_derivatives(
        Gram(),[np.stack([np.zeros((6,3)),np.ones((6,3))]),
                np.stack([2*np.ones((9,3)),3*np.ones((9,3))])])
    assert volume[0].shape==(2,2,2)
    assert np.all(volume[0][1]==1)
    assert faces[0].shape==(2,3,3) and faces[1].shape==(2,3,3)
    assert np.all(faces[1][0]==3)


def test_production_wedge_full_gram_python_boundary_preserves_mode_order():
    class Gram:
        def wedge_charge_gram_directional_derivative(self,cells,faces):
            return np.eye(2)*(np.asarray(cells)[0,0,0]+np.asarray(faces)[0,0,0])
    cells=np.zeros((2,1,18,3));faces=np.zeros((2,3,9,3))
    cells[1,0,0,0]=2;faces[1,0,0,0]=3
    out=production_wedge_charge_gram_derivatives(Gram(),cells,faces)
    assert out.shape==(2,2,2) and np.array_equal(out[1],5*np.eye(2))


def test_native_production_wedge_self_block_derivative_invariants():
    import pytest
    ng=pytest.importorskip("ngsolve")
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _build_charge_gram_wedge, _charge_basis_wedge
    try:
        mesh=MakeStructured3DMesh(nx=1,ny=1,nz=1,prism=True)
    except TypeError:
        pytest.skip("this NGSolve build cannot generate a prism mesh")
    with ng.TaskManager():
        fes=ng.HDiv(mesh,order=1)
        cb=_charge_basis_wedge(fes)
        _,gram,_,_=_build_charge_gram_wedge(fes,eps=1e-14,leafsize=256)
    def value_block(kind,host):
        ids=np.flatnonzero((np.asarray(cb["kind"])==kind)&(np.asarray(cb["host"])==host))
        return np.asarray([[gram.entry(int(i),int(j)) for j in ids] for i in ids])
    cell=np.asarray(cb["cell_nodes"]).reshape(-1,18,3)[0]
    dc=gram.wedge_volume_self_block_directional_derivative(0,cell)
    tc=gram.wedge_volume_self_block_directional_derivative(0,np.ones((18,3)))
    assert np.linalg.norm(dc+value_block(0,0))/np.linalg.norm(value_block(0,0))<2e-10
    assert np.linalg.norm(tc)<2e-12 and np.array_equal(dc,dc.T)
    face_nodes=np.asarray(cb["face_nodes"]).reshape(-1,9,3)
    for host,ft in enumerate(np.asarray(cb["face_type"])):
        nn=6 if ft==0 else 9; nodes=face_nodes[host,:nn]
        df=gram.wedge_face_self_block_directional_derivative(host,nodes)
        tf=gram.wedge_face_self_block_directional_derivative(host,np.ones((nn,3)))
        block=value_block(1,host)
        assert np.linalg.norm(df+block)/np.linalg.norm(block)<2e-10
        assert np.linalg.norm(tf)<2e-12 and np.array_equal(df,df.T)


def test_native_production_wedge_self_block_derivative_matches_general_affine_fd():
    import pytest
    ng=pytest.importorskip("ngsolve")
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _build_charge_gram_wedge, _charge_basis_wedge
    A=np.array([[.17,-.09,.04],[.06,-.13,.08],[-.03,.05,.11]])
    base=np.array([[1.0,.12,-.04],[.03,.91,.08],[-.06,.05,1.13]])
    shift=np.array([.02,-.01,.03])
    def make(eps):
        def mapping(x,y,z):
            q=np.array([x,y,z],dtype=float)
            r=base@q+eps*(A@q+shift)
            return tuple(r)
        try:
            mesh=MakeStructured3DMesh(nx=1,ny=1,nz=1,prism=True,mapping=mapping)
        except TypeError:
            pytest.skip("this NGSolve build cannot generate a prism mesh")
        with ng.TaskManager():
            fes=ng.HDiv(mesh,order=1)
            cb=_charge_basis_wedge(fes)
            _,gram,_,_=_build_charge_gram_wedge(fes,eps=1e-14,leafsize=256)
        return cb,gram
    step=2e-6
    cb0,g0=make(0.0); cbp,gp=make(step); cbm,gm=make(-step)
    kind=np.asarray(cb0["kind"]); hosts=np.asarray(cb0["host"])
    def block(gram,k,h):
        ids=np.flatnonzero((kind==k)&(hosts==h))
        return np.asarray([[gram.entry(int(i),int(j)) for j in ids] for i in ids])
    cells=np.asarray(cb0["cell_nodes"]).reshape(-1,18,3)
    faces=np.asarray(cb0["face_nodes"]).reshape(-1,9,3)
    types=np.asarray(cb0["face_type"])
    cases=[("volume",0,0,cells[0])]
    for ft,name in ((0,"tri"),(1,"quad")):
        h=int(np.flatnonzero(types==ft)[0]); nn=6 if ft==0 else 9
        cases.append((name,1,h,faces[h,:nn]))
    for name,k,h,nodes in cases:
        reference=nodes@np.linalg.inv(base).T
        velocity=reference@A.T+shift
        analytic=(g0.wedge_volume_self_block_directional_derivative(h,velocity) if k==0
                  else g0.wedge_face_self_block_directional_derivative(h,velocity))
        fd=(block(gp,k,h)-block(gm,k,h))/(2*step)
        relative=np.linalg.norm(analytic-fd)/np.linalg.norm(fd)
        assert relative<2e-7,(name,relative)


def test_native_production_wedge_full_gram_derivative_translation_scale_and_local_fd():
    import pytest
    ng=pytest.importorskip("ngsolve")
    from ngsolve.meshes import MakeStructured3DMesh
    from radia import _radia_pybind as rp
    from radia.vim._vim import (_charge_basis_wedge,_g01,_SYM5_TET,_SYM5_TRI,
                                _f64_buffer,_i32_buffer)
    try:mesh=MakeStructured3DMesh(nx=2,ny=1,nz=1,prism=True)
    except TypeError:pytest.skip("this NGSolve build cannot generate a prism mesh")
    with ng.TaskManager():cb=_charge_basis_wedge(ng.HDiv(mesh,order=1))
    cells=np.asarray(cb["cell_nodes"],dtype=float).reshape(-1,18,3)
    faces=np.asarray(cb["face_nodes"],dtype=float).reshape(-1,9,3).copy()
    types=np.asarray(cb["face_type"],dtype=np.int32)
    glo,gwo=_g01(6);gli,gwi=_g01(5)
    tri=ng.IntegrationRule(ng.ET.TRIG,5)
    ftp=np.asarray([(ip.point[0],ip.point[1]) for ip in tri]);ftw=np.asarray([ip.weight for ip in tri])
    def make(c,f,near=.6,far=1.5):
        return rp._ChargeGramHMatrix(
            wedge_cell_nodes=_f64_buffer(c),face_nodes=_f64_buffer(f),face_type=_i32_buffer(types),
            n_el=int(cb["n_el"]),n_bf=int(cb["n_bf"]),charge_host=_i32_buffer(cb["host"]),
            charge_kind=_i32_buffer(cb["kind"]),charge_expo=_i32_buffer(cb["expo"]),
            sym_tet_pts=_f64_buffer(_SYM5_TET[0]),sym_tet_w=_f64_buffer(_SYM5_TET[1]),
            sym_tri_pts=_f64_buffer(_SYM5_TRI[0]),sym_tri_w=_f64_buffer(_SYM5_TRI[1]),
            field_tri_pts=_f64_buffer(ftp),field_tri_w=_f64_buffer(ftw),
            gl_out=_f64_buffer(glo),gw_out=_f64_buffer(gwo),gl_in=_f64_buffer(gli),gw_in=_f64_buffer(gwi),
            far_tet_pts=_f64_buffer(_SYM5_TET[0]),far_tet_w=_f64_buffer(_SYM5_TET[1]),
            far_tri_pts=_f64_buffer(_SYM5_TRI[0]),far_tri_w=_f64_buffer(_SYM5_TRI[1]),
            near_grade=near,far_inner_factor=far,image_masks=np.empty(0,np.int32),
            image_signs=np.empty(0),eps=1e-14,leaf=256,eta=2.,build=False)
    n=len(cb["kind"])
    def dense(g):return np.asarray([[g.entry(i,j) for j in range(n)] for i in range(n)])
    g0=make(cells,faces);G=dense(g0)
    ones_c=np.ones_like(cells);ones_f=np.ones_like(faces)
    dt=np.asarray(g0.wedge_charge_gram_directional_derivative(ones_c,ones_f))
    ds=np.asarray(g0.wedge_charge_gram_directional_derivative(cells,faces))
    dop=g0.directional_derivative_operator("wedge",cells,faces,
        eps=1e-12,leaf=256,eta=2.0)
    probe=np.linspace(-.4,.7,n)
    np.testing.assert_allclose(dop.matvec_sym(probe),ds@probe,rtol=2e-12,atol=2e-12)
    assert np.linalg.norm(dt)<3e-11
    assert np.linalg.norm(ds+G)/np.linalg.norm(G)<3e-10
    assert np.array_equal(dt,dt.T) and np.array_equal(ds,ds.T)
    # A continuous piecewise-affine hat field is localized to the interior
    # x-plane.  It agrees exactly on duplicate cell/face nodes and preserves
    # affine hosts, just like an H1/GetTrafo P1 deformation mode.
    def local_velocity(x):
        h=np.maximum(0.,1.-2.*np.abs(x[...,0]-.5))
        return h[...,None]*np.array([.031,-.027,.023])
    vc,vf=local_velocity(cells),local_velocity(faces)
    gf=make(cells,faces,near=0.,far=0.)
    analytic=np.asarray(gf.wedge_charge_gram_directional_derivative(vc,vf))
    step=2e-6
    Gp=dense(make(cells+step*vc,faces+step*vf,near=0.,far=0.));Gm=dense(make(cells-step*vc,faces-step*vf,near=0.,far=0.))
    fd=(Gp-Gm)/(2*step)
    # Touching-but-separately-integrated host pairs have a discrete nearest-site
    # quadrature dispatch and are not a differentiable FD oracle.  Compare the
    # localized non-self derivative on the farthest cell pair; all self blocks
    # are independently locked above, while translation/scale cover full dG.
    centers=cells.mean(axis=1);dist=np.linalg.norm(centers[:,None]-centers[None,:],axis=2)
    ha,hb=np.unravel_index(np.argmax(dist),dist.shape)
    kind=np.asarray(cb["kind"]);host=np.asarray(cb["host"])
    ia=np.flatnonzero((kind==0)&(host==ha));ib=np.flatnonzero((kind==0)&(host==hb))
    aa=analytic[np.ix_(ia,ib)];ff=fd[np.ix_(ia,ib)]
    relative=np.linalg.norm(aa-ff)/np.linalg.norm(ff)
    assert relative<3e-6,(relative,ha,hb,np.linalg.norm(aa),np.linalg.norm(ff))


def test_affine_self_term_derivative_closes_tet_hex_wedge_diagonal():
    cells={
        "tet":np.array([[0.,0,0],[1.1,0,0],[.1,.9,0],[.2,.1,1.2]]),
        "hex":np.array([[0.,0,0],[1.,0,0],[1.,1,0],[0,1,0],
                        [0,0,1],[1.,0,1],[1.,1,1],[0,1,1.]]),
        "wedge":np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1],[1.,0,1],[0,1.,1]]),
    }
    rng=np.random.default_rng(741)
    for kind,nodes in cells.items():
        direction=rng.normal(size=nodes.shape)
        value,jac=affine_cell_self_energy_shape_derivative(
            kind,nodes,np.stack([np.ones_like(nodes),nodes,direction]))
        assert value>0 and np.all(np.isfinite(jac))
        assert jac[0]==0.0  # rigid translation is removed analytically
        np.testing.assert_allclose(jac[1],5*value,rtol=3e-4,atol=2e-8)
        eps=2e-6
        plus=affine_cell_self_energy_shape_derivative(
            kind,nodes+eps*direction,np.zeros((0,*nodes.shape)))[0]
        minus=affine_cell_self_energy_shape_derivative(
            kind,nodes-eps*direction,np.zeros((0,*nodes.shape)))[0]
        np.testing.assert_allclose((plus-minus)/(2*eps),jac[2],rtol=2e-3,atol=2e-6)


def test_full_charge_gram_combines_nonself_and_analytic_self_tangents():
    cells=[np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]]),
           np.array([[2.,.1,0],[3.1,.1,0],[2.,1.2,0],[2.,.1,.8]])]
    velocity=[np.array([[0,0,0],[.1,0,0],[0,.02,0],[0,0,-.03]]),
              np.array([[.03,0,0],[-.02,0,0],[0,.04,0],[0,0,.01]])]
    def geometry(cells_now):
        points=np.array([c.mean(axis=0) for c in cells_now])
        weights=[]
        for c in cells_now:
            J=np.column_stack((c[1]-c[0],c[2]-c[0],c[3]-c[0]))
            weights.append(abs(np.linalg.det(J))/6)
        return points,np.array(weights)
    points,weights=geometry(cells)
    point_velocity=np.array([[v.mean(axis=0) for v in velocity]])
    eps=2e-6
    plus_cells=[c+eps*v for c,v in zip(cells,velocity)]
    minus_cells=[c-eps*v for c,v in zip(cells,velocity)]
    plus_points,plus_weights=geometry(plus_cells); minus_points,minus_weights=geometry(minus_cells)
    rel=((plus_weights-minus_weights)/(2*eps))/weights
    model=linearize_laplace_charge_gram(points,weights,point_velocity,
        relative_weight_derivatives=rel[None,:],self_cell_types=["tet","tet"],
        self_nodes=cells,self_node_displacements=[np.array([v]) for v in velocity])
    plus=linearize_laplace_charge_gram(plus_points,plus_weights,np.zeros((0,2,3)),
        self_cell_types=["tet","tet"],self_nodes=plus_cells,
        self_node_displacements=[np.zeros((0,4,3)),np.zeros((0,4,3))]).matrix
    minus=linearize_laplace_charge_gram(minus_points,minus_weights,np.zeros((0,2,3)),
        self_cell_types=["tet","tet"],self_nodes=minus_cells,
        self_node_displacements=[np.zeros((0,4,3)),np.zeros((0,4,3))]).matrix
    assert np.all(np.diag(model.matrix)>0)
    np.testing.assert_allclose((plus-minus)/(2*eps),model.jacobian[0],rtol=2e-3,atol=2e-6)


def test_ngsolve_hdiv_mass_shape_tangent_uses_piola_weak_form():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _csr
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1)
    fes=ng.HDiv(mesh,order=1); vf=ng.VectorH1(mesh,order=1)
    velocity=ng.GridFunction(vf); velocity.Set(ng.CF((.07*ng.x,-.03*ng.y,.02*ng.z)))
    with ng.TaskManager():
        mass,dmass,dB=assemble_ngsolve_hdiv_shape_tangents(fes,[velocity],np.eye(fes.ndof))
        eps=2e-6; trial=ng.GridFunction(vf); trial.vec.data=eps*velocity.vec
        mesh.SetDeformation(trial)
        u,v=fes.TnT(); shifted=ng.BilinearForm(fes); shifted+=u*v*ng.dx; shifted.Assemble()
        plus=_csr(shifted).toarray(); mesh.UnsetDeformation()
        trial.vec.data=-eps*velocity.vec; mesh.SetDeformation(trial)
        shifted=ng.BilinearForm(fes); shifted+=u*v*ng.dx; shifted.Assemble()
        minus=_csr(shifted).toarray(); mesh.UnsetDeformation()
    np.testing.assert_allclose((plus-minus)/(2*eps),dmass[0],rtol=2e-6,atol=2e-8)
    assert mass.shape==(fes.ndof,fes.ndof) and np.count_nonzero(dB)==0


def test_production_hex_self_block_python_boundary_preserves_host_mode_order():
    class Gram:
        def hex_volume_self_block_directional_derivative(self,host,velocity):
            return np.full((2,2),host+np.sum(velocity))
    modes=[np.zeros((2,27,3)),np.ones((2,27,3))]
    blocks=production_hex_volume_self_block_derivatives(Gram(),modes)
    assert len(blocks)==2 and blocks[0].shape==(2,2,2)
    np.testing.assert_allclose(blocks[0],0)
    np.testing.assert_allclose(blocks[1],82)


def test_production_hex_face_self_block_python_boundary_preserves_host_mode_order():
    class Gram:
        def hex_face_self_block_directional_derivative(self,host,velocity):
            return np.full((3,3),host+np.sum(velocity))
    modes=[np.zeros((2,9,3)),np.ones((2,9,3))]
    blocks=production_hex_face_self_block_derivatives(Gram(),modes)
    assert len(blocks)==2 and blocks[0].shape==(2,3,3)
    np.testing.assert_allclose(blocks[0],0)
    np.testing.assert_allclose(blocks[1],28)


def test_production_tet_self_block_python_boundaries_preserve_host_mode_order():
    class Gram:
        def tet_volume_self_block_directional_derivative(self,host,velocity): return np.full((2,2),host+np.sum(velocity))
        def tet_face_self_block_directional_derivative(self,host,velocity): return np.full((3,3),host+np.sum(velocity))
    volume=production_tet_volume_self_block_derivatives(Gram(),[np.ones((2,4,3))])
    face=production_tet_face_self_block_derivatives(Gram(),[np.ones((2,3,3))])
    assert volume[0].shape==(2,2,2) and face[0].shape==(2,3,3)
    np.testing.assert_allclose(volume[0],12);np.testing.assert_allclose(face[0],9)


def test_native_production_hex_volume_self_block_derivative_invariants():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
        mapping=lambda x,y,z:(x,y,z+.12*x*y))
    fes=ng.HDiv(mesh,order=1)
    with ng.TaskManager():
        cb=_charge_basis_hex(fes,cob_quad=3)
        # ACA recompression at +/-epsilon is not differentiable; this
        # regression compares against the production dense block path.
        _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
    nodes=np.asarray(cb["cell_nodes"]).reshape(-1,27,3)[0]
    translation=gram.hex_volume_self_block_directional_derivative(0,np.ones((27,3)))
    scaling=gram.hex_volume_self_block_directional_derivative(0,nodes)
    n=scaling.shape[0]
    value=np.array([[gram.entry(i,j) for j in range(n)] for i in range(n)])
    np.testing.assert_allclose(translation,0,atol=2e-17)
    np.testing.assert_allclose(scaling,scaling.T,rtol=0,atol=0)
    # Reference Piola charges have fixed measure; uniform physical scaling
    # therefore differentiates the Laplace kernel with homogeneity -1.
    np.testing.assert_allclose(scaling,-value,rtol=3e-13,atol=3e-15)


def test_native_production_tet_self_block_derivatives_match_scaling_and_fd():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis,build_charge_gram

    gradient=np.array([[.11,-.04,.03],[.02,.08,-.05],[-.01,.06,.09]])
    offset=np.array([.013,-.021,.017]); epsilon=2e-6
    def build(step):
        mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1,
            mapping=lambda x,y,z:tuple(np.array([x,y,z])+step*epsilon*(gradient@np.array([x,y,z])+offset)))
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            cb=_charge_basis(fes,4)
            _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
        return cb,gram
    cb,gram=build(0); _,plus=build(1); _,minus=build(-1)
    kinds=np.asarray(cb["kind"]);hosts=np.asarray(cb["host"])
    for kind,nodes,key in ((0,np.asarray(cb["vV"]),"tet_volume_self_block_directional_derivative"),
                           (1,np.asarray(cb["bV"]),"tet_face_self_block_directional_derivative")):
        method=getattr(gram,key)
        for host,vertices in enumerate(nodes):
            ids=np.flatnonzero((kinds==kind)&(hosts==host))
            value=np.array([[gram.entry(int(i),int(j)) for j in ids] for i in ids])
            value_plus=np.array([[plus.entry(int(i),int(j)) for j in ids] for i in ids])
            value_minus=np.array([[minus.entry(int(i),int(j)) for j in ids] for i in ids])
            scaling=method(host,vertices)
            translation=method(host,np.ones_like(vertices))
            velocity=vertices@gradient.T+offset
            derivative=method(host,velocity)
            fd=(value_plus-value_minus)/(2*epsilon)
            np.testing.assert_allclose(translation,0,atol=8e-17)
            np.testing.assert_allclose(scaling,scaling.T,rtol=0,atol=0)
            # Flat TET ChargeGram stores physical charge monomials and both
            # physical measures: volume-volume and surface-surface blocks are
            # homogeneous of degree five and three respectively.  The Piola
            # B-map derivatives close the final B.T@G@B product separately.
            np.testing.assert_allclose(scaling,(5 if kind==0 else 3)*value,rtol=8e-12,atol=4e-15)
            np.testing.assert_allclose(derivative,fd,rtol=2e-7,atol=3e-11)

    # Reverse the physical orientation so the same regression covers both
    # signs of det(E); the analytic moments are physical (unsigned), while
    # the affine coordinate map retains the production orientation convention.
    reflection=np.diag([-1.,1.,1.])
    def build_reflected(step):
        mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1,
            mapping=lambda x,y,z:tuple((np.eye(3)+step*epsilon*gradient)@(reflection@np.array([x,y,z]))+step*epsilon*offset))
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            local_cb=_charge_basis(fes,4)
            _,local_gram,_=build_charge_gram(fes,eps=1e-10,leafsize=16,eta=2.0)
        return local_cb,local_gram
    rcb,rgram=build_reflected(0);_,rplus=build_reflected(1);_,rminus=build_reflected(-1)
    rk=np.asarray(rcb["kind"]);rh=np.asarray(rcb["host"]);vertices=np.asarray(rcb["vV"])[0]
    ids=np.flatnonzero((rk==0)&(rh==0));value=np.array([[rgram.entry(int(i),int(j)) for j in ids] for i in ids])
    derivative=rgram.tet_volume_self_block_directional_derivative(0,vertices@gradient.T+offset)
    fd=np.array([[(rplus.entry(int(i),int(j))-rminus.entry(int(i),int(j)))/(2*epsilon) for j in ids] for i in ids])
    np.testing.assert_allclose(rgram.tet_volume_self_block_directional_derivative(0,vertices),5*value,rtol=8e-12,atol=4e-15)
    np.testing.assert_allclose(derivative,fd,rtol=2e-7,atol=3e-11)


def test_native_production_tet_complete_gram_and_piola_product_derivative():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.topology_optimization import production_tet_charge_gram_derivatives
    from radia.vim._vim import _charge_basis,build_charge_gram

    gradient=np.array([[.073,-.031,.019],[.014,.052,-.027],[-.022,.041,.064]])
    offset=np.array([.009,-.015,.011]); epsilon=1e-6
    def build(step, scaling=False, translation=False):
        def mapping(x,y,z):
            point=np.array([x,y,z])
            if scaling: velocity=point
            elif translation: velocity=np.ones(3)
            else: velocity=gradient@point+offset
            return tuple(point+step*epsilon*velocity)
        mesh=MakeStructured3DMesh(hexes=False,nx=1,ny=1,nz=1,mapping=mapping)
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            cb=_charge_basis(fes,4)
            # Keep this derivative regression on one dense leaf.  Rebuilding
            # ACA factors at +/-epsilon is not a differentiable reference.
            _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
        n=len(cb["host"])
        G=np.array([[gram.entry(i,j) for j in range(n)] for i in range(n)])
        return cb,gram,G

    cb,gram,G=build(0); cbp,_,Gp=build(1); cbm,_,Gm=build(-1)
    cells=np.asarray(cb["vV"]); faces=np.asarray(cb["bV"])
    cell_v=cells@gradient.T+offset; face_v=faces@gradient.T+offset
    dG,dB=production_tet_charge_gram_derivatives(
        gram,cell_v[None,...],face_v[None,...],cb["B"])
    dG=dG[0]; dB=dB[0].toarray(); B=cb["B"].toarray()
    dop=gram.directional_derivative_operator("tet",cell_v,face_v,
        eps=1e-12,leaf=256,eta=2.0)
    probe=np.linspace(-.6,.8,dG.shape[0])
    np.testing.assert_allclose(dop.matvec_sym(probe),dG@probe,rtol=2e-12,atol=2e-12)
    fdG=(Gp-Gm)/(2*epsilon)
    np.testing.assert_allclose(dG,dG.T,rtol=0,atol=0)
    np.testing.assert_allclose(dG,fdG,rtol=4e-7,atol=8e-11)
    N=B.T@G@B
    dN=dB.T@G@B+B.T@dG@B+B.T@G@dB
    Np=cbp["B"].toarray().T@Gp@cbp["B"].toarray()
    Nm=cbm["B"].toarray().T@Gm@cbm["B"].toarray()
    np.testing.assert_allclose(dN,(Np-Nm)/(2*epsilon),rtol=7e-7,atol=2e-10)

    # Translation leaves every raw block and Piola row unchanged.
    _,gt,_=build(0,translation=True)
    dGt,rt=production_tet_charge_gram_derivatives(
        gt,np.ones_like(cells)[None,...],np.ones_like(faces)[None,...])
    np.testing.assert_allclose(dGt,0,atol=2e-16)
    np.testing.assert_allclose(rt,0,atol=2e-16)

    # Under uniform scaling raw VV/VF/FF blocks have degrees 5/4/3,
    # while dB/B is -3/-2.  Every block of B.T@G@B is therefore degree -1.
    dGs,dBs=production_tet_charge_gram_derivatives(
        gram,cells[None,...],faces[None,...],cb["B"])
    dBs=dBs[0].toarray()
    np.testing.assert_allclose(dBs.T@G@B+B.T@dGs[0]@B+B.T@G@dBs,-N,
                               rtol=2e-11,atol=2e-13)


def test_native_production_hex_face_self_block_derivative_invariants():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
        mapping=lambda x,y,z:(x,y,z+.12*x*y))
    fes=ng.HDiv(mesh,order=1)
    with ng.TaskManager():
        cb=_charge_basis_hex(fes,cob_quad=3)
        _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=16,eta=2.0)
    face_nodes=np.asarray(cb["face_nodes"]).reshape(-1,9,3)
    host=4  # warped z-face: production value and derivative share radial-Duffy
    translation=gram.hex_face_self_block_directional_derivative(host,np.ones((9,3)))
    scaling=gram.hex_face_self_block_directional_derivative(host,face_nodes[host])
    offset=8+4*host
    value=np.array([[gram.entry(offset+i,offset+j) for j in range(4)] for i in range(4)])
    np.testing.assert_allclose(translation,0,atol=5e-17)
    np.testing.assert_allclose(scaling,scaling.T,rtol=0,atol=0)
    np.testing.assert_allclose(scaling,-value,rtol=3e-13,atol=3e-15)


def test_native_affine_hex_face_self_block_derivative_matches_production_value_and_fd():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram

    velocity_gradient=np.array([[.13,-.07,.04],[.02,.09,-.05],[-.03,.06,.11]])
    velocity_offset=np.array([.01,-.02,.03])
    epsilon=2e-6

    def build(scale):
        mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
            mapping=lambda x,y,z:tuple(np.array([x,y,z])+scale*epsilon*(
                velocity_gradient@np.array([x,y,z])+velocity_offset)))
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            cb=_charge_basis_hex(fes,cob_quad=3)
            # Keep +/-epsilon on the same dense numerical path; ACA factor
            # recompression is not a differentiable finite-difference oracle.
            _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
        return cb,gram

    cb,gram=build(0); _,plus=build(1); _,minus=build(-1)
    face_nodes=np.asarray(cb["face_nodes"]).reshape(-1,9,3)
    for host,nodes in enumerate(face_nodes):
        offset=8+4*host
        value=np.array([[gram.entry(offset+i,offset+j) for j in range(4)] for i in range(4)])
        scaling=gram.hex_face_self_block_directional_derivative(host,nodes)
        translation=gram.hex_face_self_block_directional_derivative(host,np.ones((9,3)))
        velocity=nodes@velocity_gradient.T+velocity_offset
        derivative=gram.hex_face_self_block_directional_derivative(host,velocity)
        value_plus=np.array([[plus.entry(offset+i,offset+j) for j in range(4)] for i in range(4)])
        value_minus=np.array([[minus.entry(offset+i,offset+j) for j in range(4)] for i in range(4)])
        validation_difference=(value_plus-value_minus)/(2*epsilon)
        np.testing.assert_allclose(scaling,-value,rtol=2e-13,atol=3e-15)
        np.testing.assert_allclose(translation,0,atol=5e-17)
        np.testing.assert_allclose(derivative,derivative.T,rtol=0,atol=0)
        np.testing.assert_allclose(derivative,validation_difference,rtol=2e-8,atol=1e-10)


def test_native_complete_hex_charge_gram_directional_derivative():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim._vim import _charge_basis_hex,build_charge_gram
    gradient=np.array([[.13,-.07,.04],[.02,.09,-.05],[-.03,.06,.11]])
    offset=np.array([.01,-.02,.03]); epsilon=2e-6
    # Fully generic Q2 warp keeps every closest-corner/near decision away
    # from symmetry ties, so the validation difference follows one fixed
    # production quadrature branch on both sides.
    def base(x,y,z): return np.array([
        x+.071*y*z+.013*x*y,
        y+.037*x*z+.009*y*z,
        z+.053*x*y+.011*x*z,
    ])
    def build(scale):
        mesh=MakeStructured3DMesh(hexes=True,nx=1,ny=1,nz=1,
            mapping=lambda x,y,z:tuple(base(x,y,z)+scale*epsilon*(gradient@base(x,y,z)+offset)))
        fes=ng.HDiv(mesh,order=1)
        with ng.TaskManager():
            cb=_charge_basis_hex(fes,cob_quad=3)
            # Keep this regression on a single dense leaf.  Differentiating
            # ACA recompression factors is outside the ChargeGram kernel API.
            _,gram,_=build_charge_gram(fes,eps=1e-10,leafsize=256,eta=2.0)
        return cb,gram
    cb,gram=build(0); _,plus=build(1); _,minus=build(-1)
    cells=np.asarray(cb["cell_nodes"]).reshape(-1,27,3)
    faces=np.asarray(cb["face_nodes"]).reshape(-1,9,3)
    zeros_c=np.zeros_like(cells); zeros_f=np.zeros_like(faces)
    translation=gram.hex_charge_gram_directional_derivative(
        zeros_c+np.array([.2,-.1,.3]),zeros_f+np.array([.2,-.1,.3]))
    scaling=gram.hex_charge_gram_directional_derivative(cells,faces)
    derivative=gram.hex_charge_gram_directional_derivative(cells@gradient.T+offset,faces@gradient.T+offset)
    derivative_operator=gram.directional_derivative_operator(
        "hex",cells@gradient.T+offset,faces@gradient.T+offset,
        eps=1e-12,leaf=256,eta=2.0)
    n=derivative.shape[0]
    value=np.array([[gram.entry(i,j) for j in range(n)] for i in range(n)])
    value_plus=np.array([[plus.entry(i,j) for j in range(n)] for i in range(n)])
    value_minus=np.array([[minus.entry(i,j) for j in range(n)] for i in range(n)])
    np.testing.assert_allclose(translation,0,atol=2e-15)
    np.testing.assert_allclose(scaling,-value,rtol=3e-10,atol=3e-13)
    np.testing.assert_allclose(derivative,derivative.T,rtol=0,atol=0)
    operator_dense=np.array([[derivative_operator.entry(i,j) for j in range(n)] for i in range(n)])
    np.testing.assert_allclose(operator_dense,derivative,rtol=2e-13,atol=2e-15)
    probe=np.linspace(-.7,.9,n)
    np.testing.assert_allclose(derivative_operator.matvec_sym(probe),derivative@probe,
        rtol=2e-12,atol=2e-13)
    np.testing.assert_allclose(derivative,(value_plus-value_minus)/(2*epsilon),rtol=3e-7,atol=2e-10)


def test_vim_linearization_matches_analytic_two_cell_system():
    A=np.array([[3.0,-1.0],[-1.0,2.0]])
    b=np.array([1.0,0.5]); C=np.array([[1.0,2.0]])
    dA=np.array([[[0.4,0.0],[0.0,0.0]],[[0.0,0.0],[0.0,0.3]]])
    result=linearize_vim_system(A,b,C,dA)
    epsilon=1e-7
    for cell in range(2):
        shifted=np.linalg.solve(A+epsilon*dA[cell],b)
        observed=(C@shifted-result.response)/epsilon
        np.testing.assert_allclose(observed,result.response_jacobian[:,cell],rtol=2e-6,atol=2e-8)


def test_laplace_pair_shape_derivative_matches_validation_difference():
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,2.,0.]])
    weights=np.array([1.,1.5,.75]); velocity=np.array([[[.1,0,0],[0,.2,0],[-.1,0,0]]])
    rel=np.array([[.03,-.02,.01]])
    gram,derivative=linearize_laplace_pair_gram(points,weights,velocity,rel)
    eps=1e-7
    moved,_=linearize_laplace_pair_gram(points+eps*velocity[0],weights*(1+eps*rel[0]),velocity*0)
    np.testing.assert_allclose((moved-gram)/eps,derivative[0],rtol=2e-6,atol=2e-9)


def test_vim_operator_product_rule_matches_validation_difference():
    M=np.array([[2.,.2],[.2,1.5]]); B=np.array([[1.,-1.],[.5,.25]])
    G=np.array([[.8,.1],[.1,.6]]); h=np.array([2.,-1.]); inv_chi=.1
    dM=np.array([[[.1,.02],[.02,-.03]]]); dB=np.array([[[.03,0],[-.01,.02]]]); dG=np.array([[[.04,.01],[.01,-.02]]])
    lin=linearize_vim_operator(M,B,G,h,inv_chi=inv_chi,dmass=dM,dcharge_map=dB,dcharge_gram=dG)
    eps=1e-7; shifted=linearize_vim_operator(M+eps*dM[0],B+eps*dB[0],G+eps*dG[0],h,
        inv_chi=inv_chi,dmass=np.zeros_like(dM),dcharge_map=np.zeros_like(dB),dcharge_gram=np.zeros_like(dG))
    np.testing.assert_allclose((shifted.matrix-lin.matrix)/eps,lin.matrix_jacobian[0],rtol=2e-6,atol=2e-9)
    np.testing.assert_allclose((shifted.rhs-lin.rhs)/eps,lin.rhs_jacobian[0],rtol=2e-6,atol=2e-9)


def test_matrix_free_vim_directional_action_matches_dense_product_rule():
    from radia.topology_optimization import linearize_vim_operator_matrix_free
    class SymOperator:
        def __init__(self,matrix): self.matrix=np.asarray(matrix)
        def matvec_sym(self,x): return self.matrix@x
    M=np.array([[2.,.2],[.2,1.3]])
    B=np.array([[1.,.3],[-.2,.8],[.4,-.1]])
    G=np.array([[1.2,.1,.05],[.1,.9,-.03],[.05,-.03,.7]])
    dM=np.array([[[.1,.02],[.02,-.03]]])
    dB=np.array([[[.03,0],[-.01,.02],[.04,-.02]]])
    dG=np.array([[[.04,.01,0],[.01,-.02,.03],[0,.03,.01]]])
    op=linearize_vim_operator_matrix_free(M,B,SymOperator(G),inv_chi=2.5,
        dmass=dM,dcharge_map=dB,dcharge_gram=(SymOperator(dG[0]),))
    x=np.array([.7,-.4])
    A=2.5*M+B.T@G@B
    dA=2.5*dM[0]+dB[0].T@G@B+B.T@dG[0]@B+B.T@G@dB[0]
    np.testing.assert_allclose(op.matvec(x),A@x,rtol=2e-15,atol=2e-15)
    np.testing.assert_allclose(op.directional_matvec(0,x),dA@x,rtol=2e-15,atol=2e-15)
    assert op.as_scipy_linear_operator(0).shape==(2,2)


def test_lp_update_obeys_volume_and_move_limit():
    result=solve_lp_update([0.5,0.5,0.5],[-3.0,-1.0,2.0],[1.0,1.0,1.0],1.5,move_limit=0.2)
    assert np.all(np.abs(result.delta)<=0.2+1e-12)
    assert np.sum(result.density)<=1.5+1e-12
    assert result.density[0]>=result.density[1]>=result.density[2]


def test_cubit_density_journal_has_deterministic_blocks(tmp_path):
    path=tmp_path/"density.jou"
    info=write_cubit_density_journal(path,[11,12,13,14],[0.9,0.1,0.7,0.2])
    text=path.read_text(encoding="utf-8")
    assert info["solid_count"]==2 and info["void_count"]==2
    assert "group 'radia_topopt_solid' add hex 11 13" in text
    assert "block 1002 hex in group 'radia_topopt_void'" in text


def test_sequential_vim_lp_reaches_volume_constrained_material_layout():
    def linearize(density):
        A=np.eye(2); b=np.array([1.0,1.0]); C=np.eye(2)
        result=linearize_vim_system(A,b,C,np.zeros((3,2,2)))
        return VIMLinearization(
            result.state,
            np.array([density[0]-density[2],density[1]]),
            result.state_jacobian,
            np.array([[-1.0,0.0,1.0],[0.0,1.0,0.0]]),
        )
    result=optimize_vim_lp([0.5]*3,[1.0]*3,0.5,linearize,objective_weights=[1.0,0.0],move_limit=0.25,max_iterations=5)
    assert np.sum(result.density)<=1.5+1e-12
    assert result.density[0]>=result.density[2]
