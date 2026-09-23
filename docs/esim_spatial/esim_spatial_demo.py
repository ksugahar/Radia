"""Notebook-local reconstruction of P1 ESIM fields, not a solver/workbench."""
from pathlib import Path
import json
import hashlib
import numpy as np
import ngsolve as ng
from cubit_mesh_export.check import check_consistency


def load_fields(data_dir, sample_dir):
    """Use the solver's sorted boundary-vertex contract and NGSolve DOF API."""
    data_dir, sample_dir = Path(data_dir), Path(sample_dir)
    record = json.loads((data_dir / 'compute_record.json').read_text(encoding='utf-8'))
    for name, digest in record['inputs'].items():
        assert hashlib.sha256((sample_dir / name).read_bytes()).hexdigest() == digest, name
    vol = sample_dir / 'ih_bem_sample_p1.vol'
    # The BEM route uses only sibc; unused coarse air/coil regions have up to
    # 8.45% CAD discrepancy. Check them at 10%, but workpiece/sibc at 1%.
    report = check_consistency(str(vol), threshold=10.0)
    assert report['passed'], report
    for family, name in (('materials','workpiece'),('boundaries','sibc')):
        row=next(r for r in report[family] if r['name']==name)
        assert abs(row['error_pct']) < 1.0, row
    parent = ng.Mesh(str(vol))  # Preserve exported geometry; never Curve().
    import sys
    import radia
    panels = str(Path(radia.__file__).parent / 'panels')
    sys.path.insert(0, panels)
    try:
        from calc_inductance import _extract_bnd_only_inline
        mesh = _extract_bnd_only_inline(parent, 'sibc')
    finally:
        sys.path.remove(panels)
    vertices = list(range(mesh.nv))  # Exactly the solver's extracted ordering.
    fes = ng.H1(mesh, order=1)
    vertex_dofs = []
    for vnr in vertices:
        dofs = fes.GetDofNrs(ng.NodeId(ng.VERTEX, vnr))
        assert len(dofs) == 1
        vertex_dofs.append(dofs[0])
    vertex_dofs = np.asarray(vertex_dofs, dtype=int)
    cases = []
    for item in record['cases']:
        stem = f"f{item['frequency_hz']}"
        path = data_dir / (stem + '.json')
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item['result_sha256']
        result = json.loads(path.read_text(encoding='utf-8'))
        assert result['esim_converged'] and result['wp_basis_order'] == 1
        fields = []
        for key in ('esim_per_panel_H_t', 'esim_per_panel_Z_s_real'):
            values = np.asarray(result[key], dtype=float)
            assert len(values) == len(vertices) and np.isfinite(values).all()
            assert np.min(values) >= 0
            gf = ng.GridFunction(fes)
            gf.vec.FV().NumPy()[vertex_dofs] = values
            fields.append(gf)
        ht, resistance = fields
        # A local diagnostic from reconstructed accepted ESIM fields. It is
        # deliberately NOT relabelled as the solver's calibrated heat artifact.
        q_local = 0.5 * resistance * ht**2
        parent_fes = ng.H1(parent, order=1)
        parent_heat = ng.GridFunction(parent_fes)
        parent_heat.Load(str(data_dir / (stem + '_qsurf.sol')))
        heat = ng.GridFunction(fes)
        lookup = {tuple(round(float(x), 12) for x in v.point): v.nr for v in parent.vertices}
        origin_dofs = []
        for v in mesh.vertices:
            old = lookup[tuple(round(float(x), 12) for x in v.point)]
            origin_dofs.append(parent_fes.GetDofNrs(ng.NodeId(ng.VERTEX, old))[0])
        heat.vec.FV().NumPy()[vertex_dofs] = parent_heat.vec.FV().NumPy()[origin_dofs]
        p_heat = float(ng.Integrate(heat, mesh, ng.BND))
        p_local = float(ng.Integrate(q_local, mesh, ng.BND))
        assert np.isfinite([p_heat, p_local]).all() and min(p_heat, p_local) > 0
        cases.append(dict(frequency_hz=item['frequency_hz'], result=result,
                          Ht=ht, resistance=resistance, q_local=q_local, heat=heat,
                          P_heat_W=p_heat, P_local_W=p_local))
    return mesh, cases, report
