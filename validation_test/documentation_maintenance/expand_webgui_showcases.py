"""Build independent visualization supplements; never rerun historical campaigns.

prepare writes scratch notebooks under an explicit staging directory. execute
runs those notebooks on mdx. integrate preserves old cells, merges widget state,
and refreshes checksum sidecars with separate historical/new runtime provenance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import socket
import sys
from datetime import datetime, timezone

SPECS = {
    'complex_coil': {
        'path': 'docs/complex_coil_geometry/complex_coil.ipynb',
        'heading': '''## See a bent beam-steering coil and its magnetic field

Rotate the actual eight-segment CoilBuilder solid, then inspect the field beside
it. This is Radia's geometry-to-source-field capability; its MCP coil and
NGSolve tools own the operating workflow. The new field is the **finite-wire
centreline approximation** of that geometry at 1265 A, not a rerun of the
historical thick-conductor calculation below. The mesh is a sampling slab
200 mm from the coil plane along x, not an air-boundary FEM solve.
The scalar scene shows strength, and the vector scene shows direction.
''',
        'setup': '''import numpy as np
from coil_model import create_beam_steering_builder
from radia.biot_savart import h_segments_cf, h_segments_batch, MU0
from netgen.occ import Box, Pnt, OCCGeometry
builder=create_beam_steering_builder()
solid=builder.to_occ()
segments,current=builder.to_wire_segments(n_arc=60)
points=np.asarray(segments).reshape(-1,3)
lo,hi=points.min(axis=0),points.max(axis=0)
xview=float(hi[0]+0.2)
with ng.TaskManager():
    mesh=ng.Mesh(OCCGeometry(Box(Pnt(xview,lo[1],lo[2]), Pnt(xview+0.05,hi[1],hi[2]))).GenerateMesh(maxh=0.20))
    B=MU0*h_segments_cf(segments,current=current)
    probes=np.array([[xview+0.025,float(lo[1]+a*(hi[1]-lo[1])),float(lo[2]+b*(hi[2]-lo[2]))] for a in (0.1,0.5,0.9) for b in (0.1,0.5,0.9)])
    actual=np.array([B(mesh(*p)) for p in probes])
    reference=MU0*h_segments_batch(segments,probes,current=current)
    difference=float(np.linalg.norm(actual-reference)/np.linalg.norm(reference))
    assert np.isfinite(actual).all() and difference < 0.05, difference
metrics=dict(current_A=current,segments=len(segments),sampling_x_m=xview,
             relative_cf_vs_finite_segment_gap=difference,field_route='centreline Biot-Savart, not thick conductor',
             geometry_source_sha256=hashlib.sha256(Path('coil_model.py').read_bytes()).hexdigest(),
             geometry_source=Path('coil_model.py').read_text(encoding='utf-8'))
print(json.dumps({k:v for k,v in metrics.items() if k!='geometry_source'},indent=2))
''',
        'scenes': [
            ('Bent coil: real CoilBuilder CAD', "from netgen.webgui import Draw as DrawGeometry\nscene=DrawGeometry(solid,width='100%',height='480px')"),
            ('Field sampling slab mesh', "scene=Draw(mesh,name='Field_sampling_slab',draw_vol=False,draw_surf=True,width='100%',height='480px')"),
            ('Magnetic flux density magnitude (tesla)', "with ng.TaskManager():\n    scene=Draw(ng.Norm(B),mesh,name='Centreline_B_magnitude_T',draw_vol=False,draw_surf=True,autoscale=True,clipping=None,width='100%',height='480px')"),
            ('Magnetic flux density direction', "with ng.TaskManager():\n    scene=Draw(B,mesh,name='Centreline_B_vectors_T',draw_vol=False,draw_surf=True,vectors={'grid_size':12},autoscale=True,clipping=None,width='100%',height='480px')"),
        ],
        'field': True,
    },
    'mesh_fusion': {
        'path': 'docs/mesh_fusion/mesh_fusion.ipynb',
        'heading': '''## See the two-domain coupling, solution and error

This fresh small solve reuses the notebook's Nitsche interface formulation.
Inspect the two materials, the combined solution and its absolute error against
the sine reference. Separate H1 spaces are coupled on a **conforming geometric
interface**; this scene does not demonstrate arbitrary nonmatching meshes.
The source and coupling study below remain unchanged. Radia MCP's NGSolve
method guidance owns operating instructions; these views make the coupling
result visible rather than supplying a second manual.
''',
        'setup': '''import ast
source_notebook=Path('mesh_fusion.ipynb')
original=json.loads(source_notebook.read_bytes())
source=next(''.join(c['source']) for c in original['cells'] if c['cell_type']=='code' and not c.get('id','').startswith('webgui-') and 'def solve_nitsche_mortar(' in ''.join(c['source']))
tree=ast.parse(source)
# Keep the existing definitions, not the historical full sweep's main call.
tree.body=[node for node in tree.body if not isinstance(node,ast.If)]
namespace={'__name__':'mesh_fusion_preview'}
exec(compile(tree,str(source_notebook),'exec'),namespace)
with ng.TaskManager():
    mesh=namespace['build_split_mesh'](0.1,0.1)
    fes,solution=namespace['solve_nitsche_mortar'](mesh,order=2)
    left,right=solution.components
    combined=mesh.MaterialCF({'left':left,'right':right})
    exact=ng.sin(math.pi*ng.x)*ng.sin(math.pi*ng.y)/(2*math.pi**2)
    absolute_error=ng.sqrt((combined-exact)**2)
    errors=namespace['l2_error_combined'](mesh,solution)
    jump=namespace['interface_jump'](mesh,solution)
    assert errors[2] < 2e-5 and jump < 2e-6,(errors,jump)
metrics=dict(ndof=fes.ndof,l2_error=errors[2],interface_jump=jump,
             source_cell_sha256=hashlib.sha256(source.encode()).hexdigest(),
             source_cell=source,geometric_interface='conforming',order=2,maxh=0.1)
print(json.dumps({k:v for k,v in metrics.items() if k!='source_cell'},indent=2))
''',
        'scenes': [
            ('Mesh and two coupled subdomains', "scene=Draw(mesh,name='Nitsche_two_subdomains',draw_vol=True,draw_surf=True,width='100%',height='480px')"),
            ('Coupled finite-element solution', "with ng.TaskManager():\n    scene=Draw(combined,mesh,name='Nitsche_solution',draw_vol=True,draw_surf=True,autoscale=True,width='100%',height='480px')"),
            ('Absolute error against the sine solution', "with ng.TaskManager():\n    scene=Draw(absolute_error,mesh,name='Nitsche_absolute_error',draw_vol=True,draw_surf=True,autoscale=True,width='100%',height='480px')"),
        ],
        'field': True,
    },
    'electrostatics': {
        'path': 'docs/electrostatics/electrostatics.ipynb',
        'heading': r'''## See the ideal coaxial electric field

Radia's electrostatic reference formulas describe capacitance and force;
this view makes their underlying field visible. Radia MCP's electrostatics
and force tools in `radia-analysis` own the live workflow. The figures below
are **analytical reference fields**, not a new FEM or force validation run.

For inner radius $a=0.005$ m, outer radius $b=0.015$ m, and voltage
$V_0=100$ V, the infinite-coax solution [@jackson1998classical] is

$$V(r)=V_0\frac{\ln(b/r)}{\ln(b/a)},\qquad
\mathbf E=\frac{V_0}{\ln(b/a)}\frac{(x,y,0)}{x^2+y^2}.$$

The annular viewing mesh shows a 0.04 m axial slice. Its end faces are
visualization cuts: no finite-length fringing field is modeled. The mesh is
generated with Netgen OCC, **not Cubit Mesh Export**. Compare the voltage
drop with radial arrows and the stronger field near the inner electrode.
Coordinates are metres; voltage and electric field use V and V/m.
The existing numerical-validation cells below retain their original outputs.
''',
        'setup': '''from netgen.occ import Cylinder, Pnt, Dir, OCCGeometry
a, b, length, voltage = 0.005, 0.015, 0.04, 100.0
outer = Cylinder(Pnt(0,0,-length/2), Dir(0,0,1), b, length)
inner = Cylinder(Pnt(0,0,-length/2), Dir(0,0,1), a, length)
domain = outer-inner
domain.mat('dielectric_view')
mesh = ng.Mesh(OCCGeometry(domain).GenerateMesh(maxh=0.004))
mesh.Curve(3)  # Here the mesh has an in-memory CAD geometry; not a loaded .vol.
r2 = ng.x**2+ng.y**2
potential = voltage*ng.log(b/ng.sqrt(r2))/math.log(b/a)
electric = voltage/math.log(b/a)*ng.CF((ng.x/r2, ng.y/r2, 0))
checks=[]
for radius in [0.006,0.010,0.014]:
    actual=float(potential(mesh(radius,0,0)))
    expected=voltage*math.log(b/radius)/math.log(b/a)
    assert math.isclose(actual,expected,rel_tol=1e-12)
    checks.append(dict(radius_m=radius,potential_V=actual))
metrics=dict(reference='infinite coax, analytical; no fringing', samples=checks,
             n_elements=mesh.ne, materials=list(mesh.GetMaterials()))
print(json.dumps(metrics,indent=2))
''',
        'scenes': [
            ('Dielectric viewing mesh', "scene = Draw(mesh, name='Coax_dielectric_mesh', order=3, draw_vol=True, draw_surf=True, width='100%', height='480px')"),
            ('Potential: 100 V at the inner electrode, 0 V at the outer', "scene = Draw(potential, mesh, name='Coax_potential_V', order=3, draw_vol=True, draw_surf=True, autoscale=False, min=0, max=100, clipping={'x':0,'y':0,'z':1,'dist':0}, vectors=False, width='100%', height='480px')"),
            ('Electric field: outward radial direction and near-electrode concentration', "scene = Draw(electric, mesh, name='Coax_E_V_per_m', order=3, draw_vol=True, draw_surf=True, autoscale=True, clipping={'x':0,'y':0,'z':1,'dist':0}, vectors={'grid_size':20}, width='100%', height='480px')"),
        ],
        'field': True,
    },
    'background_fields': {
        'path': 'docs/background_fields/background_fields.ipynb',
        'heading': r'''## See the prescribed quadrupole before adding iron

The background-field callback is easier to understand when its direction is
visible. Radia MCP's background-field workflow in `radia-analysis` owns the
operating contract. This supplement visualizes **only the prescribed source**
$\mathbf B_s=(gy,gx,0)$ with $g=10$ T/m, matching the callback below.
It does not represent the nonlinear iron response or a newly solved magnet.

The source is both divergence-free and curl-free in this current-free viewing
region, and $|\mathbf B_s|=g\sqrt{x^2+y^2}$. The arrows reveal the quadrupole
orientation; the scalar view reveals the zero on the axis. A 20 mm cube
is only a sampling domain, not a physical magnetic outer boundary or iron
body. Netgen OCC generates this viewing mesh; it is not a Cubit export.
The original material/relaxation campaign outputs below remain unchanged.
''',
        'setup': '''from netgen.occ import Box, Pnt, OCCGeometry
domain=Box(Pnt(-0.01,-0.01,-0.01), Pnt(0.01,0.01,0.01))
domain.mat('source_view_not_iron')
mesh=ng.Mesh(OCCGeometry(domain).GenerateMesh(maxh=0.005))
gradient=10.0
source=ng.CF((gradient*ng.y, gradient*ng.x, 0))
magnitude=ng.Norm(source)
checks=[]
for point in [(0,0,0),(0.005,0,0),(0,0.005,0),(0.005,0.005,0)]:
    actual=list(source(mesh(*point)))
    expected=[gradient*point[1],gradient*point[0],0]
    assert all(math.isclose(x,y,abs_tol=1e-12) for x,y in zip(actual,expected))
    checks.append(dict(point_m=point,B_T=actual))
metrics=dict(reference='prescribed source only, no material response',
             samples=checks,n_elements=mesh.ne,gradient_T_per_m=gradient)
print(json.dumps(metrics,indent=2))
''',
        'scenes': [
            ('Sampling mesh, not an iron body', "scene = Draw(mesh, name='Quadrupole_sampling_mesh', draw_vol=True, draw_surf=True, width='100%', height='480px')"),
            ('Source-field magnitude and its axial zero', "scene = Draw(magnitude, mesh, name='Prescribed_B_magnitude_T', order=2, draw_vol=True, draw_surf=True, autoscale=False, min=0, max=0.142, clipping={'x':0,'y':0,'z':1,'dist':0}, vectors=False, width='100%', height='480px')"),
            ('Source-field direction', "scene = Draw(source, mesh, name='Prescribed_B_vector_T', order=2, draw_vol=True, draw_surf=True, autoscale=True, clipping={'x':0,'y':0,'z':1,'dist':0}, vectors={'grid_size':15}, width='100%', height='480px')"),
        ],
        'field': True,
    },
    'gmsh_animation': {
        'path': 'docs/gmsh_animation/gmsh_animation.ipynb',
        'heading': '''## Inspect the input meshes before reading displacement frames

What body is being animated? These saved WebGUI scenes reload the adjacent
`rotor.vol` and `stator.vol`, expose their element and label inventories, and
show their separate geometries. Radia MCP's GMSH workflow owns animation
and post-processing instructions; `radia-cubit` owns export guidance.

These are the existing reference meshes, **not a new Cubit export**. The
structural/quality `check-vol` gate runs before display and the reports are
saved. A missing CAD sidecar or application-specific label contract is not
silently replaced by an invented one. Each view is independently centered:
they are not a registered rotor/stator assembly or an animation. The GMSH
artifact below remains the time-dependent displacement evidence. No magnetic
field or structural-dynamics solution is implied by these mesh views.
''',
        'setup': '''from collections import Counter
from cubit_mesh_export.check import check_consistency
meshes={}
rows=[]
for name in ['rotor','stator']:
    path=Path(name+'.vol')
    report=check_consistency(str(path))
    assert report['passed'], report
    mesh=ng.Mesh(str(path))
    meshes[name]=mesh
    rows.append(dict(name=name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                     elements=dict(Counter(str(e.type) for e in mesh.Elements(ng.VOL))),
                     materials=list(mesh.GetMaterials()),boundaries=list(mesh.GetBoundaries()),
                     check_vol=report))
metrics=dict(meshes=rows)
print(json.dumps(metrics,indent=2))
''',
        'scenes': [
            ('Rotor input mesh', "scene = Draw(meshes['rotor'], name='Rotor_input_mesh', order=3, draw_vol=True, draw_surf=True, width='100%', height='480px')"),
            ('Stator reference mesh', "scene = Draw(meshes['stator'], name='Stator_reference_mesh', order=3, draw_vol=True, draw_surf=True, width='100%', height='480px')"),
        ],
        'field': False,
    },
}


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def prepare(root, stage):
    import shutil
    import nbformat as nbf
    for key, spec in SPECS.items():
        directory=stage/key
        directory.mkdir(parents=True,exist_ok=True)
        base=root/spec['path']
        cells=[nbf.v4.new_markdown_cell(spec['heading'])]
        common="""from pathlib import Path
import hashlib, json, math, platform, socket, sys
from datetime import datetime, timezone
from importlib.metadata import version
import ngsolve as ng
from ngsolve.webgui import Draw
ng.SetNumThreads(2)
"""
        cells.append(nbf.v4.new_code_cell(common+spec['setup']))
        for i,(heading,source) in enumerate(spec['scenes']):
            cells.append(nbf.v4.new_markdown_cell('### '+heading))
            cells.append(nbf.v4.new_code_cell(source+f"\nscene.GenerateHTML(filename='scene_{i}.html')\n"))
        cells.append(nbf.v4.new_code_cell("""evidence=dict(schema='radia.docs.webgui_supplement.v1',
    generated_at_utc=datetime.now(timezone.utc).isoformat(),host=socket.gethostname(),
    python_executable=sys.executable,python_version=platform.python_version(),
    versions={p:version(p) for p in ['ngsolve','netgen-mesher','anywidget','cubit-mesh-export']},
    metrics=metrics,threads=2)
Path('evidence.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
"""))
        for i,c in enumerate(cells):
            c.id=f'webgui-20260915-{key}-{i}'
        n=nbf.v4.new_notebook(cells=cells,metadata=dict(kernelspec=dict(name='python3',display_name='Python 3')))
        nbf.write(n,directory/'supplement.ipynb')
        if key=='complex_coil':
            shutil.copy2(base.parent/'coil_model.py',directory/'coil_model.py')
        if key=='mesh_fusion':
            shutil.copy2(base,directory/base.name)
        if key=='gmsh_animation':
            for name in ['rotor.vol','stator.vol']:
                shutil.copy2(base.parent/name,directory/name)
                cad=base.parent/(name+'.json')
                if cad.exists(): shutil.copy2(cad,directory/cad.name)


def execute(stage):
    import nbformat
    from nbclient import NotebookClient
    os.environ['JUPYTER_PATH']=str(Path(sys.prefix)/'share/jupyter')
    for key in SPECS:
        p=stage/key/'supplement.ipynb'
        n=nbformat.read(p,as_version=4)
        print('Execute',key,flush=True)
        NotebookClient(n,timeout=180,kernel_name='python3',resources={'metadata':{'path':str(p.parent)}}).execute()
        nbformat.write(n,p)
        print('Saved',key,flush=True)


def integrate(root,stage):
    sys.path.insert(0,str(root/'packages/radia-mcp/src'))
    from radia_mcp.document_meta.tools import document_meta_write_notebook_result_json
    for key,spec in SPECS.items():
        p=root/spec['path']
        n=json.loads(p.read_bytes())
        baseline=digest(p)
        original_codes=[c for c in n['cells'] if c['cell_type']=='code']
        supplement=json.loads((stage/key/'supplement.ipynb').read_bytes())
        evidence=json.loads((stage/key/'evidence.json').read_bytes())
        # Display saved mesh/field data first, retaining all historical cells unchanged.
        assert not any(c.get('id','').startswith('webgui-20260915-') for c in n['cells']),p
        n['cells'][1:1]=supplement['cells']
        widgets=n['metadata'].setdefault('widgets',{})
        for mime,payload in supplement['metadata']['widgets'].items():
            if mime not in widgets:
                widgets[mime]=payload
            else:
                incoming=payload.get('state',{})
                existing=widgets[mime].setdefault('state',{})
                assert not (incoming.keys() & existing.keys()), 'Widget ID collision'
                existing.update(incoming)
        n['metadata'].setdefault('radia',{}).update(notebook_role='example',webgui_required=True,
            webgui_field_required=spec['field'],visualization_supplement=dict(date='2026-09-15',
                host=evidence['host'],scope='Only added visualization cells executed; historical cells/outputs unchanged'))
        preserved=[c for c in n['cells'] if c['cell_type']=='code' and not c.get('id','').startswith('webgui-20260915-')]
        assert preserved==original_codes
        p.write_text(json.dumps(n,ensure_ascii=False,indent=1)+'\n',encoding='utf-8')
        evidence.update(notebook=spec['path'],baseline_notebook_sha256=baseline,
            notebook_sha256=digest(p),historical_code_cells_preserved=True,
            historical_code_cells_sha256=hashlib.sha256(json.dumps(original_codes,sort_keys=True).encode()).hexdigest(),
            source_cells=[dict(id=c['id'],source=''.join(c['source'])) for c in supplement['cells'] if c['cell_type']=='code'],
            supplement_sha256=digest(stage/key/'supplement.ipynb'),
            visual_qa='Saved rich output and state audited; rendered visual QA pending.')
        ep=p.with_name(p.stem+'_webgui_results.json')
        ep.write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        side=p.with_name(p.stem+'_result.json')
        old=json.loads(side.read_bytes()) if side.exists() else {}
        r=document_meta_write_notebook_result_json(str(p),overwrite=True)
        assert 'error' not in r,r
        record=json.loads(side.read_bytes())
        record['sidecar_writer_versions']=record['versions']
        record['versions']=old.get('versions',{})
        record['historical_execution']=dict(generated_at_utc=old.get('generated_at_utc'),versions=old.get('versions',{}))
        record['visualization_execution']=evidence
        record['execution_scope']='Mixed provenance: original campaigns preserved; only added visualization cells executed on mdx.'
        side.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print('Integrated',spec['path'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['prepare','execute','integrate'])
    parser.add_argument('--root',type=Path,default=Path.cwd())
    parser.add_argument('--stage',type=Path,required=True)
    parser.add_argument('--only',nargs='+',choices=list(SPECS))
    args=parser.parse_args()
    if args.only:
        SPECS={key:SPECS[key] for key in args.only}
    if args.mode=='prepare': prepare(args.root,args.stage)
    elif args.mode=='execute': execute(args.stage)
    else: integrate(args.root,args.stage)
