"""Fast saved-artifact checks; no notebook execution or solver imports."""
import hashlib
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('relative,count,field',[
    ('electrostatics/electrostatics',3,True),
    ('background_fields/background_fields',3,True),
    ('gmsh_animation/gmsh_animation',2,False),
    ('complex_coil_geometry/complex_coil',4,True),
    ('mesh_fusion/mesh_fusion',3,True),
])
def test_saved_visualization_supplement(relative,count,field):
    path=ROOT/'docs'/f'{relative}.ipynb'
    n=json.loads(path.read_bytes())
    evidence=json.loads(path.with_name(path.stem+'_webgui_results.json').read_bytes())
    assert evidence['host'].lower() in {'mdx1','mdx2'}
    assert evidence['notebook_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    historical=[c for c in n['cells'] if c['cell_type']=='code' and not c.get('id','').startswith('webgui-20260915-')]
    assert hashlib.sha256(json.dumps(historical,sort_keys=True).encode()).hexdigest()==evidence['historical_code_cells_sha256']
    cells=[c for c in n['cells'] if c['cell_type']=='code' and c.get('id','').startswith('webgui-20260915-')]
    state=n['metadata']['widgets']['application/vnd.jupyter.widget-state+json']['state']
    views=[]
    for cell in cells:
        assert cell['execution_count'] is not None
        for output in cell['outputs']:
            assert output['output_type']!='error'
            widget=output.get('data',{}).get('application/vnd.jupyter.widget-view+json')
            if widget:
                content=json.dumps(state[widget['model_id']]['state'])
                assert 'webgui' in content.lower()
                assert len(content)>1000
                views.append(widget)
    assert len(views)==count
    assert n['metadata']['radia']['webgui_required'] is True
    assert n['metadata']['radia']['webgui_field_required'] is field
    side=json.loads(path.with_name(path.stem+'_result.json').read_bytes())
    assert side['notebook_sha256']==evidence['notebook_sha256']
    assert side['visualization_execution']['host']==evidence['host']
    if relative.startswith('gmsh_animation'):
        for mesh in evidence['metrics']['meshes']:
            assert mesh['check_vol']['passed']
            assert hashlib.sha256(path.with_name(mesh['name']+'.vol').read_bytes()).hexdigest()==mesh['sha256']
