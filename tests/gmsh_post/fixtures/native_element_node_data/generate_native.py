"""Standalone Gmsh oracle for core; deliberately no MCP dependency."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import os
import sys

ROOT=Path(__file__).resolve().parent
out=ROOT/'inputs';out.mkdir(exist_ok=True)
provenance=json.loads(Path(os.environ['RADIA_GMSH_PROVENANCE']).read_text())
assert importlib.metadata.version('gmsh')==provenance['distribution_version']
assert provenance['distribution_version'].split('.')[0]=='5'
for filename,expected in provenance['runtime_payload_sha256'].items():
    assert hashlib.sha256(Path(filename).read_bytes()).hexdigest()==expected
import gmsh
result={'purpose':'Synthetic post-only core preflight fixtures; no solver input',
        'distribution_version':provenance['distribution_version'],
        'wheel_sha256':provenance['wheel_sha256'],'fixtures':[]}
for dim,etype,samples,name in [(0,15,1,'point_vector'),(1,1,2,'line_vector')]:
    gmsh.initialize(['-noconfig']);gmsh.option.setNumber('General.Verbosity',0)
    try:
        result['native_version']=gmsh.option.getString('General.Version')
        gmsh.option.setNumber('Mesh.MshFileVersion',4.2)
        gmsh.model.add(name);gmsh.model.addDiscreteEntity(dim,1)
        nodes=[17,29] if dim==0 else [17,29,41]
        gmsh.model.mesh.addNodes(dim,1,nodes,[c for i in range(len(nodes)) for c in (i,0,0)])
        gmsh.model.mesh.addElementsByType(1,etype,[11,23],nodes if dim==0 else [17,29,29,41])
        view=gmsh.view.add(name)
        data=[[10*(i+1)+j/8 for j in range(3*samples)] for i in range(2)]
        gmsh.view.addModelData(view,0,name,'ElementNodeData',[23,11],data,time=.25,numComponents=3)
        path=out/(name+'.msh');gmsh.view.write(view,str(path))
        raw=path.read_bytes();path.write_bytes(raw.replace(b'\r\n',b'\n'))
        gmsh.clear();gmsh.open(str(path));actual=gmsh.view.getModelData(gmsh.view.getTags()[-1],0)
        result['fixtures'].append(dict(file=path.name,samples_per_element=samples,
            element_type=etype,components=int(actual[4]),tags=[int(t) for t in actual[1]],
            values=[list(map(float,row)) for row in actual[2]],
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),native_sha256=hashlib.sha256(raw).hexdigest()))
    finally:gmsh.finalize()
assert not any(x=='radia_mcp' or x.startswith('radia_mcp.') for x in sys.modules)
result['mcp_imported']=False
(out/'native-oracle.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
