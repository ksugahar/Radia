import hashlib, importlib, json, platform, sys
from pathlib import Path
import numpy as np
import run_three_engine as d
root=Path(__file__).resolve().parent
module=importlib.import_module('radia.vim._solve')
build=module.build_charge_gram
original=d.vim.Solve
adapter_original=d.solve_hdiv
rows=[]
captured={}
class DiagnosticDone(Exception): pass
def gram(*args,**kwargs):
    if intorder is not None: kwargs['intorder']=intorder
    result=build(*args,**kwargs)
    gram.last_timings=getattr(build,'last_timings',{})
    return result
def solve(*args,**kwargs):
    kwargs['tol']=1e-10
    result=original(*args,**kwargs)
    captured.clear()
    captured.update({k:v for k,v in result.items() if isinstance(v,(str,int,float,bool,type(None)))})
    return result
def adapter(*args,**kwargs):
    field,diag=adapter_original(*args,**kwargs)
    rows.append(dict(intorder=intorder,field_T=field.tolist(),diagnostics=diag,solve_scalars=dict(captured)))
    raise DiagnosticDone()
module.build_charge_gram=gram
d.vim.Solve=solve
d.solve_hdiv=adapter
for intorder in (None,9,13):
    sys.argv=['run_three_engine.py','--mesh-dir',str(root/'meshes'),'--mode','linear','--hdiv-order','1','--fem-order','2','--threads','8','--output',str(root/'not_three_engine.json')]
    try: d.main()
    except DiagnosticDone: pass
baseline=np.asarray(json.loads((root/'baseline.json').read_text())['fields_T']['hdiv_mmm'])
for row in rows:
    row['relative_vector_difference_from_previous_fine']=float(np.linalg.norm(np.asarray(row['field_T'])-baseline)/np.linalg.norm(baseline))
result=dict(scope='Experimental HDiv Gram quadrature sensitivity; not three-engine pass. intorder also feeds charge-basis assembly; default inner rule tracks outer rule.',host=platform.node(),python=sys.version,rows=rows,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),solve_source_sha256=hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest())
(root/'near_quad_comparison.json').write_text(json.dumps(result,indent=2))
print(json.dumps([dict(intorder=r['intorder'],difference=r['relative_vector_difference_from_previous_fine']) for r in rows]),flush=True)
