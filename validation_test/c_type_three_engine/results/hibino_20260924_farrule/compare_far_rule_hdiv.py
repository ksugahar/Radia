import json,sys,hashlib,platform
from pathlib import Path
import numpy as np
import run_three_engine as d
root=Path(__file__).resolve().parent
rows=[]
original_solve=d.vim.Solve
original_adapter=d.solve_hdiv
class DiagnosticDone(Exception): pass
captured={}
def solve(*args,**kwargs):
    kwargs['tol']=tolerance
    if rule=='all_near_rule': kwargs['ho_far_factor']=float('inf')
    result=original_solve(*args,**kwargs)
    captured.clear()
    captured.update({k:v for k,v in result.items() if isinstance(v,(str,int,float,bool,type(None)))})
    captured['result_keys']=list(result)
    captured['field_evaluator_stats']=result.get('field_evaluator_stats')
    return result
def adapter(*args,**kwargs):
    field,diag=original_adapter(*args,**kwargs)
    rows.append(dict(rule=rule,tolerance=tolerance,field_T=field.tolist(),diagnostics=diag,solve_scalars=dict(captured)))
    raise DiagnosticDone()
d.vim.Solve=solve
d.solve_hdiv=adapter
tolerance=1e-10
for rule in ('production','all_near_rule'):
    sys.argv=['run_three_engine.py','--mesh-dir',str(root/'meshes'),'--mode','linear','--hdiv-order','1','--fem-order','2','--threads','8','--output',str(root/'not_a_three_engine_result.json')]
    try: d.main()
    except DiagnosticDone: pass
baseline=np.array(json.loads((root/'baseline.json').read_text())['fields_T']['hdiv_mmm'])
for row in rows:
    row['relative_vector_difference_from_previous_fine']=float(np.linalg.norm(np.array(row['field_T'])-baseline)/np.linalg.norm(baseline))
result=dict(scope='HDiv-only replay diagnostic; not a three-engine validation pass',host=platform.node(),python=sys.version,rows=rows,
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(root/'far_rule_comparison.json').write_text(json.dumps(result,indent=2))
print(json.dumps([{'tol':r['tolerance'],'difference':r['relative_vector_difference_from_previous_fine'],'iterations':r['diagnostics']['linear_iterations']} for r in rows]),flush=True)
