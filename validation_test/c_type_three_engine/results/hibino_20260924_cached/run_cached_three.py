"""Experimental exact-point source reuse. No solver/physics relaxation."""
import hashlib,json,sys,time,platform,faulthandler
from pathlib import Path
import numpy as np
import ngsolve as ng
import radia
import radia.kelvin_solver as kelvin
import run_three_engine as driver
root=Path(__file__).resolve().parent
output=root/'cached-three.json'
events=[]
original=kelvin.project_source_total_hodge
def event(name,**data):
    item=dict(event=name,**data); events.append(item)
    print(json.dumps(item),flush=True)
    (root/'cached-three-events.json').write_text(json.dumps(events,indent=2))
def cached_hodge(mesh,source,materials,**kwargs):
    if kwargs.get('order')!=2 or kwargs.get('bonus_intorder')!=4:
        raise ValueError('This diagnostic cache was checked only for order 2, bonus 4')
    start=time.perf_counter()
    rule=ng.IntegrationRule(ng.ET.TET,8)
    mapped=mesh.MapToAllElements(rule,mesh.Materials('|'.join(materials)))
    points=np.asarray(ng.CF((ng.x,ng.y,ng.z))(mapped)).reshape(-1,3).tolist()
    event('source_cache_begin',points=len(points),mapping_s=time.perf_counter()-start)
    preparation=time.perf_counter()
    source.PrepareCache(points)
    event('source_cache_ready',preparation_s=time.perf_counter()-preparation,cache=source.GetCacheStats())
    result=original(mesh,source,materials,**kwargs)
    event('source_hodge_complete',inclusive_s=time.perf_counter()-start,cache=source.GetCacheStats())
    return result
if __name__=='__main__':
    if output.exists() or (root/'cached-three-events.json').exists():
        raise FileExistsError('Refusing to overwrite prior run')
    provenance=dict(host=platform.node(),python=sys.version,executable=sys.executable,
        radia_path=radia.__file__,ngsolve_path=ng.__file__,
        wrapper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_cache_is_inside_engine_runtime=True,resume=False,
        warning='Experimental implementation; no speed superiority without matched error and completed gates.')
    (root/'cached-three-provenance.json').write_text(json.dumps(provenance,indent=2))
    kelvin.project_source_total_hodge=cached_hodge
    sys.argv=[str(root/'run_three_engine.py'),'--mesh-dir',str(root/'meshes'),
        '--mode','linear','--hdiv-order','1','--fem-order','2','--threads','8',
        '--relative-rms-tolerance','0.01','--output',str(output)]
    faulthandler.dump_traceback_later(2400,exit=True)
    try:
        driver.main()
        event('completed')
    except BaseException as exc:
        event('failed',error=repr(exc)); raise
    finally:
        faulthandler.cancel_dump_traceback_later()
