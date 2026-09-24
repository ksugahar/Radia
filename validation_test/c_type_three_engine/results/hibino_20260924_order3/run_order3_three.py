"""Coarse-mesh all-three order refinement diagnostic; no matched-error claim."""
import hashlib,json,sys,time,platform,faulthandler
from pathlib import Path
import numpy as np
import ngsolve as ng
import radia
import radia.kelvin_solver as kelvin
import run_three_engine as driver
from phase_memory import instrument
import psutil
import importlib.metadata
root=Path(__file__).resolve().parent
output=root/'cached-three.json'
events=[]
original=kelvin.project_source_total_hodge
original_interface=kelvin.project_source_interface_potential
original_hdiv_solve=driver.vim.Solve

def measured_hdiv_solve(*args,**kwargs):
    result=original_hdiv_solve(*args,**kwargs)
    event('hdiv_native_convergence',
          converged=result.get('last_solve_converged'),
          residual=result.get('last_solve_final_relative_residual'))
    if result.get('last_solve_converged') != 1:
        raise RuntimeError('HDiv native convergence was not confirmed')
    return result
def event(name,**data):
    item=dict(event=name,**data); events.append(item)
    print(json.dumps(item),flush=True)
    (root/'cached-three-events.json').write_text(json.dumps(events,indent=2))
def cached_hodge(mesh,source,materials,**kwargs):
    if kwargs.get('order')!=3 or kwargs.get('bonus_intorder')!=4:
        raise ValueError('This diagnostic requires source order 3, bonus 4')
    start=time.perf_counter()
    rule=ng.IntegrationRule(ng.ET.TET,10)
    mapped=mesh.MapToAllElements(rule,mesh.Materials('|'.join(materials)))
    points=np.asarray(ng.CF((ng.x,ng.y,ng.z))(mapped)).reshape(-1,3).tolist()
    event('source_cache_begin',points=len(points),mapping_s=time.perf_counter()-start)
    preparation=time.perf_counter()
    source.PrepareCache(points)
    event('source_cache_ready',preparation_s=time.perf_counter()-preparation,cache=source.GetCacheStats())
    result=original(mesh,source,materials,**kwargs)
    event('source_hodge_complete',inclusive_s=time.perf_counter()-start,cache=source.GetCacheStats(),
        auxiliary_space_role='iron_source_total_hodge',
        auxiliary_ndof=int(result['fes'].ndof),
        auxiliary_free_ndof=sum(bool(flag) for flag in result['fes'].FreeDofs()),
        auxiliary_count_scope='this projection only; not all auxiliary spaces or simultaneous system size')
    return result
def measured_interface(mesh,source,boundary,**kwargs):
    started=time.perf_counter()
    result=original_interface(mesh,source,boundary,**kwargs)
    event('source_interface_complete',runtime_s=time.perf_counter()-started,
        auxiliary_space_role='source_interface_potential',boundary=str(boundary),
        auxiliary_ndof=int(result['fes'].ndof),
        auxiliary_free_ndof=sum(bool(flag) for flag in result['fes'].FreeDofs()),
        auxiliary_count_scope='this projection only; not simultaneous system size')
    return result
if __name__=='__main__':
    if output.exists() or (root/'cached-three-events.json').exists():
        raise FileExistsError('Refusing to overwrite prior run')
    provenance=dict(host=platform.node(),python=sys.version,executable=sys.executable,
        versions={name:importlib.metadata.version(name) for name in ('radia','ngsolve','numpy')},
        logical_cpus=psutil.cpu_count(),memory=psutil.virtual_memory()._asdict(),
        radia_path=radia.__file__,ngsolve_path=ng.__file__,
        wrapper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_cache_is_inside_engine_runtime=True,resume=False,
        warning='Experimental implementation; no speed superiority without matched error and completed gates.')
    (root/'cached-three-provenance.json').write_text(json.dumps(provenance,indent=2))
    kelvin.project_source_total_hodge=cached_hodge
    driver.vim.Solve=measured_hdiv_solve
    kelvin.project_source_interface_potential=measured_interface
    process=psutil.Process()
    for name in ('solve_hdiv','solve_reduced_a','solve_omega'):
        setattr(driver,name,instrument(getattr(driver,name),event,lambda:process.memory_info().rss))
    sys.argv=[str(root/'run_three_engine.py'),'--mesh-dir',str(root/'meshes'),
        '--mode','linear','--hdiv-order','2','--fem-order','3','--threads','8',
        '--relative-rms-tolerance','0.01','--output',str(output)]
    faulthandler.dump_traceback_later(5400,exit=True)
    try:
        driver.main()
        event('completed')
    except BaseException as exc:
        event('failed',error=repr(exc)); raise
    finally:
        faulthandler.cancel_dump_traceback_later()
