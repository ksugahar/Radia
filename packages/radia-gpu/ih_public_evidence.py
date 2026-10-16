"""Export bounded IH qualification scalars; operational evidence stays private."""
import math
import re

TOP = {'schema','accepted','source_commit','wheel_sha256','manifest_sha256','binary_sha256',
       'matlab_version','toolkit_commit','toolkit_accepted','cpu','lifecycle','power'}
PRIVATE = {'loaded_mex','python','hostname','error','traceback','commands','input_path','model_path'}
CPU = {'accepted','relative_errors','true_relative_residual','scope','max_field_A_per_m'}
ERRORS = {'temperature','heat','impedance','electromagnetic'}
LIFECYCLE = {'repeat_identical','save_close_reopen','contexts_released','failed_candidate_preserved','reset'}


def fields(value, allowed):
    if not isinstance(value, dict) or set(value) - allowed:
        raise ValueError('Unknown qualification fields')


def scalar(value, limit=None):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or (limit is not None and value>limit):
        raise ValueError('Invalid qualification scalar')
    return value


def ih_public(record):
    fields(record, TOP | PRIVATE)
    if record.get('schema')!='radia.gpu.ih-installed-qualification.v1' or record.get('accepted') is not True:
        raise ValueError('Incomplete IH qualification')
    result={k:record[k] for k in ('schema','accepted')}
    for key,width in [('source_commit',40),('toolkit_commit',40),('wheel_sha256',64),('manifest_sha256',64),('binary_sha256',64)]:
        value=record.get(key)
        if not isinstance(value,str) or not re.fullmatch('[0-9a-f]{'+str(width)+'}',value):
            raise ValueError('Invalid qualification identity')
        result[key]=value
    version=record.get('matlab_version')
    if not isinstance(version,str) or not re.fullmatch(r'[0-9.]+ \(R[0-9]{4}[ab]\)( Update [0-9]+)?',version):
        raise ValueError('Invalid public MATLAB version')
    if record.get('toolkit_accepted') is not True:
        raise ValueError('Missing Toolkit qualification')
    result.update(matlab_version=version,toolkit_accepted=True)
    cpu=record.get('cpu');fields(cpu,CPU)
    errors=cpu.get('relative_errors');fields(errors,ERRORS)
    if cpu.get('accepted') is not True or set(errors)!=ERRORS:
        raise ValueError('Incomplete CPU qualification')
    result['cpu']={'accepted':True,'relative_errors':{k:scalar(errors[k],1e-8) for k in sorted(ERRORS)},'true_relative_residual':scalar(cpu.get('true_relative_residual'),1e-6)}
    life=record.get('lifecycle');fields(life,LIFECYCLE)
    if set(life)!=LIFECYCLE or any(value is not True for value in life.values()):
        raise ValueError('Incomplete lifecycle qualification')
    result['lifecycle']=dict(life)
    # B1 schema is intentionally not guessed: no public passed record until
    # the reviewed independent power/sign gate is connected by its owner.
    if record.get('power') is not None:
        raise ValueError('Power qualification schema pending reviewed B1 integration')
    result.update(power_qualification_pending=True,release_qualified=False,
                  timing_classification='numerical and lifecycle diagnostics; not performance acceptance')
    return result
