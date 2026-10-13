"""Select public acceptance records; retain raw operational evidence on the runner."""
import argparse
import json
from pathlib import Path


def poisson_public(record):
    keys = ('schema', 'passed', 'platform', 'versions', 'ngsolve_upstream_cuda_enabled',
            'cpu_threads', 'taskmanager_used', 'matlab_startup_seconds', 'gateway_checks',
            'python_peak_rss_bytes', 'source_commit', 'limitations', 'cases')
    result = {key: record[key] for key in keys if key in record}
    build_keys = ('schema', 'source_sha256', 'binary_sha256', 'source_commit',
                  'source_dirty', 'matlab', 'release', 'gpu', 'python_version',
                  'no_cpu_solver_fallback', 'cuda_targets')
    if 'build' in record:
        result['build'] = {key: record['build'][key] for key in build_keys if key in record['build']}
    result['timing_classification'] = 'installation diagnostics; not performance acceptance'
    return result


def ams_public(record):
    keys = ('schema', 'accepted', 'platform', 'versions', 'threads', 'precision',
            'cases', 'source_commit', 'source_dirty', 'taskmanager_used', 'matlab_start_seconds',
            'gpu_binary_sha256', 'execution_mode')
    result = {key: record[key] for key in keys if key in record}
    if 'build' in record:
        result['build'] = {key: record['build'][key] for key in (
            'schema', 'source_commit', 'source_dirty', 'matlab', 'sha256',
            'source_sha256', 'no_cpu_solver_fallback') if key in record['build']}
    result['timing_classification'] = 'numerical acceptance diagnostics; not idle performance acceptance'
    return result


def setup_public(record):
    keys = ('schema','accepted','timing_classification','formal_performance_acceptance',
            'taskmanager_used','ngsolve_threads','source_commit','source_dirty',
            'versions','cases','matlab_start_seconds','parallel_contract_passed','memory_measurement')
    result = {key: record[key] for key in keys if key in record}
    result['build'] = ams_public({'build':record.get('build',{})})['build']
    return result


def nonlinear_public(record):
    keys = ('schema', 'accepted', 'source_commit', 'versions', 'gpu_binary_sha256',
            'dofs', 'free_dofs', 'elements', 'csr_nonzeros', 'mesh_sha256',
            'map_assembly_relative_error', 'cpu_threads', 'gpu_host_setup_threads',
            'nonlinear_tolerance', 'linear_true_tolerance', 'independent_parity_limit',
            'cpu_samples', 'gpu', 'gpu_device', 'resident_contract', 'assembly_storage',
            'gpu_residency', 'timing_boundary', 'near_roundoff_line_search',
            'shared_fe_setup_seconds', 'matlab_startup_seconds')
    result = {key: record[key] for key in keys if key in record}
    result['build'] = ams_public({'build': record.get('build', {})})['build']
    result['timing_classification'] = 'numerical installation diagnostics; not performance acceptance'
    return result


def fmm_public(record):
    result = {key: record[key] for key in ('schema', 'accepted', 'scope', 'metrics', 'report') if key in record}
    result['build'] = {key: record.get('build', {})[key] for key in
        ('source_commit', 'source_dirty', 'source_sha256', 'binary_sha256', 'taskmanager_used')
        if key in record.get('build', {})}
    result['timing_classification'] = 'numerical and lifecycle CI diagnostics; not performance acceptance'
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--poisson', required=True, type=Path)
    parser.add_argument('--ams', type=Path)
    parser.add_argument('--setup', type=Path)
    parser.add_argument('--nonlinear', type=Path)
    parser.add_argument('--fmm', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    for name in ('result.json', 'failed-result.json'):
        path = args.poisson / name
        if path.is_file():
            record = json.loads(path.read_text(encoding='utf-8'))
            if record.get('schema') != 'radia.ngsolve-cuda-poisson.v1':
                raise ValueError('Unexpected Poisson acceptance schema')
            (args.output / name).write_text(json.dumps(poisson_public(record), indent=2)+'\n', encoding='utf-8')
    if args.ams and (args.ams / 'result.json').is_file():
        record = json.loads((args.ams / 'result.json').read_text(encoding='utf-8'))
        if record.get('schema') != 'radia.ams-cuda-validation.v1':
            raise ValueError('Unexpected AMS acceptance schema')
        (args.output / 'ams-result.json').write_text(json.dumps(ams_public(record), indent=2)+'\n', encoding='utf-8')
    if args.setup and (args.setup / 'result.json').is_file():
        record = json.loads((args.setup / 'result.json').read_text(encoding='utf-8'))
        if record.get('schema') != 'radia.ams-setup-comparison.v1':
            raise ValueError('Unexpected AMS setup comparison schema')
        (args.output / 'ams-setup-result.json').write_text(json.dumps(setup_public(record),indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if args.nonlinear and (args.nonlinear / 'result.json').is_file():
        record = json.loads((args.nonlinear / 'result.json').read_text(encoding='utf-8'))
        if record.get('schema') != 'radia.ams-nonlinear-sparse-resident.v2':
            raise ValueError('Unexpected nonlinear AMS schema')
        (args.output / 'ams-nonlinear-result.json').write_text(json.dumps(nonlinear_public(record), indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if args.fmm and (args.fmm / 'result.json').is_file():
        record = json.loads((args.fmm / 'result.json').read_text(encoding='utf-8'))
        if record.get('schema') != 'radia.fmm-physical-adapters.v1':
            raise ValueError('Unexpected FMM acceptance schema')
        (args.output / 'fmm-result.json').write_text(json.dumps(fmm_public(record), indent=2, allow_nan=False)+'\n', encoding='utf-8')
