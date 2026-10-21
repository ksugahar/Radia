"""Qualify an exact main CI artifact and reviewed stable GPU evidence."""
import argparse
import json
import math
from pathlib import Path
import re
from verify_wheel import verify


def _repeated_runs_pass(value):
    return isinstance(value, int) and not isinstance(value, bool) and value >= 2


def _qualify_stable(wheel, qualification_directory):
    if qualification_directory is None:
        raise ValueError('Stable publication needs independent review and 100 acceptance evidence')
    path = qualification_directory / (wheel['version'] + '.json')
    if not path.is_file():
        raise ValueError('Stable publication needs independent review and 100 acceptance evidence')
    record = json.loads(path.read_text(encoding='utf-8'))
    if (record.get('schema') != 'radia.gpu.stable-qualification.v1'
            or any(record.get(key) != wheel[key]
                   for key in ('version', 'source_commit', 'wheel_sha256'))):
        raise ValueError('Stable qualification differs from the exact release wheel')
    review = record.get('independent_review', {})
    if (review.get('source_commit') != wheel['source_commit']
            or review.get('outcome') != 'accepted'
            or review.get('unresolved_findings') != []
            or not review.get('reviewer') or not review.get('report')):
        raise ValueError('Stable publication needs a completed independent review')
    acceptance = record.get('build_host_acceptance', {})
    if (acceptance.get('schema') != 'radia.gpu.installed-acceptance.v1'
            or acceptance.get('passed') is not True
            or acceptance.get('source_commit') != wheel['source_commit']
            or acceptance.get('wheel_sha256') != wheel['wheel_sha256']
            or acceptance.get('manifest_sha256') != wheel['manifest_sha256']
            or acceptance.get('distribution_version') != wheel['version']
            or acceptance.get('execution_role') != 'cuda-build-host'
            or acceptance.get('device') != 'Quadro RTX 5000'
            or acceptance.get('ltspice_checked') is not True
            or not _repeated_runs_pass(acceptance.get('simulink_repeated_runs'))
            or acceptance.get('simulink_failure_cleanup') is not True
            or acceptance.get('cpu_parity_checked') is not True):
        raise ValueError('Stable publication needs exact-wheel 100 and CPU parity acceptance')
    comparisons = acceptance.get('cpu_parity', {})
    solvers = wheel.get('cuda_solvers')
    if (not solvers or not isinstance(comparisons, dict)
            or set(comparisons) != set(solvers)):
        raise ValueError('Stable CPU parity must cover every published CUDA solver')
    for name, result in comparisons.items():
        # These modules apply operators; they do not solve Ax=b. Require an
        # explicit applicability declaration instead of a fabricated residual.
        operator_only = name in ('fmm', 'projection')
        if operator_only and (result.get('operation_kind') != 'operator'
                or result.get('residual_applicability') != 'not_applicable_operator'
                or 'true_relative_residual' not in result
                or result['true_relative_residual'] is not None):
            raise ValueError('Operator parity needs explicit residual applicability')
        if not operator_only and result.get('operation_kind', 'solve') != 'solve':
            raise ValueError('Solver parity requires an original-equation true residual')
        limits = [('relative_l2', 1e-8)]
        if not operator_only:
            limits.append(('true_relative_residual', 1e-6))
        for key, limit in limits:
            value = result.get(key)
            if (isinstance(value, bool) or not isinstance(value, (int, float))
                    or not math.isfinite(value) or value < 0 or value > limit):
                raise ValueError('Stable CPU parity or true residual is invalid')


def qualify(run, repository, run_id, artifact, expected_hash,
            qualification_directory=None):
    if (str(run.get('id')) != str(run_id) or run.get('status') != 'completed'
            or run.get('conclusion') != 'success' or run.get('event') not in ('push', 'workflow_dispatch')
            or run.get('head_branch') != 'main'
            or run.get('path') != '.github/workflows/radia-simulink.yml'
            or run.get('head_repository', {}).get('full_name') != repository):
        raise ValueError('Release requires the successful repository main GPU CI run')
    wheels = list(artifact.rglob('radia_simulink-*.whl'))
    if len(wheels) != 1:
        raise ValueError('Expected exactly one GPU wheel')
    wheel = verify(wheels[0])
    if wheel['source_commit'] != run['head_sha'] or wheel['wheel_sha256'] != expected_hash:
        raise ValueError('Release wheel hash or source differs from accepted CI')
    prerelease = re.fullmatch(r'\d+\.\d+\.\d+(?:a|b|rc)\d+', wheel['version']) is not None
    if not prerelease:
        if not re.fullmatch(r'\d+\.\d+\.\d+', wheel['version']):
            raise ValueError('Unsupported GPU release version')
        _qualify_stable(wheel, qualification_directory)
    acceptance = []
    for path in artifact.rglob('result.json'):
        record = json.loads(path.read_text())
        if record.get('schema') == 'radia.gpu.installed-acceptance.v1':
            acceptance.append(record)
    if len(acceptance) != 1:
        raise ValueError('Expected one installed GPU/Simulink/LTspice acceptance record')
    record = acceptance[0]
    if (record.get('passed') is not True or record.get('ltspice_checked') is not True
            or record.get('source_commit') != wheel['source_commit']
            or record.get('wheel_sha256') != expected_hash
            or record.get('manifest_sha256') != wheel['manifest_sha256']
            or record.get('distribution_version') != wheel['version']
            or not _repeated_runs_pass(record.get('simulink_repeated_runs'))
            or record.get('simulink_failure_cleanup') is not True):
        raise ValueError('Installed acceptance does not certify the exact release wheel')
    return dict(**wheel, tag='radia-simulink-v'+wheel['version'], wheel=str(wheels[0]),
                prerelease=prerelease)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-json',required=True,type=Path)
    parser.add_argument('--repository',required=True)
    parser.add_argument('--run-id',required=True)
    parser.add_argument('--artifact',required=True,type=Path)
    parser.add_argument('--wheel-sha256',required=True)
    parser.add_argument('--qualification-directory',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    result = qualify(json.loads(args.run_json.read_text()),args.repository,args.run_id,
                     args.artifact,args.wheel_sha256,args.qualification_directory)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
