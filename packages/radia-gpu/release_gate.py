"""Qualify one successful main CI artifact for GPU alpha publication."""
import argparse
import json
from pathlib import Path
import re
from verify_wheel import verify


def qualify(run, repository, run_id, artifact, expected_hash):
    if (str(run.get('id')) != str(run_id) or run.get('status') != 'completed'
            or run.get('conclusion') != 'success' or run.get('event') != 'push'
            or run.get('head_branch') != 'main'
            or run.get('path') != '.github/workflows/radia-gpu.yml'
            or run.get('head_repository', {}).get('full_name') != repository):
        raise ValueError('Release requires the successful repository main GPU CI run')
    wheels = list(artifact.rglob('radia_gpu-*.whl'))
    if len(wheels) != 1:
        raise ValueError('Expected exactly one GPU wheel')
    wheel = verify(wheels[0])
    if wheel['source_commit'] != run['head_sha'] or wheel['wheel_sha256'] != expected_hash:
        raise ValueError('Release wheel hash or source differs from accepted CI')
    # Stable publication requires the separately requested independent review.
    if not re.fullmatch(r'\d+\.\d+\.\d+(?:a|b|rc)\d+', wheel['version']):
        raise ValueError('This publisher accepts prereleases only; stable needs independent review')
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
            or record.get('distribution_version') != wheel['version']
            or record.get('simulink_repeated_runs',0) < 2
            or record.get('simulink_failure_cleanup') is not True):
        raise ValueError('Installed acceptance does not certify the exact release wheel')
    return dict(**wheel, tag='radia-gpu-v'+wheel['version'], wheel=str(wheels[0]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-json',required=True,type=Path)
    parser.add_argument('--repository',required=True)
    parser.add_argument('--run-id',required=True)
    parser.add_argument('--artifact',required=True,type=Path)
    parser.add_argument('--wheel-sha256',required=True)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    result = qualify(json.loads(args.run_json.read_text()),args.repository,args.run_id,
                     args.artifact,args.wheel_sha256)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
