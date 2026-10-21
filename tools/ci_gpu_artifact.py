"""Seal and recover exact-SHA GPU wheels; compute roles never compile."""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import urllib.request
import zipfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seal(directory, sha):
    wheels = list(directory.glob('*.whl'))
    if len(wheels) != 1:
        raise ValueError('Expected exactly one wheel')
    (directory / 'gpu-artifact.json').write_text(json.dumps({
        'schema': 1, 'source_commit': sha,
        'files': {p.name: digest(p) for p in wheels}}, indent=2) + '\n')


def verify(directory, sha):
    receipt = json.loads((directory / 'gpu-artifact.json').read_text())
    if receipt['schema'] != 1 or receipt['source_commit'] != sha:
        raise ValueError('GPU artifact source mismatch')
    if set(receipt['files']) != {p.name for p in directory.glob('*.whl')} or len(receipt['files']) != 1:
        raise ValueError('GPU wheel inventory mismatch')
    for name, expected in receipt['files'].items():
        path = directory / name
        if Path(name).name != name or digest(path) != expected:
            raise ValueError('GPU artifact hash mismatch')


class SafeRedirect(urllib.request.HTTPRedirectHandler):
    """Storage redirects must never receive the repository credential."""
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        redirected = super().redirect_request(req, fp, code, msg, headers, newurl)
        if redirected is not None:
            redirected.remove_header('Authorization')
        return redirected


def request(route):
    repository = os.environ['GITHUB_REPOSITORY']
    req = urllib.request.Request('https://api.github.com/repos/' + repository + '/' + route,
        headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'],
                 'Accept': 'application/vnd.github+json'})
    with urllib.request.build_opener(SafeRedirect()).open(req, timeout=60) as response:
        return response.read()


def newest_main_sha():
    # Artifact listing is newest-first and paginated; only unexpired push/main
    # artifacts qualify. Commit order below determines newest main ancestry.
    available = {}
    page = 1
    while True:
        artifacts = json.loads(request(f'actions/artifacts?per_page=100&page={page}'))['artifacts']
        if not artifacts:
            break
        for artifact in artifacts:
            run = artifact.get('workflow_run', {})
            if artifact['name'] in {'gpu-simulink-wheel', 'gpu-optuna-wheel'} and not artifact['expired'] and run.get('head_branch') == 'main':
                available.setdefault(run['head_sha'], set()).add(artifact['name'])
        page += 1
    candidates = {sha for sha, names in available.items() if len(names) == 2}
    page = 1
    while candidates:
        commits = json.loads(request(f'commits?sha=main&per_page=100&page={page}'))
        if not commits:
            break
        for commit in commits:
            if commit['sha'] in candidates:
                try:
                    resolve(commit['sha'])
                except ValueError:
                    candidates.remove(commit['sha'])
                else:
                    return commit['sha']
        page += 1
    raise ValueError('No main SHA has both successful sealed GPU artifacts')


def resolve(sha):
    selected = []
    for workflow, artifact in [('radia-simulink.yml', 'gpu-simulink-wheel'),
                               ('radia-optuna.yml', 'gpu-optuna-wheel')]:
        runs = json.loads(request(f'actions/workflows/{workflow}/runs?head_sha={sha}&event=push&status=success&per_page=100'))['workflow_runs']
        found = None
        for run in runs:
            if run['head_sha'] != sha:
                continue
            artifacts = json.loads(request(f"actions/runs/{run['id']}/artifacts?per_page=100"))['artifacts']
            found = next((a for a in artifacts if a['name'] == artifact and not a['expired']), None)
            if found:
                break
        if found is None:
            raise ValueError(f'Missing successful exact-SHA artifact: {artifact}')
        selected.append((run, found, artifact))
    return selected


def fetch(output, sha):
    for run, found, artifact in resolve(sha):
        archive = request(f"actions/artifacts/{found['id']}/zip")
        target = output / artifact
        target.mkdir(parents=True, exist_ok=False)
        with zipfile.ZipFile(io.BytesIO(archive)) as bundle:
            if any(Path(n).name != n or '\\' in n for n in bundle.namelist()):
                raise ValueError('Unexpected artifact member')
            bundle.extractall(target)
        verify(target, sha)
        (target / 'download.json').write_text(json.dumps({'run_id': run['id'],
            'artifact_id': found['id'], 'source_commit': sha,
            'archive_sha256': hashlib.sha256(archive).hexdigest()}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['seal', 'fetch', 'select'])
    parser.add_argument('directory', type=Path)
    parser.add_argument('--sha')
    args = parser.parse_args()
    if args.mode == 'select':
        sha = newest_main_sha()
        args.directory.write_text(sha + '\n')
        raise SystemExit(0)
    if not args.sha or len(args.sha) != 40 or any(c not in '0123456789abcdef' for c in args.sha):
        parser.error('Expected a full commit SHA')
    (seal if args.mode == 'seal' else fetch)(args.directory, args.sha)
