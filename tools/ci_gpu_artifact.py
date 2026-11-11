"""Seal and recover exact-SHA GPU wheels; compute roles never compile."""
import argparse
import hashlib
import io
import json
import os
import shutil
from pathlib import Path
import urllib.request
import zipfile


# Compute validation of the MATLAB payloads consumes wheels built on main pushes.
COMPUTE_WHEELS = [('radia-simulink.yml', 'gpu-simulink-wheel', '&event=push')]
# The compute bundle of one radia-simulink push run: its wheel and every native build output
# the GPU acceptance executes. A pull-request run builds and seals the merge commit, which is
# not its head SHA, so only push runs are fetched. The build host compiles; the A100 only runs it.
SIMULINK_BUNDLE = [('radia-simulink.yml', 'gpu-simulink-validation-bundle', '&event=push')]
BUNDLE = 'radia-simulink-gpu-bundle.zip'
# Native Release runs only for a release tag or a dispatch at the exact SHA, never for a PR.
RADIA_WHEELS = [('build-test.yml', 'gpu-radia-wheel', ''), ('build-test.yml', 'gpu-netgen-overlay', '')]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload(directory):
    """One wheel, the patched Netgen library of tools/netgen_overlay.py, or one compute bundle."""
    return list(directory.glob('*.whl')) + list(directory.glob('nglib.dll')) + list(directory.glob(BUNDLE))


def bundle(directory, root, temp, entries):
    """Pack checkout-relative ('checkout:<path>') and RUNNER_TEMP-relative ('temp:<name>') trees."""
    directory.mkdir(parents=True, exist_ok=False)
    bases = {'checkout': root.resolve(), 'temp': temp.resolve()}
    members = {}
    for entry in entries:
        kind, _, name = entry.partition(':')
        if kind not in bases or not name or Path(name).is_absolute() or '..' in Path(name).parts:
            raise ValueError(f'Invalid bundle entry: {entry}')
        source = (bases[kind] / name).resolve()
        if bases[kind] not in source.parents or not source.exists():
            raise ValueError(f'Missing or escaping bundle entry: {entry}')
        files = [source] if source.is_file() else sorted(p for p in source.rglob('*') if p.is_file())
        if not files:
            raise ValueError(f'Empty bundle entry: {entry}')
        for path in files:
            if path.is_symlink() or bases[kind] not in path.resolve().parents:
                raise ValueError(f'Linked or escaping bundle member: {path}')
            members[kind + '/' + path.relative_to(bases[kind]).as_posix()] = path
    with zipfile.ZipFile(directory / BUNDLE, 'w', zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(members):
            archive.write(members[name], name)


def unbundle(directory, sha, root, temp):
    """Verify the sealed bundle, then restore its trees; no member may escape its base."""
    verify(directory, sha)
    bases = {'checkout': root.resolve(), 'temp': temp.resolve()}
    with zipfile.ZipFile(directory / BUNDLE) as archive:
        for name in archive.namelist():
            kind, _, rest = name.partition('/')
            if kind not in bases or not rest or chr(92) in name or bases[kind] not in (bases[kind] / rest).resolve().parents:
                raise ValueError(f'Unexpected bundle member: {name}')
            target = (bases[kind] / rest).resolve()
            if target.exists():
                raise ValueError(f'Bundle member would overwrite: {name}')
        archive.extractall(directory / 'unpacked')
    for kind, base in bases.items():
        source = directory / 'unpacked' / kind
        if source.is_dir():
            for path in sorted(p for p in source.rglob('*') if p.is_file()):
                target = base / path.relative_to(source)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(path), str(target))


def seal(directory, sha):
    files = payload(directory)
    if len(files) != 1:
        raise ValueError('Expected exactly one wheel')
    (directory / 'gpu-artifact.json').write_text(json.dumps({
        'schema': 1, 'source_commit': sha,
        'files': {p.name: digest(p) for p in files}}, indent=2) + '\n')


def verify(directory, sha):
    receipt = json.loads((directory / 'gpu-artifact.json').read_text())
    if receipt['schema'] != 1 or receipt['source_commit'] != sha:
        raise ValueError('GPU artifact source mismatch')
    if set(receipt['files']) != {p.name for p in payload(directory)} or len(receipt['files']) != 1:
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
            if artifact['name'] in {name for _, name, _ in COMPUTE_WHEELS} and not artifact['expired'] and run.get('head_branch') == 'main':
                available.setdefault(run['head_sha'], set()).add(artifact['name'])
        page += 1
    candidates = {sha for sha, names in available.items() if len(names) == len(COMPUTE_WHEELS)}
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
    raise ValueError('No main SHA has the required successful sealed GPU artifacts')


def resolve(sha, wheels=COMPUTE_WHEELS):
    selected = []
    for workflow, artifact, event in wheels:
        runs = json.loads(request(f'actions/workflows/{workflow}/runs?head_sha={sha}{event}&status=success&per_page=100'))['workflow_runs']
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


def fetch(output, sha, wheels=COMPUTE_WHEELS):
    for run, found, artifact in resolve(sha, wheels):
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
    parser.add_argument('mode', choices=['seal', 'verify', 'fetch', 'select', 'bundle', 'unbundle'])
    parser.add_argument('directory', type=Path)
    parser.add_argument('--sha')
    parser.add_argument('--radia', action='store_true', help='Fetch the radia wheel instead of the MATLAB payloads')
    parser.add_argument('--simulink-bundle', action='store_true', help='Fetch the radia-simulink compute bundle')
    parser.add_argument('--entry', action='append', default=[], help="bundle: 'checkout:<path>' or 'temp:<name>'")
    args = parser.parse_args()
    if args.mode == 'select':
        sha = newest_main_sha()
        args.directory.write_text(sha + '\n')
        raise SystemExit(0)
    if not args.sha or len(args.sha) != 40 or any(c not in '0123456789abcdef' for c in args.sha):
        parser.error('Expected a full commit SHA')
    if args.mode == 'seal':
        seal(args.directory, args.sha)
    elif args.mode == 'verify':
        verify(args.directory, args.sha)
    elif args.mode == 'bundle':
        bundle(args.directory, Path.cwd(), Path(os.environ['RUNNER_TEMP']), args.entry)
        seal(args.directory, args.sha)
    elif args.mode == 'unbundle':
        unbundle(args.directory, args.sha, Path.cwd(), Path(os.environ['RUNNER_TEMP']))
    else:
        if args.radia and args.simulink_bundle:
            parser.error('Choose one artifact set')
        fetch(args.directory, args.sha, RADIA_WHEELS if args.radia else SIMULINK_BUNDLE if args.simulink_bundle else COMPUTE_WHEELS)
