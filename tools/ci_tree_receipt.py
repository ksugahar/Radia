"""Fail-closed reuse of successful PR runs for an identical merged tree.

Receipts save runner time, never authorize release artifacts. Manual dispatch
always executes the lane and produces fresh exact-commit release artifacts.
"""
import io
import json
import os
from pathlib import Path
import urllib.request
import urllib.error
import urllib.parse
import zipfile

LANES = {
    'sparsesolv': ['ams-regression'],
    'radia-optuna': ['build-test', 'radia-optuna installed-wheel MATLAB/Simulink E2E'],
    'radia-simulink': ['cuda-mex'],
}


class ArtifactRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, message, headers, new_url):
        redirected = super().redirect_request(request, fp, code, message, headers, new_url)
        if urllib.parse.urlsplit(request.full_url).netloc != urllib.parse.urlsplit(new_url).netloc:
            redirected.remove_header('Authorization')
        return redirected


def api(path, binary=False):
    request = urllib.request.Request('https://api.github.com/repos/' + os.environ['GITHUB_REPOSITORY'] + '/' + path,
        headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'],
                 'Accept': 'application/vnd.github+json', 'X-GitHub-Api-Version': '2022-11-28'})
    with urllib.request.build_opener(ArtifactRedirect).open(request, timeout=20) as response:
        body = response.read(1024 * 1024 + 1) if binary else response.read()
    if binary and len(body) > 1024 * 1024:
        raise ValueError('oversized identity artifact')
    return body if binary else json.loads(body)


def pages(path, key=None):
    result = []
    for page in range(1, 21):
        value = api(path + ('&' if '?' in path else '?') + f'per_page=100&page={page}')
        batch = value[key] if key else value
        result.extend(batch)
        if len(batch) < 100:
            return result
    raise ValueError('pagination limit; execute instead')


def numerical_test_set(root):
    profiles = json.loads((root / 'tests/test_tier_manifest.json').read_text(encoding='utf-8'))['profiles']
    selected = {}
    for name in ('native-build-contracts', 'sparsesolv', 'solver-numerics', 'solver-heavy'):
        paths = set()
        current = name
        visited = set()
        while current:
            if current in visited:
                raise ValueError('cyclic numerical profile inheritance')
            visited.add(current)
            paths.update(profiles[current].get('paths', []))
            current = profiles[current].get('extends')
        selected[name] = sorted(paths)
    return selected


def profile_matches(lane, prior, current):
    if lane != 'sparsesolv':
        return True
    # A skipped or compact PR cannot certify the full post-merge lane.
    return (prior.get('qualification_scope') == 'full-lane'
            and prior.get('impact_required') is True
            and prior.get('test_set') == current.get('test_set')
            and bool(current.get('test_set')))


def identity(lane, git):
    result = dict(schema='radia.ci-tree.v1', lane=lane, repository=os.environ['GITHUB_REPOSITORY'],
                commit=git('rev-parse', 'HEAD').strip(), tree=git('rev-parse', 'HEAD^{tree}').strip(),
                run_id=int(os.environ['GITHUB_RUN_ID']))
    if lane == 'sparsesolv':
        result['test_set'] = numerical_test_set(Path(__file__).resolve().parents[1])
    return result


def reuse(lane, event, current, get=api, listing=pages):
    """Only a merged PR, successful original run and API-verified merge tree qualify."""
    sha = current['commit']; repository = current['repository']
    for pr in listing(f'commits/{sha}/pulls'):
        if (not pr.get('merged_at') or pr.get('merge_commit_sha') != sha
                or pr['base']['ref'] != 'main' or pr['base']['repo']['full_name'] != repository
                or pr['head']['repo']['full_name'] != repository):
            continue
        # Branch protection is authoritative. Missing permissions/configuration
        # fails closed; an unavailable API must never create a reuse exemption.
        try:
            required = get('branches/main/protection/required_status_checks')
        except urllib.error.HTTPError as error:
            if error.code != 404:
                raise
            required = {}
        # Rulesets can supplement classic protection; verify their union.
        rules = get('rules/branches/main')
        required = dict(required, checks=required.get('checks', []) + [
            dict(context=c['context'], app_id=c.get('integration_id'))
            for rule in rules if rule['type'] == 'required_status_checks'
            for c in rule['parameters']['required_status_checks']])
        contexts = set(required.get('contexts', [])) | {c['context'] for c in required.get('checks', [])}
        runs = listing(f'actions/workflows/{lane}.yml/runs?event=pull_request&head_sha={pr["head"]["sha"]}', 'workflow_runs')
        for run in runs[:20]:
            if (run.get('status') != 'completed' or run.get('conclusion') != 'success'
                    or run.get('event') != 'pull_request' or run.get('head_sha') != pr['head']['sha']
                    or run.get('path') != f'.github/workflows/{lane}.yml'
                    or run.get('head_repository', {}).get('full_name') != repository
                    or not any(p['number'] == pr['number'] for p in run.get('pull_requests', []))):
                continue
            artifacts = listing(f'actions/runs/{run["id"]}/artifacts', 'artifacts')
            candidates = [a for a in artifacts if a['name'] == f'{lane}-tree-identity' and not a['expired']]
            if len(candidates) != 1:
                continue
            with zipfile.ZipFile(io.BytesIO(get(f'actions/artifacts/{candidates[0]["id"]}/zip', True))) as archive:
                prior = json.loads(archive.read('ci-tree-receipt.json'))
            if any(prior.get(k) != current[k] for k in ('schema', 'lane', 'repository', 'tree')) or prior.get('run_id') != run['id']:
                continue
            if not profile_matches(lane, prior, current):
                continue
            commit = get('git/commits/' + prior['commit'])
            if (commit['tree']['sha'] != current['tree'] or len(commit['parents']) != 2
                    or commit['parents'][1]['sha'] != pr['head']['sha']):
                continue
            # Actions check-runs attach to the PR head, although checkout tests
            # the merge ref. The archived identity above proves that tested tree.
            checks = listing(f'commits/{pr["head"]["sha"]}/check-runs?filter=latest', 'check_runs')
            statuses = listing(f'commits/{pr["head"]["sha"]}/statuses')
            latest = {}
            for status in statuses:
                latest.setdefault(status['context'], status['state'])
            passed = {c['name'] for c in checks if c['status'] == 'completed' and c['conclusion'] == 'success'}
            passed |= {name for name, state in latest.items() if state == 'success'}
            for required_check in required.get('checks', []):
                app_id = required_check.get('app_id')
                if app_id not in (None, -1) and not any(
                    c['name'] == required_check['context'] and c.get('app', {}).get('id') == app_id
                    and c['status'] == 'completed' and c['conclusion'] == 'success' for c in checks
                ):
                    passed.discard(required_check['context'])
            if not checks or not contexts <= passed or any(c['status'] != 'completed' or c['conclusion'] not in ('success', 'skipped', 'neutral') for c in checks) or any(state != 'success' for state in latest.values()):
                continue
            jobs = listing(f'actions/runs/{run["id"]}/jobs?filter=latest', 'jobs')
            expected = 'skipped' if prior.get('impact_required') is False else 'success'
            if not set(LANES[lane]) <= {j['name'] for j in jobs if j['conclusion'] == expected}:
                continue
            return dict(current, mode='tree already validated', pr=pr['number'],
                        validated_commit=prior['commit'], run_ids=sorted({run['id']} | {
                            int(c['details_url'].split('/actions/runs/')[1].split('/')[0])
                            for c in checks if '/actions/runs/' in c.get('details_url', '')}),
                        required_checks=sorted(contexts), artifact_id=candidates[0]['id'],
                        validated_impact_required=prior.get('impact_required', True))
    return None


def route(lane, event_name, event, git, impact_required=True):
    current = identity(lane, git)
    receipt = dict(current, mode='execute', impact_required=impact_required)
    if lane == 'sparsesolv':
        # A successful partially selected PR job cannot certify the full main lane.
        receipt['qualification_scope'] = 'case-plan-not-full-lane'
        receipt['reuse_policy'] = 'disabled: use scoped heavy source receipts instead'
    if lane != 'sparsesolv' and event_name == 'push' and event.get('ref') == 'refs/heads/main':
        try:
            receipt = reuse(lane, event, current) or receipt
        except Exception as error:
            receipt['fallback'] = type(error).__name__ + ': execute full lane'
    Path('ci-tree-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    if os.environ.get('GITHUB_STEP_SUMMARY'):
        with open(os.environ['GITHUB_STEP_SUMMARY'], 'a', encoding='utf-8') as summary:
            summary.write('```json\n' + json.dumps(receipt, indent=2) + '\n```\n')
    return receipt['mode'] != 'execute'
