"""Collect private CI observations through gh api GET; never execute workflows."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import urlencode

API_VERSION = '2022-11-28'


def timestamp(value):
    result = datetime.fromisoformat(value.replace('Z', '+00:00'))
    if result.utcoffset() is None or result.utcoffset().total_seconds() != 0:
        raise ValueError('timestamps must be explicit UTC')
    return result


def elapsed(start, end):
    if start is None or end is None:
        return None
    seconds = (timestamp(end) - timestamp(start)).total_seconds()
    if seconds < 0:
        raise ValueError('negative timestamp interval')
    return seconds


def gh_get(endpoint):
    """Credentials remain inside gh. Never read tokens or echo failed payloads."""
    result = subprocess.run(['gh', 'api', '--method', 'GET', '-H',
                             f'X-GitHub-Api-Version: {API_VERSION}', endpoint],
                            capture_output=True, encoding='utf-8', errors='replace', timeout=90)
    if result.returncode:
        raise RuntimeError(f'gh GET failed (exit {result.returncode}); response not logged')
    return json.loads(result.stdout)


class Snapshot:
    def __init__(self, out, get=gh_get):
        self.out, self.get = Path(out), get
        if not self.out.is_absolute():
            raise ValueError('output must be an absolute private path')
        for parent in (self.out, *self.out.parents):
            if parent.is_symlink() or getattr(parent, 'is_junction', lambda: False)():
                raise ValueError('linked output path refused')
            if (parent/'.git').exists() or (parent/'.git').is_file():
                raise ValueError('output inside Git worktree refused')
        self.out.mkdir(parents=True, exist_ok=False)
        (self.out/'raw-api').mkdir()
        self.responses, self.errors = [], []

    def fetch(self, endpoint):
        try:
            data = self.get(endpoint)
            encoded = (json.dumps(data, sort_keys=True, ensure_ascii=True)+'\n').encode()
            path = f'raw-api/{len(self.responses):06d}.json'
            (self.out/path).write_bytes(encoded)
            self.responses.append(dict(endpoint=endpoint, path=path,
                                       sha256=hashlib.sha256(encoded).hexdigest(),
                                       retrieved_at=datetime.now(timezone.utc).isoformat()))
            return data
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError):
            # Do not persist an exception string that might contain credentials.
            self.errors.append(dict(endpoint=endpoint, reason='GET or response decoding failed'))
            return None

    def pages(self, endpoint, key=None, query=None):
        rows, page, seen = [], 1, set()
        while True:
            url = endpoint+'?'+urlencode(dict(query or {}, per_page=100, page=page))
            data = self.fetch(url)
            if data is None:
                return rows
            values = data.get(key) if key and isinstance(data, dict) else data if not key else None
            if not isinstance(values, list) or any(not isinstance(x, dict) or 'id' not in x for x in values):
                self.errors.append(dict(endpoint=url, reason='invalid page shape'))
                return rows
            ids = [x['id'] for x in values]
            if len(ids) != len(set(ids)) or seen.intersection(ids):
                self.errors.append(dict(endpoint=url, reason='duplicate or repeated page'))
                return rows
            seen.update(ids)
            rows.extend(values)
            if len(values) < 100:
                if key and data.get('total_count', len(rows)) != len(rows):
                    self.errors.append(dict(endpoint=url, reason='total_count mismatch or truncated listing'))
                return rows
            page += 1

    def save(self, name, data):
        (self.out/name).write_text(json.dumps(data, indent=2, sort_keys=True)+'\n', encoding='utf-8')


def association(pr, run):
    """Keep PR API links and exact SHA evidence; never join on branch names."""
    if any(x.get('number') == pr['number'] for x in run.get('pull_requests', [])):
        return 'api_pull_request_link'
    if run.get('head_sha') and run['head_sha'] == pr.get('merge_commit_sha'):
        return 'exact_merge_sha'
    if run.get('head_sha') and run['head_sha'] == pr.get('head', {}).get('sha'):
        return 'exact_final_head_sha'
    return None


def observed_elapsed(snapshot, start, end):
    try:
        return elapsed(start, end)
    except (ValueError, TypeError, AttributeError):
        snapshot.errors.append(dict(reason='invalid timestamp interval'))
        return None


def job_record(job, run_id, attempt, snapshot):
    return dict(run_id=run_id, attempt=attempt, job_id=job['id'], name=job.get('name'),
                status=job.get('status'), conclusion=job.get('conclusion'),
                started_at=job.get('started_at'), completed_at=job.get('completed_at'),
                duration_s=observed_elapsed(snapshot, job.get('started_at'), job.get('completed_at')))


def collect(snapshot, repo, start, until):
    if not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+', repo):
        raise ValueError('expected explicit OWNER/REPO')
    lo, hi = timestamp(start), timestamp(until)
    if lo >= hi:
        raise ValueError('empty or reversed merge window')
    prefix = f'repos/{repo}'
    # No search API or filtered runs query: avoids the search 1000-result cap.
    pulls = snapshot.pages(prefix+'/pulls', query={'state':'closed'})
    prs = []
    for pr in pulls:
        if not pr.get('merged_at'): continue
        try:
            if lo <= timestamp(pr['merged_at']) < hi:
                if type(pr.get('number')) is not int: raise ValueError('missing PR number')
                prs.append(pr)
        except (ValueError, TypeError, AttributeError):
            snapshot.errors.append(dict(reason='invalid PR identity or merge timestamp'))
    runs = snapshot.pages(prefix+'/actions/runs', 'workflow_runs') if prs else []
    selected, jobs, links = [], [], []
    for run in runs:
        matches = [(p['number'], association(p, run)) for p in prs]
        matches = [(number, method) for number, method in matches if method]
        if not matches:
            continue
        links.extend(dict(pr=number, run_id=run['id'], method=method,
                          event=run.get('event')) for number, method in matches)
        count = run.get('run_attempt')
        if type(count) is not int or count < 1:
            snapshot.errors.append(dict(run_id=run['id'], reason='missing attempt count'))
            continue
        for attempt in range(1, count+1):
            detail = snapshot.fetch(prefix+f'/actions/runs/{run["id"]}/attempts/{attempt}')
            if not isinstance(detail, dict) or detail.get('id') != run['id'] or detail.get('run_attempt') != attempt:
                snapshot.errors.append(dict(run_id=run['id'], attempt=attempt, reason='missing/mismatched attempt'))
                continue
            selected.append(dict(run_id=run['id'], attempt=attempt, event=detail.get('event'),
                                 head_sha=detail.get('head_sha'), workflow_path=detail.get('path'),
                                 status=detail.get('status'), conclusion=detail.get('conclusion'),
                                 created_at=detail.get('created_at'), started_at=detail.get('run_started_at'),
                                 updated_at=detail.get('updated_at'),
                                 queue_s=observed_elapsed(snapshot, detail.get('created_at'), detail.get('run_started_at')) if attempt == 1 else None,
                                 wall_s=None, wall_reason='API updated_at is not an exact completion timestamp'))
            rows = snapshot.pages(prefix+f'/actions/runs/{run["id"]}/attempts/{attempt}/jobs', 'jobs')
            for job in rows:
                if job.get('run_id') != run['id'] or job.get('run_attempt') != attempt:
                    snapshot.errors.append(dict(job_id=job['id'], reason='job attempt identity mismatch'))
                    continue
                jobs.append(job_record(job, run['id'], attempt, snapshot))
    summary_prs = []
    for pr in prs:
        ids = {x['run_id'] for x in links if x['pr'] == pr['number']}
        related = [j for j in jobs if j['run_id'] in ids]
        summary_prs.append(dict(number=pr['number'], created_at=pr.get('created_at'),
                               merged_at=pr['merged_at'], final_head_sha=pr.get('head',{}).get('sha'),
                               merge_commit_sha=pr.get('merge_commit_sha'),
                               open_to_merge_s=observed_elapsed(snapshot, pr.get('created_at'), pr['merged_at']),
                               associated_run_count=len(ids),
                               executed_job_count=sum(bool(j['started_at']) and j['conclusion'] != 'skipped' for j in related),
                               skipped_job_count=sum(j['conclusion'] == 'skipped' for j in related),
                               cancelled_job_count=sum(j['conclusion'] == 'cancelled' for j in related),
                               unknown_duration_job_count=sum(j['duration_s'] is None for j in related),
                               unstarted_job_count=sum(not j['started_at'] for j in related),
                               observed_sum_job_execution_s=sum(j['duration_s'] for j in related if j['duration_s'] is not None),
                               linkage_status='observed_links_not_exhaustive_history' if ids else 'unresolved'))
    result = dict(schema='radia.tests-design-metrics.v1', repository=repo,
                  merged_from=start, merged_until_exclusive=until,
                  collected_at=datetime.now(timezone.utc).isoformat(), api_version=API_VERSION,
                  collector_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  collection_complete=not snapshot.errors, errors=snapshot.errors,
                  prs=summary_prs, runs=selected, jobs=jobs, associations=links,
                  limitations=['PR API links and exact final/merge SHA only; earlier heads with no API link unresolved',
                               'No causal attribution, critical path, adoption or release qualification',
                               'Job sums are observed resource time, not user wait; missing durations not zero',
                               'Independent GET pages are not a transactional historical snapshot'],
                  responses=snapshot.responses)
    snapshot.save('summary.json', result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('repo','merged-from','merged-until','out'):
        parser.add_argument('--'+name, required=True)
    args = parser.parse_args(argv)
    # Validate input before creating snapshot or invoking gh.
    if not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+', args.repo):
        parser.error('expected OWNER/REPO')
    if timestamp(args.merged_from) >= timestamp(args.merged_until):
        parser.error('invalid UTC window')
    result = collect(Snapshot(args.out), args.repo, args.merged_from, args.merged_until)
    return 0 if result['collection_complete'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
