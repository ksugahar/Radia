"""Checkpoint V-N execution completeness; never certify a release or physics oracle.

Load explicitly with ``pytest -p validation_execution_receipt`` and put tools on
PYTHONPATH. A reviewed plan binds expected pytest nodeids to one clean Git tree.
The scheduler owns resource admission, exact native import checks, timeout and
artifact/private-evidence transport. This plugin does not supply those claims.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time

import pytest


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def git(root, *args):
    return subprocess.check_output(['git', '-c', f'safe.directory={root}', *args],
                                   cwd=root, timeout=30).decode('utf-8').strip()


def unlinked(path):
    """Do not follow output/input junctions or symlinks outside the job boundary."""
    path = Path(path).absolute()
    for part in (path, *path.parents):
        if part.is_symlink() or getattr(part, 'is_junction', lambda: False)():
            raise ValueError(f'linked evidence path: {path.name}')
    return path


def read_plan(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result: raise ValueError(f'duplicate plan key: {key}')
            result[key] = value
        return result
    plan = json.loads(unlinked(path).read_text(encoding='utf-8'), object_pairs_hook=pairs)
    if not isinstance(plan, dict): raise ValueError('execution plan must be an object')
    if set(plan) != {'schema', 'source_sha', 'source_tree', 'run_id', 'attempt',
                     'profile', 'expected_nodeids'} or plan['schema'] != 'radia.validation-execution-plan.v1':
        raise ValueError('invalid execution plan fields/schema')
    for key in ('source_sha', 'source_tree'):
        if not isinstance(plan[key], str) or not re.fullmatch('[0-9a-f]{40}', plan[key]):
            raise ValueError(f'invalid {key}')
    if not isinstance(plan['run_id'], str) or not plan['run_id']:
        raise ValueError('missing run identity')
    if type(plan['attempt']) is not int or plan['attempt'] < 1:
        raise ValueError('invalid attempt')
    profile = plan['profile']
    if (not isinstance(profile, dict) or set(profile) != {'name', 'sha256'}
            or not isinstance(profile['name'], str) or not profile['name']
            or not isinstance(profile['sha256'], str) or not re.fullmatch('[0-9a-f]{64}', profile['sha256'])):
        raise ValueError('invalid profile identity')
    expected = plan['expected_nodeids']
    if not isinstance(expected, list) or not expected or any(not isinstance(x, str) for x in expected):
        raise ValueError('expected cases must be a nonempty explicit list')
    if len(set(expected)) != len(expected): raise ValueError('duplicate expected cases')
    for nodeid in expected:
        file = nodeid.split('::', 1)[0]
        if (not file.startswith('validation_test/') or not file.endswith('.py')
                or '::' not in nodeid or '\\' in file or ':' in file or '..' in Path(file).parts):
            raise ValueError(f'unsafe/non-validation case: {nodeid}')
    return plan


class Receipt:
    def __init__(self, plan_path, source_root, receipt_path, output_root):
        self.root = unlinked(source_root).resolve()
        self.path = unlinked(receipt_path).resolve()
        self.output = unlinked(output_root).resolve()
        if self.path.is_relative_to(self.root) or self.output.is_relative_to(self.root):
            raise ValueError('evidence/output must be outside source tree')
        if self.root.is_relative_to(self.output) or self.path.is_relative_to(self.output):
            raise ValueError('numerical output scope must exclude source and receipt')
        if self.path.exists(): raise ValueError('receipt already exists; use a fresh attempt')
        self.plan = read_plan(plan_path)
        if Path(git(self.root, 'rev-parse', '--show-toplevel')).resolve() != self.root:
            raise ValueError('source root is not a Git checkout root')
        if (git(self.root, 'rev-parse', 'HEAD') != self.plan['source_sha']
                or git(self.root, 'rev-parse', 'HEAD^{tree}') != self.plan['source_tree']):
            raise ValueError('plan/source identity mismatch')
        if git(self.root, 'status', '--porcelain', '--untracked-files=all'):
            raise ValueError('source tree is dirty before validation')
        self.started = time.monotonic()
        self.data = dict(schema='radia.validation-execution.v1', purpose='V-N',
                         qualification_scope='source-nightly-only; never exact release artifact qualification',
                         plan_sha256=sha256(plan_path), identity=self.plan,
                         collected=[], phases={}, phase_durations_s={}, collection_reports=[],
                         unexpected=[], duplicate_collected=[], numerical_outputs=[],
                         workspace_changed=False, state='initializing', status='incomplete')
        self.write()

    def write(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        unlinked(self.path.parent)
        self.data['runtime_s'] = time.monotonic() - self.started
        expected = set(self.plan['expected_nodeids'])
        executed = {key for key, phases in self.data['phases'].items() if 'call' in phases}
        self.data['executed'] = sorted(executed)
        self.data['missing_collection'] = sorted(expected - set(self.data['collected']))
        self.data['missing_execution'] = sorted(expected - executed)
        temporary = self.path.with_name(self.path.name + '.pending')
        unlinked(temporary)
        with temporary.open('x', encoding='utf-8') as stream:
            json.dump(self.data, stream, sort_keys=True, indent=2, allow_nan=False)
            stream.write('\n')
        os.replace(temporary, self.path)

    def collection(self, nodeids):
        selected = [key for key in nodeids if key.split('::', 1)[0].startswith('validation_test/')]
        self.data['collected'] = sorted(set(selected))
        self.data['duplicate_collected'] = sorted({key for key in selected if selected.count(key) > 1})
        self.data['unexpected'] = sorted(set(selected) - set(self.plan['expected_nodeids']))
        self.data['state'] = 'collected'
        self.write()
        return bool(self.data['unexpected'] or self.data['duplicate_collected'])

    def phase(self, nodeid, when, outcome, wasxfail=False, duration_s=None):
        if nodeid not in self.plan['expected_nodeids']: return
        if when not in {'setup','call','teardown'} or outcome not in {'passed','failed','skipped'}:
            raise ValueError('invalid phase/outcome')
        if duration_s is not None and (type(duration_s) not in (int,float) or not math.isfinite(duration_s) or duration_s < 0):
            raise ValueError('invalid duration observation')
        phases = self.data['phases'].setdefault(nodeid, {})
        if when in phases:
            self.data.setdefault('duplicate_phases', []).append([nodeid, when])
        self.data['phase_durations_s'].setdefault(nodeid, {})[when] = duration_s
        phases[when] = ('xpassed' if outcome == 'passed' else 'xfailed') if wasxfail else outcome
        self.data['state'] = 'running'
        self.write()  # A killed child leaves an explicitly incomplete checkpoint.

    def pytest_collectreport(self, report):
        if report.outcome != 'passed':
            self.data['collection_reports'].append(dict(nodeid=report.nodeid, outcome=report.outcome))
            self.write()

    def finish(self, pytest_exit_code):
        self.data['pytest_exit_code'] = int(pytest_exit_code)
        self.data['workspace_changed'] = bool(git(self.root, 'status', '--porcelain', '--untracked-files=all'))
        self.data['final_source_sha'] = git(self.root, 'rev-parse', 'HEAD')
        self.data['final_source_tree'] = git(self.root, 'rev-parse', 'HEAD^{tree}')
        self.data['source_identity_changed'] = (self.data['final_source_sha'] != self.plan['source_sha']
                                                or self.data['final_source_tree'] != self.plan['source_tree'])
        # Numerical files are observations, not independently verified oracles.
        if self.output.exists():
            for path in sorted(self.output.rglob('*')):
                unlinked(path)
                if path.is_file() and path.suffix == '.json':
                    before = path.stat()
                    digest = sha256(path)
                    after = path.stat()
                    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                        raise ValueError('numerical output changed while hashing')
                    self.data['numerical_outputs'].append(dict(path=path.relative_to(self.output).as_posix(),
                                                               sha256=digest, bytes=after.st_size))
        self.data['numerical_output_status'] = 'observed' if self.data['numerical_outputs'] else 'none_observed'
        values = [outcome for phases in self.data['phases'].values() for outcome in phases.values()]
        failed = (int(pytest_exit_code) in (1,3) or 'failed' in values or self.data['workspace_changed'] or self.data['source_identity_changed']
                  or any(x['outcome'] == 'failed' for x in self.data['collection_reports']))
        self.write()
        incomplete = (self.data['missing_collection'] or self.data['missing_execution']
                      or self.data['unexpected'] or self.data['duplicate_collected']
                      or self.data.get('duplicate_phases')
                      or any(set(self.data['phases'].get(nodeid, {})) != {'setup','call','teardown'} for nodeid in self.plan['expected_nodeids'])
                      or any(x in {'skipped', 'xfailed', 'xpassed'} for x in values)
                      or any(x['outcome'] == 'skipped' for x in self.data['collection_reports']))
        self.data['status'] = 'failed' if failed else 'incomplete' if incomplete else 'passed' if pytest_exit_code == 0 else 'incomplete'
        self.data['state'] = 'finished'
        self.data['effective_exit_code'] = int(pytest_exit_code) or (0 if self.data['status'] == 'passed' else 1)
        self.write()
        return self.data['effective_exit_code']


def pytest_addoption(parser):
    for option in ('validation-plan', 'validation-source-root', 'validation-receipt', 'validation-output-root'):
        parser.addoption('--' + option, default=None)


def pytest_configure(config):
    options = [config.getoption(key) for key in ('validation_plan', 'validation_source_root', 'validation_receipt', 'validation_output_root')]
    if not any(options): return
    if not all(options): raise pytest.UsageError('all four validation receipt options are required')
    try:
        config._validation_receipt = Receipt(*options)
        config.pluginmanager.register(config._validation_receipt, 'validation-receipt-collection')
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        raise pytest.UsageError(str(error)) from error


def pytest_collection_finish(session):
    receipt = getattr(session.config, '_validation_receipt', None)
    if receipt and receipt.collection([item.nodeid for item in session.items]):
        raise pytest.UsageError('unapproved/duplicate validation cases; nothing may execute')


@pytest.hookimpl(trylast=True)
def pytest_sessionfinish(session, exitstatus):
    receipt = getattr(session.config, '_validation_receipt', None)
    if receipt: session.exitstatus = receipt.finish(exitstatus)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    result = yield
    report = result.get_result()
    receipt = getattr(item.config, '_validation_receipt', None)
    if receipt:
        receipt.phase(report.nodeid, report.when, report.outcome, hasattr(report, 'wasxfail'), report.duration)
