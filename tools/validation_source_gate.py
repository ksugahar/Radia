"""Execute solver-heavy validation or verify an unchanged-source trusted receipt.

This reuses source validation, NEVER certification of a newly built artifact.
Ordinary tests, native builds and release gates remain separate obligations.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import urllib.request
import urllib.parse
import zipfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone, timedelta

PROFILE = 'solver-heavy'
WORKFLOW = '.github/workflows/sparsesolv.yml'
ARTIFACT = 'solver-heavy-validation-v1'
SCHEMA = 'radia.validated-source.v1'
MAX_AGE = timedelta(days=30)
ROOT_INPUTS = {'CMakeLists.txt','Build.ps1','setup.py','pyproject.toml','MANIFEST.in',
               'tools/run_test_tier.py','tools/native_build_provenance.ps1',
               'tools/validation_source_gate.py','tests/conftest.py','conftest.py'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(data):
    return json.dumps(data, sort_keys=True, separators=(',',':')).encode()


def git(root, *args):
    return subprocess.check_output(['git','-c',f'safe.directory={root}',*args],cwd=root).decode('utf-8').strip()


def environment_contract(root):
    compiler = []
    for file in (root/'build-msvc/CMakeFiles').glob('*/CMakeCXXCompiler.cmake'):
        text=file.read_text(encoding='utf-8')
        compiler.append({key:re.findall(r'set\('+key+r' "([^"\n]+)"\)',text)
                         for key in ('CMAKE_CXX_COMPILER_ID','CMAKE_CXX_COMPILER_VERSION')})
    if len(compiler)!=1 or any(len(v)!=1 for v in compiler[0].values()):
        raise ValueError('compiler identity unknown')
    cache=(root/'build-msvc/CMakeCache.txt').read_text(encoding='utf-8')
    sdk_versions=sorted(set(re.findall(r'10\.0\.\d+\.\d+',cache)))
    if not sdk_versions: raise ValueError('Windows SDK version unknown')
    if os.environ.get('RADIA_TESTS_ALLOW_PARTIAL')!='0': raise ValueError('partial validation context')
    threads={key:os.environ.get(key) for key in ('RADIA_TEST_NGSOLVE_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')}
    if set(threads.values())!={'1'}: raise ValueError('unbounded or changed thread contract')
    return {'python':sys.version,'platform':platform.platform(),'machine':platform.machine(),
            'processor':platform.processor(),'compiler':compiler,'windows_sdk_versions':sdk_versions,'threads':threads,
            'packages':sorted((d.metadata['Name'].lower(),d.version) for d in importlib.metadata.distributions())}


def source_fingerprint(root, environment):
    manifest=json.loads((root/'tests/test_tier_manifest.json').read_text(encoding='utf-8'))
    profiles={}
    name=PROFILE
    while name:
        if name in profiles: raise ValueError('cyclic validation profile')
        profiles[name]=manifest['profiles'][name]
        name=profiles[name].get('extends')
    tracked=git(root,'ls-files','-z').split('\0')
    selected=[]
    def test_input(path):
        return path.startswith('tests/') and not (path.startswith(('tests/policy/','tests/ci_routing/','tests/validation_gate/','tests/test_ci_'))
            or path in {'tests/test_local_push_checks.py','tests/test_pre_push_check.py','tests/test_sparsesolv_ci_contract.py','tests/test_tier_manifest.json'})
    for path in tracked:
        if (path in ROOT_INPUTS or test_input(path) or path.startswith(('src/','validation_test/'))
                or (path.startswith('packages/') and not any(x in path.split('/') for x in ('tests','docs')))):
            selected.append(path)
    # All production, validation inputs/helpers/fixtures and relevant build files.
    # CI workflow routing and unrelated compact-test changes are not numerical inputs.
    if not selected or not any(p.startswith('src/') for p in selected): raise ValueError('source closure unavailable')
    for profile in profiles.values():
        for path in profile['paths']:
            if path not in selected: selected.append(path)
    rows={}
    for path in sorted(selected):
        file=root/path
        if not file.is_file() or file.is_symlink(): raise ValueError('missing/linked validation input: '+path)
        rows[path]=digest(file.read_bytes())
    changes=git(root,'diff','--name-only','HEAD','--').splitlines()
    if set(changes).intersection(rows): raise ValueError('validation inputs dirty')
    payload={'schema':SCHEMA,'profile':PROFILE,'sources':rows,'profiles':profiles,'environment':environment}
    return digest(canonical(payload)),len(rows)


def check_junit(payload):
    cases=list(ET.fromstring(payload).iter('testcase'))
    if not cases or any(c.find(tag) is not None for c in cases for tag in ('failure','error','skipped')):
        raise ValueError('validation must execute every required case successfully')
    return len(cases)


def verified_receipt(root, artifact, run, archive, fingerprint, now):
    """Authoritative GitHub API metadata is separate from artifact contents."""
    if artifact.get('expired') is not False or artifact.get('name')!=ARTIFACT: raise ValueError('expired/wrong artifact')
    if (run.get('path')!=WORKFLOW or run.get('head_branch')!='main' or run.get('status')!='completed'
            or run.get('conclusion')!='success' or run.get('event') not in ('push','workflow_dispatch')):
        raise ValueError('untrusted run provenance')
    if artifact.get('workflow_run',{}).get('id')!=run['id']: raise ValueError('artifact run mismatch')
    if artifact.get('digest')!='sha256:'+digest(archive): raise ValueError('archive digest absent/mismatch')
    expires=datetime.fromisoformat(artifact['expires_at'].replace('Z','+00:00'))
    if expires<=now: raise ValueError('expired artifact')
    updated=datetime.fromisoformat(run['updated_at'].replace('Z','+00:00'))
    if updated>now or now-updated>MAX_AGE: raise ValueError('receipt too old/future')
    # A successful PR/fork artifact can never become a reusable main receipt.
    subprocess.run(['git','-c',f'safe.directory={root}','merge-base','--is-ancestor',run['head_sha'],'origin/main'],
                   cwd=root,check=True,capture_output=True)
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        if set(z.namelist())!={'receipt.json','solver-heavy.xml'}: raise ValueError('unexpected evidence archive')
        if any(i.file_size>16*1024*1024 for i in z.infolist()): raise ValueError('oversized evidence')
        original_receipt=z.read('receipt.json')
        receipt=json.loads(original_receipt)
        xml=z.read('solver-heavy.xml')
    if (receipt.get('schema')!=SCHEMA or receipt.get('profile')!=PROFILE or receipt.get('accepted') is not True
            or receipt.get('reusable') is not True or receipt.get('fingerprint')!=fingerprint
            or receipt.get('source_sha')!=run['head_sha'] or receipt.get('run_id')!=run['id']
            or receipt.get('run_attempt')!=run['run_attempt'] or receipt.get('junit_sha256')!=digest(xml)
            or type(receipt.get('case_count')) is not int or receipt.get('case_count')!=check_junit(xml)):
        raise ValueError('receipt/input/evidence mismatch')
    return receipt,xml,original_receipt


class _SafeRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, message, headers, newurl):
        if urllib.parse.urlsplit(newurl).scheme != 'https':
            raise ValueError('non-HTTPS artifact redirect')
        redirected=super().redirect_request(request,fp,code,message,headers,newurl)
        if urllib.parse.urlsplit(newurl).hostname != 'api.github.com':
            redirected.remove_header('Authorization')
        return redirected


class GitHub:
    def __init__(self, repository, token):
        self.base='https://api.github.com/repos/'+repository
        self.token=token
    def get(self, path, binary=False):
        request=urllib.request.Request(self.base+path,headers={'Authorization':'Bearer '+self.token,
            'Accept':'application/vnd.github+json','X-GitHub-Api-Version':'2022-11-28'})
        with urllib.request.build_opener(_SafeRedirect()).open(request,timeout=30) as response:
            data=response.read(32*1024*1024+1)
        if len(data)>32*1024*1024: raise ValueError('oversized API response')
        return data if binary else json.loads(data)
    def find(self, root, fingerprint):
        # Bounded search; pagination/cache miss safely executes validation.
        artifacts=self.get('/actions/artifacts?name='+ARTIFACT+'&per_page=30')['artifacts']
        for artifact in artifacts:
            try:
                run=self.get('/actions/runs/'+str(artifact['workflow_run']['id']))
                if run.get('head_branch')!='main' or run.get('conclusion')!='success': continue
                archive=self.get('/actions/artifacts/'+str(artifact['id'])+'/zip',binary=True)
                receipt,xml,original_receipt=verified_receipt(root,artifact,run,archive,fingerprint,datetime.now(timezone.utc))
                return artifact,receipt,xml,original_receipt
            except Exception:
                continue  # Untrusted/unknown/expired is a miss, never PASS.
        return None


def run(root, out, force=False, client=None):
    out.mkdir(parents=True,exist_ok=True)
    identity={'source_sha':git(root,'rev-parse','HEAD'),'run_id':int(os.environ.get('GITHUB_RUN_ID','0')),
              'run_attempt':int(os.environ.get('GITHUB_RUN_ATTEMPT','0'))}
    result={'schema':SCHEMA,'profile':PROFILE,'current':identity,'accepted':False,
            'current_artifact_certified':False,'mode':'execute','reuse_reason':None,
            'verified_at':datetime.now(timezone.utc).isoformat()}
    fingerprint=None
    try:
        fingerprint,count=source_fingerprint(root,environment_contract(root))
        result.update(fingerprint=fingerprint,source_count=count)
    except Exception as error:
        result['reuse_reason']='unknown fingerprint: '+type(error).__name__
    found=None
    if not force and fingerprint and client:
        try: found=client.find(root,fingerprint)
        except Exception as error: result['reuse_reason']='receipt lookup unavailable: '+type(error).__name__
    if found:
        artifact,receipt,xml,original_receipt=found
        # Preserve original proof verbatim. Verification has its own current identity.
        (out/'receipt.json').write_bytes(original_receipt)
        (out/'solver-heavy.xml').write_bytes(xml)
        result.update(accepted=True,mode='verified-reuse',original=receipt,artifact_id=artifact['id'])
    else:
        evidence=out/'solver-heavy.xml'
        if evidence.is_symlink(): raise ValueError('linked output evidence')
        evidence.unlink(missing_ok=True)  # Never accept a previous attempt's XML.
        command=[sys.executable,'tools/run_test_tier.py','--profile',PROFILE,'--junitxml',str(out/'solver-heavy.xml')]
        process=subprocess.run(command,cwd=root,check=False)
        if process.returncode==0:
            try:
                xml=(out/'solver-heavy.xml').read_bytes()
                count=check_junit(xml)
                receipt={'schema':SCHEMA,'profile':PROFILE,**identity,'fingerprint':fingerprint,
                         'accepted':True,'reusable':bool(fingerprint),'case_count':count,'junit_sha256':digest(xml)}
                (out/'receipt.json').write_bytes(canonical(receipt))
                result.update(accepted=True,case_count=count)
            except Exception as error:
                result['failure']='invalid validation evidence: '+type(error).__name__
        else: result['exit_code']=process.returncode
    (out/'verification.json').write_bytes(canonical(result))
    return 0 if result['accepted'] else 1


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--force',action='store_true')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    token=os.environ.get('GH_TOKEN')
    repo=os.environ.get('GITHUB_REPOSITORY')
    client=GitHub(repo,token) if token and repo else None
    sys.exit(run(root,args.output,args.force,client))
