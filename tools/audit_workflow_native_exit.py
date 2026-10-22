"""Find discarded native exit codes in PowerShell workflow run blocks.

Conservative line-level audit, not a PowerShell interpreter. Recognizes common
native executables, literal .exe paths and call-operator invocations; here-string
bodies and backtick continuations are not treated as separate invocations.
Explicit non-PowerShell shells are outside this audit. Unknown runner defaults
are scanned as PowerShell so a matrix/custom runner cannot hide the gate.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import yaml

NATIVE = re.compile(r'^\s*(?:[\'\"]@\s*\|\s*)?(?:\$[\w:]+\s*=\s*)?(?:&\s+|(?:python(?:3(?:\.\d+)?)?|git|cmake|ctest|ninja|pip|dotnet|npm|node|gh|pwsh|powershell|signtool|robocopy|curl|tar|matlab|octave|cmd|ssh|scp|wsl|bash|7z|nuget|msbuild|cl|link|vswhere)(?:\.exe)?\s|[\w./\\:-]+\.exe(?:\s|$))', re.I)
GUARD = re.compile(r'^\s*if\s*\(\s*\$LASTEXITCODE\s*-ne\s*0\s*\)\s*\{\s*(?:exit\s+\$LASTEXITCODE|throw\b[^}]*)\s*\}\s*$', re.I)

def native_lines(script):
    lines = script.splitlines()
    result = []
    here = None
    continuation = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if here:
            if stripped.startswith(here):
                here = None
                if NATIVE.match(line): result.append((i, i))
            continue
        if stripped.startswith('#'): continue
        if continuation:
            continuation = stripped.endswith('`')
            continue
        match = NATIVE.match(line)
        if match:
            end = i
            while end + 1 < len(lines) and lines[end].rstrip().endswith('`'):
                end += 1
            result.append((i, end))
        if re.search(r'@[\'\"]\s*$', line):
            here = "'@" if line.rstrip().endswith("@'") else '"@'
        continuation = stripped.endswith('`')
    return result

def propagated(text):
    if GUARD.fullmatch(text.splitlines()[0]): return True
    # Multiline form may emit a diagnostic before the unconditional exit.
    return bool(re.match(r'^if\s*\(\s*\$LASTEXITCODE\s*-ne\s*0\s*\)\s*\{\s*(?:Write-Host[^\n]*\n\s*)*(?:exit\s+\$LASTEXITCODE|throw[^\n]*)\s*\}', text, re.I))

def segments(line):
    """Split only unquoted statement separators; quoted Python code stays whole."""
    result=[]; start=0; quote=None; escaped=False
    for i,char in enumerate(line):
        if escaped: escaped=False; continue
        if char == '`': escaped=True; continue
        if quote:
            if char == quote: quote=None
        elif char in {"'", '"'}: quote=char
        elif char == ';': result.append(line[start:i]); start=i+1
    result.append(line[start:]); return result

def unchecked(script):
    lines = script.splitlines()
    calls = native_lines(script)
    result = []
    for i,line in enumerate(lines):
        if i in {start for start, _ in calls} and sum(bool(NATIVE.match(part)) for part in segments(line)) > 1:
            result.append({'line':i+1,'end_line':i+1,'command':line.strip(),'reason':'split inline native statements into individually guarded lines'})
    for (start, end), (following, _) in zip(calls, calls[1:]):
        between = [s.strip() for s in lines[end + 1:following] if s.strip() and not s.lstrip().startswith('#')]
        # Do not count merely printing/saving LASTEXITCODE as propagation.
        if not between or not propagated('\n'.join(between)):
            result.append({'line':start + 1, 'end_line':end + 1, 'command':lines[start].strip()})
    return result

def powershell(doc, job, step):
    shell = step.get('shell') or job.get('defaults', {}).get('run', {}).get('shell') or doc.get('defaults', {}).get('run', {}).get('shell')
    if shell:
        return str(shell).split()[0].lower() in {'pwsh','powershell'}
    runner = str(job.get('runs-on','')).lower()
    return not any(x in runner for x in ('ubuntu-','macos-'))

def audit(root):
    findings = []
    for path in sorted((Path(root) / '.github/workflows').glob('*.y*ml')):
        doc = yaml.safe_load(path.read_text(encoding='utf-8'))
        for jid, job in doc.get('jobs', {}).items():
            for index, step in enumerate(job.get('steps', [])):
                if not isinstance(step.get('run'), str) or not powershell(doc, job, step): continue
                for finding in unchecked(step['run']):
                    findings.append(dict(workflow=path.relative_to(root).as_posix(),job=jid,step=index,name=step.get('name',''),**finding))
    return findings

def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--report',type=Path)
    args=parser.parse_args(); findings=audit(args.root)
    report={'schema':'radia.workflow-native-exit-audit.v1','findings':findings,'scope':'recognized line-level PowerShell native calls; not a proof for arbitrary dynamic shell programs'}
    if args.report: args.report.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
    return int(bool(findings))
if __name__=='__main__':raise SystemExit(main())
