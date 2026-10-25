"""Source-only bounded PowerShell builder write audit; never executes a builder.

Guards: ordinary builders cannot acquire a package payload write capability.
Assumes: PowerShell AST, finite commands/helpers and cmd dialect below. Helpers
are exact allow-listed boundaries, not transitively analyzed programs.
"""
import json
import re
import subprocess


def audit(path):
    quoted = str(path).replace("'", "''")
    script = r"""
$errors=$null;$tokens=$null
$a=[System.Management.Automation.Language.Parser]::ParseFile('__PATH__',[ref]$tokens,[ref]$errors)
if($errors.Count){throw $errors}
$rows=@($a.FindAll({param($n)
$n -is [System.Management.Automation.Language.CommandAst] -or
$n -is [System.Management.Automation.Language.AssignmentStatementAst] -or
$n -is [System.Management.Automation.Language.ForEachStatementAst] -or
$n -is [System.Management.Automation.Language.InvokeMemberExpressionAst] -or
$n -is [System.Management.Automation.Language.FileRedirectionAst] -or
($n -is [System.Management.Automation.Language.ExpandableStringExpressionAst] -and $n.Extent.Text.StartsWith('@"'))
},$true)|ForEach-Object {
$r=@{kind=$_.GetType().Name;text=$_.Extent.Text;line=$_.Extent.StartLineNumber}
if($_ -is [System.Management.Automation.Language.CommandAst]){$r.elements=@($_.CommandElements|ForEach-Object{$_.Extent.Text})}
if($_ -is [System.Management.Automation.Language.AssignmentStatementAst]){$r.left=$_.Left.Extent.Text;$r.right=$_.Right.Extent.Text}
if($_ -is [System.Management.Automation.Language.ForEachStatementAst]){$r.left=$_.Variable.Extent.Text;$r.right=$_.Condition.Extent.Text}
if($_ -is [System.Management.Automation.Language.InvokeMemberExpressionAst]){$r.member=$_.Member.Extent.Text;$r.args=@($_.Arguments|ForEach-Object{$_.Extent.Text})}
if($_ -is [System.Management.Automation.Language.FileRedirectionAst]){$r.target=$_.Location.Extent.Text}
$r})
ConvertTo-Json -InputObject $rows -Depth 8 -Compress
""".replace('__PATH__', quoted)
    result = subprocess.run(['pwsh', '-NoProfile', '-NonInteractive', '-Command', script],
                            capture_output=True, encoding='utf-8', errors='replace', timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    rows = json.loads(result.stdout)
    definitions = {}
    for row in rows:
        if row['kind'] in ('AssignmentStatementAst', 'ForEachStatementAst'):
            definitions.setdefault(row['left'].lower(), []).append(row['right'])
    # Roots are capabilities, not arbitrary variable names. Check all their
    # definitions when following dataflow; changing one assignment is visible.
    def expand(expr, seen=()):
        expr = expr.replace('`\n', '').strip()
        if any(bad in expr.lower().replace('\\', '/') for bad in ('/../', 'packages/cubit-mesh-export/src')):
            return '!FORBIDDEN!'
        def sub(match):
            key = match[0].lower()
            if key in ('$project_dir', '$reporoot'): return '<REPO>'
            if key == '$psscriptroot': return '<SCRIPT>'
            if key in seen: return '!CYCLE!'
            values = definitions.get(key)
            if values is None: return key
            return '(' + '|'.join(expand(v, (*seen, key)) for v in values) + ')'
        return re.sub(r'\$(?:env:)?[A-Za-z_]\w*', sub, expr)

    def safe(expr):
        value = expand(expr).lower().replace('\\', '/')
        if '!forbidden!' in value: return False
        # Prefix before any opaque filename interpolation is decisive. No
        # unresolved destination root is accepted. File names in module/install
        # tables cannot contain separators or parent traversal.
        value = value.replace('"', '').replace("'", '')
        roots = ('<repo>/build-msvc', '<repo>/src/radia/', 'src/radia',
                 '<repo> src/cubit_plugin', '<repo>/src/cubit_plugin')
        if any(root in value for root in roots):
            if 'src/cubit_plugin' in value and not any(x in value for x in ('build-pyd', 'build-ccm')):
                # Reset's Path parameter is bound only at checked calls below.
                return False
            return not re.search(r'(?:\.\.|packages/|!cycle!)', value)
        # Reset-BuildDirectory parameter: caller arguments checked as sinks.
        if '$path' in value and 'getfullpath' in value: return True
        # Explicit opt-in Radia install only, no Cubit package output.
        if 'sysconfig.get_paths' in value and 'radia' in value: return True
        return False

    writers = {'copy-item': (1, 2), 'remove-item': (1,), 'new-item': ('-Path',),
               'out-file': ('-FilePath',), 'set-content': ('-LiteralPath',),
               'add-content': ('-LiteralPath',), 'move-item': (1, 2)}
    commands = {'write-host','test-path','join-path','split-path','get-command','get-item',
                'get-content','get-filehash','get-process','stop-process','where.exe','select-object',
                'foreach-object','convertfrom-json','out-null','format-list','set-item',
                'get-nativebuildsourceidentity','clear-nativebuildprovenance','write-nativebuildprovenance',
                'get-nativebuildresult','find-cubitvisualstudio2022','assert-cubitcompilercache',
                'reset-builddirectory','git','cmd','cmd.exe','python','matlab','push-location','pop-location', *writers}
    dynamic = {'$v swherepath'.replace(' ',''), '$vswhere', '$pythonexecutable', '$cmake'}
    helpers = {'tools/native_build_provenance.ps1','tools/find_cubit.ps1',
               'tools/check_native_module_abi.py',
               *('packages/radia-simulink/tools/build_' + x + '.py' for x in
                 ('fmm_cuda','hmatrix_cuda','projection_cuda','ih_cuda','field_cuda','ams_cuda','cuda_mex'))}
    for row in rows:
        kind = row['kind']
        if kind == 'CommandAst':
            parts = row['elements']; name = parts[0].strip('"\'').lower()
            if name == '$ctest': name = '$cmake'
            if name not in commands and name not in dynamic:
                # dot-source exact local helpers, otherwise reject new commands.
                literal = name.replace('\\','/').replace('$project_dir/','').replace('$reporoot/','')
                assert literal in helpers, (path, row, 'unreviewed command/helper')
            if name in writers:
                for selector in writers[name]:
                    if isinstance(selector, str):
                        i = next((i for i,x in enumerate(parts) if x.lower() == selector.lower()), None)
                        assert i is not None, row
                        dest = parts[i+1]
                    else:
                        # Copy/Remove positional forms; parameterized native
                        # builder Remove uses LiteralPath explicitly.
                        args = [x for x in parts[1:] if not x.startswith('-')]
                        assert len(args) >= selector, row
                        dest = args[selector-1]
                    assert safe(dest), (path, row, 'unresolved/non-build write destination', expand(dest))
            if name in ('cmd','cmd.exe'):
                if '$BatchFile' in row['text']:
                    assert safe('$BatchFile') and safe('$BuildLog'), row
                    assert re.fullmatch(r'& cmd.exe /c .*\$BatchFile.*\$BuildLog.*',row['text']), row
                else:
                    assert '$vcvarsBat' in row['text'] and '&& set' in row['text'], row
            if name == 'set-item':
                assert '"env:$($matches[1])"' in parts, row
            if name == 'reset-builddirectory':
                assert safe(parts[1]), (row, expand(parts[1]))
            if name in ('python', '$pythonexecutable'):
                for token in parts:
                    for candidate in re.findall(r'(?:tools|packages)/[\w/.-]+\.py',token.replace('\\','/')):
                        assert candidate in helpers, (row, candidate)
                assert 'packages/cubit-mesh-export/src' not in expand(' '.join(parts)).lower().replace('\\','/'), row
            if name == '$cmake':
                # CMake runs only at explicit build-dir cwd or --test-dir;
                # outputs are not changed by a package destination override.
                assert 'packages/cubit-mesh-export/src' not in row['text'].lower().replace('\\','/'), row
        elif kind == 'AssignmentStatementAst' and '.' in row['left']:
            # Property assignment (mtime etc.) is a write too, independent of
            # which statement block contains it. Scalar variables are dataflow.
            if not row['left'].startswith('$') or ':' in row['left']: continue
            base = row['left'].rsplit('.',1)[0]
            assert safe(base), (path, row, 'unresolved property write')
        elif kind == 'InvokeMemberExpressionAst':
            member = row['member'].lower()
            assert member in {'trim','trimend','tolowerinvariant','startswith','getfullpath',
                              'where','count','round','getenvironmentvariable','join','contains'}, (path, row, 'unreviewed member operation')
        elif kind == 'FileRedirectionAst':
            assert row['target'].lower() == '$null' or safe(row['target']), (path, row)
        elif kind == 'ExpandableStringExpressionAst':
            if row['text'].replace('\r\n','\n').startswith('@\"\nimport pathlib'):
                assert 'write' not in row['text'].lower()
            else:
                _batch(row['text'])
    # Filename tables must remain filenames: no path escape via interpolation.
    for row in rows:
        if row['kind'] == 'AssignmentStatementAst' and row['left'].lower() in ('$modules','$installlist'):
            values = re.findall(r'["\']([^"\']+)["\']',row['right'])
            assert all('/' not in v and '\\' not in v and '..' not in v for v in values), row
    return rows


def _batch(text):
    allowed = {'@echo','echo','rem','set','setlocal','endlocal','call','if','exit',
               'mkdir','cd','for','copy','move','xcopy'}
    native = {'"$cmake_exe"','"$ctest_exe"','cl','link','cl.exe','link.exe','"$ninja_exe"'}
    text = re.sub(r'\^\r?\n\s*', ' ', text)
    for line in text.splitlines()[1:-1]:
        line=line.strip(); lower=line.lower()
        if not line or line in ('(',')',') else (') or lower.startswith('rem '): continue
        if line == '$MklBatchEnvironment': continue  # only reviewed nested environment
        for part in re.split(r'(?<![>^])&(?![&])', lower):
            if part.strip() and part != lower:
                assert part.strip().split()[0] in {'echo','exit'}, ('cmd chained token',line)
        token = lower.split()[0]
        assert token in allowed | native or re.fullmatch(r'"\$[a-z_]+"',token) and token in native, ('cmd token',line)
        # call only the selected vcvars, not arbitrary script invocation.
        if token=='call': assert 'vcvars64.bat' in lower and '$vs_path' in lower, line
        if token in ('copy','move','xcopy'):
            assert re.search(r'"\$PROJECT_DIR\\src\\radia\\[^"/]+"',line), line
        for target in re.findall(r'(?<![\d>])>(?![>&])\s*("[^"]+"|[^\s&]+)',line):
            assert target.lower() in ('nul','$null') or target.startswith('"$BUILD_DIR\\'), ('cmd redirection',line)
