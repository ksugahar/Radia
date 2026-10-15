"""Manifest-driven distribution CI impact, including conservative root content diffs."""

from fnmatch import fnmatchcase
import json
import os
from pathlib import Path
import subprocess
import tomllib
import argparse
from collections import defaultdict
import importlib.util
import math
import re
import difflib

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("ci_contract", ROOT / "tools/ci_distribution_contract.py")
try:
    CONTRACT = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(CONTRACT)
except (OSError, SyntaxError, ValueError):
    CONTRACT = None  # main emits required=true even when its configuration cannot load.


def manifest():
    return CONTRACT.load(ROOT)


# Checked against both workflow trigger lists by the lightweight contract suite.
try:
    NATIVE_PATHS = tuple(p for p in manifest()['lanes']['sparsesolv']['events']['push']['paths']
                         if not p.startswith('!') and p not in ('pyproject.toml','tests/test_tier_manifest.json'))
except (OSError, ValueError, KeyError, TypeError, AttributeError):
    NATIVE_PATHS = ()  # compatibility export; never used to exempt runtime work


def test_dependencies_only(before: str, after: str) -> bool:
    """All configuration outside the two development extras must be identical."""
    old, new = tomllib.loads(before), tomllib.loads(after)
    if _canonical(old) == _canonical(new):
        return False
    for document in (old, new):
        extras = document.get("project", {}).get("optional-dependencies", {})
        for name in ("test", "dev"):
            values = extras.pop(name, [])
            if not isinstance(values, list) or not all(isinstance(x, str) for x in values):
                return False
    return _canonical(old) == _canonical(new)


def _valid_budget(value) -> bool:
    """A time budget is a finite positive number; JSON true/false are not numbers."""
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    try:
        return math.isfinite(value) and value > 0
    except (OverflowError, ValueError):
        return False


def _typed(value):
    """Tag every leaf with its type so no string can impersonate another type."""
    if isinstance(value, dict):
        return ["dict", sorted([str(k), _typed(v)] for k, v in value.items())]
    if isinstance(value, (list, tuple)):
        return ["list", [_typed(v) for v in value]]
    return [type(value).__name__, repr(value)]


def _canonical(value) -> str:
    """Type-strict comparison form: 1, true, a date and a look-alike string all differ."""
    return json.dumps(_typed(value), ensure_ascii=False)


def ams_manifest_unchanged(before: str, after: str) -> bool:
    def relevant(text):
        data = json.loads(text)
        if data["schema"] != "radia.test-tier-manifest.v1":
            raise ValueError("unknown manifest schema")
        profiles = data.pop("profiles")
        data.pop("impact_rules", None)  # AMS invokes the profile without --since.
        selected = {}
        name = "sparsesolv"
        while name:
            if name in selected:
                raise ValueError("cyclic AMS profile inheritance")
            profile = profiles[name]
            paths = profile.get("paths", [])
            if not isinstance(paths, list) or not all(isinstance(p, str) for p in paths):
                raise ValueError("invalid AMS paths")
            budget = profile.get("max_elapsed_seconds")
            if budget is not None and not _valid_budget(budget):
                raise ValueError("invalid AMS time budget")
            selected[name] = profile
            name = profile.get("extends")
            if name is not None and not isinstance(name, str):
                raise ValueError("invalid AMS parent")
        for profile_name in ("solver-numerics", "solver-heavy"):
            if profile_name in profiles:
                selected[profile_name] = profiles[profile_name]
        data["profiles"] = selected
        paths = [path for profile in selected.values() for path in profile.get("paths", [])]
        if not paths or len(paths) != len(set(paths)):
            raise ValueError("empty or duplicate AMS test paths")
        return data
    return _canonical(relevant(before)) == _canonical(relevant(after))


NEWLINE = chr(10)


def code_tokens(text, language):
    """Conservative lexer: preserve strings, discard comments/spacing, fail closed."""
    tokens = []; i = 0
    while i < len(text):
        if text[i].isspace():
            # PowerShell ends statements at newlines; keep one boundary token so joining two
            # statements is a change, while indentation and blank/comment-only lines are not.
            if language == 'powershell' and text[i] == NEWLINE and tokens and tokens[-1] != NEWLINE:
                tokens.append(NEWLINE)
            i += 1; continue
        if language == 'powershell' and text.startswith('<#', i):
            end = text.find('#>', i+2)
            if end < 0: raise ValueError('unterminated block comment')
            i = end+2; continue
        if text[i] == '#':
            bracket = re.match(r'#\[(=*)\[', text[i:])
            if language == 'cmake' and bracket:
                end = text.find(']'+bracket[1]+']', i+len(bracket[0]))
                if end < 0: raise ValueError('unterminated CMake comment')
                i = end+len(bracket[1])+2
            else:
                end = text.find('\n', i); i = len(text) if end < 0 else end  # keep the newline as a boundary
            continue
        start = i
        here = re.match(r'@([\'"])\r?\n', text[i:]) if language == 'powershell' else None
        bracket = re.match(r'\[(=*)\[', text[i:]) if language == 'cmake' else None
        if here or bracket:
            closing = here[1]+'@' if here else ']'+bracket[1]+']'
            end = text.find(('\n' if here else '')+closing, i+len((here or bracket)[0]))
            if end < 0: raise ValueError('unterminated multiline string')
            i = end+len(closing)+(1 if here else 0)
        elif text[i] in "\"'":
            quote = text[i]; i += 1
            while i < len(text):
                if text[i] == ('`' if language == 'powershell' else '\\'):
                    i += 2; continue
                if text[i] == quote:
                    if language == 'powershell' and i+1 < len(text) and text[i+1] == quote:
                        i += 2; continue
                    i += 1; break
                i += 1
            else: raise ValueError('unterminated string')
        else:
            word = re.match(r'[\w$./:-]+', text[i:])
            i += len(word[0]) if word else 1
        tokens.append(text[start:i])
    return tokens


def _leaves(value, prefix=()):
    if isinstance(value, dict):
        if not value: yield prefix, {}
        for key, child in value.items(): yield from _leaves(child, prefix+(key,))
    else: yield prefix, value


def _toml_owners(old, new, config):
    before, after = tomllib.loads(old), tomllib.loads(new)
    a, b = dict(_leaves(before)), dict(_leaves(after)); result = set()
    for key in a.keys() | b.keys():
        if key in a and key in b and _canonical(a[key]) == _canonical(b[key]): continue
        if key in (('project','optional-dependencies','test'), ('project','optional-dependencies','dev')):
            for value in (a.get(key, []), b.get(key, [])):
                if not isinstance(value, list) or not all(isinstance(v,str) for v in value):
                    raise ValueError('invalid dependency extra')
        if key[:2] == ('project','optional-dependencies') and len(key) == 3:
            result.update(config['extras'][key[2]])
        elif key[:2] in (('tool','setuptools'),) and len(key) >= 3:
            if key[2] in ('package-data','exclude-package-data') and len(key) == 4:
                result.update(config['packages'][key[3]])
            else: raise ValueError('unknown setuptools setting')
        elif len(key) >= 2 and key[0] == 'project' and key[1] in config['project_keys']:
            result.add('radia')
        elif key == ('build-system','requires') or key == ('build-system','build-backend'):
            result.update(config['build_system'])
        elif key[:2] in tuple(tuple(p.split('.')) for p in config['tooling_sections']):
            result.add('repo')
        else: raise ValueError(f'unknown TOML key: {key}')
    return result


def _ps_parts(text, switches):
    tokens = code_tokens(text, 'powershell'); stack = []; pairs = {}
    for i,t in enumerate(tokens):
        if t == '{': stack.append(i)
        elif t == '}':
            if not stack: raise ValueError('unbalanced PowerShell block')
            pairs[stack.pop()] = i
    if stack: raise ValueError('unbalanced PowerShell block')
    delimiters=[]
    for token in tokens:
        if token in ('(','['): delimiters.append(token)
        elif token in (')',']'):
            if not delimiters or delimiters.pop()!= {')':'(',']':'['}[token]:
                raise ValueError('unbalanced PowerShell expression')
    if delimiters: raise ValueError('unbalanced PowerShell expression')
    parts = defaultdict(list); i = 0
    while i < len(tokens):
        if (i+4 < len(tokens) and tokens[i].lower() == 'if' and tokens[i+1] == '('
                and tokens[i+3:i+5] == [')','{'] and tokens[i+2].lstrip('$') in switches):
            switch = tokens[i+2].lstrip('$'); end = pairs[i+4]
            parts[switch].append(tuple(tokens[i+5:end]))
            parts['unknown'].append(tuple(tokens[i:i+5])+('SCOPED_BODY',))
            i = end+1
        else:
            parts['unknown'].append(tokens[i]); i += 1
    return dict(parts)


def _ps_owners(old,new,config):
    a,b = _ps_parts(old,config),_ps_parts(new,config)
    if a.get('unknown') != b.get('unknown'): raise ValueError('unscoped PowerShell edit')
    result=set()
    for switch in (a.keys() | b.keys())-{'unknown'}:
        if a.get(switch)!=b.get(switch): result.update(config[switch])
    # PowerShell whitespace is not provably harmless (backtick continuation, statement
    # boundaries, here-strings): a textual change with identical tokens selects everyone.
    if not result and old != new:
        raise ValueError('unproven PowerShell formatting change')
    return result


def _cmake_parts(text, config):
    tokens=code_tokens(text,'cmake'); i=0; scopes=[]; parts=defaultdict(list)
    while i<len(tokens):
        if i+1>=len(tokens) or tokens[i+1]!='(' or not re.fullmatch(r'\w+',tokens[i]):
            raise ValueError('invalid CMake command')
        name=tokens[i].lower(); start=i; i+=2; depth=1; args=[]
        while i<len(tokens) and depth:
            token=tokens[i]; depth += (token=='(')-(token==')')
            if depth: args.append(token)
            i+=1
        if depth: raise ValueError('unbalanced CMake command')
        owner=None
        if name=='if':
            owner=config['options'].get(args[0]) if len(args)==1 else None
            scopes.append(owner)
        elif name=='endif':
            if not scopes: raise ValueError('unbalanced CMake if')
            owner=scopes.pop()
        elif name in ('else','elseif'):
            if not scopes: raise ValueError('unbalanced CMake branch')
            # Conditional alternate branches must be reviewed globally.
            if scopes[-1]: raise ValueError('scoped CMake alternate branch')
        elif name=='option' and args: owner=config['options'].get(args[0])
        elif (name.startswith('target_') or name in ('add_library','add_executable','pybind11_add_module','set_target_properties')) and args:
            owner=config['targets'].get(args[0].strip('"'))
        if owner is None:
            owner=next((o for o in reversed(scopes) if o is not None),None)
        parts[tuple(owner) if owner else ('unknown',)].append(tuple(tokens[start:i]))
    if scopes: raise ValueError('unclosed CMake if')
    return dict(parts)


def _cmake_owners(old,new,config):
    a,b=_cmake_parts(old,config),_cmake_parts(new,config)
    if a.get(('unknown',))!=b.get(('unknown',)): raise ValueError('unscoped CMake edit')
    result=set()
    for owner in (a.keys()|b.keys())-{('unknown',)}:
        if a.get(owner)!=b.get(owner): result.update(owner)
    if not result and code_tokens(old,'cmake') != code_tokens(new,'cmake'):
        raise ValueError('CMake command ordering changed')
    return result


def _manifest_parts(text):
    result=[]
    for raw in text.splitlines():
        if '\\#' in raw: raise ValueError('escaped sdist comment requires review')
        line=raw.split('#',1)[0].strip()
        if not line: continue
        fields=line.split()
        arity={'include':2,'exclude':2,'graft':2,'prune':2,'recursive-include':3,'recursive-exclude':3,
               'global-include':2,'global-exclude':2}
        if fields[0] not in arity or len(fields)<arity[fields[0]]:
            raise ValueError('unknown sdist directive')
        if fields[0] in ('graft','prune') and len(fields)!=2:
            raise ValueError('invalid sdist directory directive')
        result.append(tuple(fields))
    return result


def _sdist_owners(old,new,data):
    result=set(); a,b=_manifest_parts(old),_manifest_parts(new); changed=[]
    for op,i,j,k,l in difflib.SequenceMatcher(a=a,b=b,autojunk=False).get_opcodes():
        if op!='equal': changed.extend(a[i:j]+b[k:l])
    for fields in changed:
        if fields[0].startswith('global-'): raise ValueError('global sdist directive changed')
        paths=fields[1:2] if fields[0].startswith('recursive-') else fields[1:]
        for path in paths:
            if path in ('LICENSE','COPYRIGHT.txt','README.md','BUILD.md','setup.py'):
                result.add('radia'); continue
            if path in ('docs','tests','src','packages') or path.startswith(('..','/')):
                raise ValueError('cross-owner sdist glob')
            probe=path.split('*',1)[0].rstrip('/')
            if not probe: raise ValueError('global sdist glob')
            owners=CONTRACT.owner_rule(data,probe+'/__input__' if '.' not in probe.rsplit('/',1)[-1] else probe)['owners']
            if owners==['repo']: raise ValueError('unattributed sdist input')
            result.update(owners)
            # A directory/glob can contain narrower shared ownership rules.
            for rule in data['rules']:
                if rule['pattern'].startswith(probe+'/') or ('*' in path and rule['pattern'].startswith(probe)):
                    result.update(rule['owners'])
    return result


def _tier_owners(old,new,data,config):
    def unique(pairs):
        d={}
        for key,value in pairs:
            if key in d: raise ValueError('duplicate JSON key')
            d[key]=value
        return d
    a,b=json.loads(old,object_pairs_hook=unique),json.loads(new,object_pairs_hook=unique); result=set()
    for d in (a,b):
        if d.get('schema')!='radia.test-tier-manifest.v1' or not set(d)<= {'schema','profiles','impact_rules'}:
            raise ValueError('unknown tier schema')
        profiles=d['profiles']
        for name,profile in profiles.items():
            if not set(profile)<={'description','paths','extends','max_elapsed_seconds'}: raise ValueError('unknown profile key')
            if not isinstance(profile.get('paths',[]),list) or not all(isinstance(p,str) for p in profile.get('paths',[])): raise ValueError('invalid profile paths')
            if 'max_elapsed_seconds' in profile and not _valid_budget(profile['max_elapsed_seconds']): raise ValueError('invalid budget')
            seen=set(); current=name
            while current is not None:
                if current in seen: raise ValueError('tier inheritance cycle')
                seen.add(current); current=profiles[current].get('extends')
    for name in a['profiles'].keys()|b['profiles'].keys():
        present=(name in a['profiles'], name in b['profiles'])
        pa,pb=a['profiles'].get(name,{}).copy(),b['profiles'].get(name,{}).copy()
        pa.pop('description',None);pb.pop('description',None)
        if all(present) and _canonical(pa)==_canonical(pb): continue
        if name in config: result.update(config[name])
        else:
            paths=pa.get('paths',[])+pb.get('paths',[])
            if not paths or pa.get('extends') or pb.get('extends'): raise ValueError('unknown profile inheritance')
            for path in paths: result.update(CONTRACT.owner_rule(data,path)['owners'])
    ia,ib=a.get('impact_rules',{}),b.get('impact_rules',{})
    for rules in (ia,ib):
        if not isinstance(rules,dict): raise ValueError('invalid impact_rules')
        for key,value in rules.items():
            if not isinstance(key,str) or not isinstance(value,list) or not all(isinstance(v,str) for v in value):
                raise ValueError('invalid impact rule')
    for path in ia.keys()|ib.keys():
        if path in ia and path in ib and _canonical(ia[path])==_canonical(ib[path]): continue
        for p in [path]+ia.get(path,[])+ib.get(path,[]): result.update(CONTRACT.owner_rule(data,p)['owners'])
    return result


def root_owners(path,old,new,data):
    """Return owners and an auditable reason; every uncertainty selects all owners."""
    all_owners=set(CONTRACT.owner_rule(data,path)['owners'])
    try:
        if old is None or new is None: raise ValueError('missing/deleted/renamed blob')
        rules=data['root_content_diff']['attribution']
        if path=='pyproject.toml': owners=_toml_owners(old,new,rules[path])
        elif path=='Build.ps1': owners=_ps_owners(old,new,rules[path])
        elif path=='CMakeLists.txt': owners=_cmake_owners(old,new,rules[path])
        elif path=='MANIFEST.in': owners=_sdist_owners(old,new,data)
        elif path=='tests/test_tier_manifest.json': owners=_tier_owners(old,new,data,rules[path])
        else: raise ValueError('unknown root config')
        if not owners <= all_owners|{'repo'}: raise ValueError('unknown attributed owner')
        return owners,'verified content diff' if owners else 'comment/whitespace/prose-only change'
    except (ValueError,KeyError,TypeError,AttributeError,IndexError) as exc:
        return all_owners,f'{type(exc).__name__}: {exc}; all file owners (fail closed)'


def lane_required(lane,changed,blobs,data):
    if not changed: return True,'empty/unknown change set; fail closed'
    selected_owners=set(); reasons=[]; definition=data['lanes'][lane]
    for path in sorted(changed):
        if path in CONTRACT.ROOT_CONFIGS:
            try: old,new=blobs(path)
            except (OSError,ValueError,KeyError,subprocess.SubprocessError): old=new=None
            owners,reason=root_owners(path,old,new,data)
            # Preserve the existing AMS-only profile/dev-extra exemptions.
            if lane=='sparsesolv' and old is not None and new is not None and 'fail closed' not in reason:
                try:
                    if path=='pyproject.toml' and test_dependencies_only(old,new): owners=set()
                    if path=='tests/test_tier_manifest.json' and ams_manifest_unchanged(old,new): owners=set()
                except (ValueError,KeyError,TypeError,AttributeError): pass
            selected_owners.update(owners); reasons.append(path+': '+reason)
        else:
            try: CONTRACT.owner_rule(data,path)
            except ValueError: return True,'unowned input; fail closed'
            events=[p for e,p in definition['events'].items() if e in ('push','pull_request')]
            if any(CONTRACT.selected(path,p) for p in events):
                if lane=='sparsesolv' and path.endswith('.md'): continue
                return True,'existing non-root lane input: '+path
    required=bool(selected_owners.intersection(definition['owners']))
    return required,'; '.join(reasons) or 'no lane inputs'


def changed_paths(raw):
    """Parse --name-status -z, preserving both sides of renames/copies/deletions."""
    fields=raw.rstrip('\0').split('\0'); result=set(); i=0
    while i<len(fields):
        status=fields[i];i+=1
        if not re.fullmatch(r'[ACDMRTUXB][0-9]*',status): raise ValueError('invalid diff status')
        count=2 if status[0] in 'RC' else 1
        if i+count>len(fields): raise ValueError('truncated diff status')
        result.update(fields[i:i+count]);i+=count
    return result


def required_for_event(lane,event_name,event,git,data):
    try:
        return _required_for_event(lane,event_name,event,git,data)
    except Exception as exc:  # malformed event payloads must also select the lane
        return True,f'impact unknown ({type(exc).__name__}); fail closed'


def _required_for_event(lane,event_name,event,git,data):
    if event_name not in ('push','pull_request') or event.get('ref','').startswith('refs/tags/'):
        return True,'release/manual/unknown event: preserve full lane'
    base=event.get('before') if event_name=='push' else event.get('pull_request',{}).get('base',{}).get('sha')
    if not base or set(base)=={'0'}: return True,'missing comparison base; fail closed'
    try:
        if event_name=='pull_request': base=git('merge-base',base,'HEAD').strip()
        if not base: raise ValueError('missing merge base')
        changed=changed_paths(git('diff','--name-status','-z','--find-renames',base,'HEAD','--'))
        return lane_required(lane,changed,lambda p:(git('show',f'{base}:{p}'),git('show',f'HEAD:{p}')),data)
    except Exception as exc:  # any unanticipated failure must select the lane
        return True,f'impact unknown ({type(exc).__name__}); fail closed'


def native_required(event_name, event, git):
    """Compatibility API for callers supplying a path-only diff; one decision engine."""
    if event_name not in ('push','pull_request'): return True,'manual or unknown event'
    base='HEAD^1' if event_name=='pull_request' else event.get('before')
    if not base or set(base)=={'0'}: return True,'no comparison base'
    changed=set(filter(None,git('diff','--name-only','-z',base,'HEAD','--').split('\0')))
    if changed and all(p.endswith('.md') for p in changed): return False,'documentation-only change'
    return lane_required('sparsesolv',changed,lambda p:(git('show',f'{base}:{p}'),git('show',f'HEAD:{p}')),manifest())


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lane",default="sparsesolv")
    args=parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]

    def git(*args):
        return subprocess.check_output(
            ["git", "-c", f"safe.directory={root.as_posix()}", "-C", str(root), *args],
            encoding="utf-8", stderr=subprocess.PIPE,
        )

    try:
        event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text(encoding="utf-8"))
        data=manifest()
        CONTRACT.validate_schema(data)
        if args.lane not in data["lanes"]: raise ValueError("unknown lane")
        required, reason = required_for_event(args.lane, os.environ["GITHUB_EVENT_NAME"], event, git, data)
    except Exception as exc:  # any unanticipated failure must select the lane
        required, reason = True, f"impact unknown ({type(exc).__name__}); fail closed"
    safe_reason = reason.replace('\r',' ').replace('\n',' ')
    print(f"{args.lane} required={required}: {safe_reason}")
    with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as output:
        output.write(f"required={str(required).lower()}\n")
        output.write("reason=" + safe_reason + "\n")


if __name__ == "__main__":
    main()
