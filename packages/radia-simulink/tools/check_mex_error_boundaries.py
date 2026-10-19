"""Conservative syntactic fence for explicit MATLAB error exits.

This is deliberately not a general C++ type checker. Every explicit mexErrMsg*
call must be in a narrowly recognized POD-only gateway; unknown syntax fails.
It does not certify MATLAB internals or transitive errors raised by other APIs.
"""
from pathlib import Path
import re

TOKENS = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'', re.S)
FUNCTION = re.compile(r'\b(?:void|bool)\s+(\w+)\s*\([^;{}]*\)\s*(?:noexcept)?\s*\{')
ERROR = re.compile(r'\bmexErrMsg\w*\s*\(')
GATEWAY = re.compile(
    r'(?:if\(nr!=7\|\|nl!=2\)mexErrMsgIdAndTxt\(STRING,STRING\);)?'
    r'charerrorMessage\[2048\]=\{0\};'
    r'if\(!(?:(?:radia_mex::checked\(\w+,)|dispatchChecked\()'
    r'nl,out,nr,in,errorMessage,sizeof\(errorMessage\)\)\)\{'
    r'(?:diagnosticTrace\(STRING\);)?'
    r'mexErrMsgIdAndTxt\(STRING,STRING,errorMessage\);\}')

def inspect_text(source: str) -> list[str]:
    code = TOKENS.sub(lambda m: ' ' if m[0].startswith(('/',)) else 'STRING', source)
    covered = set()
    failures = []
    for definition in FUNCTION.finditer(code):
        start = definition.end(); depth = 1; end = start
        while end < len(code) and depth:
            depth += (code[end] == '{') - (code[end] == '}'); end += 1
        body = code[start:end-1]
        calls = list(ERROR.finditer(body))
        if not calls:
            continue
        covered.update(start + call.start() for call in calls)
        if definition[1] not in ('mexFunction', 'RadiaProjectionGateway') or not GATEWAY.fullmatch(re.sub(r'\s+', '', body)):
            failures.append(f'{definition[1]}: error exit outside recognized POD-only gateway')
    if any(call.start() not in covered for call in ERROR.finditer(code)):
        failures.append('Unrecognized error exit outside a checked function definition')
    return failures

def audit(native: Path) -> list[str]:
    failures = []
    for source in sorted(native.rglob('*')):
        if source.suffix in ('.cu', '.cpp', '.c', '.h', '.hpp', '.cuh'):
            failures.extend(f'{source.name}: {error}' for error in inspect_text(source.read_text(encoding='utf-8-sig')))
    return failures

def require_boundaries() -> None:
    failures = audit(Path(__file__).resolve().parents[1] / 'native')
    if failures:
        raise RuntimeError('Unsafe explicit MEX error gateway:\n' + '\n'.join(failures))

if __name__ == '__main__':
    native = Path(__file__).resolve().parents[1] / 'native'
    failures = audit(native)
    if failures:
        raise SystemExit('\n'.join(failures))
    print('Explicit MEX error exits satisfy the conservative POD-only gateway fence')
