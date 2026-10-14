"""Verify the narrowly audited privacy redactions of historical native evidence.

Execution hashes remain historical. Published byte hashes are a separate
binding; numerical C++ tokens and the MATLAB generator outside its one
runtime metadata pair must remain identical.
"""
import hashlib
import ast
import json
import re

AUDITED_PATHS = {
    "validation_test/maglev/team28_axisym_fem_height_sweep.py": "python-runtime-pair-v1",
    "validation_test/maglev/team28_arnoldi_height_sweep.py": "python-runtime-pair-v1",
    "docs/maglev/demos/team28/team28_axisym_fem.py": "python-module-docstring-v1",
    "validation_test/mixed_omega_trace_reuse/run.py": "python-keyword-runtime-pair-v1",
    "validation_test/maglev/team28_axisym_fem_height_sweep.json": "json-root-runtime-redaction-v1",
    "validation_test/maglev/ecb_foster_lorentz_3d_reference.py": "python-runtime-pair-v1",
    "src/core/rad_hacapk_hdiv.cpp": "cpp-tokens-v1",
    "src/core/rad_hacapk_hdiv.h": "cpp-tokens-v1",
    "validation_test/radia_mcp/generate_motor_angle_family_mex_artifact.m": "matlab-runtime-pair-v1",
}


def text_sha256(text):
    return hashlib.sha256(text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")).hexdigest()


def cpp_tokens(text):
    # C++ phase-2 line splicing precedes comment recognition. Preserve raw and
    # quoted literals, maximal-munch punctuators and directive line boundaries.
    text = re.sub(r"\\\r?\n", "", text)
    tokens = []
    i = 0
    line_start = True
    directive = False
    operators = ("%:%:", "<<=", ">>=", "->*", "...", "<=>", "##", "::", "->", "++", "--", "<<", ">>", "<=", ">=", "==", "!=", "&&", "||", "+=", "-=", "*=", "/=", "%=", "&=", "|=", "^=", ".*", "<:", ":>", "<%", "%>", "%:")
    while i < len(text):
        char = text[i]
        if char.isspace():
            if char == "\n":
                if directive:
                    tokens.append("<directive-end>")
                directive = False
                line_start = True
            i += 1
            continue
        if text.startswith("//", i):
            stop = text.find("\n", i)
            i = len(text) if stop < 0 else stop
            continue
        if text.startswith("/*", i):
            stop = text.find("*/", i + 2)
            if stop < 0:
                raise ValueError("unterminated C++ comment")
            i = stop + 2
            continue
        if line_start and char == "#":
            directive = True
        line_start = False
        raw = re.match(r'(?:u8|u|U|L)?R"([^ ()\\\t\r\n]{0,16})\(', text[i:])
        if raw:
            ending = ")" + raw.group(1) + '"'
            stop = text.find(ending, i + raw.end())
            if stop < 0:
                raise ValueError("unterminated C++ raw literal")
            stop += len(ending)
            tokens.append(text[i:stop]); i = stop
            continue
        if char in ("\"", "'"):
            start = i; i += 1
            while i < len(text):
                if text[i] == "\\":
                    i += 2
                elif text[i] == char:
                    i += 1
                    break
                else:
                    i += 1
            else:
                raise ValueError("unterminated C++ literal")
            tokens.append(text[start:i])
            continue
        word = re.match(r"[A-Za-z_][A-Za-z_0-9]*|(?:[0-9]|\.[0-9])(?:[A-Za-z_0-9.]|(?<=[eEpP])[+-])*", text[i:])
        if word:
            tokens.append(word.group()); i += word.end()
            continue
        operator = next((op for op in operators if text.startswith(op, i)), char)
        tokens.append(operator); i += len(operator)
    return tokens


def semantic_sha256(path, text):
    method = AUDITED_PATHS[path]
    if method == "cpp-tokens-v1":
        canonical = json.dumps(cpp_tokens(text), ensure_ascii=True, separators=(",", ":"))
    elif method == "python-module-docstring-v1":
        tree = ast.parse(text)
        if not (tree.body and isinstance(tree.body[0], ast.Expr) and isinstance(tree.body[0].value, ast.Constant) and isinstance(tree.body[0].value.value, str)):
            raise ValueError("expected one leading module docstring")
        tree.body.pop(0)
        canonical = ast.dump(tree, include_attributes=False)
    elif method == "json-root-runtime-redaction-v1":
        value = json.loads(text)
        if "host" in value and not isinstance(value["host"], str):
            raise ValueError("root runtime host must be a string")
        value.pop("host", None)
        value.pop("source_privacy_redactions", None)
        canonical = json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))
    elif method == "python-keyword-runtime-pair-v1":
        old = "host=platform.node()"
        new = "platform_class=platform.system()"
        if text.count(old) + text.count(new) != 1:
            raise ValueError("expected exactly one audited Python keyword runtime pair")
        canonical = text.replace(old, "<runtime-metadata-pair>").replace(new, "<runtime-metadata-pair>")
    elif method == "python-runtime-pair-v1":
        old = '"host": socket.gethostname()'
        new = '"platform_class": platform.system()'
        if text.count(old) + text.count(new) != 1:
            raise ValueError("expected exactly one audited Python runtime pair")
        canonical = text.replace(old, "<runtime-metadata-pair>").replace(new, "<runtime-metadata-pair>")
    else:
        old = '"hostname", hostName'
        new = '"platform_class", string(computer(\'arch\'))'
        if text.count(old) + text.count(new) != 1:
            raise ValueError("expected exactly one audited MATLAB runtime pair")
        canonical = text.replace(old, "<runtime-metadata-pair>").replace(new, "<runtime-metadata-pair>")
    return text_sha256(canonical)


def verifies_redaction(path, historical_sha, published_text, record):
    return (
        path in AUDITED_PATHS
        and record.get("method") == AUDITED_PATHS[path]
        and record.get("historical_text_sha256") == historical_sha
        and record.get("published_text_sha256") == text_sha256(published_text)
        and record.get("semantic_sha256") == semantic_sha256(path, published_text)
    )
