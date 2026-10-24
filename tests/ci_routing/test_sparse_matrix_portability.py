"""Guards: no Python source selects sparse storage by compiler-specific names.
Assumes: a Git checkout with UTF-8 tracked Python sources; no native imports.
"""
from pathlib import Path
import re
import subprocess


def test_no_compiler_specific_sparse_matrix_references():
    root = Path(__file__).resolve().parents[2]
    paths = subprocess.check_output(
        ["git", "ls-files", "--", "src", "packages", "tests", "validation_test"],
        cwd=root, text=True).splitlines()
    forbidden = re.compile(r"\bSparseMatrix(?:Symmetric)?(?:double|d)\b")
    references = []
    for relative in paths:
        if not relative.endswith(".py"):
            continue
        for number, line in enumerate((root / relative).read_text(encoding="utf-8").splitlines(), 1):
            if forbidden.search(line):
                references.append(f"{relative}:{number}: {line.strip()}")
    assert not references, "compiler-specific sparse matrix references:\n" + "\n".join(references)
