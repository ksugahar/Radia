#!/usr/bin/env python
"""Run debug_quick.py and capture all output."""

import subprocess
import sys

result = subprocess.run(
    [sys.executable, 'debug_quick.py'],
    cwd=r's:\Radia\01_GitHub',
    capture_output=True,
    text=True
)

print("=== STDOUT ===")
print(result.stdout)
print("=== STDERR ===")
print(result.stderr)
print("=== RETURN CODE ===")
print(result.returncode)
