#!/usr/bin/env python3
"""
Fix import paths in test files after reorganization
"""

import os
import re

def fix_test_file(filepath):
    """Fix the build directory path in a test file"""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # Pattern 1: build_dir = os.path.join(os.path.dirname(__file__), 'build', 'lib', 'Release')
    old_pattern1 = r"build_dir = os\.path\.join\(os\.path\.dirname\(__file__\), 'build', 'lib', 'Release'\)"
    new_pattern1 = """# Add build output directory to path
# Tests are now in tests/ subdirectory, so go up one level to find build/
project_root = os.path.dirname(os.path.dirname(__file__))
build_dir = os.path.join(project_root, 'build', 'lib', 'Release')"""

    if re.search(old_pattern1, content):
        # Remove the comment line before it if it exists
        content = re.sub(r"# Add build output directory to path\s*\n" + old_pattern1,
                        new_pattern1, content)
        if "project_root" not in content:
            content = re.sub(old_pattern1, new_pattern1, content)

    # Pattern 2: sys.path.insert(0, 'build/lib/Release')
    old_pattern2 = r"sys\.path\.insert\(0, 'build/lib/Release'\)"
    new_pattern2 = """project_root = os.path.dirname(os.path.dirname(__file__))
build_dir = os.path.join(project_root, 'build', 'lib', 'Release')
sys.path.insert(0, build_dir)"""

    if re.search(old_pattern2, content):
        content = re.sub(old_pattern2, new_pattern2, content)

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

    print(f"[OK] Fixed: {filepath}")

# Fix all test files
test_dir = r'S:\Visual_Studio\02_Visual_Studio_2022_コマンドライン_コンパイル\04_Radia\tests'

for filename in ['test_simple.py', 'test_radia.py', 'test_advanced.py', 'test_parallel_performance.py']:
    filepath = os.path.join(test_dir, filename)
    if os.path.exists(filepath):
        fix_test_file(filepath)

# Fix benchmark files
benchmark_dir = os.path.join(test_dir, 'benchmarks')
for filename in os.listdir(benchmark_dir):
    if filename.endswith('.py'):
        filepath = os.path.join(benchmark_dir, filename)
        fix_test_file(filepath)

print("\n[OK] All test imports fixed!")
