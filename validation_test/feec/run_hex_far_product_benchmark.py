"""Run the separately compiled native microbenchmark and retain its provenance.

Compile benchmark_hex_far_product.cpp as a shared library with optimized C++17,
the selected Python's pip mkl-devel include/lib directories, and mkl_rt.
This measures a contraction kernel, NOT end-to-end magnet performance.
"""
import argparse
import ctypes
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from datetime import datetime, timezone


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compiler-flags', required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    code = '''import ctypes,os,sys
from pathlib import Path
if os.name == 'nt':
    directory = os.add_dll_directory(str(Path(sys.prefix)/'Library'/'bin'))
lib = ctypes.CDLL(sys.argv[1])
lib.run_benchmark.restype = ctypes.c_int
raise SystemExit(lib.run_benchmark())
'''
    run = subprocess.run([sys.executable, '-c', code, str(args.library.resolve())],
                         text=True, capture_output=True, timeout=120, check=True)
    data = json.loads(run.stdout)
    data.update(schema='radia.hex-far-product.microbenchmark.v1',
                timestamp=datetime.now(timezone.utc).isoformat(),
                machine=platform.node(), platform=platform.platform(),
                python=sys.version, mkl=importlib.metadata.version('mkl'),
                compiler_flags=args.compiler_flags,
                scope='synthetic kernel only; not full Gram or solve acceptance')
    data['source_sha256'] = {}
    for path in (root/'src/core/rad_hex_far_product.h',
                 Path(__file__).with_name('benchmark_hex_far_product.cpp'), args.library):
        data['source_sha256'][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, indent=2))


if __name__ == '__main__':
    main()
