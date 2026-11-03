"""Pybind11 and plain MEX on the same analytic Gram kernel: bitwise action on one input; HDF5 and JSON records."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess

import numpy as np

from hdiv_hdf5 import read_hdiv_h5, write_hdiv_h5


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matlab", type=Path, required=True)
    parser.add_argument("--mex-dir", type=Path, required=True)
    parser.add_argument("--matlab-package", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    from radia.gpu import _native
    tetra = np.array([[0., 0., 0.], [.4, 0., 0.], [0., .4, 0.], [0., 0., .4]])
    second = .5*tetra + np.array([.6, .1, .1])
    cells = np.concatenate([tetra.ravel(), second.ravel()])
    faces = np.concatenate([tetra[[0, 1, 2]].ravel(), second[[0, 2, 3]].ravel()])
    x = np.random.default_rng(20261103).standard_normal(4)
    gram = _native("_hdiv_cuda")._ResidentHDiv.from_analytic_constant_block(cells, faces, backend="gpu")
    try:
        python = gram.apply(x)
    finally:
        gram.close()
    source, output = args.out/"input.h5", args.out/"mex.h5"
    write_hdiv_h5(source, dict(cells=cells, faces=faces, x=x), dict(units="SI", labels=["analytic-gram"]))
    call = (f"addpath('{args.mex_dir.as_posix()}','{args.matlab_package.as_posix()}');"
            f"[f,m]=radia.gpu.readHDivH5('{source.as_posix()}');g=radia.gpu.AnalyticChargeGram(f.cells(:),f.faces(:));"
            f"y=g.apply(f.x(:));clear g;assert(radia_hdiv_device('live_count')==0);"
            f"radia.gpu.writeHDivH5('{output.as_posix()}',struct('y',y),m);")
    with (args.out/"matlab.log").open("wb") as stream:
        subprocess.run([str(args.matlab), "-batch", call], check=True, stdout=stream, stderr=subprocess.STDOUT)
    mex = read_hdiv_h5(output)[0]["y"].ravel()
    result = dict(schema="radia.hdiv-analytic-gram-bitwise.v1", source_sha=os.environ.get("RADIA_SOURCE_SHA"), charges=x.size,
                  accepted=python.tobytes() == mex.tobytes(), python_sha256=hashlib.sha256(python.tobytes()).hexdigest(),
                  mex_sha256=hashlib.sha256(mex.tobytes()).hexdigest(),
                  maximum_absolute_difference=float(np.max(np.abs(python - mex))))
    (args.out/"result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not result["accepted"]:
        raise SystemExit(f"Analytic Gram frontends differ: {result}")


if __name__ == "__main__":
    main()
