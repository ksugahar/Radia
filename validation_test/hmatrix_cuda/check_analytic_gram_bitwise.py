"""Pybind11 and plain MEX on the same analytic Gram and inverse-BH kernels: bitwise on one input; HDF5 and JSON records."""
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
    knots, coefficients, base = np.array([0., 1., 2.]), np.array([1., 0., 1., 0., 1., 3., 4., 2.]), np.array([0., .75, 6.])
    scalars = np.array([1.5, 7., 4.875, 1.])
    amplitudes = np.array([.1, 1., 1.4, 3., 100.])
    samples = _native("_hdiv_cuda")._inverse_bh_samples(knots, coefficients, base, amplitudes, *scalars, backend="gpu")
    python_samples = np.concatenate([samples[name] for name in ("secant", "tangent", "energy")])
    source, output = args.out/"input.h5", args.out/"mex.h5"
    write_hdiv_h5(source, dict(cells=cells, faces=faces, x=x, knots=knots, coefficients=coefficients, base=base, scalars=scalars,
                               amplitudes=amplitudes), dict(units="SI", labels=["analytic-gram"]))
    call = (f"addpath('{args.mex_dir.as_posix()}','{args.matlab_package.as_posix()}');"
            f"[f,m]=radia.gpu.readHDivH5('{source.as_posix()}');g=radia.gpu.AnalyticChargeGram(f.cells(:),f.faces(:));"
            f"y=g.apply(f.x(:));clear g;assert(radia_hdiv_device('live_count')==0);"
            f"law=struct('knots',f.knots(:),'coefficients',reshape(f.coefficients(:),4,[]).','integral_base',f.base(:),'cap',f.scalars(1),"
            f"'barrier_slope',f.scalars(2),'last_field',f.scalars(3),'initial_reluctivity',f.scalars(4));"
            f"s=radia.gpu.HDivBHMaterial.inverseSamples(law,f.amplitudes(:));"
            f"radia.gpu.writeHDivH5('{output.as_posix()}',struct('y',y,'samples',[s.secant;s.tangent;s.energy]),m);")
    with (args.out/"matlab.log").open("wb") as stream:
        subprocess.run([str(args.matlab), "-batch", call], check=True, stdout=stream, stderr=subprocess.STDOUT)
    written = read_hdiv_h5(output)[0]
    mex, mex_samples = written["y"].ravel(), written["samples"].ravel()
    result = dict(schema="radia.hdiv-analytic-gram-bitwise.v1", source_sha=os.environ.get("RADIA_SOURCE_SHA"), charges=x.size,
                  accepted=python.tobytes() == mex.tobytes() and python_samples.tobytes() == mex_samples.tobytes(),
                  samples_python_sha256=hashlib.sha256(python_samples.tobytes()).hexdigest(),
                  samples_mex_sha256=hashlib.sha256(mex_samples.tobytes()).hexdigest(), python_sha256=hashlib.sha256(python.tobytes()).hexdigest(),
                  mex_sha256=hashlib.sha256(mex.tobytes()).hexdigest(),
                  maximum_absolute_difference=float(np.max(np.abs(python - mex))))
    (args.out/"result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not result["accepted"]:
        raise SystemExit(f"Analytic Gram frontends differ: {result}")


if __name__ == "__main__":
    main()
