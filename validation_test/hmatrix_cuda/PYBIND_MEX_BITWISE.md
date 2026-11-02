# HDiv device same-kernel frontend parity

Gate (1) of the CPU/GPU execution policy: pybind11 and the plain MEX on the same
kernel, bitwise. It is not CPU/GPU accuracy (separate tests compare the device
route with the CPU-assembled equation) and not a performance measurement.

## Build

`tools/build_python_gpu.py` compiles the device HDiv kernels once with nvcc
(`--fmad=false`, no fast math) into one object. The pybind11 module and
`packages/radia-simulink/tools/build_hdiv_mex.py` (CMake `matlab_add_mex`, no
mexcuda) both link that object; the MEX build receipt
(`radia.hdiv-mex-build.v1`) records the object hash and the pybind11 binary hash
of the same source SHA. Identical compiler controls therefore hold by
construction; the comparator checks the object hash instead of a flag list.

## Case `device-flat-tet`

`hdiv_device_parity.py make-input` writes the canonical input HDF5: the NGSolve
device-construction tables of a Netgen 6:1:1 bar (BDM2), device ACA with leaf
size 32 so low-rank leaves are exercised, the production inverse-BH payload of a
fixed soft-iron table, the input vectors and every solver control. Both
producers read only this file and hash its bytes.

- `hdiv_device_parity.py pybind` runs `_hdiv_cuda`.
- `produceHDivDeviceMexRecord.m` runs `radia_hdiv_device` (index base and
  layout: HDF5 row-major tables are read flat, the inverse-BH coefficients as
  the MATLAB `(knots-1)x4` matrix the gateway expects).

Outputs, FP64 in `/outputs`: `apply` (B^T G B on a fixed vector),
`linear_solution` (mass-diagonal PCG), `nonlinear_solution` and
`nonlinear_residual_history` (reference-diagonal frozen-tangent Newton).
Both producers fail if ACA built no low-rank leaf; the MEX producer also fails
if a device handle remains live.

## Receipts and comparison

Each record carries a UTF-8 `receipt_json` attribute: `frontend` (`pybind11` or
`mex`), the shared build identity (`source_sha`, `kernel_object_sha256`,
`cuda_compiler`, `python_binary_sha256`, `mex_binary_sha256`), `input_sha256`,
`device_uuid`, `architecture`, `case`, `output_names`, and the hash of the binary
actually loaded (`artifact_sha256`). The comparator requires equal shared fields,
loaded binaries equal to their build receipts, and exact shape, dtype and bytes
(signed zero included); no tolerance.

```sh
python validation_test/hmatrix_cuda/check_pybind_mex_bitwise.py \
  --python-record python.h5 --mex-record mex.h5 --report parity.json
```

Record the input, both records, the report and the build receipts with their
hashes in private evidence. Windows only (the plain MEX targets the
radia-simulink platform); Linux pybind11 builds are outside this gate.
