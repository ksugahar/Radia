# Bounded HDiv interchange

`hdiv_hdf5.py` and `radia.gpu.readHDivH5` / `writeHDivH5` implement
`radia.hdiv-hdf5.v1`. This interchange is for independent CPU validation;
it is not a Python GPU adapter. No MAT fallback is available.

Files contain SI metadata and ordered UTF-8 label strings. MATLAB readers
return labels as a canonical row cell array, independently of JSON decoder
orientation. Label order and string values are preserved.

Numerical datasets use explicit little-endian types and file axes in
row/column order. Rank-one arrays are returned as MATLAB column vectors;
validators normalize RHS and expected vectors with `(:)` and assert matching
sizes before subtraction, preventing implicit expansion. Non-square rank-two
arrays preserve their row/column dimensions.
Complex values have separate real and imaginary components.
Each component has a SHA256 of its C-order bytes. Canonical CSR uses uint64
zero-based indices and row pointers, with explicit matrix dimensions.
MATLAB currently writes dense outputs; the independent CPU producer writes
CSR fixtures. Empty arrays and empty sparse operators are explicitly unsupported.
Readers check types, dimensions, indices, hashes and a host allocation budget.

`checkHDivHdf5` exercises the public 2-by-3 CSR fixture with three nonzeros
(different row count and value count), dense and complex output, and UTF-8
metadata. Its MATLAB execution is a required separate runtime check;
Python contract tests alone do not certify the MATLAB adapter.
