# LTspice double RAW regression

`ltspice_numdgt7.cir` generated both RAW fixtures on LAB on 2026-09-10.
LTspice executable file version: 26.0.2.1 (product version 26.0.1.1).
Two copies of the same netlist were run with `-Run -b`, adding `-ascii` for
the text reference. Both owned processes exited with code zero before files
were collected. No Radia RAW parser was used to generate the expectations.

The MATLAB test compares all names and values in binary and ASCII output with
relative tolerance 1e-12 and absolute tolerance 1e-15. It also requires the
binary header's `double` flag. The fixture is small and does not require an
LTspice installation when running the regression.

Run `runtests('tests/matlab/test_ltspice_data_safety.m')` in MATLAB.
The scoped `MATLAB LTspice data safety` CI runs this test through a dedicated
MATLAB Engine on mdx without installing LTspice or building MEX binaries.
The initial 13-test run passed through Python's MATLAB Engine. Layout tests
also exercise mixed float/double and complex records, CRLF, unsupported flags,
truncation, extra bytes, and fail-loud transient state injection.

Known boundaries: FastAccess and unrecognized RAW flags are rejected, not
guessed. Offset is preserved as a header property, not automatically added to
the time axis. Hierarchical inductor state reinjection is unsupported and
raises instead of silently dropping currents. This change does not repair
process lifecycle, dependency staging, SimRunner, or the per-step S-Function
architecture; those remain separate review items.
