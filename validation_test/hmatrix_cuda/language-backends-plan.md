# Separate language bindings from execution backends

Status: proposed migration following the user's 2026-10-18 request. The new
MATLAB CPU and Python GPU routes below are **not implemented or accepted**.

| Language | CPU | CUDA GPU |
| --- | --- | --- |
| Python | Existing `radia` / pybind11 / NGSolve | Proposed optional native extension owned by `radia-simulink` |
| MATLAB | Proposed thin native CPU entry owned by `radia-simulink` | Existing `radia-simulink` / mexcuda |

Keep the CPU `radia` distribution independent of CUDA and MATLAB. Keep all
MATLAB integration in `radia-simulink`. Python GPU users should call a native CUDA
binding, rather than require a MATLAB Engine session. This requires a separate
supported build and redistribution check; existing mexcuda artifacts do not
establish standalone Python portability or licensing clearance.

The first useful vertical slice is the bounded HDiv geometry/material workflow:
shared HACApK, constitutive energy/gradient/tangent, row-major storage, residual
certification and native CUDA kernels; thin language-specific gateways. NGSolve
continues to own FE spaces, orientation, mapped quadrature and charge maps.
Do not duplicate these semantics in a Python or MATLAB implementation.

Use explicit backend selection and fail on an unavailable backend. Preserve
Python caller-owned TaskManager and MATLAB-owned host scheduling; never link
NGSolve TaskManager into a mexcuda module. Contexts own their device, stream,
buffers and lifetime; checked handles or native language-owned objects replace
raw pointer APIs. A Python CUDA route needs explicit interaction with external
CUDA contexts and stream ownership before accepting other libraries' arrays.

Implement one thin slice first, then compare all four routes on the same public
mesh and BH law. Independent CPU numerical comparison, original-equation true
residuals, unsupported-input errors, reload/destruction, transferred bytes,
setup/solve/release times and memory budgets remain mandatory. A second binding
to the same kernel is not a second independent numerical reference.

The existing GPU prepared-context optimization remains useful: immutable FE and
material maps plus geometry H-matrix can be shared behind either language
binding. Distinguish reusable preparation from per-solve and complete one-shot
timing in both APIs. No release or new-backend availability is claimed here.
