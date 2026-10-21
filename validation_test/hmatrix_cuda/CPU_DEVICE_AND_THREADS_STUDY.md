# CPU/device execution and MATLAB thread control: study notes

Status: design study; no new CPU/GPU route or numerical acceptance is claimed.

## NGSolve's execution model

NGSolve constructs FE spaces and operators with its ordinary host APIs. Its
Poisson tutorial constructs device counterparts explicitly with
`CreateDeviceMatrix()` and `CreateDeviceVector(copy=True)`, then applies a
`DeviceCGSolver`. This makes operator/vector placement visible rather than
claiming that importing a module accelerates every FE operation.
[Official device solver tutorial](https://docu.ngsolve.org/latest/i-tutorials/unit-5.5-cuda/poisson_cuda.html).

The newer common GPU layer provides `GetGPUDevice()`, `GetCPUDevice()`, buffers,
queues and explicit H2D/D2H transfer. `Finish()` synchronizes launches.
`ngsolve.gpu` can select a host reference device, so inspect the actual backend;
a non-null device does not establish CUDA execution. Its shared-source CPU
reference is useful for implementation differential tests but is not Radia's
independent physics oracle. Runtime compilation in this tutorial is study-only;
it does not authorize NVRTC in Radia's P1 scope.
[Common GPU tutorial](https://docu.ngsolve.org/nightly/i-tutorials/unit-5.6-gpu/commonGPU.html),
[pinned selection source](https://github.com/NGSolve/ngsolve/blob/736417e0624fa555b5db82281382860ef4a8c9bc/python/gpu.py).

## Three distinct MATLAB parallel controls

| Control | Meaning | Radia qualification needed |
|---|---|---|
| `maxNumCompThreads(n)` | Maximum computational threads used by MATLAB; returns the previous maximum; setting is session-local | Record requested/reported limit. It is not proof that an external MEX OpenMP runtime obeys the limit. |
| `parpool("Threads", n)` | Creates a pool of thread workers for supported parallel MATLAB work | Worker count is distinct from BLAS or MEX thread count. Do not create a pool implicitly or reuse unqualified persistent handles concurrently. |
| Native OpenMP team size | Threads executing a MEX-owned host parallel region | Bound explicitly at each supported region; validate exception containment, MATLAB API restrictions and absence of nested oversubscription. |

[MathWorks computational-thread control](https://www.mathworks.com/help/matlab/ref/maxnumcompthreads.html),
[MathWorks parallel pools](https://www.mathworks.com/help/parallel-computing/parpool.html).

An intended MATLAB usage pattern is a scoped computational-thread setting:

```matlab
previous = maxNumCompThreads(2);
restoreThreads = onCleanup(@() maxNumCompThreads(previous));
% Call an explicitly selected and qualified CPU or GPU operator here.
```

This is an illustrative pattern, not a benchmark or an implemented Radia API.
Avoid changing a user's session persistently. Thread-worker concurrency needs
separate qualification before using a session-wide setter. A job CPU envelope
and actual measured thread usage remain separate evidence.

## Proposed application of the lessons

1. Keep NGSolve FE assembly, orientation and coefficient semantics authoritative.
2. Select host/device operator and vector placement explicitly; preserve a
   resident operator over repeated solves and transfer only at declared boundaries.
3. For MATLAB, compile CPU and GPU paths through mexcuda. The CPU branch must
   not initialize a CUDA context. Do not change execution after allocation by
   silently copying or falling back.
4. Expose CPU resource limits separately from device choice once the native
   thread-control capabilities are qualified. CPU BLAS, native OpenMP and pool
   workers must not each consume the full process budget simultaneously.
5. Record selected route, requested/effective thread settings, transfer/setup/
   solve/recovery time, true residual, memory and object lifetime. Benchmark
   timing requires synchronization at the relevant device boundary.

P1 remains the bounded memory/failure-lifecycle qualification. This study does
not expand its execution scope, add production Python GPU or enable NVRTC.


## CPU context trace acceptance protocol

Status: proposed protocol; no CPU MEX no-context acceptance is claimed.
The current CSR gateway calls `mxInitGPU()` at the start of `dispatchCsr`.
A CPU selector inside that initialized GPU path cannot satisfy the contract.
Separate host argument checking and CPU computation from GPU initialization,
array conversion, device allocation and launches before qualifying either route.

Observe fresh owned MATLAB processes from before MEX loading through explicit
CPU execution, invalid-input handling, repeated calls, `clear mex` and shutdown.
Run separately with a device present and genuinely absent. Hiding a device with
an environment variable is not evidence of physical absence; changing a shared
device or driver is outside this protocol.

The observer must report its version, driver/runtime API coverage, dynamically
resolved entrypoint coverage, dropped events and the complete observation
interval. Capture context creation and primary-context retention, including
transitive calls. An empty or incomplete trace is inconclusive. CUDA DLL loading
alone is not evidence that a context exists, and device-monitor absence alone
is not proof that no context was created.

Use a separate explicit GPU process as a positive control: the same observer
must detect its context. Observer-only controls must establish that tracing
does not itself initialize CUDA inside the CPU target. Do not call `gpuDevice`,
`gpuArray` or GPU initialization/probing routines merely to test the CPU route.
No compatible observer has yet been qualified for this protocol.

Bind raw traces, commands, observer configuration and numerical results to the
same source and binary hashes. Require independent CPU agreement and the true
original-system residual in addition to the context result. Distinguish PASS,
FAIL and INCONCLUSIVE; unsupported tracing never becomes a pass. Keep requested
MATLAB computational threads, native effective threads and pool-worker counts
separate, as above. These checks do not establish runtime redistribution rights.
