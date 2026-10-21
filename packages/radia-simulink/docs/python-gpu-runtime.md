# Python GPU runtime dependencies

The Windows Python GPU route uses the NVIDIA CUDA 12.8 runtime dependencies selected by `radia-simulink[gpu]`. MATLAB and MATLAB Engine are not runtime prerequisites for this route. Kernel compilation remains a separate mexcuda build operation.

After explicit GPU selection, the Python adapter validates the installed runtime package versions and registers their DLL directories. CPU dispatch does not import this loader. Missing or different versions raise an error; there is no MATLAB runtime search or CPU fallback. Before registering paths, the adapter also rejects a CUDA dependency already loaded from another location: registering a new directory cannot replace a library already selected by Windows. Directory handles remain alive until the Python process exits. This path registration does not itself create a CUDA context.

The qualified component versions follow NVIDIA's [CUDA 12.8 GA component table](https://docs.nvidia.com/cuda/archive/12.8.0/cuda-toolkit-release-notes/index.html): cudart 12.8.57, cuBLAS 12.8.3.14, cuSOLVER 11.7.2.55, cuSPARSE 12.5.7.53 and nvJitLink 12.8.61. A compatible NVIDIA driver and the qualified Python/NGSolve/native ABI are also required. NVRTC and an external CUDA compiler are not dependencies.

NVIDIA runtime libraries are supplied by their separate packages, with their own license terms; they are not copied into the Radia wheel. Technical runtime verification does not grant distribution rights. Inclusion of the optional HDiv binding remains subject to its artifact, dependency, license and final release gates.

## Selected resident AMG payload

For a wheel containing the qualified resident AMG component, import
`ResidentAMG` from `radia_simulink.resident_amg` and select `backend="gpu"`
explicitly. The adapter takes an already constructed, checked CSR hierarchy and
its coarse inverse; hierarchy construction and physical assembly remain with
the caller. CPU solver workflows remain separate in `radia`.

The adapter checks the installed component receipt and binary hash before
registering NVIDIA runtime directories. Importing the Python module alone does
not import NGSolve, MATLAB Engine or CUDA. Missing qualification raises; it does
not select a different backend. Repeated `solve` calls reuse the owner. Use a
context manager or call `close()` to release it; repeated close is harmless and
solve after close raises.

Packaging selects this component with `RADIA_RESIDENT_AMG_BUILD`. Staging
requires matching clean-source build and runtime qualification records,
independent CPU numerical comparisons, lifecycle checks and distribution
evidence. Successful compilation alone does not qualify a resident AMG wheel.
