# Third-party notices

Radia's own component licenses and the HACApK notices are in `LICENSE`.

Optional Python native bindings incorporate pybind11 headers. The BSD 3-Clause notice is reproduced in `licenses/PYBIND11-LICENSE.txt` from pybind11 3.1.0.

The HDiv host staging interface uses NGSolve and Netgen headers and dynamically links their libraries. Their LGPL 2.1 license texts, taken from the qualified 6.2.2608 dependency wheels, are reproduced in `licenses/NGSOLVE-LICENSE.txt` and `licenses/NETGEN-LICENSE.txt`. These dependencies are installed separately. Corresponding upstream sources are available from [NGSolve](https://github.com/NGSolve/ngsolve) and [Netgen](https://github.com/NGSolve/netgen); Radia binding and build sources are identified by the package's source-commit provenance. These notices do not restrict replacement or modification of those dependencies under their licenses; native ABI compatibility remains a technical requirement.

Python GPU execution obtains CUDA runtime libraries from the separately installed NVIDIA packages. Their license files remain in those packages. No NVIDIA, MATLAB, NGSolve or Netgen DLL is copied into this wheel. See the [CUDA 12.8 license agreement](https://docs.nvidia.com/cuda/archive/12.8.0/eula/index.html) and each installed package's license. MATLAB and Parallel Computing Toolbox remain execution requirements for the MATLAB/MEX interface, not for the qualified Python GPU interface.

This notice records dependency provenance and attribution. It is not a substitute for the distribution acceptance required by the release process.
