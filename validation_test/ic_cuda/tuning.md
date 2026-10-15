# GPU IC synchronization and triangular-solve tuning

Natural-order fixed-shift IC(0) retains the original free-DOF factor action.
Validated constrained zero couplings are omitted from the factor graph; the
original free-DOF CSR equations remain available for the final true residual.
The native check uses the projected original operator (constrained diagonals
= 1, constrained RHS = 0). An explicit upper
CSR retains the transpose values for non-transpose SpSV. Device scalars and
cached eight-step CUDA graphs batch iteration control; device state freezes
updates at first convergence or the exact iteration limit. There is no fallback.
The context owns the graph, stream and buffers and retains the existing device
restoration, typed-handle and unloadability contracts. cuSPARSE 12.8 documents
[SpSV graph capture support](https://docs.nvidia.com/cuda/archive/12.8.0/cusparse/index.html#cusparsespsv).

Historical overlay accepted 2026-10-15: 360 independent original free-DOF
solution/residual checks across
10,460 / 32,233 / 106,272 DoFs, CPU IC/AMS 1/16 threads and GPU IC/AMS.
Maximum vector error 1.122e-11; true residual 5.621e-12.
Independent CPU factor application agrees within 6.1e-16 at every size.
Both GPUs passed 25 rejection gates and lifecycle/unload checks. The added
iteration-limit boundary test checks convergence at the exact required count
and rejection one step earlier. 33 focused package/parity tests passed.

The base 0.1.5rc1 wheel remains hash-verified; the optimized IC is an explicit
experimental overlay built from clean source, not a newly published wheel.
Full six-route [baseline](results/gpu-ic-baseline-scaling.json),
[tuned](results/gpu-ic-tuned-scaling.json), [phase/comparison](results/ic-tuning-comparison.json)
and [independent validation](results/ic-tuning-validation.json) JSON records are retained.

| Total DoFs | Old GPU IC reuse (s) | Tuned GPU IC reuse (s) | CPU16 IC reuse (s) | Tuned GPU IC fresh (s) | CPU16 IC fresh (s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10,460 | 0.2647 | 0.2217 | 0.2554 | 0.4269 | 0.3465 |
| 32,233 | 0.4326 | 0.2569 | 0.4216 | 0.7320 | 0.6378 |
| 106,272 | 2.0862 | 1.6665 | 1.9803 | 3.2338 | 2.6733 |

At 106,272 DoFs the microprofile triangular pair changes from 2.688 to 2.157 ms.
Factor levels drop from 658 to 514; numerical factor construction
changes from 0.1394 to 0.1100 s.
Triangular dependency depth is the main remaining bottleneck; host scalar
readbacks are secondary on this GPU. The profiles are separate from full
acceptance medians. CPU1/16, GPU AMS, all setup/transfer phases, reused RHS
scales and sampled process/device-endpoint memory remain in the full JSON.

Retained explicit IC buffers at this size grow from 70.88 to 76.73 MB;
setup owned-buffer bounds grow from 81.90 to 89.57 MB.
These IC-only bounds exclude caller GPU inputs/outputs and opaque CUDA/library/graph
allocations. Assembly bounds and process/device endpoints are separate JSON
fields; these are not whole-device transient peaks. At the largest size, median
host FE basis setup is 1.382 s within 1.553 s GPU setup, versus 0.671 s CPU16 setup.
Fresh GPU IC remains slower than CPU16 IC because host FE
setup is serial and graph/factor construction must be paid once. GPU AMS
remains substantially faster for these HCurl systems.

Review fix: the context now retains an explicit 4 MiB cuBLAS workspace through
graph/handle teardown. The historical binary hashes, times and memory figures
above precede that fix; they do not certify the rebuilt binary. Revalidate the
new binary and exact formal wheel separately. The native `profile` command is
diagnostic-only, excluded from acceptance timing medians.
