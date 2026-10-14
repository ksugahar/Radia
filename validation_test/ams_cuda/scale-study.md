# Installed 0.1.4 FEM/AMS scaling

The independently generated affine unit-cube series has 10,460 / 32,233 /
106,272 total HCurl DOFs and 7,571 / 26,839 / 91,533 free DOFs. The formulation
is `curl(u)*curl(v)+u*v`, unit reluctivity and positive unit gauge mass,
constant SI current density `(1,0,0)`, and homogeneous tangential boundaries.
These elementary cases do not certify a winding, material contrast, nonlinear
device or frequency-domain HX problem.

[Machine-readable results](results/scale014-a100-20261013.json) contain every
sample, construction and solve phase, transfer counters, memory samples,
independent errors, original-equation residuals and actual loaded MEX hashes.
The exact released `radia-gpu==0.1.4` wheel SHA-256 remains
`6b6cc5034365a1aebd6dc26c1d17b423ac509e874e021f70549905ed6f8ec154`.
No GPU kernel or package was rebuilt for this evaluation.

The idle A100 host also ran the independent Python CPU routes on its Xeon
Platinum 8368 allocation (18 reported cores). CPU NGSolve threads were 1 or
16, CPU BLAS threads 1, CPU AMS coarse factorization threads 1, and GPU host
setup threads 1. Python built a separate NGSolve `sparsecholesky` accuracy
reference for each mesh. The GPU consumed the checked `.vol` directly through
MATLAB/mexcuda; CPU-generated operators or solver snapshots were not inputs.

## Results

Warm medians in seconds, including construction and transfers for the first
solve of each newly constructed system:

| Total DOFs | CPU AMS, 1 thread | CPU AMS, 16 threads | GPU AMS | CPU IC, 16 threads |
| ---: | ---: | ---: | ---: | ---: |
| 10,460 | 0.2584 | 0.2767 | 0.2323 | 0.3457 |
| 32,233 | 0.8159 | 0.5157 | 0.5881 | 0.6530 |
| 106,272 | 2.9079 | 1.3152 | 1.8632 | 2.6938 |

Matrix and preconditioner retained, zero initial guesses, changing only current
scale; medians across scales `2, 0.5, 1`, including GPU result download:

| Total DOFs | CPU AMS, 1 thread | CPU AMS, 16 threads | GPU AMS | CPU IC, 16 threads |
| ---: | ---: | ---: | ---: | ---: |
| 10,460 | 0.04676 | 0.02133 | 0.01811 | 0.25711 |
| 32,233 | 0.14711 | 0.03514 | 0.03986 | 0.42811 |
| 106,272 | 0.64653 | 0.08311 | 0.07462 | 2.01482 |

| Total DOFs | CPU AMS iterations | GPU AMS iterations | CPU IC iterations |
| ---: | ---: | ---: | ---: |
| 10,460 | 20 | 20 | 232 |
| 32,233 | 21 | 22 | 364 |
| 106,272 | 26 | 27 | 728 |

CPU AMS iteration counts match at 1 and 16 threads. IC uses the same native
fixed diagonal shift 1.05, natural ordering, no diagonal scaling and no ABMC
at both thread settings. Thus this is a comparison against that declared IC
configuration, not a claim about an optimally tuned parallel IC. GPU IC is
unavailable in the selected wheel and was not simulated or substituted.

The series demonstrates a modest increase in AMS iterations and a much larger
increase for this IC configuration. It also shows why GPU solve throughput
and first-solve latency must be distinguished: at 106k DOFs, GPU AMS reuse is
1.11 times faster than CPU16 AMS, but GPU fresh construction plus solve is
1.42 times slower. In one representative 106k GPU sample, host FE basis work
took 1.391 s, CUDA upload/assembly 0.0112 s, and the subsequent first solve
0.0763 s. Host preparation is the main first-solve cost in this configuration.

## Measurement and memory boundaries

Each route has one separately recorded first initialized call, followed by
two forward/reverse balanced cycles (four measured system constructions per
route and mesh). Every construction performs four zero-guess solves at
current scales `1, 2, 0.5, 1`. Independent vector parity and original-equation
residual checks run outside the timers for every result. Admission limits
remain relative solution error `1e-7` and true residual `1e-6`; observed
measured maxima are approximately `1.12e-11` and `5.62e-12` respectively.

MATLAB startup and GPU initialization are separate. Initialized first calls
and warmed medians are not OS/device cold-start measurements. Each timed
sample requires CPU utilization at most 10%, GPU at most 5%, and declared
free-memory reserves. A short utilization tail stopped an initial run; its
results and original runner were retained. The final runner keeps those
thresholds and records up to five bounded settling checks before admission.

At 106k, sampled absolute process peak RSS was approximately 435--438 MiB for
CPU AMS, 416--417 MiB for CPU IC and 1,412 MiB for MATLAB/GPU AMS. These are
process footprints, not equal-scope kernel allocation comparisons: the CPU
process retains the independent reference operator, and MATLAB sampling also
includes serialization around the API. Before/after values and sampled RSS
increases are in the JSON. A representative GPU AMS setup reports a native
phase-owned peak bound of 198,863,784 bytes, and separate PCG workspace and
coarse-action cache counters. Phase bounds and device-free endpoints do not
certify an aggregate transient whole-device peak.

## Reproduction

`prepare_scale_meshes.py --out <new-directory> --checker <checked-check.py>`
generates the public meshes and versioned label receipts. Re-running this
generator reproduced all three evaluated `.vol` files byte for byte with
NGSolve/Netgen 6.2.2607. Preserve the exact checker, contract and native
reference identities when staging a validation job.

`measure_scale_installed.py --job <verified-job> --wheel <exact-0.1.4-wheel>
--out <new-results-directory> --cycles 2` runs the study on the selected idle
GPU evaluation host. The job manifest binds meshes, receipts, wheel verifier,
MATLAB helper, DLLs and independent CPU modules. The old CPU AMS bridge's
verified build/source receipt and the IC module's actual SHA are recorded.
These CPU modules are validation dependencies, not GPU runtime prerequisites.
Only owned MATLAB Engines are closed. All successful and stopped-run raw
evidence was recovered and hash-verified privately before scratch cleanup.
