# Common IH and IC candidate acceptance plan

This is integration preparation, not release qualification. The current published versions remain unchanged. Claude review and findings are pending.

The isolated integration combines current main, IH temperature-BH source 266580d23fb7b799c943020a285efdc8d8cc97eb and IC source 63ef045b120edd36260c1921d99f6ba885e479ce. All native modules must be rebuilt from the final clean integration commit after review fixes. Historical binaries and overlay measurements cannot qualify that wheel.

## Required acceptance sequence

1. Record each actual Claude finding as accepted-fix, rebutted with evidence, or deferred with owner and release impact. Review the final integrated source, including packaging and CI.
2. Build csr_pcg, ic, ams, hmatrix, projection, fmm and ih sequentially in the isolated CI environment. Verify every source and binary receipt, including RADIA_IH_BUILD, against the same final commit. Freeze one wheel and its SHA256.
3. Install that wheel in isolated environments on both declared acceptance machines. Verify archive manifest, installed manifest and actual MATLAB `which` paths and loaded MEX hashes. Reject a source, binary or installed-manifest mismatch even if the version matches.
4. Run independent CPU parity and original-equation residual checks for IC/AMS/HDiv and IH. IH uses checkTemperatureBH.m plus check_temperature_bh_cpu.py: five updates, accepted temperature/heat/EM/Zs retention after rejected candidates, reset, invalid tables, temperature-range errors, repeat simulation and failure cleanup. Shared geometry assembly is explicitly disclosed; this is not independent FEM certification.
5. Use the installed IH builder and samples for a separate Simulink check: typed ports, official Toolkit model_check, save/close/reopen, repeated run, actual installed handles and zero live contexts. Run existing common SLX and LTspice gates as well. CI build inclusion alone does not satisfy this step.
6. FMM uses explicit order 8 and separation 2, independent direct CPU reference, declared operator applicability and relative L2 <=1e-8 for the declared numerical case. Do not extend this tolerance claim to every analytical physical case.
7. After checking queues, process CPU deltas, GPU and host memory, run one job at a time for fresh same-host transfer-inclusive measurements. Record fixed CPU threads, first-call and warm ABBA samples, setup/upload/evaluate/download/destroy boundaries, exact owned buffers separately from sampled and opaque peaks. No cold-process claim without fresh processes. Retain GPU reserve max(2 GiB,20%).
8. Recover raw success/failure evidence and hashes privately. Public evidence uses neutral machine roles and independently authored public inputs. Validate exact-source Native Release, final CI, machine gates and independent review before publication; publish the accepted wheel bytes without rebuilding.

## Remaining preparation gaps

The CI now selects the IH exact-commit artifact for packaging, but the new run_installed.py IH numerical/CPU/Toolkit runner and pinned official Toolkit CI invocation still need execution against the final wheel before release. The reviewed B1 power/sign gate is not yet integrated; B3 remains open. The new final integrated build and both-machine installed-wheel runs are pending. No qualification records may be synthesized from older standalone IH results.

IH research limitations: flat P1 closed genus-zero body, bounded complete reaction assembly, user-supplied coil impedance, temperature-interpolated single-valued BH with no hysteresis history/minor-loop loss, staggered thermal coupling and bounded surface-layer thermal model. Independent FEM and experimental validation remain unverified. The synthetic high-current test is not a laboratory operating instruction. IC tuned reuse may improve speed while fresh setup remains slower and owned memory increases; disclose the measured scope rather than claiming universal acceleration.
