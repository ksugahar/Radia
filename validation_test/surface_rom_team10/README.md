# Surface-trace POD for nonlinear transients with a changing excitation pattern

This lane asks whether a fixed-basis reduced model of a nonlinear transient
eddy-current problem gains anything from building its basis on the steel
surface trace.  The answer here is yes, but only when the excitation pattern
changes.  On the two-coil TEAM 10 variant below, the surface-trace basis meets
a 10 % target with 12 columns and runs one transient in 0.2-0.3 s.  Plain
snapshot POD needs 30 columns and about 5 s.  A residual-driven greedy
selection finds the missing excitation pattern without solving any candidate
full model.  With it, two cheap training runs (about 550 s together) are
enough; one exact full-model run takes about 2200 s.

The lane is a research study run by hand on a compute host. It is not a
production accuracy guarantee or a pytest lane. Current direct solves use
SparseCholesky without a fallback. The archived 2026-09-27 results were
produced with the earlier direct solver; their source hashes and measured
numbers are preserved. They do not validate the migrated implementation or
its performance. The separate migration smoke checks constrained true
residuals and a short nonlinear trajectory, not the full p3 campaign.

## Model

- **Geometry and material:** TEAM Problem 10 (Nakata et al.; specification
  as in Yan and Jin, PIER 153, 2015).
  - steel channels and plate, sigma = 7.505e6 S/m;
  - TEAM 13 B-H table;
  - quarter model: z = 0 natural symmetry, y = 0 tangential-A Dirichlet;
  - Kelvin open boundary through `radia.kelvin_geometry` /
    `radia.kelvin_material`.
- **Two coils** (`TEAM_COIL2=1`):
  - the TEAM 10 main coil;
  - a small independently driven web coil (centre x = 123.7 mm, 400 AT).
  - Because the steel saturates, a combined excitation is not the sum of the
    single-coil responses.
- **Full model (FOM):**
  - A-phi formulation (HCurl without gradients on all space x H1 on the
    steel, phi scaled by dt);
  - backward Euler, energy-minimising Newton with line search;
  - p3 at 6 mm, 329 183 DOF;
  - 60 steps over 120 ms.
  - The p3 / 6 mm FOM was checked against p3 / 3 mm (1-2 %) and against a 1D
    slab reference before this study.
- **Reduced model:** u = Phi q + E i(t).
  - **E:** one exact air solution per coil.
  - **Phi:** split into an A block and a phi block, each orthonormalised in
    the energy metric.  This split is required for stable time stepping.
  - **Nonlinear steel term:** ECSW hyper-reduction (non-negative NNLS
    weights, numpy/BLAS online kernel).
- **Bases compared at equal column count:**
  - `surface_lift`: POD of the steel-surface A trace, extended exactly into
    the air and by a nu_ref/sigma lifting operator into the steel, plus
    interior POD;
  - `surface_zero`: the same, with a zero interior continuation;
  - `global`: plain snapshot POD;
  - `global_src`: plain POD plus the exact source columns E.
- **Error measure:** relative steel-B error over all time steps,
  sqrt(sum |B_rom - B_fom|^2 / sum |B_fom|^2), against the exact 60-step FOM.

## Waveforms

The two-coil waveforms are in `transient.waveform`.  The rise waveform is
162 turns x 5.64 A x (1 - exp(-t/50 ms)); the pulse waveform decays after
15 ms.

| role | names |
|---|---|
| tests (exact FOM, 60 steps) | `c2_rise` (web coil only), `mix_add` (both, same sign), `mix_oppose` (both, opposite sign) |
| training pool (cheap FOM: 30 steps, one Newton step) | `c1_rise`, `c1_pulse`, `c2_pulse`, `c2_pulse_x2`, `mixp_add`, `mixp_oppose`, `c1_rise_half` |
| trap candidate | `c1_rise_neg` = -`c1_rise` exactly (odd B-H law), so a sound indicator must never pick it |

## Results (hibino / mdx2, 2026-09-27)

### 1. The excitation pattern must be in the training set, and the surface basis uses it best

Source: `results/team10_pattern_study_p3.json`.  The table shows the
worst-case steel-B error over the three tests with ECSW 3e-3.

| training | basis | k10 r5 | k16 r8 |
|---|---|---|---|
| S1 = `c1_rise` + `c1_pulse` (web coil never seen) | all bases | 29-76 % | 29-76 % |
| S2 = `c1_rise` + `c2_pulse` | `surface_lift` | **7.3 %** | **5.3 %** |
| | `surface_zero` | 12.8 % | 5.7 % |
| | `global` | 9.6 % | 8.0 % |
| | `global_src` | 24.2 % | 12.1 % |

- The advantage is not the exact air columns alone: `global_src` is worse
  than `global`.
- At small rank, the nu_ref continuation into the steel matters:
  `surface_zero` is worse than `surface_lift`.

### 2. Smallest basis at the 10 % target

Source: `results/team10_basis_10pct_p3.json`.  Training is the
automatically selected `c1_rise` + `c2_pulse_x2`.

| (k, r) | columns | `surface_lift` worst / elements / online | `global` worst / elements / online |
|---|---|---|---|
| (4, 2) | 12 | **8.7 %** / 20 / 0.3 s | 15.8 % / 27 / 0.4 s |
| (6, 3) | 18 | 7.1 % / 26 / 0.9 s | 14.9 % / 36 / 1.3 s |
| (8, 4) | 24 | 7.0 % / 29 / 1.1 s | 10.2 % / 43 / 1.7 s |
| (10, 5) | 30 | 6.8 % / 37 / 2.1 s | **9.6 %** / 58 / 5.1 s |
| (16, 8) | 48 | 6.2 % / 46 / 4.5 s | 8.7 % / 96 / 7.9 s |

ECSW tolerance is 1e-2 in this table.  At 3e-2, `surface_lift` (4, 2) still
passes: 8.8 %, 15 elements, 0.2 s.

### 3. Automatic selection of the training pattern

Source: `results/team10_pattern_select_p3.json`.  The selection starts from
`c1_rise` and adds one pool member per round.

- **residual:** RB-greedy.  Run the current reduced model on every candidate
  and pick the largest dual-norm FOM residual eta.  No candidate FOM is read
  before its selection.
- **gp:** BayPOD-like.  Pick the largest Gaussian-process posterior variance
  over the pattern features (a1, a2, shape).

Test errors below are for the `surface_lift` basis at k16 r8.

| round | residual picks | c2_rise / mix_add / mix_oppose | gp picks | c2_rise / mix_add / mix_oppose |
|---|---|---|---|---|
| 0 | - | 61.2 / 64.7 / 27.2 % | - | same |
| 1 | `c2_pulse_x2` | 2.6 / 6.2 / 3.4 % | `c2_pulse_x2` | same |
| 2 | `mixp_oppose` | 2.5 / 6.3 / 2.2 % | **`c1_rise_neg`** (trap) | 2.6 / 6.3 / 3.5 % |
| 3 | `c2_pulse` | 1.8 / 5.8 / 2.0 % | `mixp_oppose` | 2.4 / 6.5 / 2.4 % |

- **Round 0:** eta is 0.39-0.92 for the web-coil candidates and at most
  0.015 for the main-coil ones.  The trap has the lowest eta in every round.
- **GP rule:** it spends round 2 on the trap, because it sees only feature
  distances.

### 4. Cheaper indicator

Source: `results/team10_indicator_check_p3.json`, rebuilt from the run log;
see the provenance section.

- **Why only the steel rows are needed:** for the `surface_lift` basis the
  air rows of the FOM residual vanish by construction.  The measured
  air/total ratio at step 30 is 0.0 for every candidate.
- **The steel indicator:**
  - assembles only the steel Variations plus K_rest u on the gamma rows;
  - measures the residual with the 57k-DOF steel lifting operator;
  - is sampled every 5th step.
- **Cost:** about 4 s per candidate, against about 110 s for the full
  dual-norm indicator; one selection round falls from about 40 min to about
  110 s.
- **Ranking:**
  - both indicators separate the four web-coil candidates from the three
    main-coil candidates by a factor of 20 or more, and both rank the trap
    last;
  - they order `c2_pulse` and `c2_pulse_x2` differently (2.53 vs 2.47, steel;
    0.77 vs 0.92, dual).
- **Scale:** the steel indicator is about 3x the dual one.  Its 10 % stop
  threshold (provisionally eta < 0.1) is calibrated on these same data only.

### 5. Small selection basis

Source: `results/team10_pattern_select_p3_sel42.json`, hibino.

- **Selection basis:** (4, 2), 12 columns, ECSW 1e-2.
- **Evaluation basis:** unchanged at k16 r8.

| | selection basis k16 r8 | selection basis k4 r2 |
|---|---|---|
| one round (build / candidate ROMs / indicator) | 38-46 / 30-38 / 27-30 s, about 100-110 s | 13-18 / 4 / 27-32 s, **about 49 s** |
| round 0 ranking (steel indicator) | `c2_pulse` 2.53 > `c2_pulse_x2` 2.47 > `mixp_add` 1.23 > `mixp_oppose` 1.14 >> main-coil group <= 0.07 | `c2_pulse` 1.01 > `c2_pulse_x2` 1.00 > `mixp_add` 0.51 > `mixp_oppose` 0.49 >> main-coil group <= 0.021 |
| picks, rounds 1 and 2 | `c2_pulse`, `c2_pulse_x2` | the same |
| tests after round 1 (k16 r8) | 2.2 / 5.3 / 3.8 % | the same |

- **Same selection:** the small basis makes the same selections, so the test
  errors are identical.  Every round passes the 10 % target.
- **Where the time goes now:** the indicator (about 4.5 s per candidate).
  The seven candidate reduced runs together take about 4 s.
- **Trap candidate:** with the k4 r2 basis the mirror `c1_rise_neg` is no
  longer last.  Its eta (0.017) equals the truncation floor of the training
  pattern itself; by odd symmetry the mirror is represented exactly as well
  as `c1_rise`.  That eta is still about 60 times below the web-coil group,
  so it is never picked.
- **Stop rule:** eta scales with the selection basis.  After round 1 the pool
  maximum is 0.24 at k4 r2 against 0.12 at k16 r8.  A stop rule should
  therefore compare the pool against this floor rather than against a fixed
  number.  The following section evaluates that rule on the recorded calibration set.

### 6. Stop-rule calibration

Source: `results/team10_stop_calibration_p3.json`, `run_stop_calibration.py`,
hibino, 4766 s.

- **Training sets:** `c1_rise` plus every subset of 0-2 pool members, 29 sets.
- **Bases:** each set is scored with both the (4, 2) and the (8, 4) basis.
  The same basis serves selection and production.
- **Stop quantity:** the largest eta over the pool members that are not
  trained.
- **Truth:** the worst of the three test errors.

| basis | worst < 10 % | 10-20 % | >= 20 % | corr(log eta, log err), all / below 20 % |
|---|---|---|---|---|
| (4, 2) | 14 | 8 | 7 | 0.67 / 0.51 |
| (8, 4) | 18 | 4 | 7 | 0.82 / 0.22 |

Rule A stops when the largest pool eta is below tau.  "Stops" counts the
sets on which the rule stops; the other columns describe those sets.

| tau | (4, 2): stops / at >= 10 % / worst stopped | (8, 4): stops / at >= 10 % / worst stopped |
|---|---|---|
| 0.10 | 0 / 0 / - | 5 / 0 / 7.3 % |
| 0.20 | 2 / 1 / 10.1 % | 14 / 2 / 10.5 % |
| 0.25 | 13 / 3 / 16.9 % | 14 / 2 / 10.5 % |
| 0.50 | 16 / 3 / 16.9 % | 18 / 3 / 11.6 % |
| 1.00 | 19 / 5 / 16.9 % | 22 / 4 / 13.5 % |

- **Gross failures in this calibration set:** eta separates them.  Every set with a test
  error of 20 % or more has pool eta >= 1.0 (the web coil is untrained), so no
  tau <= 1.0 ever stops on one.
- **The 10 % line:** eta does not resolve it.  Below 20 % the correlation is
  weak, and near 10 % the error depends on which patterns share the small
  basis.  At (4, 2), every set that holds `c1_rise_half` ends at 13-17 %: a
  low-amplitude pattern spends part of a 12-column basis.
- **Empirical thresholds at 10 %:** no observed false stop with tau = 0.195 at (4, 2), which
  stops 1 of 14 passing sets, and tau = 0.11 at (8, 4), which stops 5 of 18.
- **Rule B:** the ratio to the training floor is not better.  Its values with no observed false stop
  (1.79 and 0.93) stop 4 and 2 passing sets.
- **Resulting rule:**
  - eta < 0.5 is a candidate selection heuristic: no missing-coil failures
    were observed in this finite calibration set, and the worst stopped set
    was 17 %. It is not a coverage or error bound for unseen excitations;
  - the (8, 4) basis with tau = 0.1 had no observed 10 % failures here.
    A required 10 % accuracy needs independent held-out full-model checks
    or a justified error bound; this threshold alone does not establish it.
  - These thresholds come from one geometry and pool.  The premise that the
    tests resemble the pool is untested beyond it.

## Limits

- **Where the advantage holds:** only with a changing excitation pattern.
  With one coil and a changing saturation state (amplitude x0.3 or x2,
  switch-off), plain POD was as good or better in 7 of 8 cases.
- **Size of the advantage:** it shrinks with rank, to 2-3 points at k16 r8.
- **Coverage:** one geometry and one mesh (p3, 6 mm).
- **Unseen patterns:** training without the web coil fails for every basis.
  The residual selection is the remedy.
- **Agreement with measurement:** this lane does not validate the FOM against
  the TEAM 10 measurements.  The converged FOM rises faster than the
  digitised curves (S1 at 20 ms: 0.80 T vs 0.52 T), and that mismatch is
  open.

## Running

The SparseCholesky migration smoke is reproducible with
`python check_sparsecholesky.py --output results/sparsecholesky_smoke.json`.
The recorded coarse p1 case has 22,415 DOFs. Its five constrained relative
true residuals are at most 2.18e-11 (limit 1e-7), and the three-step nonlinear
trajectory finishes with finite fields. See
`results/sparsecholesky_smoke_20260929.json`. This does not establish p3
memory capacity or reproduce the historical ROM accuracy/timing campaign.

Heavy.  Run on hibino (one job at a time) or on an idle mdx host.  The
selection run committed about 27 GB of private memory on mdx2.  That host has
no pagefile (commit limit 57.4 GB), and a concurrent MATLAB session holding
about 10 GB was enough to make it fail.

```powershell
$env:TEAM_CASE = 'team10'; $env:TEAM_QUARTER = '1'; $env:TEAM_COIL2 = '1'; $env:NGS_THREADS = '24'
# cheap training pool (about 270 s each)
foreach ($w in 'c1_rise','c1_pulse','c2_pulse','mixp_add','mixp_oppose','c2_pulse_x2','c1_rise_half') { python fom_aphi.py 0.006 3 $w 30 1 }
# exact tests (about 2200 s each)
foreach ($w in 'c2_rise','mix_add','mix_oppose') { python fom_aphi.py 0.006 3 $w 60 }
python run_pattern_study.py 0.006 3
python run_basis_10pct.py 0.006 3
python run_pattern_select.py 0.006 3 3      # SELECT_ONLY_RESIDUAL, SELECT_CHECK, SELECT_BASIS=k,r,tol, SELECT_ETA_STOP, SELECT_TAG
```

- **Snapshot files:** `fom_aphi.py` writes `results/<tag>/<wave>.npy`
  (about 160 MB each, git-ignored) and `<wave>.json`.  The JSONs are kept
  here.  The `.npy` files are regenerated by the commands above.
- **Coarse smoke run:** use `0.01 1` for `maxh` and `order`.

## Provenance

- **Recorded source hashes:** every result JSON records the SHA-256 of the
  sources that produced it.
- **`team10_pattern_select_p3.json`:** produced by the earlier dual-norm
  version of `run_pattern_select.py`.  The current file adds the steel
  indicator, the small selection basis and the stop rule.
- **`team10_indicator_check_p3.json`:** parsed from the recovered `run.log`.
  That run reached the mdx2 commit limit after printing rounds 0-2 and before
  writing its own JSON.  The JSON records the log's SHA-256.
