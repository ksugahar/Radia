# Complex cases for the BDDC wirebasket AMS

`../bddc_coarse_cases.py` solves three time-harmonic eddy-current problems on
`HCurl(order=2, nograds=True, complex=True)` with NGSolve BDDC + COCR (tol
1e-8), once with the direct wirebasket inverse (sparse Cholesky) and once
with `coarsetype="sparsesolv_ams"` (4 cycles, edge wirebasket, lean coarse
levels).  Each uses the reduced field A = A_s + A_r, A_r = 0 on the outer
boundary, and an `eps*nu` mass with eps = 1e-6.  The Hiruma coil is in
`../hiruma/`.

- `sphere`: Cu sphere, a = 5 mm, sigma = 5.8e7, in a uniform 1 T AC field
  inside a 60 mm air box.  Reference: the analytic (Smythe) induced dipole,
  loss = (w/2) |Im m_z| B0.
- `disk`: disk r = 10 mm, t = 0.5 mm, mu_r = 100, sigma = 1e7, in a uniform
  axial H0 inside an 80 x 80 x 40 mm air box.  Reference: stored
  axisymmetric volume-averaged Bz / (mu0 H0) from
  `../../vim_coupled/results_magnetic_conductor_disk_adjudication.json`.
- `plate`: Cu plate 17.72 x 17.72 x 2 mm, conductor-only mesh, A_r = 0 on its
  surface (the A1 sweep setting).  Reference: the direct run.

mdx1 (disk) and mdx2 (sphere, plate), 8 threads, commit 41a044187, one binary
(`summary_20260928.json`); each host also carried one job of another session.

| case | f | dofs | coarse | iterations | setup s | solve s | peak GB | vs reference |
|---|---|---|---|---|---|---|---|---|
| sphere | 700 Hz (a/delta 2.0) | 349,648 | direct | 46 | 76.6 | 15.7 | 10.7 | loss 0.14 % |
| sphere | 700 Hz | 349,648 | AMS | 44 | 1.0 | 1.8 | 0.71 | loss 0.14 % |
| sphere | 10 kHz (a/delta 7.6) | 349,648 | direct | 50 | 77.7 | 16.6 | 10.2 | loss 4.0 % |
| sphere | 10 kHz | 349,648 | AMS | 41 | 1.0 | 1.6 | 0.72 | loss 4.0 % |
| disk | 100 Hz | 213,733 | direct | 43 | 33.1 | 8.2 | 5.7 | Bz 0.38 % |
| disk | 100 Hz | 213,733 | AMS | 55 | 0.8 | 1.6 | 0.49 | Bz 0.38 % |
| disk | 10 kHz | 213,733 | direct | 71 | 34.8 | 13.5 | 6.0 | Bz 1.97 % |
| disk | 10 kHz | 213,733 | AMS | 61 | 0.7 | 1.6 | 0.49 | Bz 1.97 % |
| disk | 10 kHz, finer | 702,468 | direct | 65 | 243.5 | 57.1 | 24.9 | Bz 1.37 % |
| disk | 10 kHz, finer | 702,468 | AMS | 85 | 1.9 | 8.1 | 1.30 | Bz 1.37 % |
| plate | 1 kHz | 742,783 | direct | 39 | 60.4 | 17.7 | 10.1 | - |
| plate | 1 kHz | 742,783 | AMS | 34 | 2.1 | 5.0 | 1.30 | loss = direct to 2e-13 |
| plate | 100 kHz | 742,783 | direct | 50 | 62.7 | 22.7 | 10.0 | - |
| plate | 100 kHz | 742,783 | AMS | 31 | 2.1 | 4.9 | 1.31 | loss = direct to 1e-11 |

The two coarse solvers give the same field in every case (loss to 2e-9 or
better, the solver tolerance being 1e-8), so the reference errors are the
discretization's: the 10 kHz sphere has 0.5 mm conductor elements for a
0.66 mm skin depth, and the disk error falls from 1.97 % to 1.37 % with the
finer mesh.  Unlike Hiruma (where the direct coarse solve needs about a third
of the AMS iterations) the iteration counts are close here, and AMS is 11-36x
faster (setup + solve) with 8-19x less memory; the direct cost is again the
wirebasket factorization.  The finer direct disk run started with the host
at 89 % CPU, so its time may be inflated.

Not run: the A-V round wire (a compound HCurl x H1 space) and the Kelvin
sphere (a periodic space).  The wirebasket AMS requires an HCurl space whose
lowest-order dofs are numbered by edge and fails loudly when they are not.

## Order 3 (`summary_20260928_p3.json`)

At order 3 the wirebasket is still the lowest-order edge block (the new
cell dofs are local), so the same `sparsesolv_ams` applies.  Same host split
and binary; `--order 3`:

| case | dofs | coarse | iterations | setup s | solve s | peak GB | vs reference |
|---|---|---|---|---|---|---|---|
| sphere 700 Hz | 293,071 | direct | 103 | 70.9 | 28.7 | 9.3 | loss 0.16 % |
| sphere 700 Hz | 293,071 | direct, edge wirebasket | 85 | 71.7 | 23.4 | 9.3 | loss 0.16 % |
| sphere 700 Hz | 293,071 | AMS | 62 | 0.9 | 2.9 | 0.78 | loss 0.16 % |
| sphere 700 Hz | 949,723 | direct (both wirebaskets) | fails: bad allocation during Assemble | | | | |
| sphere 700 Hz | 949,723 | AMS | 67 | 2.6 | 10.8 | 2.2 | loss 0.14 % |
| disk 10 kHz | 582,178 | direct | 174 | 213.3 | 117.0 | 19.8 | Bz 0.88 % |
| disk 10 kHz | 582,178 | direct, edge wirebasket | 171 | 210.0 | 108.8 | 21.4 | Bz 0.88 % |
| disk 10 kHz | 582,178 | AMS | 101 | 1.7 | 9.6 | 1.4 | Bz 0.88 % |
| plate 1 kHz | 318,700 | direct | 93 | 9.8 | 10.6 | 2.9 | - |
| plate 1 kHz | 318,700 | direct, edge wirebasket | 77 | 9.6 | 8.8 | 2.8 | - |
| plate 1 kHz | 318,700 | AMS | 52 | 1.1 | 2.7 | 0.77 | loss = direct to 5e-13 |
| plate 1 kHz | 1,435,724 | direct (both wirebaskets) | fails: bad allocation during Assemble | | | | |
| plate 1 kHz | 1,435,724 | AMS | 53 | 4.6 | 12.5 | 3.2 | - |

The coarse solvers again give the same field.  At order 3 AMS also needs
fewer outer iterations than the direct wirebasket inverse (62 vs 85-103, 101
vs 171-174, 52 vs 77-93), and it is 5-29x faster where the direct solver
runs at all; the two largest meshes (0.95M and 1.44M dofs) exceed the direct
solver's memory on mdx2 (about 45 GB free) and solve with AMS in 13-17 s.
The disk error against the axisymmetric reference drops from 1.37-1.97 %
(order 2) to 0.88 %.