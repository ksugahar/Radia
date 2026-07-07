# HDiv-VIM Cube Benchmark

This directory contains the HDiv-VIM counterpart of the soft-iron cube
benchmark used in the 2026-08-25 static/rotating-machine MMMM manuscript.

The driver is intentionally separate from the older hexahedron/tetrahedron
benchmark scripts because those results came from the retired six-face MMMM
path.  Current Radia soft-iron demag validation should use the RT1 HDiv-VIM
entry:

```powershell
python validation_test\cube_uniform_field\hdiv\bench_hdiv_cube.py `
  --sizes 2 4 6 `
  --material both `
  --bh-table "W:\02_学会資料\2026年度\2026_08_25_静止器・回転機@八戸\MMMM@菅原・伊田・矢野\原稿\figures\BH.txt" `
  --output "W:\02_学会資料\2026年度\2027_01_静止器・回転機@\HDiv要素MMM@菅原\results\hdiv_cube_initial_lab_20260707.json"
```

Problem definition:

- 1 m cube centered at the origin.
- Structured pure-hex `N x N x N` mesh.
- Applied field `H0 = 200 kA/m` in `+z`.
- Linear `mu_r=1000` or nonlinear BH table.
- Solver: `radia.vim.Solve`, RT1 HDiv-VIM, analytic HACApK charge Gram.

Recorded fields include `n_el`, `ndof`, `n_charge`, `iters`, `M_avg`,
`demag`, wall time, memory snapshot, and the solver's `hmat_stats`.

The LAB initial run on 2026-07-07 reached `N=6` only.  This is a smoke /
setup result, not a publication timing claim.  Larger scaling points should be
run on mdx and labelled as mdx validation before being used in the 2027-01
presentation.
