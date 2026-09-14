# Analytical branch-source preservation

Recovered without modification from the lab source directory:
`W:\999_菅原賢悟\19_磁気ヒステリシス\2024_IGTE_共同研究\2024_03_08_H-input_B-input\Potter_Schmulian`.

`CASE_02.m` constructs analytical branches, optimizes endpoints with
`fminsearch`, and resamples with `makima` before saving `B_input.mat`.
These are model-generated values, not a measured specimen. The branch spacing
is 0.095 T, with 20 levels up to 1.9 T. Original source names and bytes are kept.
The label `Ms=1.91` in the generator is used as a flux-density-scale quantity;
it must not be reinterpreted as SI magnetization in A/m.

`../../verify_fixture_lineage.py` reidentifies the branches with the current
Python implementation and compares all 128 fields with the committed NPZ.
`../../fixture_lineage.json` stores the successful comparison and source hashes.
It does not rerun the MATLAB generator or the coupled rod/cube solve.
