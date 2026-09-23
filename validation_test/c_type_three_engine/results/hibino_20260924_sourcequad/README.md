# Source-load quadrature sensitivity (HDiv only)

This is an experimental diagnostic, NOT a three-engine validation pass or
matched-accuracy performance comparison. Checked fine mesh, mu_r=1000, BDM1,
eight threads, Gram eps=1e-14, CG tolerance=1e-10 and direct field evaluation.
Only the source LinearForm integral gets bonus_intorder=0,4,8; Gram and mass
assembly are untouched. The script patches one exact expression in a disposable
venv and restores it in finally. Original/experimental source hashes are in JSON.
Source snapshots, wheel, inputs and logs are in the recovered archive.

Bonus 8 changes the full observation vector by 2.3637805988164055e-5 relative
to the previous fine baseline (approximately 0.00236%). This is much smaller
than the roughly 1% medium-to-fine HDiv change. It does not prove continuum
accuracy or rule out every integration/discretization issue. All three native
solves reported convergence. No production source or manuscript was changed.

Recovery: S:/Radia/validation_artifacts/hdiv_sourcequad_20260924/recovery.zip
SHA256: d0a28016dbc0d5a96621fbcdb2e2cf332dd35aaa5ce43194dc0062c2c0e459e8
Identical hash verified on hibino and LAB before committing these results.
Remote task directory and archive are to be removed after this commit;
the shared C:/temp/radia-ctype-family is not owned by this diagnostic.
