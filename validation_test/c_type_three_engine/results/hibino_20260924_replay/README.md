# HDiv fine-mesh replay diagnostic

Not a three-engine validation pass. The unchanged driver's HDiv adapter is invoked with the same checked fine meshes and coil setup; a task-owned exception exits after collecting its result, before other engines. The driver resets native objects for each run. No comparison result is fabricated for omitted engines.

Two tol=1e-8 runs reproduce the previous fine field exactly (relative vector difference 0). Tightening to 1e-10 changes the full observation vector by 2.2859223789941344e-10 relative to the previous fine result; iterations rise from 40 to 46. These effects cannot explain the approximately 1% medium-to-fine change. Root cause is not established; geometry/discretization and native field evaluation remain candidates. No algebraic residual is inferred from field stability.

Runtime verified after execution: Python 3.12.10 at C:/temp/hdiv-replay-20260924/venv/Scripts/python.exe; Radia 5.0.0 from this venv's Lib/site-packages/radia; NGSolve 6.2.2606 from C:/Program Files/Python312/Lib/site-packages/ngsolve. Same CI wheel as the baseline, archived with scripts, inputs, baseline JSON and log. Foreground SSH python -u replay_hdiv.py; OMP/MKL/OPENBLAS each 8; driver threads 8. No other compute process was present at preflight.

Recovery archive S:/Radia/validation_artifacts/hdiv_replay_20260924/recovery.zip has matching remote/LAB SHA256 3ebbc0652763e04e4296e47e47c87b4f6ec3ef000f22abb635550bce6fab0865. Remote task-owned scratch is removed after this evidence commit; shared mesh family remains.
