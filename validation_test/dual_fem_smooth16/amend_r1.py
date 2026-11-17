"""Amendment R1 (2026-11-08): finer 3D reference levels.

Reason: none of the eight 3D cases met the predeclared reference-convergence
test at REFERENCE = (n=3, p=4).  The test is the last-step drift of Q_A and of
Q_Omega, each at most 0.1 x |coarse midpoint - coarse 2:1|.  This amendment
changes only the reference side.  The coarse level (n=3, p=2), the tested
weights and the decision rule are unchanged.  The coarse results were not used
to choose these levels.

New 3D reference: n=6, p=4.  The drift step is n=6, p=3.  The 2D results stay
as they are.
Usage: python -I amend_r1.py --out <results dir>
"""

import argparse, json, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REFERENCE_R1 = {"n": 6, "p": 4}
LEVELS_R1 = [(6, 3), (6, 4)]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--cases", default="")
    ap.add_argument("--logical-cpu", type=int, default=0)
    args = ap.parse_args()
    sys.path.insert(0, str(HERE))
    import run as R
    if not 0 <= args.logical_cpu < 64:
        raise ValueError("--logical-cpu must be in [0, 64)")
    R.LOGICAL, R.MASK = args.logical_cpu, 1 << args.logical_cpu
    resources = R.apply_limits()
    import ngsolve as ng
    import cases as C, smooth16 as S
    out = Path(args.out)
    sources = {n: R.sha(HERE / n) for n in ("cases.py", "smooth16.py", "run.py", "amend_r1.py")}
    (out / "owner-r1.json").write_text(json.dumps(dict(resources=resources, sources=sources,
                                                       ngsolve=ng.__version__, reference=REFERENCE_R1), indent=2))
    wanted = set(filter(None, args.cases.split(",")))
    for case in C.cases(ng):
        cid, dim = case[0], case[1]
        if dim != 3 or (wanted and cid not in wanted):
            continue
        base = json.loads((out / f"{cid}.json").read_text())
        started = time.perf_counter()
        rows = [S.run_level(ng, case, n, p, C.ASSEMBLY_BONUS, C.READOUT_EXTRA) for n, p in LEVELS_R1]
        blend = S.reference_blend_test(base["p_ladder"], base["h_ladder"] + rows, C.COARSE[3], REFERENCE_R1)
        rec = dict(case=cid, amendment="R1", reference_levels=rows, blend_test=blend,
                   seconds=time.perf_counter() - started, sources=sources)
        (out / f"{cid}-r1.json").write_text(json.dumps(rec, indent=2, default=float))
        print(cid, "R1 done", round(rec["seconds"], 1), "s", blend["closer"], flush=True)


if __name__ == "__main__":
    main()
