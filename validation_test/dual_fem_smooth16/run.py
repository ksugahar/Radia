"""Job-limited runner for the smooth periodic A/Omega study.

Usage (venv with NGSolve 6.2.2608, Windows):
  python -I run.py --out <results dir> [--cases F1-2D,F1-3D] [--smoke] [--logical-cpu N]

Limits: one logical CPU (group 0, ``--logical-cpu``), BelowNormal priority,
a job object with KILL_ON_JOB_CLOSE and a 16 GiB job-memory cap, and one
thread for BLAS/OpenMP.  The limits are read back and recorded in owner.json.
"""

import argparse, ctypes, ctypes.wintypes as w, datetime, hashlib, json, os, platform, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOGICAL = 0
MASK = 1 << LOGICAL
JOB_MEMORY = 16 << 30
_HANDLES = []


def apply_limits():
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[name] = "1"
    k = ctypes.WinDLL("kernel32", use_last_error=True)
    k.GetCurrentProcess.restype = w.HANDLE
    k.GetCurrentThread.restype = w.HANDLE

    class GroupAffinity(ctypes.Structure):
        _fields_ = [("Mask", ctypes.c_size_t), ("Group", w.WORD), ("Reserved", w.WORD * 3)]

    class Basic(ctypes.Structure):
        _fields_ = [("times", ctypes.c_int64 * 2), ("flags", w.DWORD), ("minimum", ctypes.c_size_t),
                    ("maximum", ctypes.c_size_t), ("active", w.DWORD), ("affinity", ctypes.c_size_t),
                    ("priority", w.DWORD), ("scheduling", w.DWORD)]

    class Extended(ctypes.Structure):
        _fields_ = [("basic", Basic), ("io", ctypes.c_uint64 * 6), ("process_memory", ctypes.c_size_t),
                    ("job_memory", ctypes.c_size_t), ("peak_process", ctypes.c_size_t), ("peak_job", ctypes.c_size_t)]

    k.SetThreadGroupAffinity.argtypes = [w.HANDLE, ctypes.POINTER(GroupAffinity), ctypes.POINTER(GroupAffinity)]
    k.GetThreadGroupAffinity.argtypes = [w.HANDLE, ctypes.POINTER(GroupAffinity)]
    k.SetPriorityClass.argtypes = [w.HANDLE, w.DWORD]
    k.GetPriorityClass.argtypes = [w.HANDLE]; k.GetPriorityClass.restype = w.DWORD
    k.CreateJobObjectW.restype = w.HANDLE
    k.SetInformationJobObject.argtypes = [w.HANDLE, ctypes.c_int, ctypes.c_void_p, w.DWORD]
    k.QueryInformationJobObject.argtypes = [w.HANDLE, ctypes.c_int, ctypes.c_void_p, w.DWORD, ctypes.c_void_p]
    k.AssignProcessToJobObject.argtypes = [w.HANDLE, w.HANDLE]
    proc = k.GetCurrentProcess()
    target = GroupAffinity(); target.Group = 0; target.Mask = MASK
    if not k.SetThreadGroupAffinity(k.GetCurrentThread(), ctypes.byref(target), ctypes.byref(GroupAffinity())):
        raise ctypes.WinError(ctypes.get_last_error())
    if not k.SetPriorityClass(proc, 0x4000):
        raise ctypes.WinError(ctypes.get_last_error())
    job = k.CreateJobObjectW(None, None)
    info = Extended(); info.basic.flags = 0x2000 | 0x200; info.job_memory = JOB_MEMORY
    if not job or not k.SetInformationJobObject(job, 9, ctypes.byref(info), ctypes.sizeof(info)):
        raise ctypes.WinError(ctypes.get_last_error())
    if not k.SetInformationJobObject(job, 14, ctypes.byref(target), ctypes.sizeof(target)):
        raise ctypes.WinError(ctypes.get_last_error())
    if not k.AssignProcessToJobObject(job, proc):
        raise ctypes.WinError(ctypes.get_last_error())
    _HANDLES.append(job)
    back = GroupAffinity()
    k.GetThreadGroupAffinity(k.GetCurrentThread(), ctypes.byref(back))
    rb = Extended()
    k.QueryInformationJobObject(job, 9, ctypes.byref(rb), ctypes.sizeof(rb), None)
    if back.Group != 0 or back.Mask != MASK or k.GetPriorityClass(proc) != 0x4000 \
            or not rb.basic.flags & 0x200 or rb.job_memory != JOB_MEMORY:
        raise RuntimeError("resource readback mismatch")
    return dict(group=0, logical=[LOGICAL], mask=hex(MASK), threads=1, priority="BelowNormal",
                job_memory_bytes=JOB_MEMORY, kill_on_close=True, readback=True)


def _free_memory():
    class M(ctypes.Structure):
        _fields_ = [("len", w.DWORD), ("load", w.DWORD), ("total", ctypes.c_uint64), ("avail", ctypes.c_uint64),
                    ("tp", ctypes.c_uint64), ("ap", ctypes.c_uint64), ("tv", ctypes.c_uint64), ("av", ctypes.c_uint64),
                    ("ae", ctypes.c_uint64)]
    m = M(); m.len = ctypes.sizeof(M)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    return int(m.avail)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--cases", default="")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--logical-cpu", type=int, default=0,
                    help="logical CPU in processor group 0 that the job is pinned to")
    args = ap.parse_args()
    global LOGICAL, MASK
    if not 0 <= args.logical_cpu < 64:
        raise ValueError("--logical-cpu must be in [0, 64)")
    LOGICAL, MASK = args.logical_cpu, 1 << args.logical_cpu
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    resources = apply_limits()
    import ngsolve as ng
    sys.path.insert(0, str(HERE))
    import cases as C
    import smooth16 as S
    sources = {name: sha(HERE / name) for name in ("cases.py", "smooth16.py", "run.py")}
    owner = dict(study="dual FEM smooth periodic A/Omega",
                 start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 resources=resources, sources=sources, ngsolve=ng.__version__, python=sys.version,
                 platform=platform.platform(), smoke=args.smoke,
                 free_memory_bytes=_free_memory())
    (out / ("owner-smoke.json" if args.smoke else "owner.json")).write_text(json.dumps(owner, indent=2))
    wanted = set(filter(None, args.cases.split(",")))
    for case in C.cases(ng):
        cid, dim = case[0], case[1]
        if wanted and cid not in wanted:
            continue
        started = time.perf_counter()
        pl, hl = C.P_LADDER[dim], C.H_LADDER[dim]
        p_orders = pl["orders"][:3] if args.smoke else pl["orders"]
        h_ns = hl["n"][:2] if args.smoke else hl["n"]
        mesh0 = S.torus(ng, dim, pl["n"])
        check = S.verify_truth(ng, mesh0, dim, S.truth(ng, dim, case[3]))
        p_rows = [S.run_level(ng, case, pl["n"], p, C.ASSEMBLY_BONUS, C.READOUT_EXTRA) for p in p_orders]
        h_rows = [S.run_level(ng, case, n, hl["order"], C.ASSEMBLY_BONUS, C.READOUT_EXTRA) for n in h_ns]
        summary = S.summarize(cid, p_rows, h_rows)
        ref = C.REFERENCE[dim] if not args.smoke else {"n": pl["n"], "p": p_orders[-1]}
        try:
            blend = S.reference_blend_test(p_rows, h_rows, C.COARSE[dim], ref)
        except KeyError as e:
            blend = dict(error=f"missing level {e}")
        record = dict(case=cid, dim=dim, family=case[2], truth_check=check, p_ladder=p_rows, h_ladder=h_rows,
                      rates=summary, blend_test=blend, seconds=time.perf_counter() - started, sources=sources)
        name = f"{'smoke-' if args.smoke else ''}{cid}.json"
        (out / name).write_text(json.dumps(record, indent=2, default=float))
        print(cid, "done", round(record["seconds"], 1), "s", blend.get("closer"), summary["verdict"], flush=True)


if __name__ == "__main__":
    main()
