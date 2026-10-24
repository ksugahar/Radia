"""Audit saved F1 gate arithmetic, without assembling or solving any model.

This derived audit is not the executed solver runner. The exact executed runner
and original artifact hashes are preserved in result.json provenance.
"""
import argparse
import json
import math
from pathlib import Path
import sys

def atomic_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf8")

def check():
    cases = {(c["m"],c["p"]):c for c in RESULT["cases"]}
    tests = []
    def gate(name, actual, bound, ok):
        tests.append({"name":name,"actual":actual,"bound":bound,"status":"PASS" if ok else "FAIL"})
    def roundoff_plateau(ei, ej):
        # Diagnostic only, never a PASS waiver or a relaxation of an accuracy gate.
        # Relative-error quantities within 64 eps, or identical within 16 ulps,
        # cannot support a credible observed rate and are UNRESOLVED.
        return min(ei, ej) <= 0 or max(ei, ej) <= 64*sys.float_info.epsilon or abs(ei-ej) <= 16*max(math.ulp(ei),math.ulp(ej))
    for c in cases.values():
        for key in ("true_residual","direct_true_residual"):
            gate(f"{c['m']}/{c['p']}/{key}",c[key],1e-6,c[key]<=1e-6)
        for key,e in c["parity"].items():gate(f"parity/{c['m']}/{c['p']}/{key}",e,1e-6,e<=1e-6)
    keys = list(next(iter(cases.values()))["errors"])
    for key in keys:
        for p in (1,2,3):
            e=cases[4,p]["errors"][key]
            if p == 1:
                tests.append({"name":f"finest/{p}/{key}","actual":e,"status":"RECORDED_ONLY",
                              "scope":"APPROVAL-20261024-IH-AMS-F1-FROZEN-03: finest <=1% applies to p2/3 only"})
            else:
                gate(f"finest/{p}/{key}",e,.01,e<=.01)
            for m,n in ((2,4),):
                a,b=cases[m,p],cases[n,p]
                ei,ej=a["errors"][key],b["errors"][key]
                ratioh=a["hmax_over_a"]/b["hmax_over_a"]
                if roundoff_plateau(ei,ej) or ratioh<=1:
                    tests.append({"name":f"h/{m}/{n}/{p}/{key}","actual_errors":[ei,ej],"status":"UNRESOLVED","reason":"unchanged h, zero error or roundoff plateau"});continue
                rate=math.log(ei/ej)/math.log(ratioh)
                gate(f"h/{m}/{n}/{p}/{key}",rate,p-.5,ej<ei and rate>=p-.5)
        for m in (2,4):
            for p in (1,2):
                ei,ej=cases[m,p]["errors"][key],cases[m,p+1]["errors"][key]
                if roundoff_plateau(ei,ej):
                    tests.append({"name":f"p_monotone/{m}/{p}/{key}","actual_errors":[ei,ej],"status":"UNRESOLVED","reason":"zero error or roundoff plateau"})
                    if m==4: tests.append({"name":f"p_ratio/{p}/{key}","status":"UNRESOLVED" if key in ("H","B","J") else "RECORDED_ONLY","actual_errors":[ei,ej],"reason":"no credible rate"})
                    continue
                gate(f"p_monotone/{m}/{p}/{key}",[ei,ej],"strict decreasing",0<ej<ei)
                if m==4:
                    ratio=ej/ei if ei>0 else None
                    if key in ("H","B","J"):
                        gate(f"p_ratio/{p}/{key}",ratio,2*math.sqrt(2)/16,ratio is not None and ratio<=2*math.sqrt(2)/16)
                    else:tests.append({"name":f"p_ratio/{p}/{key}","actual":ratio,"status":"RECORDED_ONLY"})
    atomic_json(ROOT/"gates.json", {"checks":tests,"decision":"PASS" if all(x["status"] in ("PASS","RECORDED_ONLY") for x in tests) else "FAIL_OR_UNRESOLVED"})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    ROOT = args.report.parent
    RESULT = json.loads(Path(__file__).with_name("result.json").read_text(encoding="utf8"))
    check()
    generated = ROOT / "gates.json"
    if generated != args.report:
        generated.replace(args.report)
    report = json.loads(args.report.read_text(encoding="utf8"))
    print(json.dumps({"decision": report["decision"], "checks": len(report["checks"])}))
    raise SystemExit(0 if report["decision"] == "PASS" else 1)
