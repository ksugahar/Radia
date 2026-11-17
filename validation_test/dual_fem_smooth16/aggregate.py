"""Aggregate the 16 per-case JSON files into summary.json and a Markdown table.

Usage: python -I aggregate.py <results dir>
"""

import hashlib, json, sys
from pathlib import Path

ORDER = [f"F{i}-{d}D" for i in range(1, 9) for d in (2, 3)]


def main(root):
    root = Path(root)
    rows, files = [], {}
    for cid in ORDER:
        path = root / f"{cid}.json"
        if not path.is_file():
            rows.append(dict(case=cid, status="missing"))
            continue
        files[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        d = json.loads(path.read_text())
        r = d["rates"]
        r1 = root / f"{cid}-r1.json"
        if r1.is_file():  # amendment R1: finer 3D reference, coarse level unchanged
            files[r1.name] = hashlib.sha256(r1.read_bytes()).hexdigest()
            b = json.loads(r1.read_text())["blend_test"]
            reference_source = "R1 (n=6, p=4)"
        else:
            b = d["blend_test"]
            reference_source = "predeclared"
        ref = b["reference"]
        e = b["coarse_errors"]
        rows.append(dict(
            case=cid, status="completed", family=d["family"], reference_source=reference_source,
            q_star=ref["analytic_check"]["q_star"], q_ref=ref["q_ref"],
            ref_rel_half_width=ref["half_width"] / ref["q_ref"],
            q_ref_minus_q_star_rel=ref["analytic_check"]["q_ref_minus_q_star"] / ref["analytic_check"]["q_star"],
            q_star_inside=ref["analytic_check"]["q_star_inside"],
            reference_converged=b["reference_converged"],
            coarse_rel={k: v["relative"] for k, v in e.items()},
            closer=b["closer"], r_coarse=b["r_coarse"],
            h_orders_omega=r["h_orders_omega"], h_orders_a=r["h_orders_a"],
            p_decay_omega=r["p_decay_omega"], p_decay_a=r["p_decay_a"],
            rate_verdict=r["verdict"], excluded=r["excluded_levels"], seconds=d["seconds"]))
    done = [x for x in rows if x["status"] == "completed"]
    agg = dict(
        attempted=len(ORDER), completed=len(done),
        closer_counts_decided={k: sum(1 for x in done if x["reference_converged"] and x["closer"] == k)
                               for k in ("midpoint", "2:1")},
        pending=[x["case"] for x in done if not x["reference_converged"] or x["closer"] == "inconclusive"],
        reference_converged=sum(1 for x in done if x["reference_converged"]),
        rate_verdict_counts={k: sum(1 for x in done if x["rate_verdict"] == k)
                             for k in ("equal", "mixed", "different", "undetermined")},
        q_star_inside_reference=sum(1 for x in done if x["q_star_inside"]))
    (root / "summary.json").write_text(json.dumps(dict(aggregate=agg, cases=rows, files=files), indent=2))
    f = lambda v, s="{:+.3e}": "-" if v is None else s.format(v)
    lines = ["| case | Q_ref half-width (rel) | Q_ref−Q* (rel) | coarse A | coarse Ω | coarse 1:1 | coarse 2:1 | closer | r(coarse) | h order Ω/A (finest pair) | p decay Ω/A | rate verdict |",
             "|---|---:|---:|---:|---:|---:|---:|---|---:|---|---|---|"]
    for x in rows:
        if x["status"] != "completed":
            lines.append(f"| {x['case']} | missing |" + " |" * 11)
            continue
        c = x["coarse_rel"]
        ho = (f"{x['h_orders_omega'][-1]:.2f}/{x['h_orders_a'][-1]:.2f}" if x["h_orders_omega"] else "-")
        pd = (f"{x['p_decay_omega']:.2f}/{x['p_decay_a']:.2f}" if x["p_decay_omega"] is not None else "-")
        lines.append(f"| {x['case']} | {f(x['ref_rel_half_width'], '{:.1e}')} | {f(x['q_ref_minus_q_star_rel'])} | "
                     f"{f(c['A'])} | {f(c['Omega'])} | {f(c['midpoint'])} | {f(c['two_to_one'])} | {x['closer']} | "
                     f"{f(x['r_coarse'], '{:.3f}')} | {ho} | {pd} | {x['rate_verdict']} |")
    (root / "summary.md").write_text("\n".join(lines) + "\n\n" + json.dumps(agg, indent=2) + "\n")
    print("\n".join(lines)); print(json.dumps(agg, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
