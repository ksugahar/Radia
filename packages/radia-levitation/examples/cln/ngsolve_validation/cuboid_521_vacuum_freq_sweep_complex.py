"""NGSolve complex frequency sweep on cuboid 5x2x1 vacuum (air-box AIR_SCALE=5).

Solve at complex s = i*2*pi*f for each of 30 ELF frequencies and compare:
  - Y_NGSolve(i omega)  (this is the "true" answer for the air-box-truncated problem)
  - Y_ELF(i omega)      (the BEM answer, for the true vacuum problem)

The DIFFERENCE between these two reveals the air-box truncation error.

Then Foster-fit Y_NGSolve(i omega) and compare its tau_lead with:
  - Kameari Stage 0 tau (104.4 us at AIR_SCALE=5)
  - ELF Foster tau_lead (11.51 us)

If Foster-fit-NG gives ~104 us  -> air-box problem genuinely has tau ~104 us (Kameari is faithful, but to wrong physics)
If Foster-fit-NG gives ~11 us   -> Kameari is broken (air-box and BEM agree but Kameari doesn't recover it)

(Bigger picture: this disambiguates "Kameari iteration is broken" vs
"air-box problem is fundamentally different from radiation problem").
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from netgen.occ import Box, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from collections import deque
from math import pi
import csv, json, time
import numpy as np
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

AIR_SCALE = 5
H_COND = 0.30e-3
H_AIR = 1.5e-3
ORDER = 2


def build_spanning_tree(mesh):
    nv = mesh.nv
    visited = [False] * nv
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        adj[v0].append((v1, ed.nr))
        adj[v1].append((v0, ed.nr))
    visited[0] = True
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for vn, edn in adj[v]:
            if not visited[vn]:
                visited[vn] = True
                tree_edges.append(edn)
                queue.append(vn)
    return tree_edges


def build_geometry():
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = H_COND
    Bx, By, Bz = AIR_SCALE * ax / 2, AIR_SCALE * ay / 2, AIR_SCALE * az / 2
    air = Box(Pnt(-Bx, -By, -Bz), Pnt(Bx, By, Bz))
    air.mat("air").bc("outer_box")
    air.maxh = H_AIR
    return OCCGeometry(Glue([air - cuboid, cuboid]))


def main():
    ngsglobals.msg_level = 0
    print("=== NGSolve complex frequency sweep on air-box mesh (AIR_SCALE=5) ===")
    print(f"  Conductor: 5x2x1 mm Cu")
    print(f"  Air box: {AIR_SCALE}x conductor")
    print(f"  h_cond = {H_COND*1000} mm, h_air = {H_AIR*1000} mm, order = {ORDER}\n")

    # Load ELF freq list and m_z values
    elf_data = []
    csv_path = Path(__file__).parent.parent / "elf_validation" / "Yzz_elf_30f.csv"
    with open(csv_path) as fp:
        r = csv.reader(fp); next(r)
        for row in r:
            elf_data.append({
                "freq": float(row[0]),
                "mz_re_elf": float(row[1]),
                "mz_im_elf": float(row[2]),
            })
    print(f"Loaded {len(elf_data)} ELF freqs from {csv_path.name}\n")

    # Build mesh
    geo = build_geometry()
    print("Generating mesh...")
    mesh = Mesh(geo.GenerateMesh(maxh=H_AIR))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}")

    # Complex HCurl with tree-cotree gauge
    fes = HCurl(mesh, order=ORDER, dirichlet="outer_box", complex=True, nograds=True)
    print(f"  HCurl ndof = {fes.ndof}")

    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  Tree-cotree masked {masked} edges, active = {sum(fd)}\n")

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0})
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5  # B_0 = 1, z-dir

    u, v = fes.TnT()
    results = []

    print(f"Sweeping 30 frequencies at complex s = i*2*pi*f...\n")
    for k, d in enumerate(elf_data):
        f_hz = d["freq"]
        omega = 2 * pi * f_hz
        s_complex = 1j * omega
        t0 = time.time()

        a = BilinearForm(fes)
        a += (1.0/mu0) * curl(u) * curl(v) * dx
        a += s_complex * sigma_cf * u * v * dx

        rhs = LinearForm(fes)
        rhs += -s_complex * sigma_cf * A_ext * v * dx("conductor")

        with TaskManager():
            a.Assemble()
            rhs.Assemble()
            try:
                inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
            except Exception as e:
                print(f"  f={f_hz:.2e}: solver failed: {e}")
                continue
        gf = GridFunction(fes)
        with TaskManager():
            gf.vec.data = inv * rhs.vec

        # m_z = (1/2) integral_cond (r x J_total)_z dV
        # J_total = -s sigma (A_ext + A_ind) in conductor (induced eddy current)
        J_total = -s_complex * sigma_cf * (A_ext + gf)
        # m_z = (1/2) integral (x J_y - y J_x) dV
        m_int = Integrate((x * J_total[1] - y * J_total[0]) * 0.5 * dx("conductor"), mesh)
        m_re = float(m_int.real)
        m_im = float(m_int.imag)

        elapsed = time.time() - t0
        rel_re = (m_re / d["mz_re_elf"] - 1) * 100 if d["mz_re_elf"] != 0 else float('nan')
        print(f"  [{k+1:2d}/30] f={f_hz:.2e} Hz: m_NG=({m_re:.4e}+{m_im:.4e}j), "
              f"ELF=({d['mz_re_elf']:.4e}+{d['mz_im_elf']:.4e}j)  ({elapsed:.1f}s)")

        results.append({
            "freq": f_hz,
            "omega": omega,
            "m_re_NG": m_re,
            "m_im_NG": m_im,
            "m_re_ELF": d["mz_re_elf"],
            "m_im_ELF": d["mz_im_elf"],
        })

    # Save
    out_path = Path(__file__).parent / "cuboid_521_vacuum_NG_30f.csv"
    with open(out_path, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["freq_Hz", "m_re_NG", "m_im_NG", "m_re_ELF", "m_im_ELF",
                    "ratio_re", "ratio_im", "ratio_abs"])
        for r in results:
            mNG = r["m_re_NG"] + 1j * r["m_im_NG"]
            mELF = r["m_re_ELF"] + 1j * r["m_im_ELF"]
            w.writerow([r["freq"], r["m_re_NG"], r["m_im_NG"],
                        r["m_re_ELF"], r["m_im_ELF"],
                        r["m_re_NG"] / r["m_re_ELF"] if r["m_re_ELF"] != 0 else float('nan'),
                        r["m_im_NG"] / r["m_im_ELF"] if r["m_im_ELF"] != 0 else float('nan'),
                        abs(mNG) / abs(mELF) if abs(mELF) > 0 else float('nan')])
    print(f"\nSaved to {out_path}")

    # Foster fit on NGSolve data
    from scipy.optimize import least_squares
    s_arr = np.array([1j * 2*pi*r["freq"] for r in results])
    m_NG_arr = np.array([r["m_re_NG"] + 1j*r["m_im_NG"] for r in results])
    m_ELF_arr = np.array([r["m_re_ELF"] + 1j*r["m_im_ELF"] for r in results])
    M_NG = -m_NG_arr / (mu0 * s_arr)
    M_ELF = -m_ELF_arr / (mu0 * s_arr)

    def fit_freq(M_arr_target, N):
        init = np.logspace(np.log10(0.5e-6), np.log10(200e-6), N)
        def residual(x):
            taus = np.exp(x[:N])
            a = np.exp(x[N:2*N]) * np.sign(M_arr_target.real[0])
            fit = np.zeros_like(s_arr, dtype=complex)
            for tk, ak in zip(taus, a):
                fit += ak / (1 + s_arr*tk)
            return np.concatenate([(fit-M_arr_target).real/np.abs(M_arr_target),
                                   (fit-M_arr_target).imag/np.abs(M_arr_target)])
        M0 = abs(M_arr_target[0])
        x0 = np.concatenate([np.log(init), np.log(np.ones(N)*M0/N)])
        res = least_squares(residual, x0, method='trf', max_nfev=5000)
        taus = np.exp(res.x[:N])
        a = np.exp(res.x[N:2*N]) * np.sign(M_arr_target.real[0])
        order = np.argsort(-taus)
        # rms
        M_fit = np.zeros_like(s_arr, dtype=complex)
        for tk, ak in zip(taus[order], a[order]):
            M_fit += ak / (1 + s_arr*tk)
        rms = np.sqrt(np.mean(np.abs(M_fit-M_arr_target)**2 / np.abs(M_arr_target)**2)) * 100
        return taus[order], a[order], rms

    print("\n" + "=" * 78)
    print("FOSTER FIT COMPARISON")
    print("=" * 78)
    for N in [4, 6, 8]:
        taus_ng, a_ng, rms_ng = fit_freq(M_NG, N)
        taus_elf, a_elf, rms_elf = fit_freq(M_ELF, N)
        print(f"\n--- N = {N} stage ---")
        print(f"  NGSolve air-box: tau = {[f'{t*1e6:.2f}' for t in taus_ng]} us, rms = {rms_ng:.3f}%")
        print(f"  ELF (vacuum):    tau = {[f'{t*1e6:.2f}' for t in taus_elf]} us, rms = {rms_elf:.3f}%")
        print(f"  ratio tau_lead (NG/ELF) = {taus_ng[0]/taus_elf[0]:.3f}")

    print("\n=== Critical comparison ===")
    print(f"  Kameari Stage 0 tau     = 104.37 us (AIR_SCALE=5, from breakdown.json)")
    taus_ng4, _, _ = fit_freq(M_NG, 4)
    print(f"  Foster-fit NG tau_lead  = {taus_ng4[0]*1e6:.2f} us (this script, AIR_SCALE=5)")
    print(f"  Foster-fit ELF tau_lead = 11.51 us")
    print()
    if abs(taus_ng4[0]*1e6 - 104.37) / 104.37 < 0.10:
        print("  [Kameari faithful to air-box problem; air-box != BEM]")
    elif abs(taus_ng4[0]*1e6 - 11.51) / 11.51 < 0.10:
        print("  [Kameari iteration is broken: NG freq sweep agrees with ELF, but Kameari ladder differs]")
    else:
        print(f"  [Different from both: NG = {taus_ng4[0]*1e6:.2f} us]")

    out_json = Path(__file__).parent / "cuboid_521_vacuum_NG_freq_sweep.json"
    out_data = {
        "method": "NGSolve complex freq sweep on air-box mesh",
        "AIR_SCALE": AIR_SCALE,
        "H_COND_mm": H_COND*1000, "H_AIR_mm": H_AIR*1000, "order": ORDER,
        "ne": mesh.ne, "ndof": fes.ndof,
        "frequencies": [r["freq"] for r in results],
        "m_NG_complex": [(r["m_re_NG"], r["m_im_NG"]) for r in results],
        "m_ELF_complex": [(r["m_re_ELF"], r["m_im_ELF"]) for r in results],
        "Foster_fit": {
            "NG_4stage": {
                "taus_us": [t*1e6 for t in fit_freq(M_NG, 4)[0].tolist()],
                "a": fit_freq(M_NG, 4)[1].tolist(),
                "rms_pct": fit_freq(M_NG, 4)[2],
            },
            "ELF_4stage": {
                "taus_us": [t*1e6 for t in fit_freq(M_ELF, 4)[0].tolist()],
                "a": fit_freq(M_ELF, 4)[1].tolist(),
                "rms_pct": fit_freq(M_ELF, 4)[2],
            },
        },
        "Kameari_stage0_tau_us": 104.37,
    }
    out_json.write_text(json.dumps(out_data, indent=2))
    print(f"\nResults: {out_json}")


if __name__ == "__main__":
    main()
