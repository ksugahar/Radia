#!/usr/bin/env python
"""Loop-mode deflation on the C-type electromagnet (paper headline geometry).

The eigenvalue study's deflation was proven on TOY models (cubes <= 8^3, rings).
HACApK only "sells" on a LARGER, real problem -- so here we run the *same*
matrix-free loop-mode deflation (`SetDeflateNullspace`) on the full C-type yoke
at the paper's three mesh densities (6x6x6 / 10x10x10 / 20x20x20 = 1200 / 3150 /
27600 elements = 7200 / 18900 / 165600 DOF), full model (NO IMA, so deflation is
fully valid as validated on cubes/rings).

For each case we solve with HACApK BiCGSTAB, deflation OFF vs ON, and record the
*linear* (BiCGSTAB) iteration count, nonlinear iterations, wall time, peak
memory, installed deflation cycles, and B_z at the gap centre. The research
question: does loop deflation reduce the linear iteration count (-> time, memory
pressure) at scale on the REAL nonlinear magnet -- and how does it compare with
a high-mu_r LINEAR run where the near-null loop cluster is most active?

Mesh source: ELF reference .meg (LAB-local validation data); results are written
to a committed JSON next to this script so the figure can be regenerated.

Usage:
  python deflation_benchmark.py                 # nonlinear + linear, 6x6x6 & 10x10x10
  python deflation_benchmark.py --with-20       # also the 165600-DOF headline (SLOW)
  python deflation_benchmark.py --cases nonlinear:6x6x6 linear:10x10x10
"""
import argparse
import json
import os
import platform
import sys
import time
from datetime import datetime

import numpy as np
import psutil

work_dir = os.path.dirname(os.path.abspath(__file__))            # .../nonlinear
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(work_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, work_dir)

import radia as rad
from coil_model import create_racetrack_coil

ELF_BASE = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT"
SCALE = 0.001          # mm -> m
COIL_AT = 20000.0
ALPHA = 1.0            # deflation shift strength (O(1), as on the cubes)


def elf_dir(density):
    return os.path.join(ELF_BASE, f"ELF_MMB8T_EIEM2_{density}")


def load_geometry(density):
    nodes, elements = {}, []
    with open(os.path.join(elf_dir(density), "ELF_magic.meg")) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == 'MGR1':
                nodes[int(p[1])] = np.array([float(p[3]), float(p[4]), float(p[5])])
            elif p[0] == 'MMB8T':
                elements.append([int(p[i]) for i in range(4, 12)])
    return nodes, elements


def load_bh():
    return np.loadtxt(os.path.join(work_dir, "BH.txt"), comments='#').tolist()


def load_elf_bz(density):
    try:
        with open(os.path.join(elf_dir(density), "ELF_magic.mag")) as f:
            for line in f:
                p = line.split()
                if p and p[0] == 'M3GB' and len(p) >= 8 and int(p[1]) == 1:
                    return float(p[7])
    except OSError:
        pass
    return None


def build_model(nodes, elements, material):
    rad.UtiDelAll()
    if material[0] == 'nonlinear':
        mat = rad.MatSatIsoTab(load_bh())
    else:                                   # ('linear', mu_r)
        mat = rad.MatLin(float(material[1]))
    hexes = []
    for nid in elements:
        verts = [[nodes[n][0]*SCALE, nodes[n][1]*SCALE, nodes[n][2]*SCALE] for n in nid]
        h = rad.ObjHexahedron(verts, [0, 0, 0]); rad.MatApl(h, mat); hexes.append(h)
    yoke = rad.ObjCnt(hexes)
    coil = create_racetrack_coil(COIL_AT)
    return rad.ObjCnt([yoke, coil]), yoke, len(hexes)


def peak_mem_mb():
    m = psutil.Process(os.getpid()).memory_info()
    return (m.peak_wset if hasattr(m, 'peak_wset') else m.rss) / (1024*1024)


def solve_case(nodes, elements, material, deflate, max_nl):
    rad.SetDeflateNullspace(bool(deflate), ALPHA)
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    model, yoke, ne = build_model(nodes, elements, material)
    t0 = time.perf_counter()
    res = rad.Solve(model, 0.0001, max_nl, 2)      # method 2 = HACApK
    dt = time.perf_counter() - t0
    st = rad.GetSolveStats()
    Bz = float(np.array(rad.Fld(model, 'b', [0, 0, 0]))[2])
    return {
        'ne': ne, 'dof': ne*6,
        'nl_iter': int(res[2]) if len(res) > 2 else -1,
        'lin_iter': st.get('linear_iterations', -1),
        'defl_cycles': st.get('deflation_cycles', -1),
        't_solve': dt,
        't_hmat': st.get('t_hmatrix_build', float('nan')),
        'peak_mb': peak_mem_mb(),
        'Bz_mT': Bz*1000.0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cases', nargs='*', default=None,
                    help='material:density tokens, e.g. nonlinear:6x6x6 linear:10x10x10')
    ap.add_argument('--mu-r', type=float, default=1e4, help='mu_r for linear cases')
    ap.add_argument('--with-20', action='store_true', help='include 20x20x20 (165600 DOF, SLOW)')
    ap.add_argument('--max-nl', type=int, default=100, help='max nonlinear iterations')
    args = ap.parse_args()

    if args.cases:
        cases = []
        for tok in args.cases:
            mlab, dens = tok.split(':')
            cases.append((('nonlinear',) if mlab == 'nonlinear' else ('linear', args.mu_r), dens))
    else:
        dens_list = ['6x6x6', '10x10x10'] + (['20x20x20'] if args.with_20 else [])
        cases = [(('nonlinear',), d) for d in dens_list] + \
                [(('linear', args.mu_r), d) for d in dens_list]

    print(f"C-type loop-mode deflation benchmark (full model, COIL={COIL_AT:.0f} AT, "
          f"alpha={ALPHA}, BiCGSTAB tol=1e-8)")
    print(f"{'material':<14} {'dens':<9} {'dof':>7} | "
          f"{'lin(off)':>8} {'lin(on)':>7} {'cyc':>5} | "
          f"{'t_off':>7} {'t_on':>6} {'ratio':>5} | {'Bz off/on (mT)':>16}")

    results = []
    for material, dens in cases:
        mlabel = material[0] if material[0] == 'nonlinear' else f"linear mu={material[1]:.0e}"
        nodes, elements = load_geometry(dens)
        off = solve_case(nodes, elements, material, False, args.max_nl)
        on = solve_case(nodes, elements, material, True, args.max_nl)
        ratio = on['t_solve']/off['t_solve'] if off['t_solve'] > 0 else float('nan')
        elf_bz = load_elf_bz(dens)
        results.append({'material': mlabel, 'density': dens, 'dof': off['dof'],
                        'elf_bz_mT': (elf_bz*1000.0 if elf_bz is not None else None),
                        'off': off, 'on': on})
        print(f"{mlabel:<14} {dens:<9} {off['dof']:>7} | "
              f"{off['lin_iter']:>8} {on['lin_iter']:>7} {on['defl_cycles']:>5} | "
              f"{off['t_solve']:>7.1f} {on['t_solve']:>6.1f} {ratio:>5.2f} | "
              f"{off['Bz_mT']:>7.1f}/{on['Bz_mT']:<7.1f}")

    rad.SetDeflateNullspace(False, 0.0)
    rad.UtiDelAll()

    out = os.path.join(work_dir, "deflation_benchmark_results.json")
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({'timestamp': datetime.now().isoformat(),
                   'hostname': platform.node(),
                   'benchmark': 'ctype_loop_deflation',
                   'problem': {'coil_AT': COIL_AT, 'alpha': ALPHA,
                               'bicgstab_tol': 1e-8, 'model': 'full (no IMA)'},
                   'results': results}, f, indent=2, ensure_ascii=False)
    print(f"\nResults -> {out}")


if __name__ == '__main__':
    main()
