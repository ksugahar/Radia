"""Closed-torus magnetostatic energy with a physical, conserved loop current.

For a divergence-free surface current, c(J)=integral J.grad(phi)/(2*pi) dS
equals current across a meridional cut: coarea gives the average over phi,
and conservation makes every cut current equal. This avoids assigning
coefficient-space ones. NGSolve owns the Piola basis and surface divergence.
This is a zero-resistance surface-current model, not a finite-frequency port.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
from scipy.linalg import null_space


def constrained_energy(matrix, divergence, current):
    basis = null_space(divergence, rcond=1e-11)
    raw_reduced = basis.T @ matrix @ basis
    symmetry = np.linalg.norm(raw_reduced-raw_reduced.T)/np.linalg.norm(raw_reduced)
    if not np.isfinite(symmetry) or symmetry > 1e-6:
        raise ValueError('Reduced magnetic energy matrix is not symmetric within quadrature tolerance')
    reduced = (raw_reduced + raw_reduced.T)/2
    # Strict convexity on the divergence-free subspace certifies a minimum.
    np.linalg.cholesky(reduced)
    load = basis.T @ current
    if np.linalg.norm(load) <= 1e-6*np.linalg.norm(current):
        raise ValueError('Mesh has no resolved divergence-free toroidal current; refine/check topology')
    solution = basis @ np.linalg.solve(reduced, load)
    solution /= current @ solution
    energy = float(solution @ matrix @ solution)
    gradient = reduced @ (basis.T @ solution)
    kkt = float(np.linalg.norm(gradient-energy*load)/
                max(np.linalg.norm(gradient), abs(energy)*np.linalg.norm(load)))
    if not np.isfinite(kkt) or kkt > 1e-10:
        raise ValueError('Constrained magnetic energy KKT residual exceeds tolerance')
    return solution, energy, basis.shape[1], {
        'reduced_relative_symmetry_error': float(symmetry),
        'reduced_positive_definite': True, 'relative_kkt_residual': kkt}


def run(notebook, output, curvaturesafety, bonus):
    import ngsolve as ng
    from ngsolve.bem import LaplaceSL
    ng.SetNumThreads(2)
    raw = Path(notebook).read_bytes()
    nb = json.loads(raw)
    source = '\n'.join(''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code')
    definitions = [n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef)
                   and n.name == 'create_torus_mesh']
    namespace = {'np': np}
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(notebook), 'exec'), namespace)
    started = time.perf_counter()
    with ng.TaskManager():
        mesh = namespace['create_torus_mesh'](0.05, 0.005, curvaturesafety)
        space = ng.HDivSurface(mesh, order=0)
        if space.ndof > 2000:
            raise RuntimeError(f'Dense budget exceeded: {space.ndof}')
        scalar = ng.SurfaceL2(mesh, order=0)
        u, v = space.TnT()
        q = scalar.TestFunction()
        measure = ng.ds(bonus_intorder=bonus)
        operator = LaplaceSL(u.Trace()*measure, use_fmm=False)*v.Trace()*measure
        matrix = np.array(operator.mat.ToDense().NumPy()) * (4e-7*np.pi)
        dform = ng.BilinearForm(trialspace=space, testspace=scalar)
        dform += ng.div(u.Trace())*q*measure
        dform.Assemble()
        divergence = np.array(dform.mat.ToDense().NumPy())
        r2 = ng.x*ng.x+ng.y*ng.y
        grad_phi = ng.CF((-ng.y/r2, ng.x/r2, 0))
        load = ng.LinearForm(space)
        load += (ng.InnerProduct(v.Trace(), grad_phi)/(2*np.pi))*measure
        load.Assemble()
        current = load.vec.FV().NumPy().copy()
        area = float(ng.Integrate(1, mesh, ng.BND))
    active = sorted({d for element in mesh.Elements(ng.BND)
                     for d in space.GetDofNrs(element) if d >= 0})
    matrix = matrix[np.ix_(active, active)]
    divergence = divergence[:, active]
    current = current[active]
    solution, inductance, nullity, optimality = constrained_energy(matrix, divergence, current)
    signs = np.where(np.arange(len(active))%2, -1., 1.)
    transformed, invariant, _, _ = constrained_energy(
        signs[:, None]*matrix*signs[None, :], divergence*signs[None, :], current*signs)
    conservation = float(np.linalg.norm(divergence@solution) /
                         (np.linalg.norm(divergence)*np.linalg.norm(solution)))
    report = {
        'schema': 'radia.torus-physical-current.v1', 'host': platform.node(),
        'python': platform.python_version(), 'ngsolve': ng.__version__,
        'numpy': np.__version__, 'notebook_sha256': hashlib.sha256(raw).hexdigest(),
        'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'curvaturesafety': curvaturesafety, 'bonus_intorder': bonus, 'threads': 2,
        'ndof': space.ndof, 'active_boundary_dofs': len(active),
        'divergence_free_dimension': nullity,
        'optimality': optimality,
        'current_A': float(current@solution), 'relative_conservation_residual': conservation,
        'L_H': inductance, 'basis_transformed_L_H': invariant,
        'relative_basis_invariance_error': float(abs(invariant-inductance)/abs(inductance)),
        'relative_symmetry_error': float(np.linalg.norm(matrix-matrix.T)/np.linalg.norm(matrix)),
        'relative_surface_area_error': area/(4*np.pi**2*.05*.005)-1,
        'elapsed_s': time.perf_counter()-started,
        'scope': 'Divergence-free unit toroidal current, minimum magnetic energy; no finite-frequency, resistive, or experimental acceptance.',
    }
    Path(output).write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)
    assert inductance > 0
    assert 1e-8 < inductance < 1e-6, 'Invalid scale for this 50 mm torus'
    assert conservation < 1e-10
    np.testing.assert_allclose(current@solution, 1., atol=1e-10)
    np.testing.assert_allclose(invariant, inductance, rtol=1e-9)
    np.testing.assert_allclose(transformed, signs*solution, rtol=1e-8, atol=1e-8)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--notebook', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--curvaturesafety', type=float, default=.65)
    parser.add_argument('--bonus', type=int, default=4)
    args = parser.parse_args()
    run(args.notebook, args.output, args.curvaturesafety, args.bonus)
