"""Complete body reaction using the transpose of the incident-field maps.

Work is bilinear in phasors. Conjugation is applied only to terminal power.
The surface-current kernel and Poisson projection match the forward solvers.
"""
from __future__ import annotations

import numpy as np

MU_0 = 4e-7*np.pi


def surface_coil_reaction_rhs(emf, omega):
    """Physical EMF to the unscaled Laplace-SL coil saddle RHS."""
    if not np.isfinite(omega) or omega <= 0:
        raise ValueError("Coil reaction needs finite positive frequency")
    return -np.asarray(emf, dtype=complex)/(1j*omega*MU_0)


def _iterate_complete_current(initial, response, map_current, *, max_iter, tol,
                              relax, verbose=False):
    """Certify body and coil residuals at one unchanged returned state."""
    current = np.asarray(initial, dtype=complex).copy()
    residual = float('inf')
    for iteration in range(max_iter):
        state = response(current)
        body_residual = state['body']['residual']
        if not np.isfinite(body_residual) or body_residual > 1e-6:
            raise RuntimeError("Strong body residual exceeds 1e-6")
        mapped = map_current(state['emf'])
        residual = float(np.linalg.norm(mapped-current)/max(np.linalg.norm(current), 1e-30))
        if not np.isfinite(residual):
            raise RuntimeError("Nonfinite strong coil residual")
        if verbose:
            print(f'  iter {iteration}: coil fixed-point residual={residual:.3e}, '
                  f'body residual={body_residual:.3e}')
        if residual <= tol and iteration > 0:
            return current, state, residual, iteration+1
        current = (1-relax)*current+relax*mapped
    raise RuntimeError(f'Strong coupling did not converge: coil residual={residual:.3e}, '
                       f'tol={tol:.3e}, max_iter={max_iter}')


def electric_incident_vertex_load(poisson, impedance, total_tangential_field):
    """Transpose of incident surface-Poisson projection, including loop Ht."""
    field = np.asarray(total_tangential_field, dtype=complex)
    tri, g, area = poisson._tri, poisson._g, poisson._areas
    from .surface_impedance import PanelSurfaceImpedance
    z = np.asarray(impedance.values if isinstance(impedance, PanelSurfaceImpedance)
                   else impedance, dtype=complex)
    if z.ndim == 0:
        z = np.full(len(tri), z, dtype=complex)
    if (field.shape != (len(tri), 3) or z.shape != (len(tri),)
            or not np.all(np.isfinite(field)) or not np.all(np.isfinite(z))
            or np.any(z.real < 0)):
        raise ValueError("Complete reaction needs finite passive face impedance and total Ht")
    local = -area[:, None]*z[:, None]*np.einsum('tkd,td->tk', g, field)
    load = np.zeros(poisson._nv+1, dtype=complex)
    for k in range(3):
        np.add.at(load, tri[:, k], local[:, k])
    lam = (poisson._lu.solve(load.real, trans='T')
           + 1j*poisson._lu.solve(load.imag, trans='T'))[:poisson._nv]
    grad_lam = np.einsum('tkd,tk->td', g, lam[tri])
    vertex_load = np.zeros((poisson._nv, 3), dtype=complex)
    for k in range(3):
        np.add.at(vertex_load, tri[:, k], -area[:, None]*grad_lam/3)
    return vertex_load


def surface_source_reaction_load(source_centroids, source_areas, body_points,
                                 body_centroids, body_areas, body_current,
                                 incident_vertex_load, omega,
                                 *, pair_budget=1_000_000):
    """Return EMF load on each Cartesian source-panel average current."""
    source = np.asarray(source_centroids, dtype=float)
    source_area = np.asarray(source_areas, dtype=float)
    points = np.asarray(body_points, dtype=float)
    centers = np.asarray(body_centroids, dtype=float)
    moments = np.asarray(body_areas)[:, None]*np.asarray(body_current)
    hload = np.asarray(incident_vertex_load)
    if (source.shape != (len(source_area), 3) or centers.shape != moments.shape
            or points.shape != hload.shape or not np.isfinite(omega) or omega <= 0
            or pair_budget < 1 or np.any(source_area <= 0)
            or any(not np.all(np.isfinite(x)) for x in
                   [source, source_area, points, centers, moments, hload])):
        raise ValueError("Invalid complete-reaction geometry or field")
    result = np.zeros(source.shape, dtype=complex)
    width = max(1, int(pair_budget)//max(len(points), len(centers), 1))
    for start in range(0, len(source), width):
        stop = min(start+width, len(source))
        d = centers[:, None, :]-source[None, start:stop, :]
        distance = np.maximum(np.linalg.norm(d, axis=2), 1e-12)
        magnetic = MU_0/(4*np.pi)*np.einsum('ts,td->sd', 1/distance, moments)
        d = points[:, None, :]-source[None, start:stop, :]
        r2 = np.maximum(np.sum(d*d, axis=2), 1e-60)
        electric = np.sum(np.cross(d, hload[:, None, :])
                          / (r2*np.sqrt(r2))[:, :, None], axis=0)/(4*np.pi)
        result[start:stop] = source_area[start:stop, None]*(1j*omega*magnetic+electric)
    return result


def surface_current_average_maps(fes):
    """Sparse coefficient-to-face-current maps matching native integration."""
    from ngsolve import BND, BilinearForm, CF, Integrate, SurfaceL2, ds
    from scipy.sparse import coo_matrix, diags
    mesh = fes.mesh
    space = SurfaceL2(mesh, order=0, dual_mapping=False)
    elements = list(mesh.Elements(BND))
    areas = np.asarray(Integrate(CF(1), mesh, BND, element_wise=True))
    if np.any(areas <= 0) or len(elements) != space.ndof:
        raise ValueError("Complete reaction requires nondegenerate coil surface panels")
    dofs = [space.GetDofNrs(el)[0] for el in elements]
    rows = [el.nr for el in elements]
    current, test = fes.TrialFunction(), space.TestFunction()
    maps = []
    for component in range(3):
        form = BilinearForm(trialspace=fes, testspace=space)
        form += current.Trace()[component]*test*ds
        form.Assemble()
        r, c, v = form.mat.COO()
        matrix = coo_matrix((v, (r, c)), shape=(space.ndof, fes.ndof)).tocsr()
        maps.append(diags(1/areas[rows]) @ matrix[dofs])
    return maps


def strong_surface_impedance(impedance, triangle_count):
    """Validate tagged face values without accepting ambiguous nodal arrays."""
    from .surface_impedance import PanelSurfaceImpedance
    if isinstance(impedance, PanelSurfaceImpedance):
        if impedance.values.shape != (triangle_count,):
            raise ValueError("Strong panel Zs needs one value per body BND triangle")
        return impedance
    if np.ndim(impedance) != 0:
        raise ValueError("Strong impedance needs a scalar or tagged PanelSurfaceImpedance")
    value = complex(impedance)
    if not np.isfinite(value) or value.real < 0:
        raise ValueError("Strong impedance requires finite passive values")
    return value


def strong_impedance_metadata(impedance, areas):
    """Aggregate Z is diagnostic only; retain exact face values for exports."""
    from .surface_impedance import PanelSurfaceImpedance
    if not isinstance(impedance, PanelSurfaceImpedance):
        return dict(Z_s=complex(impedance))
    return dict(Z_s=complex(areas @ impedance.values/areas.sum()),
                Z_s_summary_kind='area-weighted-report-only',
                Z_s_per_panel=impedance.values.copy())


def solve_complete_body(solver, poisson, phi_inc, impedance, omega, a_inc,
                        *, loop_dof=False, hacapk=False, gmres=None):
    """Solve the body at this exact coil state and retain its total current."""
    if solver.order != 1 or solver.mesh.GetCurveOrder() > 1 or solver.mesh.deformation is not None:
        raise ValueError("Complete strong reaction requires undeformed flat P1 body panels")
    from .surface_impedance import PanelSurfaceImpedance
    z = strong_surface_impedance(impedance, len(poisson._tri))
    face_z = (z.values if isinstance(z, PanelSurfaceImpedance)
              else np.full(len(poisson._tri), z))
    # Exact equal-face values retain the adopted scalar arithmetic.
    solve_z = complex(face_z[0]) if np.all(face_z == face_z[0]) else z
    if hacapk and isinstance(z, PanelSurfaceImpedance):
        raise ValueError("Strong panel Zs currently requires the dense body backend")
    import hashlib
    from ngsolve import BND
    # Preparation must also bind vertex/DOF numbering, not only face coordinates.
    mesh_points = np.asarray([tuple(v.point) for v in solver.mesh.vertices], dtype='<f8')
    mesh_tri = np.asarray([[v.nr for v in t.vertices] for t in solver.mesh.Elements(BND)], dtype='<i8')
    identity = hashlib.sha256(mesh_points.tobytes()+mesh_tri.tobytes()).hexdigest()
    frozen = getattr(solver, '_strong_body_geometry_identity', identity)
    if frozen != identity:
        raise ValueError("Strong body geometry changed; rebuild the coupled solver")
    solver._strong_body_geometry_identity = identity
    value_hash = hashlib.sha256(face_z.astype('<c16').tobytes()).hexdigest()
    prepared_key = (identity, float(omega), value_hash)
    if loop_dof:
        if hacapk:
            raise ValueError("Loop strong reaction requires dense body operators")
        import time
        from .bem_loop_extension import solve_loop_extended
        start = time.perf_counter()
        out = solve_loop_extended(solver, phi_inc, solve_z, omega, a_inc, _reuse_prepared=True)
        field, heat = out['H_t_tri'], out['P_total']
        residual = out['linear_residual_rel']
        metadata = dict(wp_loop_alpha=out['alpha'], wp_loop_theta_jump=out['theta_jump'],
            wp_loop_cut_n_vertices=out['cut_n_vertices'], wp_loop_P_frozen=out['P_frozen'],
            wp_loop_H_t_frozen=out['Ht_frozen'],
            wp_loop_screening_ratio=heat/max(out['P_frozen'], 1e-300),
            t_loop_dof_s=time.perf_counter()-start,
            wp_loop_faraday_residual=out['faraday_residual_rel'],
            wp_loop_impedance_assembly=out['loop_impedance_assembly'],
            body_operator_reused=out['loop_system_reused'],
            body_geometry_reused=out['loop_geometry_reused'])
    else:
        out = (solver.solve_hacapk(phi_inc, solve_z, omega, **(gmres or {})) if hacapk
               else solver.solve(phi_inc, solve_z, omega, _prepared_key=prepared_key))
        field = -np.einsum('tkd,tk->td', poisson._g, out['phi_vec'][poisson._tri])
        heat = out['P_density']*out['area']
        residual = out['linear_residual_rel']
        metadata = {} if hacapk else dict(body_operator_reused=out['operator_reused'],
                                         body_factorization_reused=out['factorization_reused'])
    if not np.isfinite(residual) or residual > 1e-6:
        raise RuntimeError(f"Strong body residual {residual:.3e} exceeds 1e-6")
    areas, normals = poisson._areas, poisson._n_hat
    local_heat = .5*face_z.real*np.sum(abs(field)**2, axis=1)
    integrated = float(areas @ local_heat)
    if not np.isclose(integrated, heat, rtol=1e-10, atol=1e-30):
        raise RuntimeError("Strong total-field heat disagrees with the body solve")
    return dict(field=field, current=np.cross(normals, field),
                heat=float(heat), q_tri=local_heat, residual=float(residual),
                h_rms=float(np.sqrt(areas @ np.sum(abs(field)**2, axis=1)/areas.sum())),
                loop_metadata=dict(body_impedance_value_hash=value_hash, **metadata))
