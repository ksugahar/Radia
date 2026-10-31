"""Independent convolution/Duffy references for finite-volume PEEC entries.

Reference: Neumann volume energy (Ruehli, IEEE MTT 1974,
DOI 10.1109/TMTT.1974.1128204). The reference below transforms the singular
six-dimensional self integral to three smooth two-dimensional integrals;
it does not reuse the production Laplace/Gaussian factorization.
This exploratory report is not a fin adoption verdict or stored baseline.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import itertools
import math
from pathlib import Path
import subprocess
import time
import numpy as np


def duffy_self_H(widths_m, current_axis, order):
    """Integrate 8*prod(d_i-|s_i|)/|s| over the positive difference box.

    Scale differences s_i=d_i*y_i and split the unit cube according to its
    largest y coordinate. In each sector use (t,t*u,t*v), Jacobian t**2;
    integrate t*(1-t)*(1-t*u)*(1-t*v) exactly. No empirical coefficients.
    """
    d = np.asarray(widths_m, dtype=float)
    if d.shape != (3,) or not np.isfinite(d).all() or np.any(d <= 0):
        raise ValueError('Three positive finite dimensions required')
    z,w = np.polynomial.legendre.leggauss(order)
    u = (z[:,None]+1)/2;v = (z[None,:]+1)/2
    weight = w[:,None]*w[None,:]/4
    polynomial = .5-(1+u+v)/3+(u+v+u*v)/4-u*v/5
    denominator = (1/np.sqrt(d[0]**2+(d[1]*u)**2+(d[2]*v)**2)
                   +1/np.sqrt((d[0]*u)**2+d[1]**2+(d[2]*v)**2)
                   +1/np.sqrt((d[0]*u)**2+(d[1]*v)**2+d[2]**2))
    integral = 8*np.prod(d)**2*np.sum(weight*polynomial*denominator)
    area = np.prod(np.delete(d,current_axis))
    return float(1e-7*integral/area**2)


def convolution_mutual_H(box_a, axis_a, box_b, axis_b, order, sign_a=1, sign_b=1):
    """Six-volume Neumann integral via piecewise-linear difference densities.

    For s=x-y, w(s)=length(A intersect (B+s)); the 6D integral equals
    integral prod(w_i(s_i))/|s| ds. Split at all overlap breakpoints and zero.
    Each positive-octant cell containing the origin uses three Duffy sectors;
    its radial polynomial is integrated exactly. Other cells use a 3D Gauss
    rule. No production Gaussian/Laplace kernel, fitted value or self formula
    is called. Different box sizes and overlap are handled directly.
    """
    a, b = np.asarray(box_a, dtype=float), np.asarray(box_b, dtype=float)
    if (a.shape != (2, 3) or b.shape != (2, 3) or not np.isfinite([a, b]).all()
            or np.any(a[1] <= a[0]) or np.any(b[1] <= b[0])):
        raise ValueError('Two finite positive Cartesian boxes required')
    if axis_a not in (0, 1, 2) or axis_b not in (0, 1, 2) or sign_a not in (-1, 1) or sign_b not in (-1, 1):
        raise ValueError('Explicit Cartesian axes and orientations required')
    if axis_a != axis_b:
        return 0.  # Neumann J_a dot J_b, not a singular 0/0 quotient
    z, weights = np.polynomial.legendre.leggauss(order)
    uv = (z+1)/2
    u, v = uv[:, None], uv[None, :]
    weights_uv = weights[:, None]*weights[None, :]/4
    intervals = []
    for axis in range(3):
        lo, hi = a[0, axis]-b[1, axis], a[1, axis]-b[0, axis]
        knots = {lo, hi, a[0, axis]-b[0, axis], a[1, axis]-b[1, axis]}
        if lo < 0 < hi:
            knots.add(0.)
        knots = sorted(knots)
        def overlap(s):
            return max(0., min(a[1, axis], b[1, axis]+s)-max(a[0, axis], b[0, axis]+s))
        pieces = []
        for left, right in zip(knots[:-1], knots[1:]):
            if right <= 0:
                pieces.append((-right, -left, overlap(right), overlap(left)))
            else:
                pieces.append((left, right, overlap(left), overlap(right)))
        intervals.append(pieces)
    terms = []
    for cell in itertools.product(*intervals):
        lo = np.array([c[0] for c in cell]);hi = np.array([c[1] for c in cell])
        width = hi-lo
        alpha = np.array([c[2] for c in cell]);beta = np.array([c[3]-c[2] for c in cell])
        if np.all(lo == 0):
            sector_integrals = []
            for dominant in range(3):
                c = [None]*3;c[dominant] = np.ones_like(u*v)
                remaining = [i for i in range(3) if i != dominant]
                c[remaining[0]], c[remaining[1]] = u, v
                # integral_0^1 t prod(alpha_i+beta_i*c_i*t) dt
                radial = np.zeros_like(u*v)
                for selected in itertools.product((False, True), repeat=3):
                    product = np.ones_like(radial)
                    for i in range(3):
                        product = product*(beta[i]*c[i] if selected[i] else alpha[i])
                    radial += product/(2+sum(selected))
                denominator = np.sqrt(sum((width[i]*c[i])**2 for i in range(3)))
                sector_integrals.append(float(np.prod(width)*np.sum(weights_uv*radial/denominator)))
            terms.append(math.fsum(sector_integrals))
        else:
            nodes = [lo[i]+width[i]*uv for i in range(3)]
            weighted = [weights*width[i]/2*(alpha[i]+beta[i]*uv) for i in range(3)]
            denominator = np.sqrt(nodes[0][:, None, None]**2+nodes[1][None, :, None]**2+nodes[2][None, None, :]**2)
            terms.append(float(np.sum(weighted[0][:, None, None]*weighted[1][None, :, None]
                                      *weighted[2][None, None, :]/denominator)))
    area_a = np.prod(np.delete(a[1]-a[0], axis_a))
    area_b = np.prod(np.delete(b[1]-b[0], axis_b))
    return float(sign_a*sign_b*1e-7*math.fsum(terms)/(area_a*area_b))


def _add_box(builder, box, axis, sign=1):
    size = box[1]-box[0]
    direction = np.zeros(3);direction[axis] = sign
    builder.add_segment(box.mean(axis=0), direction, size[axis],
                        size[1 if axis == 2 else 2], size[1 if axis == 0 else 0], 5.8e7)


def _compression_report(peec_matrices, native):
    """Caller owns TaskManager; operator parity is not independent physics."""
    n = 193
    geometry = dict(centers=np.column_stack((np.arange(n)*.0025, np.zeros(n), np.zeros(n))),
                    directions=np.tile([0., 0., 1.], (n, 1)), lengths=np.full(n, .003),
                    widths=np.full(n, .001), heights=np.full(n, .0005), sigmas=np.full(n, 5.8e7))
    builder = peec_matrices.PEECBuilder(volume_order=64)
    for i in range(n):
        builder.add_segment(geometry['centers'][i], geometry['directions'][i],
                            geometry['lengths'][i], geometry['widths'][i], geometry['heights'][i], geometry['sigmas'][i])
    started = time.perf_counter();L, _, _, _ = builder.build(include_star=False)
    dense_s = time.perf_counter()-started
    entries = np.array([[builder.inductance_entry(i, j) for j in range(n)] for i in range(n)])
    exact = L.tobytes() == entries.tobytes()
    if not exact:
        raise RuntimeError('Dense and on-demand volume entries differ in raw bytes')
    manager = native.HACApKPEECManager(**geometry, volume_order=64)
    started = time.perf_counter()
    if not manager.BuildHMatrix(aca_eps=1e-12, leaf_size=8, eta=1.5, max_rank=64, print_level=0):
        raise RuntimeError('Explicit finite-volume HACApK build failed')
    hierarchical_s = time.perf_counter()-started
    vectors = np.random.default_rng(381).normal(size=(4, n))
    errors = []
    for vector in vectors:
        actual = np.asarray(manager.MatVec(vector));original = L@vector
        if not np.isfinite(actual).all():
            raise RuntimeError('Nonfinite HACApK action')
        errors.append(float(np.linalg.norm(actual-original)/np.linalg.norm(original)))
    stats = manager.GetStats()
    return dict(status='DIAGNOSTIC_NOT_ADOPTION', n=n, geometry=geometry['centers'].tolist(),
                dimensions_m=[.0005, .001, .003], current_axis=2, volume_order=64,
                aca=dict(eps=1e-12, leaf_size=8, eta=1.5, max_rank=64), vector_seed=381,
                dense_vs_entry_raw_bytes_equal=exact, hacapk_stats=stats,
                low_rank_observed=stats['n_lowrank'] > 0, original_dense_action_relative_errors=errors,
                dense_build_s=dense_s, hacapk_build_s=hierarchical_s,
                note='Same kernel operator parity; not independent physics, fin speedup, GPU parity or an adoption verdict')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    from radia import peec_matrices
    from radia import _radia_pybind as native
    import radia
    import ngsolve as ng
    module = Path(peec_matrices.__file__).resolve()
    root = Path(__file__).resolve().parents[2]
    source = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    for imported in (radia, peec_matrices, native):
        if not Path(imported.__file__).resolve().is_relative_to(root):
            raise RuntimeError('Wrong worktree import: '+imported.__file__)
    rows=[];started=time.perf_counter()
    result=dict(status='INCOMPLETE',source_sha=source,criterion_approval='PENDING',
                physics_verdict=None,gpu_parity='NOT_RUN',radia_path=radia.__file__,
                native_path=str(module),native_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),
                core_native_path=native.__file__,core_native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
                conductivity_S_per_m=5.8e7,backend='cpu',units='SI',
                entry_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                reference='Neumann six-volume energy; self Duffy and general difference-density convolution/Duffy, no production kernel used',
                native_orders=[16,32,64,128],reference_orders=[64,128],rows=rows,mutual_rows=[])
    def save():
        result['runtime_s'] = time.perf_counter()-started
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    save()
    with ng.TaskManager():
        for dims in ((1e-3,1e-3,1e-3),(4e-3,1e-3,1e-3),(2e-3,2e-3,.25e-3)):
            for axis in (0,1,2):
                reference = {str(n):duffy_self_H(dims,axis,n) for n in (32,64,128)}
                values={};box=np.array([-np.array(dims)/2,np.array(dims)/2])
                for n in (16,32,64,128):
                    builder=peec_matrices.PEECBuilder(volume_order=n)
                    _add_box(builder,box,axis)
                    value=float(builder.inductance_entry(0,0))
                    values[str(n)]=dict(value_H=value,
                        relative_error=abs(value-reference['128'])/abs(reference['128']))
                rows.append(dict(dimensions_m=dims,current_axis=axis,
                    independent_reference_H=reference,
                    reference_64_to_128_relative_spread=abs(reference['64']-reference['128'])/abs(reference['128']),production=values))
                save()
        box_a=np.array([[0.,0.,0.],[1.,.6,.4]])*1e-3
        cases=(
            ('separated-x',box_a,0,np.array([[1.5,.2,.1],[2.3,.9,.7]])*1e-3,0,1),
            ('face-touching-x',box_a,0,np.array([[1.,0.,0.],[1.7,.6,.4]])*1e-3,0,1),
            ('overlapping-unequal-y',box_a,1,np.array([[.2,.1,-.15],[1.3,.9,.45]])*1e-3,1,1),
            ('separated-z',box_a,2,np.array([[1.5,.2,.1],[2.3,.9,.7]])*1e-3,2,1),
            ('overlapping-orthogonal-zero',box_a,0,np.array([[.2,.1,-.15],[1.3,.9,.45]])*1e-3,1,1),
            ('overlapping-unequal-y-reversed',box_a,1,np.array([[.2,.1,-.15],[1.3,.9,.45]])*1e-3,1,-1))
        for name,a,axis_a,b,axis_b,sign_b in cases:
            reference={str(n):convolution_mutual_H(a,axis_a,b,axis_b,n,sign_b=sign_b) for n in (64,128)}
            values={}
            for n in (16,32,64,128):
                builder=peec_matrices.PEECBuilder(volume_order=n)
                _add_box(builder,a,axis_a);_add_box(builder,b,axis_b,sign_b)
                value=float(builder.inductance_entry(0,1));error=abs(value-reference['128'])
                if not math.isfinite(value):raise RuntimeError('Nonfinite native mutual entry')
                values[str(n)]=dict(value_H=value,absolute_error_H=error,
                    relative_error=error/abs(reference['128']) if reference['128'] else None)
            spread=abs(reference['64']-reference['128'])
            result['mutual_rows'].append(dict(case=name,box_a_m=a.tolist(),box_b_m=b.tolist(),axis_a=axis_a,axis_b=axis_b,sign_a=1,sign_b=sign_b,
                independent_reference_H=reference,reference_64_to_128_absolute_spread_H=spread,
                reference_64_to_128_relative_spread=spread/abs(reference['128']) if reference['128'] else None,production=values))
            save()
        result['dense_hacapk_operator'] = _compression_report(peec_matrices,native)
    result['status']='DIAGNOSTIC_NOT_ADOPTION'
    save()


if __name__=='__main__':main()
