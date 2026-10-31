"""Independent Duffy reference for finite rectangular-volume PEEC self terms.

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
from pathlib import Path
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    from radia import peec_matrices
    import radia
    module = Path(peec_matrices.__file__).resolve()
    rows=[];started=time.perf_counter()
    for dims in ((1e-3,1e-3,1e-3),(4e-3,1e-3,1e-3),(2e-3,2e-3,.25e-3)):
        reference = {str(n):duffy_self_H(dims,0,n) for n in (32,64,128)}
        values={}
        for n in (16,32,64,128):
            builder=peec_matrices.PEECBuilder(volume_order=n)
            builder.add_segment(np.zeros(3),np.array([1.,0.,0.]),
                                dims[0],dims[2],dims[1],5.8e7)
            value=float(builder.inductance_entry(0,0))
            values[str(n)]=dict(value_H=value,
                relative_error=abs(value-reference['128'])/abs(reference['128']))
        rows.append(dict(dimensions_m=dims,current_axis=0,
                         independent_reference_H=reference,production=values))
    result=dict(status='DIAGNOSTIC_NOT_ADOPTION',runtime_s=time.perf_counter()-started,
                radia_path=radia.__file__,native_path=str(module),
                native_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),rows=rows)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8',newline='\n')


if __name__=='__main__':main()
