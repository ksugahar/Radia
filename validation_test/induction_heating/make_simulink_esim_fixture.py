"""Self-authored volume/coil fixture and independent strong ESIM CLI reference."""
from pathlib import Path
import argparse
import json
import numpy as np
import ngsolve as ng
from netgen.occ import Box, Cylinder, OCCGeometry, Pnt, X, Z
from radia.panels.calc_inductance import build_argparser, run_inductance


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--genus',type=int,choices=(0,1),default=0)
    parser.add_argument('--backend',choices=('intree-dense','hacapk'),default='intree-dense')
    args=parser.parse_args(); root=args.directory.resolve(); root.mkdir(parents=True,exist_ok=True)
    if args.genus == 0:
        shape=Box(Pnt(0,0,0),Pnt(.01,.01,.01))
        maxh=.012
    else:
        shape=Cylinder(Pnt(0,0,-.002),Z,.01,.004)-Cylinder(Pnt(0,0,-.003),Z,.006,.006)
        maxh=.006
    shape.mat('workpiece'); shape.faces.name='sibc'
    OCCGeometry(shape).GenerateMesh(maxh=maxh).Save(str(root/'workpiece.vol'))
    coil=(Box(Pnt(.025,0,.1),Pnt(.035,.002,.102)) if args.genus == 0 else Box(Pnt(.025,.1,0),Pnt(.035,.102,.002))); coil.mat('coil'); coil.faces.name='body'
    coil.faces.Min(X).name='source'; coil.faces.Max(X).name='sink'
    OCCGeometry(coil).GenerateMesh(maxh=.01).Save(str(root/'coil.vol'))
    np.savetxt(root/'bh.txt',np.c_[np.array([0,1,10,100,1000,10000,100000]),4e-7*np.pi*100*np.array([0,1,10,100,1000,10000,100000])])
    argv=['--coil-solver','bem-a','--vol',str(root/'workpiece.vol'),'--coil-vol',str(root/'coil.vol'),
        '--wp-label','sibc','--sigma','5e6','--mu-r','100','--frequency','7000','--current','2',
        '--coil-sigma','5.8e7','--coil-source-name','source','--coil-sink-name','sink',
        '--coupling-mode','strong','--wp-bem-backend',args.backend,'--wp-loop-dof','auto',
        '--h1-order','1','--impedance-model','esim','--esim-per-panel','--bh-file',str(root/'bh.txt'),
        '--esim-panel-evaluator','table','--esim-tol','1e-3','--esim-max-iter','15',
        '--half-thickness','.01','--msh-output',str(root/'reference.msh')]
    ng.SetNumThreads(1)
    with ng.TaskManager():
        result=run_inductance(build_argparser().parse_args(argv))
    if result.get('error'): raise RuntimeError(result['error'])
    assert result['esim_converged'] and result['esim_panel_evaluation']['final_certification']=='direct-all-panels'
    (root/'reference.json').write_text(json.dumps(result,allow_nan=False),encoding='utf-8')
    return 0


if __name__=='__main__': raise SystemExit(main())
