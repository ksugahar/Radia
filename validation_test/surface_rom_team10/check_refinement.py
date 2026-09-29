"""Verify same-factor residual refinement and rejection of an unusable inverse."""
import argparse,os,json,weakref
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
os.environ.update(TEAM_CASE='team10',TEAM_QUARTER='1',TEAM_COIL2='1')
import ngsolve as ng
import numpy as np
from aphi_model import ModelAphi,run_fom
from transient import waveform
live_inverses=weakref.WeakSet()
class Inverse:
    def __init__(self,actual):
        self.actual=actual
        live_inverses.add(self)
    def __mul__(self,v):return self.actual*v
class Matrix:
    def __init__(self,actual,scale):self.actual=actual;self.scale=scale
    def Inverse(self,*a,**kw):
        assert not live_inverses, 'previous factor retained during next factorization'
        return Inverse(self.scale*self.actual.Inverse(*a,**kw))
    def __mul__(self,v):return self.actual*v
class Form:
    def __init__(self,actual,scale):self.actual=actual;self.scale=scale
    @property
    def mat(self):return Matrix(self.actual.mat,self.scale)
    def __getattr__(self,name):return getattr(self.actual,name)
class EquilibriumForm(Form):
    def AssembleLinearization(self,*args,**kwargs):
        raise AssertionError('equilibrium must be checked before forming a correction')
ng.SetNumThreads(2)
with ng.TaskManager():
    m=ModelAphi(.02,order=1,dt=.002);real=m.a
    m.a=Form(real,.999)
    _,hist,info=run_fom(m,waveform('mix_add'),1)
    assert not live_inverses, 'factor retained after the solve'
    assert info['refinement_solves']>0 and info['max_relative_linear_residual']<=1e-7
    m.a=Form(real,0.)
    try:run_fom(m,waveform('mix_add'),1)
    except RuntimeError as error:
        assert 'true linear residual' in str(error)
        rejection=str(error)
    else:raise AssertionError('zero inverse was accepted')
    m.a=EquilibriumForm(real,1.)
    _,equilibrium,equilibrium_info=run_fom(m,lambda t:[0.]*len(m.fs),2)
    assert equilibrium_info['linear_solves']==0
    assert all(h['newton_converged'] and h['newton_stop_reason']=='nonlinear_residual'
               and h['relative_nonlinear_residual']==0 for h in equilibrium)
    m.a=Form(real,1.)
    snaps,small,small_info=run_fom(m,lambda t:1e-6*np.asarray(waveform('mix_add')(t)),1)
    assert np.linalg.norm(snaps[-1])>0 and small_info['linear_solves']>0
    assert small[-1]['newton_stop_reason']=='nonlinear_residual'
    assert small[-1]['relative_nonlinear_residual']<=1e-10
result=dict(ngsolve=ng.__version__,ndof=m.fes.ndof,imperfect_inverse=info,zero_inverse_rejected=rejection,
            equilibrium=equilibrium,equilibrium_info=equilibrium_info,small_source=small,small_source_info=small_info)
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result))
