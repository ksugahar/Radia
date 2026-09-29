"""Verify same-factor residual refinement and rejection of an unusable inverse."""
import argparse,os,json
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
os.environ.update(TEAM_CASE='team10',TEAM_QUARTER='1',TEAM_COIL2='1')
import ngsolve as ng
from aphi_model import ModelAphi,run_fom
from transient import waveform
class Matrix:
    def __init__(self,actual,scale):self.actual=actual;self.scale=scale
    def Inverse(self,*a,**kw):return self.scale*self.actual.Inverse(*a,**kw)
    def __mul__(self,v):return self.actual*v
class Form:
    def __init__(self,actual,scale):self.actual=actual;self.scale=scale
    @property
    def mat(self):return Matrix(self.actual.mat,self.scale)
    def __getattr__(self,name):return getattr(self.actual,name)
ng.SetNumThreads(2)
with ng.TaskManager():
    m=ModelAphi(.02,order=1,dt=.002);real=m.a
    m.a=Form(real,.999)
    _,hist,info=run_fom(m,waveform('mix_add'),1)
    assert info['refinement_solves']>0 and info['max_relative_linear_residual']<=1e-7
    m.a=Form(real,0.)
    try:run_fom(m,waveform('mix_add'),1)
    except RuntimeError as error:
        assert 'true linear residual' in str(error)
        rejection=str(error)
    else:raise AssertionError('zero inverse was accepted')
result=dict(ngsolve=ng.__version__,ndof=m.fes.ndof,imperfect_inverse=info,zero_inverse_rejected=rejection)
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result))

