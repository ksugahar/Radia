"""Behavioral checks with public synthetic geometry; no customer data."""
import json,tempfile,importlib.util
from pathlib import Path
import numpy as np
from netgen.occ import WorkPlane,Axis,Axes,Pnt,Vec
from planar_face_to_coilbuilder import face_to_coilbuilder,write_builder_module

def main(output):
    outer=WorkPlane().Circle(.15).Face();inner=WorkPlane().Circle(.12).Face();face=outer-inner
    builder,cert=face_to_coilbuilder(face,current_A=-3,turns=10,extrusion_vector_m=(0,0,.02),section_width_m=.03)
    assert builder.current==-30 and cert['closure_error_m']<1e-12 and cert['segments']==4
    import radia as rad
    obj=rad.ObjCnt(builder.to_radia())
    field=np.asarray(rad.Fld(obj,'b',[[0,0,.01]]),dtype=float).reshape(-1,3)[0]
    expected_b=4e-7*np.pi*builder.current/(2*.03)*(np.arcsinh(2*.15/.02)-np.arcsinh(2*.12/.02))
    field_error=float(np.linalg.norm(field-[0,0,expected_b])/abs(expected_b))
    assert field_error<1e-6,field_error
    expected=face.Extrude(Vec(0,0,.02));actual=builder.to_occ()
    geometry_rms=float(((expected-actual).mass+(actual-expected).mass)/expected.mass)
    assert geometry_rms<1e-8,geometry_rms
    def ring(radius):
        return WorkPlane().MoveTo(-.1,-radius).LineTo(.1,-radius).Arc(radius,180).LineTo(-.1,radius).Arc(radius,180).Close().Face()
    racetrack=ring(.07)-ring(.04)
    rb,rt=face_to_coilbuilder(racetrack,current_A=5,turns=4,extrusion_vector_m=(0,0,.015),section_width_m=.03)
    assert rt['segments']==4 and rt['closure_error_m']<1e-12 and rt['boundary_sample_max_deviation_m']<1e-12
    scaled,sc=face_to_coilbuilder(face.Scale(Pnt(0,0,0),1000),current_A=3,turns=10,extrusion_vector_m=(0,0,.02),length_scale_to_m=.001,section_width_m=.03)
    assert abs(sc['generated_analytic_volume_m3']-cert['generated_analytic_volume_m3'])<1e-14
    moved=face.Rotate(Axis(Pnt(0,0,0),Vec(1,0,0)),90).Move(Vec(.7,-.2,.5))
    rotated,rc=face_to_coilbuilder(moved,current_A=1,turns=2,extrusion_vector_m=(0,-.02,0),section_width_m=.03)
    assert rc['closure_error_m']<1e-12
    for bad,kwargs in ((outer,{}),(face,{'section_width_m':.04}),(face,{'extrusion_vector_m':(.02,0,0)})):
        params={'current_A':1,'turns':1,'extrusion_vector_m':(0,0,.02)};params.update(kwargs)
        try:face_to_coilbuilder(bad,**params)
        except ValueError:pass
        else:raise AssertionError('unsupported geometry/contract accepted')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    modulepath=output/'synthetic_build_coil.py';write_builder_module(builder,modulepath)
    spec=importlib.util.spec_from_file_location('synthetic_build_coil',modulepath);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    reloaded=mod.build_coil(ampere_turns_A=45)
    assert reloaded.current==45 and np.linalg.norm(reloaded.segments[-1].end_pos-reloaded.segments[0].start_pos)<1e-12
    result={'passed':True,'checks':['circular ribbon','line-and-arc racetrack','explicit CAD unit conversion','arbitrary plane and translation','SI excitation and signed turns','exact OCC volume reconstruction','unsupported face/width/extrusion rejected','portable build_coil module reload','analytic finite-section coil-centre field and polarity'],'relative_symmetric_difference_volume':geometry_rms,'analytic_centre_field_relative_error':field_error,'circle_certificate':cert,'rotated_certificate':rc}
    (output/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);main(ap.parse_args().output)

