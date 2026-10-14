"""Publishable original numeric element cases for independent RadFld parity."""
import numpy as np

REPRESENTATIONS=('tet','polynomial_tet','cubic_tet','sextic_quartic_tet','curved_tet')

def element_case(kind,size,targets):
    if kind not in REPRESENTATIONS:
        raise ValueError(kind)
    rng=np.random.default_rng(20261012+size)
    origins=rng.uniform(-0.5,0.5,(size,3))
    tet=np.array([[0.,0.,0.],[0.03,0.,0.],[0.,0.025,0.],[0.,0.,0.02]])
    tri=tet[[0,1,2]]
    vc,sc={'tet':(16,22),'polynomial_tet':(32,22),'cubic_tet':(32,29),
           'sextic_quartic_tet':(96,44),'curved_tet':(34,24)}[kind]
    v=np.zeros((size,vc));s=np.zeros((size,sc))
    spec=dict(version=1.,representation=kind)
    if kind=='curved_tet':
        edges=((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))
        nodes=np.vstack([tet]+[(tet[a]+tet[b])[None,:]/2 for a,b in edges])
        # Non-affine P2 geometry; leave vertex positions fixed.
        nodes[4:]+=np.array([0.0002,-0.00015,0.0003])
        faces=np.vstack([tri,(tri[0]+tri[1])/2,(tri[1]+tri[2])/2,(tri[2]+tri[0])/2])
        faces[3:,2]+=0.0003
        v[:,:30]=(nodes[None,:,:]+origins[:,None,:]).reshape(size,30)
        s[:,:18]=(faces[None,:,:]+origins[:,None,:]).reshape(size,18)
        v[:,30:]=rng.normal(size=(size,4));s[:,18:]=rng.normal(size=(size,6))
        gx,gw=np.polynomial.legendre.leggauss(8)
        spec.update(gauss_points=((gx+1)/2)[:,None],gauss_weights=(gw/2)[:,None])
    else:
        v[:,:12]=(tet[None,:,:]+origins[:,None,:]).reshape(size,12)
        s[:,:9]=(tri[None,:,:]+origins[:,None,:]).reshape(size,9)
        v[:,12:]=rng.normal(size=(size,vc-12));s[:,9:]=rng.normal(size=(size,sc-9))
        # Keep high-order terms material at the element length scale.
        monomial_scales=np.array([0.03**(-d) for d in range(7)
                                  for ax in range(d+1) for ay in range(d-ax+1)])
        v[:,12:]*=monomial_scales[:vc-12]
        if sc==22:
            s[:,9:]*=np.array([1.]+[0.03**-1]*3+[0.03**-2]*9)
        else:
            s[:,9:]*=monomial_scales[:sc-9]
        spec.update(volume_origins=origins.copy(),surface_origins=origins.copy())
    spec.update(volume=v,surface=s)
    obs=rng.uniform(1.5,3.,(targets,3))
    near=np.array([[.004,.005,.006],[.012,.008,.0001],[.015,.01,0.],[.032,.002,.001]])
    count=min(targets,len(near));obs[:count]=origins[0]+near[:count]
    return spec,obs
