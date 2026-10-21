"""Fast HDF5 AMG-transfer cache integrity and invalidation protection."""
import argparse,json,tempfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.sparse import diags,csr_matrix
from mixed_omega_amg_cache import save_transfers,load_transfers,symmetric_row_bound


def main():
    parser=argparse.ArgumentParser(__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    record=dict(passed=False,rejections=[])
    try:
        a=diags([-np.ones(64),4*np.ones(65),-np.ones(64)],[-1,0,1]).tocsr()
        d=1/np.sqrt(a.diagonal());expected=float(np.max(np.asarray(abs(diags(d)@a@diags(d)).sum(axis=1))))
        np.testing.assert_allclose(symmetric_row_bound(a),expected,rtol=1e-14,atol=0)
        p=csr_matrix((np.ones(65),(np.arange(65),np.arange(65)//8)),shape=(65,9))
        levels=[SimpleNamespace(A=a,P=p),SimpleNamespace(A=p.T@a@p)]
        b=np.ones((65,1));settings=dict(strength=0.,damping=1.5)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'cache.h5';save_transfers(path,levels,b,parameters=settings)
            loaded,receipt=load_transfers(path,a,b,parameters=settings)
            np.testing.assert_allclose(loaded[-1].A.toarray(),(p.T@a@p).toarray(),rtol=0,atol=1e-14)
            assert receipt['reused'] and receipt['galerkin_operators_rebuilt']
            changed=a*(1+2e-15);rebuilt,_=load_transfers(path,changed,b,parameters=settings)
            np.testing.assert_allclose(rebuilt[-1].A.toarray(),(p.T@changed@p).toarray(),rtol=0,atol=1e-14)
            reused,details=load_transfers(path,changed,b,parameters=settings,rebuild_galerkin=False)
            assert details['coarse_blocks_reused'] and details['fine_operator_is_current']
            np.testing.assert_array_equal(reused[0].A.data,changed.data)
            np.testing.assert_array_equal(reused[-1].A.toarray(),levels[-1].A.toarray())
            fine=reused[0].A;coarse=reused[-1].A.toarray()
            # Independent dense action of the symmetric two-pre/two-post cycle.
            d=1/np.sqrt(fine.diagonal());bound=np.max(np.asarray(abs(diags(d)@fine@diags(d)).sum(axis=1)))
            weights=(1.5/bound/fine.diagonal())[:,None]
            identity=np.eye(65);action=np.zeros_like(identity)
            for _ in range(2):action+=weights*(identity-fine@action)
            action+=p@np.linalg.solve(coarse,p.T@(identity-fine@action))
            for _ in range(2):action+=weights*(identity-fine@action)
            symmetry=float(np.linalg.norm(action-action.T)/np.linalg.norm(action))
            minimum=float(np.linalg.eigvalsh((action+action.T)/2).min())
            assert symmetry<1e-13 and minimum>0
            record.update(reused_cycle_symmetry_defect=symmetry,reused_cycle_minimum_eigenvalue=minimum)
            def reject(name,operator,candidate,parameters):
                try:load_transfers(path,operator,candidate,parameters=parameters)
                except ValueError:record['rejections'].append(name)
                else:raise AssertionError('Invalid cache accepted: '+name)
            reject('changed_operator',a*1.001,b,settings)
            reject('changed_structure',diags(np.ones(65)),b,settings)
            reject('changed_candidate',a,b*1.001,settings)
            reject('changed_settings',a,b,dict(strength=.1,damping=1.5))
            try:save_transfers(path,levels,b,parameters=settings)
            except FileExistsError:record['rejections'].append('overwrite')
            else:raise AssertionError('Existing cache overwritten')
            import h5py
            with h5py.File(path,'r+') as f:f['p0/val'][0]+=.001
            reject('corrupted_payload',a,b,settings)
        record.update(passed=True,fresh_galerkin_action_verified=True,roundoff_admission_checked=True)
    except Exception as error:record['failure']=repr(error);raise
    finally:
        args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':main()
