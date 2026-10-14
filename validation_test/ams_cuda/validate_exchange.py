"""Python-to-MATLAB HDF5 acceptance, optionally with an actual >2 GB array."""
from pathlib import Path
import argparse
import json
import platform
import sys
import time
import numpy as np
import matlab.engine
import importlib.metadata as md

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'packages/radia-gpu/src'))
from radia_gpu.ams_exchange import write_exchange, file_sha256


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    record={'schema':'radia.ams-hdf5-validation.v1','accepted':False,
            'platform':platform.platform(),'a100_started':False,
            'versions':{p:md.version(p) for p in ('numpy','h5py','matlabengine')},'cases':[]}
    eng=None
    try:
        eng=matlab.engine.start_matlab('-nodesktop -nosplash')
        eng.addpath(str(ROOT/'packages/radia-gpu/matlab'),nargout=0)
        eng.addpath(str(Path(__file__).parent),nargout=0)
        record['matlab']=eng.version()
        for large in ([False,True] if args.large else [False]):
            payload={'matrix':np.arange(15,dtype=np.float64).reshape(3,5),
                     'indices':np.array([0,2,4],dtype=np.int32),'empty':np.empty(0,dtype=np.int32),
                     'snapshot':{'schema':'radia.ams-accelerator.v1','beta_zero':False,
                         'spaces':[{'levels':[{'l1':np.array([1.,2.,3.])}]}]}}
            count=268_435_457 if large else 0
            if large:
                # Disk-backed source: no extra 2 GB Python array allocation.
                source=args.out/'source.bin'
                values=np.memmap(source,dtype=np.float64,mode='w+',shape=(count,))
                for start in range(0,count,1_000_000):values[start:start+1_000_000]=1
                values.flush();payload['large']=values
            path=args.out/('large.hdf5' if large else 'small.hdf5')
            timer=time.perf_counter();receipt=write_exchange(path,payload);write_seconds=time.perf_counter()-timer
            timer=time.perf_counter();result=eng.check_exchange(str(path),float(count));read_seconds=time.perf_counter()-timer
            entry={'accepted':bool(result['accepted']),'large':large,**receipt,
                   'write_seconds':write_seconds,'matlab_read_and_check_seconds':read_seconds,
                   'sha256':file_sha256(path)}
            if large:assert receipt['file_bytes']>2**31 and count*8>2**31
            record['cases'].append(entry)
            print(json.dumps(entry),flush=True)
            if large:
                del payload['large'];del values;source.unlink()
        record['accepted']=True
    except Exception as exc:
        record['error']=repr(exc);raise
    finally:
        if eng is not None:eng.quit()
        (args.out/'result.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--large',action='store_true')
    main(parser.parse_args())
