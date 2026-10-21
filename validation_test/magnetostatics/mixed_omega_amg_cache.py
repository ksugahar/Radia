"""Explicit, checked HDF5 transfer reuse for the resident AMG investigation."""
import hashlib
import json
from types import SimpleNamespace
import numpy as np
from scipy.sparse import csr_matrix

SCHEMA = 'mixed-omega.amg-transfer-cache.v3'
ROUND_OFF_LIMIT = 1e-13


def _arrays(matrix):
    matrix = matrix.tocsr(copy=True)
    matrix.sum_duplicates(); matrix.eliminate_zeros(); matrix.sort_indices()
    if max(matrix.shape)>np.iinfo(np.int32).max or matrix.nnz>np.iinfo(np.int32).max:raise OverflowError('Cached CSR int32 capacity')
    return dict(row=np.asarray(matrix.indptr,dtype=np.int32),
                col=np.asarray(matrix.indices,dtype=np.int32),
                val=np.asarray(matrix.data,dtype=np.float64),
                shape=np.asarray(matrix.shape,dtype=np.int64))


def symmetric_row_bound(operator):
    """Bound rho(D^-1/2 A D^-1/2) without two scaled CSR copies."""
    a=operator.tocsr(copy=False);diagonal=a.diagonal()
    if np.any(diagonal<=0) or not np.isfinite(diagonal).all() or np.any(np.diff(a.indptr)<=0):raise ValueError('AMG bound requires nonempty rows and positive diagonal')
    inverse_sqrt=1/np.sqrt(diagonal)
    terms=np.abs(a.data);terms*=inverse_sqrt[a.indices]
    bound=float(np.max(inverse_sqrt*np.add.reduceat(terms,a.indptr[:-1])))
    if not np.isfinite(bound) or bound<=0:raise ValueError('AMG normalized row bound must be finite and positive')
    return bound


def _write(group, values, digest):
    for name, array in values.items():
        group.create_dataset(name,data=array)
        digest.update(np.ascontiguousarray(array).view(np.uint8))


def _read(group, names, digest):
    values={name:group[name][:] for name in names}
    for array in values.values():digest.update(np.ascontiguousarray(array).view(np.uint8))
    return values


def _matrix(values):
    if values['shape'].shape!=(2,) or np.any(values['shape']<1):raise ValueError('Cached CSR shape')
    n,m=map(int,values['shape']);row,col,val=(values[k] for k in ('row','col','val'))
    if row.shape!=(n+1,) or row[0]!=0 or row[-1]!=len(val) or len(col)!=len(val):raise ValueError('Cached CSR dimensions')
    if np.any(np.diff(row)<0) or np.any(col<0) or np.any(col>=m) or not np.isfinite(val).all():raise ValueError('Cached CSR entries')
    return csr_matrix((val,col,row),shape=(n,m))


def save_transfers(path, hierarchy, candidate, *, parameters):
    import h5py
    header=json.dumps(dict(schema=SCHEMA,parameters=parameters,levels=len(hierarchy)),sort_keys=True)
    digest=hashlib.sha256(header.encode())
    with h5py.File(path,'x') as f:
        f.attrs['header']=header
        _write(f.create_group('operator'),_arrays(hierarchy[0].A),digest)
        _write(f.create_group('candidate'),dict(values=np.asarray(candidate,dtype=np.float64)),digest)
        for i,level in enumerate(hierarchy[:-1]):_write(f.create_group(f'p{i}'),_arrays(level.P),digest)
        for i,level in enumerate(hierarchy[1:],1):_write(f.create_group(f'a{i}'),_arrays(level.A),digest)
        f.attrs['payload_sha256']=digest.hexdigest()
    return digest.hexdigest()


def load_transfers(path, operator, candidate, *, parameters, rebuild_galerkin=True):
    """Reuse interpolation, optionally retaining coarse preconditioner blocks.

    An independently assembled identical FE form can differ in summation low
    bits. The recorded round-off guard is an admission check for interpolation
    reuse, not a relaxation of either solve or field-accuracy gates. The fine
    matrix is always current. Retained coarse blocks are preconditioning only;
    smoothers and the checked SPD coarse inverse action are recomputed.
    """
    import h5py
    with h5py.File(path,'r') as f:
        header=f.attrs['header'];metadata=json.loads(header)
        if metadata.get('schema')!=SCHEMA or metadata.get('parameters')!=parameters:raise ValueError('AMG cache settings mismatch')
        levels=metadata['levels']
        if not 1<=levels<=12:raise ValueError('AMG cache hierarchy depth')
        digest=hashlib.sha256(header.encode())
        reference=_matrix(_read(f['operator'],('row','col','val','shape'),digest))
        b=_read(f['candidate'],('values',),digest)['values']
        transfers=[_matrix(_read(f[f'p{i}'],('row','col','val','shape'),digest)) for i in range(levels-1)]
        coarse=[_matrix(_read(f[f'a{i}'],('row','col','val','shape'),digest)) for i in range(1,levels)]
        if digest.hexdigest()!=f.attrs['payload_sha256']:raise ValueError('AMG cache payload hash mismatch')
    current=operator.tocsr(copy=True);current.sum_duplicates();current.eliminate_zeros();current.sort_indices()
    if reference.shape!=current.shape or not np.array_equal(reference.indptr,current.indptr) or not np.array_equal(reference.indices,current.indices):raise ValueError('AMG cache operator structure mismatch')
    delta=float(np.linalg.norm(reference.data-current.data)/np.linalg.norm(current.data))
    if not np.isfinite(delta) or delta>ROUND_OFF_LIMIT:raise ValueError('AMG cache operator values mismatch')
    if b.shape!=candidate.shape:raise ValueError('AMG cache candidate shape mismatch')
    bdelta=float(np.linalg.norm(b-candidate)/np.linalg.norm(candidate))
    if not np.isfinite(bdelta) or bdelta>ROUND_OFF_LIMIT:raise ValueError('AMG cache candidate values mismatch')
    hierarchy=[]
    for i,transfer in enumerate(transfers):
        if transfer.shape[0]!=current.shape[0]:raise ValueError('AMG cache transfer dimension mismatch')
        restriction=transfer.T.tocsr()
        hierarchy.append(SimpleNamespace(A=current,P=transfer,R=restriction))
        if rebuild_galerkin:
            current=(restriction@current@transfer).tocsr()
            current.sum_duplicates();current.sort_indices()
        else:
            current=coarse[i]
            if current.shape!=(transfer.shape[1],transfer.shape[1]):raise ValueError('Cached coarse operator dimensions')
    hierarchy.append(SimpleNamespace(A=current))
    return hierarchy,dict(reused=True,operator_relative_roundoff=delta,candidate_relative_roundoff=bdelta,
                          payload_sha256=digest.hexdigest(),galerkin_operators_rebuilt=rebuild_galerkin,
                          coarse_blocks_reused=not rebuild_galerkin,fine_operator_is_current=True)

NONLINEAR_SCHEMA = 'mixed-omega.nonlinear-amg-transfers.v1'


def nonlinear_layout(row, col):
    digest=hashlib.sha256()
    digest.update(np.ascontiguousarray(row,dtype=np.int64).view(np.uint8))
    digest.update(np.ascontiguousarray(col,dtype=np.int64).view(np.uint8))
    return digest.hexdigest()


def save_nonlinear_transfers(path, transfers, candidate, *, parameters, layout):
    """Cache only interpolation for a changing SPD operator, never old values."""
    import h5py
    header=json.dumps(dict(schema=NONLINEAR_SCHEMA,parameters=parameters,levels=len(transfers)+1,layout=layout),sort_keys=True)
    digest=hashlib.sha256(header.encode())
    with h5py.File(path,'x') as f:
        f.attrs['header']=header
        _write(f.create_group('candidate'),dict(values=np.asarray(candidate,dtype=np.float64)),digest)
        for i,transfer in enumerate(transfers):_write(f.create_group(f'p{i}'),_arrays(transfer),digest)
        f.attrs['payload_sha256']=digest.hexdigest()
    return digest.hexdigest()


def load_nonlinear_transfers(path, candidate, *, parameters, layout):
    """Admit matching layout/candidate; every Galerkin matrix is built fresh."""
    import h5py
    with h5py.File(path,'r') as f:
        header=f.attrs['header'];metadata=json.loads(header)
        if metadata.get('schema')!=NONLINEAR_SCHEMA or metadata.get('parameters')!=parameters or metadata.get('layout')!=layout:raise ValueError('Nonlinear AMG cache layout/settings mismatch')
        levels=metadata['levels']
        if not 1<=levels<=12:raise ValueError('Nonlinear AMG depth mismatch')
        digest=hashlib.sha256(header.encode())
        b=_read(f['candidate'],('values',),digest)['values']
        transfers=[_matrix(_read(f[f'p{i}'],('row','col','val','shape'),digest)) for i in range(levels-1)]
        if digest.hexdigest()!=f.attrs['payload_sha256']:raise ValueError('Nonlinear AMG payload hash mismatch')
    if b.shape!=candidate.shape:raise ValueError('Nonlinear AMG candidate mismatch')
    delta=float(np.linalg.norm(b-candidate)/max(np.linalg.norm(candidate),1e-300))
    if not np.isfinite(delta) or delta>ROUND_OFF_LIMIT:raise ValueError('Nonlinear AMG candidate mismatch')
    return transfers,dict(reused=True,transfers_only=True,payload_sha256=digest.hexdigest(),fine_operator_is_current=True,galerkin_operators_rebuilt=True)
