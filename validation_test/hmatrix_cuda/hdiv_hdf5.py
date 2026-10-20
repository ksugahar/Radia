"""Typed HDiv validation interchange; independent CPU only, never a GPU adapter."""
import hashlib
import json
from pathlib import Path
import re
import h5py
import numpy as np
from scipy import sparse

SCHEMA = 'radia.hdiv-hdf5.v1'
TYPES = {'float64': '<f8', 'float32': '<f4', 'int64': '<i8', 'uint64': '<u8',
         'int32': '<i4', 'uint32': '<u4', 'uint8': '|u1', 'bool': '|u1'}


def _path(path):
    path = Path(path)
    if path.suffix.lower() not in ('.h5', '.hdf5'):
        raise ValueError('HDF5 only; MAT and MAT-v7.3 are not supported')
    return path


def _write_component(group, name, value, dtype):
    value = np.asarray(value, dtype=np.dtype(dtype))
    ds = group.create_dataset(name, data=value)
    ds.attrs['sha256'] = hashlib.sha256(value.tobytes(order='C')).hexdigest()


def _write_values(group, value):
    value = np.asarray(value)
    if value.ndim == 0:
        value = value.reshape(1)
    if value.ndim not in (1, 2) or not value.size:
        raise ValueError('Only nonempty rank-one/two arrays are supported')
    group.attrs['shape_json'] = json.dumps(list(value.shape))
    if value.dtype.kind == 'c':
        if value.dtype.itemsize not in (8, 16):
            raise ValueError('Unsupported complex precision')
        dtype = 'float32' if value.dtype.itemsize == 8 else 'float64'
        group.attrs['dtype'] = dtype
        group.attrs['representation'] = 'complex-split'
        _write_component(group, 'real', value.real, TYPES[dtype])
        _write_component(group, 'imag', value.imag, TYPES[dtype])
    else:
        dtype = value.dtype.name
        if dtype not in TYPES:
            raise ValueError('Unsupported numerical dtype: ' + dtype)
        group.attrs['dtype'] = dtype
        group.attrs['representation'] = 'real'
        _write_component(group, 'real', value, TYPES[dtype])


def write_hdiv_h5(path, fields, metadata):
    path = _path(path)
    if (metadata.get('units') != 'SI' or not isinstance(metadata.get('labels'), list)
            or any(not isinstance(label, str) for label in metadata['labels'])):
        raise ValueError('Explicit SI units and label list required')
    if not fields or any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', n) for n in fields):
        raise ValueError('Nonempty named numerical fields required')
    text = json.dumps(metadata, ensure_ascii=False, allow_nan=False, sort_keys=True)
    with h5py.File(path, 'x') as f:
        f.attrs.update(schema=SCHEMA, units='SI', index_base=np.uint8(0),
                       layout='file_axes_row_column_C_order', text_encoding='UTF-8',
                       metadata_json=text, field_names_json=json.dumps(sorted(fields)))
        for name, value in fields.items():
            group = f.create_group('fields/' + name)
            if sparse.issparse(value):
                value = value.tocsr(copy=True)
                if not value.has_canonical_format or not value.nnz:
                    raise ValueError('Nonempty canonical CSR required; no implicit coalescing')
                group.attrs['kind'] = 'csr'
                group.attrs['matrix_shape_json'] = json.dumps(list(value.shape))
                _write_values(group.create_group('values'), value.data)
                _write_component(group, 'indptr', value.indptr, '<u8')
                _write_component(group, 'indices', value.indices, '<u8')
            else:
                group.attrs['kind'] = 'dense'
                _write_values(group, value)


def _component(group, name, shape, dtype, budget):
    ds = group[name]
    expected = np.dtype(dtype)
    if not isinstance(ds, h5py.Dataset) or ds.shape != tuple(shape) or ds.dtype.str != expected.str:
        raise ValueError('Dataset shape/dtype/byte-order mismatch')
    size = int(np.prod(shape, dtype=object)) * expected.itemsize
    if size > budget[0]:
        raise ValueError('Explicit host read budget exceeded')
    budget[0] -= size
    value = ds[...]
    if hashlib.sha256(value.tobytes(order='C')).hexdigest() != ds.attrs.get('sha256'):
        raise ValueError('Component SHA256 mismatch')
    return value


def _values(group, budget):
    shape = json.loads(group.attrs['shape_json'])
    if len(shape) not in (1, 2) or any(type(n) is not int or n < 1 for n in shape):
        raise ValueError('Invalid array dimensions')
    dtype = group.attrs['dtype']; representation = group.attrs['representation']
    if dtype not in TYPES or representation not in ('real', 'complex-split'):
        raise ValueError('Unsupported declared representation')
    expected = {'real'} if representation == 'real' else {'real', 'imag'}
    if set(group.keys()) != expected:
        raise ValueError('Unexpected or missing components')
    value = _component(group, 'real', shape, TYPES[dtype], budget)
    if representation == 'complex-split':
        if dtype not in ('float32', 'float64'):
            raise ValueError('Complex integer representation forbidden')
        value = value + 1j * _component(group, 'imag', shape, TYPES[dtype], budget)
    elif dtype == 'bool':
        if not np.isin(value, [0, 1]).all():
            raise ValueError('Invalid logical encoding')
        value = value.astype(bool)
    return value


def read_hdiv_h5(path, max_bytes=512 * 2**20):
    path = _path(path)
    if max_bytes < 1 or not h5py.is_hdf5(path):
        raise ValueError('Checked HDF5 and positive read budget required')
    budget = [max_bytes]
    with h5py.File(path, 'r') as f:
        if (f.attrs.get('schema') != SCHEMA or f.attrs.get('units') != 'SI'
                or f.attrs.get('index_base') != 0 or f.attrs.get('layout') != 'file_axes_row_column_C_order'
                or f.attrs.get('text_encoding') != 'UTF-8' or set(f.keys()) != {'fields'}):
            raise ValueError('HDiv HDF5 contract mismatch')
        metadata = json.loads(f.attrs['metadata_json'])
        if (metadata.get('units') != 'SI' or not isinstance(metadata.get('labels'), list)
                or any(not isinstance(label, str) for label in metadata['labels'])):
            raise ValueError('Invalid physical metadata')
        names = json.loads(f.attrs['field_names_json'])
        if len(set(names)) != len(names) or set(names) != set(f['fields'].keys()):
            raise ValueError('Field inventory mismatch')
        result = {}
        for name in names:
            g = f['fields'][name]
            if g.attrs['kind'] == 'dense':
                result[name] = _values(g, budget)
            elif g.attrs['kind'] == 'csr':
                shape = json.loads(g.attrs['matrix_shape_json'])
                if len(shape) != 2 or any(type(n) is not int or n < 1 or n >= 2**63 for n in shape) or set(g.keys()) != {'values', 'indptr', 'indices'}:
                    raise ValueError('Invalid CSR matrix dimensions/components')
                values = _values(g['values'], budget).reshape(-1)
                ptr = _component(g, 'indptr', [shape[0] + 1], '<u8', budget)
                idx = _component(g, 'indices', [values.size], '<u8', budget)
                if ptr[0] != 0 or ptr[-1] != values.size or np.any(ptr[1:] < ptr[:-1]) or np.any(idx >= shape[1]):
                    raise ValueError('Invalid zero-based CSR indices')
                matrix = sparse.csr_matrix((values, idx.astype(np.int64), ptr.astype(np.int64)), shape=shape)
                if not matrix.has_canonical_format:
                    raise ValueError('CSR must be canonical')
                result[name] = matrix
            else:
                raise ValueError('Unknown field kind')
    return result, metadata
