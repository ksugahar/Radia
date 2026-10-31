"""Guards: explicit volume API, unsupported-model rejection, bounded dense indices and operator consistency.

Assumes: required current-source native PEEC module, SI Cartesian bars, CPU.
No analytic accuracy, adopted baseline, fin validity or GPU claim is made.
"""
import numpy as np
import pytest
from radia import peec_matrices


@pytest.mark.parametrize('order',[7,257])
def test_native_rejects_invalid_volume_quadrature(order):
    with pytest.raises(ValueError, match='quadrature order'):
        peec_matrices.PEECBuilder(volume_order=order)


def test_volume_dispatch_is_shared_by_dense_and_entry_api():
    builder=peec_matrices.PEECBuilder(volume_order=32)
    for center in ([0.,0.,0.],[.002,0.,0.]):
        builder.add_segment(np.array(center),np.array([1.,0.,0.]),.001,.001,.001)
    L,R,_,_=builder.build(include_star=False)
    entries=np.array([[builder.inductance_entry(i,j) for j in range(2)] for i in range(2)])
    assert L.dtype==np.float64 and L.shape==(2,2), (L.dtype,L.shape)
    assert R.shape==(2,), R.shape
    assert L.tobytes()==entries.tobytes(), (L,entries)
    assert np.isfinite(L).all(), L
    with pytest.raises(IndexError):builder.inductance_entry(-1,0)


@pytest.mark.parametrize('direction,section', [([1.,0.,0.],1),([2**-.5,2**-.5,0.],0)])
def test_unsupported_volume_geometry_raises_instead_of_filament_fallback(direction,section):
    builder=peec_matrices.PEECBuilder(volume_order=32)
    builder.add_segment(np.zeros(3),np.array(direction),.001,.001,.001,
                        cross_section_type=section)
    with pytest.raises(ValueError,match='Cartesian'):
        builder.inductance_entry(0,0)



def test_dense_volume_index_budget_fails_before_matrix_allocation():
    # Public fail-loud capacity contract: no large matrix or field is evaluated.
    import math
    count=math.isqrt(np.iinfo(np.int32).max)+1
    builder=peec_matrices.PEECBuilder(volume_order=8)
    center=np.zeros(3);direction=np.array([1.,0.,0.])
    for _ in range(count):
        builder.add_segment(center,direction,.001,.001,.001)
    with pytest.raises(ValueError,match='Dense PEEC exceeds signed 32-bit index budget'):
        builder.build(include_star=False)
