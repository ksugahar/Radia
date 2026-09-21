import pytest
from radia.coil_builder import CoilBuilder

pytest.importorskip('netgen.occ')


def test_touching_sections_are_not_positive_volume_overlap():
    coil=CoilBuilder(100).set_cross_section(2,2).add_straight(10).add_straight(10)
    result=coil.audit_segment_volume_overlaps()
    assert result['passed']
    assert result['scope']=='pairwise_segment_volume_overlap'


def test_overlapping_segments_report_physical_volume():
    coil=CoilBuilder(100).set_cross_section(2,2).add_straight(10)
    coil.set_start([0,5,0]).add_straight(10)
    result=coil.audit_segment_volume_overlaps()
    assert not result['passed']
    assert result['overlaps'][0]['segments']==[0,1]
    assert result['overlaps'][0]['overlap_volume']==pytest.approx(20)


@pytest.mark.parametrize('tol',[-1,float('nan'),1])
def test_invalid_overlap_tolerance_is_rejected(tol):
    with pytest.raises(ValueError):
        CoilBuilder(100).audit_segment_volume_overlaps(tol)
