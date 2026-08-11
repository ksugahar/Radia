import numpy as np

from radia._topopt_graph import (
    best_admissible_singleton,
    binary_graph_interface_energy,
    candidate_face_adjacency,
    connected_graph_front_beam,
    minimax_driving_potential,
    ngsolve_facet_measure_graph,
    terminal_l1_curvature_energy,
    update_graph_front_trust,
)


def test_best_singleton_is_global_and_independent_of_graph_seed_order():
    representatives=[611,256,401,96,566]
    predicted=[177.5,24.8,28.6,3.75,4.86]
    admissible={96,566}
    assert best_admissible_singleton(
        representatives,predicted,current_ratio=6.02,
        is_valid=lambda value:value in admissible)==96


def test_best_singleton_skips_physically_invalid_global_minimum():
    assert best_admissible_singleton(
        [4,2,9],[0.4,0.7,0.6],current_ratio=1.0,
        is_valid=lambda value:value!=4)==9


def chain(n):
    return tuple(np.asarray([value for value in (i-1,i+1)
                             if 0<=value<n],dtype=np.int64)
                 for i in range(n))


def test_minimax_drive_uses_raw_response_columns():
    current=np.array([2.0,-2.0]);target=np.zeros(2);band=np.ones(2)
    delta=np.array([[-1.0,1.0],[1.0,-1.0]])
    drive,subgradient=minimax_driving_potential(
        current,target,band,delta)
    np.testing.assert_allclose(subgradient,[0.5,-0.5])
    np.testing.assert_allclose(drive,[1.0,-1.0])
    np.testing.assert_array_equal(delta,[[-1.0,1.0],[1.0,-1.0]])


def test_candidate_graph_excludes_cumulative_family_alternatives():
    elements=chain(5)
    adjacency=candidate_face_adjacency(
        ([0],[1],[2,3],[3,4]),elements,
        exclusion_groups=np.array([-1,-1,7,7]))
    assert 1 in adjacency[0]
    assert 0 in adjacency[1]
    assert 3 not in adjacency[2]
    assert 2 not in adjacency[3]


def test_binary_cut_and_terminal_curvature_prefer_coherent_front():
    adjacency=chain(6)
    compact=np.array([1,1,1,0,0,0],dtype=bool)
    checker=np.array([1,0,1,0,1,0],dtype=bool)
    assert binary_graph_interface_energy(compact,adjacency)==1.0
    assert binary_graph_interface_energy(checker,adjacency)==5.0

    # Three radial stations at each end.  Equal terminal depths have zero
    # curvature; a one-cell spike has a positive L1 second difference.
    radial=np.repeat([0.0,1.0,2.0],4)
    longitudinal=np.tile([1.0,4.0,6.0,9.0],3)
    designable=np.ones(12,dtype=bool)
    smooth=np.ones(12,dtype=bool)
    spike=smooth.copy();spike[[4,7]]=False
    assert terminal_l1_curvature_energy(
        smooth,radial,longitudinal,designable,total_length=10.0)==0.0
    assert terminal_l1_curvature_energy(
        spike,radial,longitudinal,designable,total_length=10.0)>0.0


def test_connected_beam_cannot_teleport_except_from_a_second_seed():
    adjacency=chain(5)
    current=np.array([4.0]);target=np.array([0.0]);band=np.array([1.0])
    # Nodes 0 and 4 are individually useful.  With one component, the beam
    # must traverse nodes 1..3 before combining them.  With two seeded
    # components, the two endpoints may collaborate immediately.
    delta=np.array([[-1.5,0.2,0.2,0.2,-1.5]])
    one=connected_graph_front_beam(
        current_response=current,response_target=target,response_band=band,
        candidate_response_delta=delta,adjacency=adjacency,
        seed_indices=np.array([0,4]),maximum_size=2,
        maximum_components=1,beam_width=32,proposal_limit=20)
    assert not any(set(value.candidate_indices)=={0,4} for value in one)
    two=connected_graph_front_beam(
        current_response=current,response_target=target,response_band=band,
        candidate_response_delta=delta,adjacency=adjacency,
        seed_indices=np.array([0,4]),maximum_size=2,
        maximum_components=2,beam_width=32,proposal_limit=20)
    assert any(set(value.candidate_indices)=={0,4} for value in two)


def test_connected_beam_respects_exclusion_and_raw_response_prediction():
    adjacency=(np.array([1,2]),np.array([0,2]),np.array([0,1]))
    delta=np.array([[-1.0,-2.0,-1.5],[0.0,0.0,0.0]])
    proposals=connected_graph_front_beam(
        current_response=np.array([4.0,0.0]),response_target=np.zeros(2),
        response_band=np.ones(2),candidate_response_delta=delta,
        adjacency=adjacency,seed_indices=np.arange(3),maximum_size=3,
        maximum_components=1,beam_width=32,proposal_limit=20,
        exclusion_groups=np.array([3,3,-1]))
    assert all(not {0,1}.issubset(set(value.candidate_indices))
               for value in proposals)
    best=proposals[0]
    np.testing.assert_allclose(
        best.predicted_response,
        np.array([4.0,0.0])+delta[:,best.candidate_indices].sum(axis=1))


def test_trust_update_shrinks_holds_and_expands_from_exact_agreement():
    shrink=update_graph_front_trust(
        budget=8,minimum_budget=1,maximum_budget=20,
        current_ratio=10.0,predicted_ratio=5.0,actual_ratio=10.0,
        selected_size=8,interface_weight=0.1)
    assert (shrink.action,shrink.budget_after)==("shrink",4)
    assert shrink.interface_weight_after>0.1
    hold=update_graph_front_trust(
        budget=8,minimum_budget=1,maximum_budget=20,
        current_ratio=10.0,predicted_ratio=8.0,actual_ratio=8.8,
        selected_size=4,interface_weight=0.1)
    assert (hold.action,hold.budget_after)==("hold",8)
    expand=update_graph_front_trust(
        budget=8,minimum_budget=1,maximum_budget=20,
        current_ratio=10.0,predicted_ratio=8.0,actual_ratio=8.1,
        selected_size=8,interface_weight=0.1)
    assert expand.action=="expand"
    assert expand.budget_after>8


def test_ngsolve_facet_measure_graph_uses_physical_hex_areas():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    mesh=MakeStructured3DMesh(hexes=True,nx=2,ny=1,nz=1,
                              mapping=lambda x,y,z:(2.0*x,3.0*y,4.0*z))
    with ng.TaskManager():
        adjacency,exterior,weights=ngsolve_facet_measure_graph(mesh)
    assert tuple(map(tuple,adjacency))==((1,),(0,))
    # Shared x-normal face: 3*4 = 12 square units.
    np.testing.assert_allclose(weights[(0,1)],12.0,rtol=1e-12,atol=1e-12)
    # Each unit-x cell has external faces 12 + 2*(1*4) + 2*(1*3) = 26.
    np.testing.assert_allclose(exterior,[26.0,26.0],rtol=1e-12,atol=1e-12)
