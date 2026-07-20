from copy import deepcopy
import math

from radia_mcp.radia_ngsolve.wave_energy_identity_v55 import ANTENNA, RESONATOR, validate_public_v55_identity


CASE_IDS = {"v55_public_resonator_loaded_unloaded_q_coupling_linewidth_energy_owner_mismatch", "v55_public_antenna_efficiency_accepted_radiated_loss_gain_directivity_owner_mismatch"}


def _payload():
    gen=lambda name,fields:{"generation":name,**{field:name for field in fields}}
    f0=2.4e9; ql=1200.0; q0=3000.0; qe=2000.0; power=1.0
    resonator={**gen("resonator-v55",("q_generation","coupling_generation","linewidth_generation","energy_generation","owner_generation","result_generation")),"resonance_frequency_hz":f0,"result_resonance_frequency_hz":f0,"loaded_q":ql,"result_loaded_q":ql,"unloaded_q":q0,"result_unloaded_q":q0,"external_q":qe,"result_external_q":qe,"coupling_beta":q0/qe,"result_coupling_beta":q0/qe,"linewidth_hz":f0/ql,"result_linewidth_hz":f0/ql,"dissipated_power_w":power,"result_dissipated_power_w":power,"stored_energy_j":ql*power/(2*math.pi*f0),"result_stored_energy_j":ql*power/(2*math.pi*f0),"monitor_owner":"monitor:v55","result_monitor_owner":"monitor:v55","result_sha256":"d"*64,"accepted_result_sha256":"d"*64}
    accepted=1.0;radiated=.75;loss=.25;eta=.75;directivity=6.0;gain=directivity+10*math.log10(eta)
    antenna={**gen("antenna-v55",("power_generation","efficiency_generation","gain_generation","directivity_generation","owner_generation","result_generation")),"accepted_power_w":accepted,"result_accepted_power_w":accepted,"radiated_power_w":radiated,"result_radiated_power_w":radiated,"loss_power_w":loss,"result_loss_power_w":loss,"radiation_efficiency":eta,"result_radiation_efficiency":eta,"directivity_dbi":directivity,"result_directivity_dbi":directivity,"gain_dbi":gain,"result_gain_dbi":gain,"farfield_owner":"farfield:v55","result_farfield_owner":"farfield:v55","result_sha256":"e"*64,"accepted_result_sha256":"e"*64}
    return {"runs":[{RESONATOR:resonator,ANTENNA:antenna}]}


def test_v55_positive_public_artifacts_are_accepted(): assert all(validate_public_v55_identity(_payload()).values())
def test_v55_frozen_counterfactuals_are_rejected():
    p=deepcopy(_payload());p["runs"][0][RESONATOR]["result_monitor_owner"]="monitor:stale";p["runs"][0][ANTENNA]["result_farfield_owner"]="farfield:stale";assert not all(validate_public_v55_identity(p).values())
def test_v55_self_consistent_nonphysical_artifacts_are_rejected():
    p=deepcopy(_payload());p["runs"][0][RESONATOR]["stored_energy_j"]=p["runs"][0][RESONATOR]["result_stored_energy_j"]=-1.0;p["runs"][0][ANTENNA]["radiation_efficiency"]=p["runs"][0][ANTENNA]["result_radiation_efficiency"]=1.5;assert not all(validate_public_v55_identity(p).values())
def test_v55_malformed_values_reject_without_raising():
    p=deepcopy(_payload());p["runs"][0][RESONATOR]["loaded_q"]=[1200.0];p["runs"][0][ANTENNA]["accepted_power_w"]=[1.0];assert not all(validate_public_v55_identity(p).values())


def test_v55_numeric_digests_are_rejected():
    p = deepcopy(_payload())
    numeric_digest = int("1" * 64)
    for identity in (RESONATOR, ANTENNA):
        p["runs"][0][identity]["result_sha256"] = numeric_digest
        p["runs"][0][identity]["accepted_result_sha256"] = numeric_digest
    assert not all(validate_public_v55_identity(p).values())
