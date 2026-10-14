function check_topology_information_gate(seed)
r=radia_mcp_matlab.topology_information_gate(seed); disp(jsonencode(struct("tool","matlab_topology_information_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:TopologyInformationGateFailed","Topology information gate failed."); end
end
