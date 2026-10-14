function check_actor_critic_gate(seed)
r=radia_mcp_matlab.actor_critic_gate(seed); disp(jsonencode(struct("tool","matlab_actor_critic_gate","ok",r.ok,"result",r)));
if ~r.ok, error("radia_mcp_matlab:ActorCriticGateFailed","Actor critic gate failed."); end
end
