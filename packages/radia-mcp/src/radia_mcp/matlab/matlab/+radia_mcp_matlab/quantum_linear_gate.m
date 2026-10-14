function report = quantum_linear_gate(seed)
arguments
    seed (1,1) double {mustBeInteger,mustBeNonnegative} = 107
end
rng(seed,"twister"); h=[1 1;1 -1]/sqrt(2); x=[0 1;1 0]; i2=eye(2);
cnot=[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0];
unitarity=[norm(h'*h-i2,"fro"),norm(cnot'*cnot-eye(4),"fro")];
psi=randn(4,1)+1i*randn(4,1); psi=psi/norm(psi); normError=abs(norm(cnot*psi)-1);
% One-bit Deutsch test: constant oracle I, balanced oracle X.
zero=[1;0]; one=[0;1]; start=kron(zero,one); prep=kron(h,h)*start;
constant=kron(i2,i2); balanced=[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0];
outConst=kron(h,i2)*constant*prep; outBal=kron(h,i2)*balanced*prep;
pFirstOneConst=sum(abs(outConst(3:4)).^2); pFirstOneBal=sum(abs(outBal(3:4)).^2);
bad=[1 1;0 1]; badRejected=norm(bad'*bad-i2,"fro")>.1;
loadCost=8; quantumCompute=1; readCost=8; endToEnd=loadCost+quantumCompute+readCost;
checks=struct("gates_unitary",max(unitarity)<1e-12,"state_norm_preserved",normError<1e-12, ...
    "deutsch_constant_zero",pFirstOneConst<1e-12,"deutsch_balanced_one",abs(pFirstOneBal-1)<1e-12, ...
    "nonunitary_rejected",badRejected,"io_cost_accounted",endToEnd>quantumCompute);
ids=["pFo1xs9PVjI","QopR9mjifN4","y-CheJv41Js","4FL-ZlC2ais","zjZS2nceAAY"];
topics=["quantum_ecosystem_io","deutsch_jozsa","two_qubit_unitary_function","quantum_hadamard_gate","qubit_probability"];
lessons=repmat(struct("video_id","","topic","","public_url",""),numel(ids),1);
for k=1:numel(ids), lessons(k)=struct("video_id",ids(k),"topic",topics(k),"public_url","https://www.youtube.com/watch?v="+ids(k)); end
report=struct("schema","radia-mcp.matlab.quantum-linear.v1","lessons",lessons, ...
    "seed",seed,"units",struct("amplitude","1","cost","arbitrary_time_unit"), ...
    "unitarity_errors",unitarity,"state_norm_error",normError, ...
    "deutsch_constant_p1",pFirstOneConst,"deutsch_balanced_p1",pFirstOneBal, ...
    "end_to_end_cost",endToEnd,"quantum_result_is_cae_ground_truth",false, ...
    "promotion_requires_forward_solver",true,"checks",checks,"ok",all(structfun(@logical,checks)));
end
