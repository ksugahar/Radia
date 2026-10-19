function result=linearizeVIM(A,b,C,dA,options)
%LINEARIZEVIM Analytically linearize A(rho)m=b(rho), y=C(rho)m.
arguments
 A double {mustBeFinite}; b double {mustBeFinite}; C double {mustBeFinite}; dA double {mustBeFinite}
 options.db double=double.empty; options.dC double=double.empty
end
n=size(A,1); if size(A,2)~=n||numel(b)~=n, error("radia:topopt:Shape","A must be square and b must match."); end
if size(C,2)~=n||ndims(dA)~=3||size(dA,2)~=n||size(dA,3)~=n, error("radia:topopt:Shape","C or dA has incompatible shape."); end
cells=size(dA,1); db=options.db; if isempty(db), db=zeros(cells,n); end
dC=options.dC; if isempty(dC), dC=zeros(cells,size(C,1),n); end
if ~isequal(size(db),[cells,n])||~isequal(size(dC),[cells,size(C,1),n]), error("radia:topopt:Shape","db or dC has incompatible shape."); end
state=A\b(:); rhs=zeros(n,cells);
for k=1:cells, rhs(:,k)=db(k,:)'-squeeze(dA(k,:,:))*state; end
stateJacobian=(A\rhs)'; response=C*state; responseJacobian=zeros(size(C,1),cells);
for k=1:cells, responseJacobian(:,k)=C*stateJacobian(k,:)'+reshape(dC(k,:,:),size(C))*state; end
result=struct("schema","radia.topopt.vim-linearization/v1","state",state, ...
 "response",response,"state_jacobian",stateJacobian,"response_jacobian",responseJacobian);
end
