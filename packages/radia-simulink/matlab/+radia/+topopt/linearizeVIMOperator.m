function result=linearizeVIMOperator(mass,chargeMap,chargeGram,appliedCoefficients,invChi,dMass,dChargeMap,dChargeGram,dAppliedCoefficients)
%LINEARIZEVIMOPERATOR Product-rule tangent of A=invChi*M+B'*G*B, b=M*h.
arguments
 mass (:,:) double; chargeMap (:,:) double; chargeGram (:,:) double
 appliedCoefficients (:,1) double; invChi (1,1) double
 dMass (:,:,:) double; dChargeMap (:,:,:) double; dChargeGram (:,:,:) double
 dAppliedCoefficients double=double.empty
end
n=size(mass,1); q=size(dMass,1);
if size(mass,2)~=n||size(chargeMap,2)~=n||~isequal(size(chargeGram),[size(chargeMap,1),size(chargeMap,1)])||numel(appliedCoefficients)~=n
 error("radia:topopt:Shape","Incompatible VIM operator shapes.");
end
if ~isequal(size(dMass),[q,n,n])||~isequal(size(dChargeMap),[q,size(chargeMap)])||~isequal(size(dChargeGram),[q,size(chargeGram)])
 error("radia:topopt:Shape","VIM derivative shape mismatch.");
end
if isempty(dAppliedCoefficients), dAppliedCoefficients=zeros(q,n); end
if ~isequal(size(dAppliedCoefficients),[q,n]), error("radia:topopt:Shape","Applied-field derivative shape mismatch."); end
A=invChi*mass+chargeMap'*chargeGram*chargeMap; rhs=mass*appliedCoefficients;
dA=zeros(q,n,n); db=zeros(q,n);
for k=1:q
 dM=squeeze(dMass(k,:,:)); dB=squeeze(dChargeMap(k,:,:)); dG=squeeze(dChargeGram(k,:,:));
 dA(k,:,:)=invChi*dM+dB'*chargeGram*chargeMap+chargeMap'*dG*chargeMap+chargeMap'*chargeGram*dB;
 db(k,:)=dM*appliedCoefficients+mass*dAppliedCoefficients(k,:)';
end
result=struct("schema","radia.topopt.vim-operator-linearization/v1","matrix",A,"rhs",rhs, ...
 "matrix_jacobian",dA,"rhs_jacobian",db);
end
