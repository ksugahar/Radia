function result=linearizeLaplacePairGram(points,weights,displacementModes,relativeWeightDerivatives)
%LINEARIZELAPLACEPAIRGRAM Analytic off-diagonal Laplace pair shape tangent.
arguments
 points (:,:) double; weights (:,1) double {mustBePositive}
 displacementModes (:,:,:) double
 relativeWeightDerivatives double=double.empty
end
n=size(points,1); dim=size(points,2); q=size(displacementModes,1);
if ~isequal(size(displacementModes),[q,n,dim]), error("radia:topopt:Shape","Displacement modes must be q-by-n-by-d."); end
if isempty(relativeWeightDerivatives), relativeWeightDerivatives=zeros(q,n); end
if ~isequal(size(relativeWeightDerivatives),[q,n]), error("radia:topopt:Shape","Weight derivative shape mismatch."); end
gram=zeros(n); derivative=zeros(q,n,n);
for i=1:n
 for j=i+1:n
  delta=points(i,:)-points(j,:); distance=norm(delta);
  if distance==0, error("radia:topopt:Points","Distinct sample points are required."); end
  value=weights(i)*weights(j)/(4*pi*distance); gram(i,j)=value; gram(j,i)=value;
  for k=1:q
   velocity=squeeze(displacementModes(k,i,:)-displacementModes(k,j,:))';
   tangent=value*(relativeWeightDerivatives(k,i)+relativeWeightDerivatives(k,j)-dot(delta,velocity)/distance^2);
   derivative(k,i,j)=tangent; derivative(k,j,i)=tangent;
  end
 end
end
result=struct("schema","radia.topopt.laplace-pair-linearization/v1","gram",gram,"jacobian",derivative, ...
 "self_term_policy","analytic self-panel value and derivative must be supplied separately");
end
