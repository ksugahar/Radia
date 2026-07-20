function gradient=magneticShieldRMSGradient(response,responseJacobian)
%MAGNETICSHIELDRMSGRADIENT RMS leakage-field gradient from analytic VIM data.
arguments
 response (:,1) double {mustBeFinite}
 responseJacobian (:,:) double {mustBeFinite}
end
if isempty(response)||size(responseJacobian,1)~=numel(response), error("radia:topopt:Shape","Response/Jacobian shape mismatch."); end
rmsValue=sqrt(mean(response.^2));
if rmsValue==0, error("radia:topopt:ZeroRMS","RMS response derivative is undefined at zero."); end
gradient=(response'*responseJacobian/(numel(response)*rmsValue))';
end
