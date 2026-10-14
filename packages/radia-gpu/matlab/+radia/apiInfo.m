function info = apiInfo()
%APIINFO Return legacy MEX API and handle-registry information.

info = radia.internal.callMex('api.info');
end
