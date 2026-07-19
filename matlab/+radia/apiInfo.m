function info = apiInfo()
%APIINFO Return MEX API, handle-registry, and TaskManager information.

info = radia.internal.callMex('api.info');
end
