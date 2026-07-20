function strategy = getClusterStrategy()
%GETCLUSTERSTRATEGY Return the active HACApK cluster strategy.
strategy = radia.internal.callMex('radia.GetClusterStrategy');
end
