function strategy = hluClusterStrategy()
%HLUCLUSTERSTRATEGY Return the active HACApK cluster strategy.
strategy = radia.internal.callMex('hlu.cluster_strategy');
end
