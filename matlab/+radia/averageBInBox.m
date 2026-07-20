function B = averageBInBox(M, sourceMin, sourceMax, targetMin, targetMax)
%AVERAGEBINBOX Closed-form average B over an axis-aligned target cuboid.
B = radia.internal.callMex('radia.AverageBInBox', double(M), ...
    double(sourceMin), double(sourceMax), double(targetMin), double(targetMax));
end
