function A = averageDemagTensor(sourceMin, sourceMax, targetMin, targetMax)
%AVERAGEDEMAGTENSOR Closed-form cuboid average demagnetizing tensor.
A = radia.internal.callMex('radia.AverageDemagTensor', double(sourceMin), ...
    double(sourceMax), double(targetMin), double(targetMax));
end
