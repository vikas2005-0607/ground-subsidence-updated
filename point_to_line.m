function d = point_to_line(pt,v1,v2)
%POINT_TO_LINE Perpendicular distance from point PT to the infinite line V1-V2.
% Works with 2-D or 3-D coordinates.

pt = pt(:).';
v1 = v1(:).';
v2 = v2(:).';

lineVector = v1-v2;
pointVector = pt-v2;

denom = norm(lineVector);

if denom == 0
    error('v1 and v2 must be different points.');
end

if numel(pt) == 3
    d = norm(cross(lineVector,pointVector))/denom;
elseif numel(pt) == 2
    d = abs(lineVector(1)*pointVector(2) - ...
        lineVector(2)*pointVector(1))/denom;
else
    error('Only 2-D and 3-D points are supported.');
end
end