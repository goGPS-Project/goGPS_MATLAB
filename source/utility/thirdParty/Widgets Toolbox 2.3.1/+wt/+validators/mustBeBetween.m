function mustBeBetween(A,minA,maxA)
% Validate that value between a min and max (inclusive)
%
% Syntax:
%     properties
%         PropertyName dataType {wt.validators.mustBeBetween(PropertyName,minVal,maxVal)}
%     end
%

% Copyright 2020-2021 The MathWorks, Inc.

if any(A(:) > maxA | A(:) < minA)
    error('wt:validators:mustBeBetween',...
        'Value must be between %g and %g.', minA, maxA);
end