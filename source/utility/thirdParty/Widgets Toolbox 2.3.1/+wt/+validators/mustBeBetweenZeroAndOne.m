function mustBeBetweenZeroAndOne(A)
% Validate that value between 0 and 1 (inclusive)
%
% Syntax:
%     properties
%         PropertyName dataType {wt.validators.mustBeBetweenZeroAndOne}
%     end

% Copyright 2020-2021 The MathWorks, Inc.

wt.validators.mustBeBetween(A,0,1);