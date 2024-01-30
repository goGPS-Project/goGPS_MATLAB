function mustBeBetweenZeroAnd100(A)
% Validate that value between 0 and 100 (inclusive)
%
% Syntax:
%     properties
%         PropertyName dataType {wt.validators.mustBeBetweenZeroAnd100}
%     end

% Copyright 2020-2021 The MathWorks, Inc.

wt.validators.mustBeBetween(A,0,100);