function [ minv, idx] = minNoNan(A)
% take the minimun exscluding nan values
A(isnan(A)) = Inf;
[ minv, idx] = min(A);
minv(minv == Inf) = nan;
end