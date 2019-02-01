function [bie, sols, p_sols] = BIE(a,Q,n_cand)
if nargin < 3
    n_cand = 8;
end
% compute best integer equivariant using fir 8 terms
if exist('decorrel') > 0
    [dQ,Z,L,D,da,iZt] = decorrel(Q,a);
else
    dQ = Q;
    da = a;
    iZt = eye(size(Q));
end
[zfixed,sqnorm] = ssearch(da,L,D,8);
p_sols = exp(-1/2*sqnorm);
p_sols = p_sols./sum(p_sols);
bie = sum(zfixed.*repmat(p_sols,size(zfixed,1),1),2);
bie = iZt*bie;
if nargout > 2
    sols = iZt*zfixed;
end
end
