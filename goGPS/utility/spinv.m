function X = spinv(A,tol,method)
%PINV  sparse Pseudoinverse. 
%   Copyright 1984-2013 The MathWorks, Inc. 
if nargin < 3 || strcmpi(method,'svd')
    [U,S,V] = svds(A,size(A,1));
    s = diag(S);
    if nargin < 2
        tol = max(size(A)) * eps(norm(s,inf));
    end
    r1 = sum(s > tol)+1;
    V(:,r1:end) = [];
    U(:,r1:end) = [];
    s(r1:end) = [];
    s = 1./s(:);
    X = bsxfun(@times,V,s.')*U';
elseif strcmpi(method,'qr')
    [H,R,E] = qr(A);
    dr = diag(R);
    if isempty(tol)
    tol = 1e4*eps(max(abs(dr)));
    end
    k = abs(dr)>tol;
    X = E*[inv(R(k,k))*H(:,k)'; sparse(sum(~k),size(H,2))];
end