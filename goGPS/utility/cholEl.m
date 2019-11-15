function [R,idx_rm] = cholEl(A,tol)
% perform cholesky eliminating rank deficient columns 
% A must be sparse
%
% SYNTAX:
% [R,idx_rm] = cholEl(A,tol)
R = zeros(size(A));
idx_rm = [];
g = 1;
idx_a = 0;
idx_r = 1;
while g ~=0 && idx_r < (size(A,1) - length(idx_rm))
    idx_a = idx_a + g;
    [Rt,g] = chol(A(idx_a:end, idx_a:end), 'lower');
    if g == 0
        diag_el = diag(Rt);
        idx_rm_r = diag_el < tol;
        idx_rm = [idx_rm; idx_a-1+find(idx_rm_r)];
        R(idx_r:(end-length(idx_rm)),idx_r:(end-length(idx_rm))) = Rt;
    else
        diag_el = diag(Rt(1:(g - 1),:));
        idx_rm_r = diag_el < tol;
        idx_rm = [idx_rm; idx_a-1+find(idx_rm_r); idx_a-1+g];
        R(idx_r + g - 1:end-1,:) = R(idx_r + g:end,:);
        R(idx_r:(end - length(idx_rm)),idx_r + (0:(g-1-1-sum(idx_rm_r)))) = Rt(1:size(Rt,1) ~= g,:);
        A(idx_a + g:end, idx_a + g:end) = A(idx_a + g:end, idx_a + g:end) - R(idx_r + g-1-sum(idx_rm_r):(end - length(idx_rm)), idx_r + (0:(g-1-1-sum(idx_rm_r)))) * R(idx_r + g-1-sum(idx_rm_r):(end - length(idx_rm)), idx_r + (0:(g-1-1-sum(idx_rm_r))))';
    end
    idx_r = idx_r +g  -sum(idx_rm_r) -1;
end
R(:,end-length(idx_rm)+1:end) = [];
R(end-length(idx_rm)+1:end,:) = [];
end

