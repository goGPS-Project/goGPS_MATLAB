function B = rowNormalize(A)
% DESCRIPTION: normalize each row of the matrix
A_n_r = sqrt(sum(A.^2,2,'omitnan'));
B = A./repmat(A_n_r,1,size(A,2));
end