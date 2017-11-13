function B = rowFirstNonZero(A)
%DESCRIPTION: get first non zero element of each column of a matrix
[~, c] = max( A ~=0, [], 2 );
d = (c-1)*size(A,2) +[1:size(A,1)]';
B=A(d);
end
