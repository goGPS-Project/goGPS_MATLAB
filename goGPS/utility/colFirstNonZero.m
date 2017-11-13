function B = colFirstNonZero(A)
%DESCRIPTION: get first non zero element of each column of a matrix
[~, c] = max( A ~=0, [], 1 );
d = (0:(size(A,2)-1))*size(A,1) +c;
B=A(d);
end
