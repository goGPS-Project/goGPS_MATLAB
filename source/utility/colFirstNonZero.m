function B = colFirstNonZero(A)
    % Get first non zero element of each column of a matrix
    % SYNTAX:
    % B = colFirstNonZero(A)
    [~, c] = max( A ~=0, [], 1 );
    d = (0 : (size(A, 2) - 1)) * size(A, 1) + c;
    B = A(d);
end
