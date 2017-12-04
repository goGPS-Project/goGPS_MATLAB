function idx = idxCharLines(mat,line)
%DESCRIPTION: given a matrix of char (mat) return the index of the lines
%equal to line
 idx = sum(mat == repmat(line,size(mat,1),1) ,2) == size(mat,2);
end