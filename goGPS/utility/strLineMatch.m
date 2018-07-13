function idx = strLineMatch(mat,line)
%DESCRIPTION: given a matrix of char (mat) return the index of the lines
%equal to line
    n_idx = size(line,1);
    if n_idx == 1
        idx= sum(mat == repmat(line,size(mat,1),1) ,2) == size(mat,2); % to be backward compatible return logical array
    else
        idx = zeros(n_idx,1);
        for i = 1: length(idx)
            idx_t = find(sum(mat == repmat(line(i,:),size(mat,1),1) ,2) == size(mat,2));
            if ~isempty(idx_t)
                idx(i) = idx_t;
            end
        end
    end
end
