function idx_l = num2LogIdx(idx_n, max_idx)
% convert a munerical index to a logical index
%
% SYNTAX
%    idx = num2LogIdx(idx, max_idx)

idx_l = false(max_idx,1);
idx_l(idx_n) = true;
end