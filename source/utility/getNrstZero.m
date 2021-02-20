function val = getNrstZero(vals)
% get the value nearest to zero
[val, idx] = minNoNan(abs(vals));
val = val .*sign(vals(idx + ([1:size(vals,2)]-1)*size(vals,1)));
end