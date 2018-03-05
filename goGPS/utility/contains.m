function [out] = contains(str, pattern)
% implementation fo constains function, not presents in all matlab
out = ~isempty(strfind(str,pattern));
end