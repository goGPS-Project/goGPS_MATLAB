function [L,D] = ldldecom_m(Qahat)
 % find ltdl decommposition using ldlt matlab decomposition 
 %
 % NOTE: to be used in matlab
 [L,D] = ldl(flipud(fliplr(Qahat)),'upper'); L = fliplr(flipud(L));
 D = flipud(diag(D));
end