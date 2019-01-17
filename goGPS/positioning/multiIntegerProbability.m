function [mass_pt] = multiIntegerProbability(Q,a)
% compute multi integer weighted probability
Z = eye(size(Q));
[dQ,Z,L,D,da,iZt] = decorrel(Q,a);
R = chol(dQ);
iQ = inv(dQ);
sd = sqrt(diag(dQ));
mass_pt = zeros(size(da));
nc = det(2*pi*Q)^(0.5);
w_tot = zeros(size(da));
for i =1:10000000
    ra = a+ R*rand(size(a));
    ra = round(ra);
    w = nc*exp(-1/2*ra'*iQ*ra);
    mass_pt = mass_pt+w.*ra;
    w_tot = w_tot + w;
end
mass_pt = mass_pt./w_tot;
mass_pt = iZt*(mass_pt);
end
