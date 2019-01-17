function [mass_pt] = multiIntegerProbability(Q,a)
% compute multi integer weighted probability
Z = eye(size(Q));
[dQ,Z,L,D,da,iZt] = decorrel(Q,a);
R = chol(dQ);
iQ = inv(dQ);
mass_pt = zeros(size(da));
nc = det(2*pi*Q)^(0.5);
w_tot = zeros(size(da));
n_samples = 100000;
for i =1:n_samples
    ra = da+ R*rand(size(a));
    ra = round(ra);
    u = ra - da;
    w = 1;%nc*exp(-1/2*u'*iQ*u);
    mass_pt = mass_pt+w.*ra;
    w_tot = w_tot + w;
end
mass_pt = mass_pt./w_tot;
mass_pt = iZt*(mass_pt);
end
