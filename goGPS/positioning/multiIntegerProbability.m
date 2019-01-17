function [mass_pt] = multiIntegerProbability(Q,a)
% compute multi integer weighted probability
Z = eye(size(Q));
[dQ,Z,L,D,da,iZt] = decorrel(Q,a);
R = chol(dQ);
n_samples = 500000;
%nc = det(2*pi*Q)^(0.5);
%nc*exp(-1/2*u'*iQ*u);
ra = da+ R*rand(numel(a),n_samples);
ra = round(ra);
mass_pt = sum(ra,2)/n_samples;
mass_pt = iZt*(mass_pt);
end
