function [mass_pt] = bayesianAmbFixing(a,Q)
% compute multi integer weighted probability 
if true
[dQ,Z,L,D,da,iZt] = decorrel(Q,a);
else
    dQ = Q;
    da = a;
    iZt = eye(size(Q));
end
R = chol(dQ);
n_samples = 500000;
%nc = det(2*pi*Q)^(0.5);
%nc*exp(-1/2*u'*iQ*u);
ra = repmat(da,1,n_samples)+ R*randn(numel(a),n_samples);
ra(abs(fracFNI(ra)) > 0.1) = nan;
ra = round(ra);
n_valid_sml = sum(~isnan(ra),2);
mass_pt = sum(nan2zero(ra),2)./n_valid_sml;
%mass_pt = mode(ra,2);
mass_pt = iZt*(mass_pt);
end
