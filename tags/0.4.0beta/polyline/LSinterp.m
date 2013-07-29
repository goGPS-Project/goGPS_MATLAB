function [m,q,s2m,s2q] = LSinterp(x0,y0,s2x,s2y,sxy);

xm = mean(x0);
ym = mean(y0);

x0m = x0 - xm;
y0m = y0 - ym;

% approximate angular coefficient
if sum(x0m.^2) > 0
    m0 = sum(x0m.*y0m) / sum(x0m.^2);
else
    m0 = 0;
end

% weights
q = (s2y + m0^2 * s2x - 2*m0 * sxy).^(-1);

% normal matrix
N = [sum(q), sum(q.*x0m); sum(q.*x0m), sum(q.*(x0m.^2))];
invN = inv(N);

% known term
Y = [sum(q.*y0); sum(q.*y0.*x0m)];

% estimated parameters
par = invN * Y;

% output
m = par(2);
q = par(1) - m*xm;

s2m = invN(2,2);
s2q = invN(1,1);   