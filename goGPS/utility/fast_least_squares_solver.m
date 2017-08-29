function [x, Cxx, sigma02_hat, v_hat, Cyy] = fast_least_squares_solver(y0, b, A, Q)

[n, m] = size(A);

%least-squares solution
P = A' / Q;
N = P * A;
Y = (y0-b);
L = P * Y;
x = N\L;

%estimation of the variance of the observation error
y_hat = A * x + b;
v_hat = y0 - y_hat;
V = v_hat';
T = Q \ v_hat;
sigma02_hat = (V * T) / (n - m);

%covariance matrix
try
    Cxx = sigma02_hat * cholinv(N);
catch
    Cxx = sigma02_hat * N^-1;
end

if nargout == 5
    Cyy = sigma02_hat * A/N*A';
end
