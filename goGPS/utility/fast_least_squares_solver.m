function [x, Cxx, s02_hat, v_hat, Cyy] = fast_least_squares_solver(y0, b, A, Q)

[n, m] = size(A);

% least-squares solution
P = A' / Q;
N = P * A;

try
    N_inv = cholinv(full(N));
catch
    N_inv = N^-1;
end

Y = (y0 - b);
L = P * Y;
x = N_inv * L;

% estimation of the variance of the observation error
y_hat = A * x + b;
v_hat = y0 - y_hat;
T = Q \ v_hat;
s02_hat = (v_hat' * T) / (n - m);

% covariance matrix
Cxx = s02_hat * N_inv;

if nargout == 5
    A = A';
    Cyy = s02_hat * A' * N_inv * A;
end
