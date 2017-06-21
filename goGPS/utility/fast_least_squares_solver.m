function [x, Cxx, sigma02_hat, v_hat] = fast_least_squares_solver(y0, b, A, Q)

[n, m] = size(A);

%least-squares solution
K = A';
P = Q\A;
N = K*P;
Y = (y0-b);
R = Q\Y;
L = K*R;
x = N\L;

%estimation of the variance of the observation error
y_hat = A*x + b;
v_hat = y0 - y_hat;
V = v_hat';
T = Q\v_hat;
sigma02_hat = (V*T)/(n-m);

%covariance matrix
Cxx = sigma02_hat*(N^-1);