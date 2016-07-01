function  [Z,residual] = plane_fitting(x, y, z, xi, yi)

%planar fitting (least squares)
X = [ones(size(x)) x y]; %design matrix
a = X\z;                    %parameters

Z0 = X*a;                    %fitted values
% plot3(x,y,Z0,'go','MarkerSize',7)
%plot fitting plane
min_x1 = min(x);
max_x1 = max(x);
min_x2 = min(y);
max_x2 = max(y);
step = 10;
residual = z - Z0;
x = min_x1 : (max_x1 - min_x1)/step : max_x1;
y = min_x2 : (max_x2 - min_x2)/step : max_x2;
[X1,X2] = meshgrid(x,y);

Z = zeros(size(X1));

%evaluate the plane at the requested coordinates
X = [ones(size(xi)) xi yi];
Z = X*a;


