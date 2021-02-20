%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scmidt quasi normalized associated functions discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 14;

% co to colatitude
theta = -1:0.001:1;
P = zeros(n,n,length(theta));
for j = 1:length(theta)
    for i = 0:n-1;
        P(1:i+1,i+1,j) = legendre(i,theta(j),'sch');
    end
end

% discretize derivative
theta = -0.999:0.001:0.999;
Pp = zeros(n,n,length(theta));
Pm = zeros(n,n,length(theta));
d_cos = 0.0001;
for j = 1:length(theta)
    for i = 0:n-1;
        Pp(1:i+1,i+1,j) = legendre(i,cos(acos(theta(j))+d_cos/2),'sch');
        Pm(1:i+1,i+1,j) = legendre(i,cos(acos(theta(j))-d_cos/2),'sch');
    end
end
dP = (Pp - Pm)/d_cos;